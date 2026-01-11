#ifndef STAGE1_H
#define STAGE1_H
using namespace std;
#include "parm.h"
#include "MurmurHash3.h"
#include <cstring>
#include <vector>
#include <iostream>
// Stage1 bucket: quickly records flow continuity
struct Stage1Bucket {
    uint8_t continuity : 4; // Continuity window count (0-15), records flow arrivals in consecutive windows
    uint8_t arrival    : 1; // Recent window number (0/1), used for 0/1 alternation detection
    uint8_t jump       : 1; // Promotion flag (0/1), marks whether promoted to Stage2

    Stage1Bucket() : continuity(0), arrival(0), jump(0) {}

    bool empty() const {
        return continuity == 0 && arrival == 0 && jump == 0;
    }

    void reset() {
        continuity = 0;
        arrival = 0;
        jump = 0;
    }
};

// Stage1: detects candidate stable flows
class Stage1Filter {
private:
    vector<vector<Stage1Bucket>> buckets; // Multi-row hash table: d rows × m buckets
    vector<uint32_t> hashSeeds;           // Hash seed for each row, ensuring different hash functions per row

public:
    // Constructor: accepts memory parameter (bytes)
    explicit Stage1Filter(size_t memoryBytes = STAGE1_MEMORY_BYTES) {
        size_t rows = STAGE1_ROWS;
        size_t bucketSize = sizeof(Stage1Bucket);
        size_t perRowBytes = (rows > 0) ? (memoryBytes / rows) : 0;
        size_t bucketsPerRow = (bucketSize > 0) ? max<size_t>(1, perRowBytes / bucketSize) : 1;

        // Initialize multi-row hash table
        buckets.resize(rows, vector<Stage1Bucket>(bucketsPerRow));
        hashSeeds.resize(rows);
        
        // Generate different hash seeds for each row
        for (size_t i = 0; i < rows; i++) {
            hashSeeds[i] = 0x100 + (uint32_t)i * 0x123;
        }
    }

    // Flow arrives, returns whether promoted to Stage2
    bool processPacket(const char* flowID, uint32_t windowSeq) {
        uint8_t cur = windowSeq % 2; // Calculate current window number (0/1)
        bool allContinuity5 = true;  // Check if all rows reached continuity threshold

        // Case 1: Check if flow is already promoted in all rows
        bool allJumped = true;
        for (size_t i = 0; i < STAGE1_ROWS; i++) {
            uint32_t h = 0;
            MurmurHash3_x86_32(flowID, KEY_LEN, hashSeeds[i], &h);
            Stage1Bucket& b = buckets[i][h % buckets[i].size()];
            if (!b.jump) {
                allJumped = false;
                break;
            }
        }

        if (allJumped) {
            // Case 1: Flow already promoted, only update arrival field, other fields unchanged
            for (size_t i = 0; i < STAGE1_ROWS; i++) {
                uint32_t h = 0;
                MurmurHash3_x86_32(flowID, KEY_LEN, hashSeeds[i], &h);
                Stage1Bucket& b = buckets[i][h % buckets[i].size()];
                b.arrival = cur;
            }
            return true; // Promoted to Stage2
        }

        // Handle cases 2, 3, 4: Check flow continuity
        for (size_t i = 0; i < STAGE1_ROWS; i++) {
            uint32_t h = 0;
            MurmurHash3_x86_32(flowID, KEY_LEN, hashSeeds[i], &h);
            Stage1Bucket& b = buckets[i][h % buckets[i].size()];

            if (b.empty()) {
                // Case 2: Bucket empty, initialize continuity=1, set arrival
                b.continuity = 1;
                b.arrival = cur;
                allContinuity5 = false;
            } else if (b.arrival == cur) {
                // Case 3: Repeated arrival in same window, no update
                if (b.continuity != 15) {
                    allContinuity5 = false;
                }
            } else {
                // Case 4: Different window arrival, increment continuity count
                if (b.continuity < 15) {
                    b.continuity++;
                }
                // If continuity reached maximum, keep unchanged
                b.arrival = cur;
                if (b.continuity != 15) {
                    allContinuity5 = false;
                }
            }
        }

        // Check if all rows reached continuity threshold
        if (allContinuity5) {
            // Flow promotion: set jump flag in all rows
            for (size_t i = 0; i < STAGE1_ROWS; i++) {
                uint32_t h = 0;
                MurmurHash3_x86_32(flowID, KEY_LEN, hashSeeds[i], &h);
                Stage1Bucket& b = buckets[i][h % buckets[i].size()];
                b.jump = 1;
            }

            return true;
        }
        return false;
    }

    // Reset buckets not present in current window
    void resetBuckets(uint32_t windowSeq) {
        uint8_t cur = windowSeq % 2;
        for (auto& row : buckets) {
            for (auto& b : row) {
                if (b.empty()) continue;
                // Reset buckets not accessed in current window
                if (b.arrival != cur) b.reset();
            }
        }
    }
};

#endif