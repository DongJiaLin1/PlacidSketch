#ifndef STAGE1_H
#define STAGE1_H

#include "parm.h"
#include "MurmurHash3.h"
#include <vector>

struct Stage1Bucket {
    uint8_t continuity : 4;
    uint8_t arrival    : 1;
    uint8_t jump       : 1;

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

class Stage1Filter {
private:
    std::vector<std::vector<Stage1Bucket>> buckets;
    std::vector<uint32_t> hashSeeds;

public:
    explicit Stage1Filter(size_t memoryBytes = STAGE1_MEMORY_BYTES) {
        size_t rows = STAGE1_ROWS;
        size_t bucketSize = sizeof(Stage1Bucket);
        size_t perRowBytes = (rows > 0) ? (memoryBytes / rows) : 0;
        size_t bucketsPerRow = (bucketSize > 0) ? std::max<size_t>(1, perRowBytes / bucketSize) : 1;

        buckets.resize(rows, std::vector<Stage1Bucket>(bucketsPerRow));
        hashSeeds.resize(rows);

        for (size_t i = 0; i < rows; i++) {
            hashSeeds[i] = 0x100 + (uint32_t)i * 0x123;
        }
    }

    bool processPacket(const char* flowID, uint32_t windowSeq) {
        uint8_t cur = windowSeq % 2;
        bool allContinuity5 = true;

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
            for (size_t i = 0; i < STAGE1_ROWS; i++) {
                uint32_t h = 0;
                MurmurHash3_x86_32(flowID, KEY_LEN, hashSeeds[i], &h);
                Stage1Bucket& b = buckets[i][h % buckets[i].size()];
                b.arrival = cur;
            }
            return true;
        }

        for (size_t i = 0; i < STAGE1_ROWS; i++) {
            uint32_t h = 0;
            MurmurHash3_x86_32(flowID, KEY_LEN, hashSeeds[i], &h);
            Stage1Bucket& b = buckets[i][h % buckets[i].size()];

            if (b.empty()) {
                b.continuity = 1;
                b.arrival = cur;
                allContinuity5 = false;
            } else if (b.arrival == cur) {
                if (b.continuity != 15) {
                    allContinuity5 = false;
                }
            } else {
                if (b.continuity < 15) {
                    b.continuity++;
                }
                b.arrival = cur;
                if (b.continuity != 15) {
                    allContinuity5 = false;
                }
            }
        }

        if (allContinuity5) {
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

    void resetBuckets(uint32_t windowSeq) {
        uint8_t cur = windowSeq % 2;
        for (auto& row : buckets) {
            for (auto& b : row) {
                if (b.empty()) continue;
                if (b.arrival != cur) b.reset();
            }
        }
    }
};

#endif
