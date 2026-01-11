#ifndef STAGE2_H
#define STAGE2_H
using namespace std;
#include "parm.h"
#include "stage3.h"
#include "MurmurHash3.h"
#include <numeric>
#include <cstring>
#include <string>
#include <array>
// Stage2 bucket structure
struct Stage2Bucket {
    array<uint8_t, SUBFLOW_WINDOWS + 1> cx{};
    uint8_t initialized_flags;
    uint8_t ck1: 3;
    uint8_t ck2: 3;
    uint8_t ck1_is_null: 1;
    uint8_t ck2_is_null: 1;

    Stage2Bucket() : ck1(1), ck2(1), ck1_is_null(false), ck2_is_null(false), initialized_flags(0) {
        cx.fill(0);
    }

    bool empty() const { return initialized_flags == 0; }

    bool isCounterNull(uint8_t index) const {
        return !(initialized_flags & (1U << index));
    }

    void reset() {
        cx.fill(0);
        initialized_flags = 0;
        ck1 = 1;
        ck2 = 1;
        ck1_is_null = false;
        ck2_is_null = false;
    }

    uint32_t countWindowNumber(){
        uint32_t num = 0;
        for(uint32_t i=0;i<SUBFLOW_WINDOWS + 1;i++){
            if(initialized_flags & (1U << i)){
                num++;
            }
        }
        return num;
    }

    void initializeNewWindow(uint8_t window, uint32_t absoluteWindow) {
        cx[window] = 1;
        initialized_flags |= (1U << window);

        bool useCk1 = (absoluteWindow % 2 == 0);
        if (useCk1) {
            ck1 = 1;
            ck1_is_null = false;
        } else {
            ck2 = 1;
            ck2_is_null = false;
        }
    }


    // Check stability using relative rebirth algorithm with alternating CK fields
    bool checkStability(uint8_t y1, uint8_t y2, uint32_t absoluteWindow) const {
        const uint32_t base = (1u << COUNTER_BITS);
        const uint32_t cx1 = cx[y1];
        const uint32_t cx2 = cx[y2];

        bool useCk1 = (absoluteWindow % 2 == 0);
        if (isCounterNull(y1) || isCounterNull(y2)) {
            return false;
        }
        if (useCk1) {
            if (ck1_is_null) {
                return false;
            }
            if (ck1 > 2) {
                return false;
            } else if (ck1 == 2) {
                return (cx1 + base - cx2) <= ALPHA_THRESHOLD;
            } else if (ck1 == 1) {
                return abs(int(cx2) - int(cx1)) <= ALPHA_THRESHOLD;
            } else {
                return (cx2 + base - cx1) <= ALPHA_THRESHOLD;
            }
        } else {
            if (ck2_is_null) {
                return false;
            }
            if (ck2 > 2) {
                return false;
            } else if (ck2 == 2) {
                return (cx1 + base - cx2) <= ALPHA_THRESHOLD;
            } else if (ck2 == 1) {
                return abs(int(cx2) - int(cx1)) <= ALPHA_THRESHOLD;
            } else {
                return (cx1 + base - cx2) <= ALPHA_THRESHOLD;
            }
        }

    }
};

// Stage2 monitor: monitors flows for stability
class Stage2Monitor {
private:
    vector<vector<Stage2Bucket>> buckets;
    vector<uint32_t> rowSeeds;
    uint32_t hashSeed;
    Stage3Merger &stage3;
    size_t rows = 0;
    size_t bucketsPerRow = 0;

    struct SelectedBucket {
        Stage2Bucket *bucket;
    };

    static inline bool isPowerOfTwo(uint32_t v) {
        return v && ((v & (v - 1)) == 0);
    }

    static inline uint32_t fastReduce(uint32_t hash, uint32_t range) {
        return static_cast<uint32_t>((static_cast<uint64_t>(hash) * static_cast<uint64_t>(range)) >> 32);
    }

    inline uint32_t computeRowHash(const char* flowID, uint32_t row) const {
        uint32_t h;
        MurmurHash3_x86_32(flowID, KEY_LEN, rowSeeds[row], &h);
        h ^= (row * 0x9e3779b9u) ^ 0x85ebca6bu;
        return h;
    }

    inline uint32_t indexForRow(const char* flowID, uint32_t row, uint32_t range) const {
        uint32_t h = computeRowHash(flowID, row);
        if (isPowerOfTwo(range)) {
            return h & (range - 1);
        } else {
            return fastReduce(h, range);
        }
    }
    static uint32_t calculateConsecutiveWindows(const Stage2Bucket &bucket, uint32_t currentWindow) {
        const uint32_t R = SUBFLOW_WINDOWS + 1;
        uint32_t count = 0;
        for (uint32_t i = 1; i <= SUBFLOW_WINDOWS; ++i) {
            if (currentWindow < i) break;
            uint32_t window = currentWindow - i;
            uint8_t y = window % R;
            if (bucket.isCounterNull(y)) break;
            count++;
        }
        return count;
    }

    static float calculateVariance(const vector<float> &data) {
        if (data.size() < 2) {
            return numeric_limits<float>::infinity();
        }

        float sum = accumulate(data.begin(), data.end(), 0.0f);
        float mean = sum / data.size();

        float variance = 0.0f;
        for (float val: data) {
            variance += (val - mean) * (val - mean);
        }
        return variance / (data.size() - 1);
    }

    static float calculateDirectVariance(const Stage2Bucket &bucket, uint32_t startWindow) {
        const uint32_t R = SUBFLOW_WINDOWS + 1;
        vector<float> directFreqs;
        directFreqs.reserve(SUBFLOW_WINDOWS);

        for (uint32_t i = 0; i < SUBFLOW_WINDOWS; ++i) {
            uint32_t windowIndex = (startWindow + i) % R;
            if (bucket.isCounterNull(static_cast<uint8_t>(windowIndex))) {
                return numeric_limits<float>::infinity();
            }
            uint8_t counter = bucket.cx[windowIndex];
            directFreqs.push_back(static_cast<float>(counter));
        }

        if (directFreqs.size() < 2) {
            return numeric_limits<float>::infinity();
        }
        return calculateVariance(directFreqs);
    }

    // Calculate offset variance to handle counter wrap-around
    static float calculateOffsetVariance(const Stage2Bucket &bucket, uint32_t startWindow) {
        const uint32_t R = SUBFLOW_WINDOWS + 1;
        const uint32_t base = (1u << COUNTER_BITS);
        const uint32_t half = base >> 1;
        vector<float> offsetFreqs;
        offsetFreqs.reserve(SUBFLOW_WINDOWS);

        for (uint32_t i = 0; i < SUBFLOW_WINDOWS; ++i) {
            uint32_t windowIndex = (startWindow + i) % R;
            if (bucket.isCounterNull(static_cast<uint8_t>(windowIndex))) break;
            uint32_t v = bucket.cx[windowIndex];
            uint32_t adj = (v + half) % base;  // Offset adjustment for wrap-around
            offsetFreqs.push_back(static_cast<float>(adj));
        }
        
        if (offsetFreqs.size() < 2) {
            return numeric_limits<float>::infinity();
        }
        return calculateVariance(offsetFreqs);
    }

    static float calculateMeanFrequency(const Stage2Bucket &bucket, uint32_t startWindow) {
        const uint32_t R = SUBFLOW_WINDOWS + 1;
        float sum = 0.0f;
        uint32_t count = 0;

        for (uint32_t i = 0; i < SUBFLOW_WINDOWS; ++i) {
            uint32_t windowIndex = (startWindow + i) % R;
            if (bucket.isCounterNull(static_cast<uint8_t>(windowIndex))) {
                return numeric_limits<float>::infinity();
            }
            sum += static_cast<float>(bucket.cx[windowIndex]);
            count++;
        }

        return count > 0 ? sum / count : 0.0f;
    }

    // Update CK fields on counter rebirth: increment current CK, decrement previous CK
    static bool updateCKOnRebirth(Stage2Bucket* bucket, uint32_t absoluteWindow) {
        bool useCk1 = (absoluteWindow % 2 == 0);

        if (useCk1) {
            if (!bucket->ck1_is_null && bucket->ck1 < 6) {
                bucket->ck1++;  // Increment ck1 on rebirth
            }
            uint32_t windowNum = bucket->countWindowNumber();
            if (windowNum != 1){
                if (!bucket->ck2_is_null) {
                    if (bucket->ck2 != 0) {
                        bucket->ck2--;  // Decrement ck2 from previous window
                    } else {
                        bucket->ck2_is_null = true;
                        return false;
                    }
                }
                else{
                    return false;
                }
            }

        } else {
            if (!bucket->ck2_is_null && bucket->ck2 < 6) {
                bucket->ck2++;
            }
            uint32_t windowNum = bucket->countWindowNumber();
            if (windowNum != 1){
                if (!bucket->ck1_is_null) {
                    if (bucket->ck1 != 0) {
                        bucket->ck1--;
                    } else {
                        bucket->ck1_is_null = true;
                        return false;
                    }
                }
                else{
                    return false;
                }
            }

        }
        return true;
    }

public:
    explicit Stage2Monitor(Stage3Merger &s3, size_t memoryBytes = STAGE2_MEMORY_BYTES) : hashSeed(0x200), stage3(s3) {
        rows = STAGE2_ROWS;
        size_t bucketSize = sizeof(Stage2Bucket);
        size_t perRowBytes = (rows > 0) ? (memoryBytes / rows) : 0;
        bucketsPerRow = (bucketSize > 0) ? max<size_t>(1, perRowBytes / bucketSize) : 1;

        buckets.resize(rows);
        for (size_t i = 0; i < rows; ++i) {
            buckets[i].resize(bucketsPerRow);
        }
        rowSeeds.resize(rows);

        for (size_t i = 0; i < rows; ++i) {
            rowSeeds[i] = hashSeed + static_cast<uint32_t>(i);
        }
    }

    void processPotentialFlow(const char* flowID, uint32_t currentWindow) {
        const uint32_t R = SUBFLOW_WINDOWS + 1;
        const uint8_t y_current = currentWindow % R;
        const uint8_t y_prev = (currentWindow - 1) % R;
        const uint8_t y_prev_prev = (currentWindow - 2) % R;

        vector<SelectedBucket> selected;

        for (uint32_t f = 0; f < rows; ++f) {
            uint32_t k = indexForRow(flowID, f, static_cast<uint32_t>(bucketsPerRow));
            Stage2Bucket* cell = &buckets[f][k];
            selected.push_back(SelectedBucket{cell});
        }

        bool hasEmpty = false;
        bool hasCurrentNull = false;

        for (auto& SelectBucket : selected) {
            if (!hasEmpty && SelectBucket.bucket->empty()) {
                hasEmpty = true;
            }
            if (!hasCurrentNull && SelectBucket.bucket->isCounterNull(y_current)) {
                hasCurrentNull = true;
            }
        }
        if (hasEmpty) {
            for (auto& SelectBucket : selected) {
                if (SelectBucket.bucket->empty()) {
                    SelectBucket.bucket->initializeNewWindow(y_current, currentWindow);
                }
            }
            return;
        }

        if (!hasCurrentNull) {
            for (auto& SelectBucket : selected) {
                Stage2Bucket* bucketPtr = SelectBucket.bucket;
                if (bucketPtr->cx[y_current] < ((1u << COUNTER_BITS) - 1)) {
                    bucketPtr->cx[y_current]++;
                } else {
                    bucketPtr->cx[y_current] = 0;
                    bool flag5 = updateCKOnRebirth(bucketPtr, currentWindow);
                    if (!flag5) {
                        bucketPtr->reset();
                        bucketPtr->initializeNewWindow(y_current, currentWindow);
                    }
                }
            }
            return;
        }
        vector<Stage2Bucket*> nullBuckets;
        for (auto &SelectedBucket : selected) {
            if (SelectedBucket.bucket->isCounterNull(y_current)) {
                nullBuckets.push_back(SelectedBucket.bucket);
            }
        }
        if (nullBuckets.empty()) return;

        bool havepassed = false;

        for (auto *b : nullBuckets) {
            uint32_t windowNum = b->countWindowNumber();
            if (windowNum > 2 && (b->isCounterNull(y_prev) || b->isCounterNull(y_prev_prev))) {
                b->reset();
                b->initializeNewWindow(y_current, currentWindow);
            }
            else if (windowNum <= 2 && b->isCounterNull(y_prev) ){
                b->reset();
                b->initializeNewWindow(y_current, currentWindow);
            }
            else if (windowNum <= 2 && !b->isCounterNull(y_prev) ){
                b->initializeNewWindow(y_current, currentWindow);
            }
            else {
                // Check stability using relative rebirth algorithm
                if (!b->checkStability(y_prev, y_prev_prev, currentWindow)) {
                    b->reset();

                    b->initializeNewWindow(y_current, currentWindow);
                }
                else {
                    uint32_t consecutiveWindows = calculateConsecutiveWindows(*b, currentWindow);
                    if (consecutiveWindows == SUBFLOW_WINDOWS) {
                        if(!havepassed){
                            // Calculate mean and variance, use minimum of direct and offset variance
                            uint32_t w = currentWindow - SUBFLOW_WINDOWS;
                            float meanFreq = calculateMeanFrequency(*b, w);
                            float varDirect = calculateDirectVariance(*b, w);
                            float varOffset = calculateOffsetVariance(*b, w);
                            float variance = min(varDirect, varOffset);
                            
                            // Pass stable subflow to Stage3 if variance is below threshold
                            if (variance <= STABLE_THRESHOLD) {
                                stage3.processSteadySubflow(flowID, w, meanFreq, variance);
                                havepassed = true;
                            }
                        }

                        b->reset();
                        b->initializeNewWindow(y_current, currentWindow);
                    }
                    else {
                        b->initializeNewWindow(y_current, currentWindow);
                    }
                }
            }
        }
        return;
    }
};
#endif