#ifndef STAGE2_H
#define STAGE2_H

#include "parm.h"
#include "stage3.h"
#include "MurmurHash3.h"
#include <algorithm>
#include <array>
#include <cstring>
#include <limits>
#include <numeric>
#include <vector>

struct Stage2Bucket {
    std::array<uint8_t, SUBFLOW_WINDOWS + 1> cx{};
    uint8_t initialized_flags;
    uint8_t ck1;
    uint8_t ck2;
    uint8_t ck1_is_null : 1;
    uint8_t ck2_is_null : 1;

    Stage2Bucket() : ck1(1), ck2(1), ck1_is_null(0), ck2_is_null(0), initialized_flags(0) {
        cx.fill(0);
    }

    bool empty() const { return initialized_flags == 0; }

    bool isCounterNull(uint8_t index) const {
        return (initialized_flags & (1U << index)) == 0;
    }

    void reset() {
        cx.fill(0);
        initialized_flags = 0;
        ck1 = 1;
        ck2 = 1;
        ck1_is_null = 0;
        ck2_is_null = 0;
    }

    uint32_t countWindowNumber() {
        uint32_t num = 0;
        for (uint32_t i = 0; i < SUBFLOW_WINDOWS + 1; i++) {
            if (initialized_flags & (1U << i)) {
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

    static uint32_t calculateConsecutiveWindows(const Stage2Bucket& bucket, uint32_t currentWindow) {
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

    static float calculateVariance(const std::vector<float>& data) {
        if (data.size() < 2) {
            return std::numeric_limits<float>::infinity();
        }
        float sum = std::accumulate(data.begin(), data.end(), 0.0f);
        float mean = sum / data.size();
        float variance = 0.0f;
        for (float val : data) {
            variance += (val - mean) * (val - mean);
        }
        return variance / (data.size() - 1);
    }

    static float calculateDirectVariance(const Stage2Bucket& bucket, uint32_t startWindow) {
        const uint32_t R = SUBFLOW_WINDOWS + 1;
        std::vector<float> directFreqs;
        directFreqs.reserve(SUBFLOW_WINDOWS);
        for (uint32_t i = 0; i < SUBFLOW_WINDOWS; ++i) {
            uint32_t windowIndex = (startWindow + i) % R;
            if (bucket.isCounterNull(static_cast<uint8_t>(windowIndex))) {
                return std::numeric_limits<float>::infinity();
            }
            directFreqs.push_back(static_cast<float>(bucket.cx[windowIndex]));
        }
        if (directFreqs.size() < 2) {
            return std::numeric_limits<float>::infinity();
        }
        return calculateVariance(directFreqs);
    }

    static float calculateMeanFrequency(const Stage2Bucket& bucket, uint32_t startWindow) {
        const uint32_t R = SUBFLOW_WINDOWS + 1;
        float sum = 0.0f;
        uint32_t count = 0;
        for (uint32_t i = 0; i < SUBFLOW_WINDOWS; ++i) {
            uint32_t windowIndex = (startWindow + i) % R;
            if (bucket.isCounterNull(static_cast<uint8_t>(windowIndex))) {
                return std::numeric_limits<float>::infinity();
            }
            sum += static_cast<float>(bucket.cx[windowIndex]);
            count++;
        }
        return count > 0 ? sum / count : 0.0f;
    }

    static bool updateCKOnRebirth(Stage2Bucket* bucket, uint32_t absoluteWindow) {
        bool useCk1 = (absoluteWindow % 2 == 0);
        if (useCk1) {
            if (!bucket->ck1_is_null && bucket->ck1 < 6) {
                bucket->ck1++;
            }
            uint32_t windowNum = bucket->countWindowNumber();
            if (windowNum != 1) {
                if (!bucket->ck2_is_null) {
                    if (bucket->ck2 != 0) {
                        bucket->ck2--;
                    } else {
                        bucket->ck2_is_null = true;
                        return false;
                    }
                } else {
                    return false;
                }
            }
        } else {
            if (!bucket->ck2_is_null && bucket->ck2 < 6) {
                bucket->ck2++;
            }
            uint32_t windowNum = bucket->countWindowNumber();
            if (windowNum != 1) {
                if (!bucket->ck1_is_null) {
                    if (bucket->ck1 != 0) {
                        bucket->ck1--;
                    } else {
                        bucket->ck1_is_null = true;
                        return false;
                    }
                } else {
                    return false;
                }
            }
        }
        return true;
    }
};

class Stage2Monitor {
private:
    std::vector<std::vector<Stage2Bucket>> buckets;
    std::vector<uint32_t> rowSeeds;
    uint32_t hashSeed;
    Stage3Merger& stage3;
    size_t rows = 0;
    size_t bucketsPerRow = 0;
    uint32_t currentReportWindow = UINT32_MAX;
    std::vector<uint32_t> reportedFlows;

    static inline bool isPowerOfTwo(uint32_t v) {
        return v && ((v & (v - 1)) == 0);
    }

    static inline uint32_t fastReduce(uint32_t hash, uint32_t range) {
        return static_cast<uint32_t>((static_cast<uint64_t>(hash) * static_cast<uint64_t>(range)) >> 32);
    }

    inline uint32_t indexForRow(const char* flowID, uint32_t row, uint32_t range) const {
        uint32_t h;
        MurmurHash3_x86_32(flowID, KEY_LEN, rowSeeds[row], &h);
        h ^= (row * 0x9e3779b9u) ^ 0x85ebca6bu;
        if (isPowerOfTwo(range)) {
            return h & (range - 1);
        } else {
            return fastReduce(h, range);
        }
    }

public:
    explicit Stage2Monitor(Stage3Merger& s3, size_t memoryBytes = STAGE2_MEMORY_BYTES)
        : hashSeed(0x200), stage3(s3) {
        rows = STAGE2_ROWS;
        size_t cellSize = sizeof(Stage2Bucket);
        size_t perRowBytes = (rows > 0) ? (memoryBytes / rows) : 0;
        bucketsPerRow = (cellSize > 0) ? std::max<size_t>(1, perRowBytes / cellSize) : 1;

        buckets.resize(rows);
        for (size_t i = 0; i < rows; ++i) {
            buckets[i].resize(bucketsPerRow);
        }

        rowSeeds.resize(rows);
        for (size_t i = 0; i < rows; ++i) {
            rowSeeds[i] = hashSeed + static_cast<uint32_t>(i * 131);
        }
    }

    void processPacket(const char* flowID, uint32_t currentWindow) {
        uint8_t curSlot = static_cast<uint8_t>(currentWindow % (SUBFLOW_WINDOWS + 1));
        for (uint32_t row = 0; row < rows; ++row) {
            uint32_t idx = indexForRow(flowID, row, static_cast<uint32_t>(bucketsPerRow));
            Stage2Bucket& bucket = buckets[row][idx];

            if (bucket.empty()) {
                bucket.initializeNewWindow(curSlot, currentWindow);
            } else {
                if (bucket.isCounterNull(curSlot)) {
                    if (Stage2Bucket::calculateConsecutiveWindows(bucket, currentWindow) == 0) {
                        Stage2Bucket::updateCKOnRebirth(&bucket, currentWindow);
                    }
                    bucket.initializeNewWindow(curSlot, currentWindow);
                } else {
                    if (bucket.cx[curSlot] < (1u << COUNTER_BITS) - 1u) {
                        bucket.cx[curSlot]++;
                    }
                }
            }
        }
    }

    bool processPotentialFlow(const char* flowID, uint32_t currentWindow) {
        if ((currentWindow + 1) % SUBFLOW_WINDOWS != 0) {
            return false;
        }

        if (currentWindow != currentReportWindow) {
            currentReportWindow = currentWindow;
            reportedFlows.clear();
        }

        uint32_t flowValue = 0;
        std::memcpy(&flowValue, flowID, sizeof(uint32_t));
        if (std::find(reportedFlows.begin(), reportedFlows.end(), flowValue) != reportedFlows.end()) {
            return false;
        }

        bool found = false;
        float finalMean = 0.0f;
        float finalVariance = 0.0f;
        Stage2Bucket* selectedBucket = nullptr;

        for (uint32_t row = 0; row < rows; ++row) {
            uint32_t idx = indexForRow(flowID, row, static_cast<uint32_t>(bucketsPerRow));
            Stage2Bucket& bucket = buckets[row][idx];

            if (!bucket.empty()) {
                uint32_t consecutiveWindows = Stage2Bucket::calculateConsecutiveWindows(bucket, currentWindow);
                if (consecutiveWindows >= SUBFLOW_WINDOWS) {
                    uint32_t startWindow = currentWindow + 1 - SUBFLOW_WINDOWS;
                    float mean = Stage2Bucket::calculateMeanFrequency(bucket, startWindow);
                    float var = Stage2Bucket::calculateDirectVariance(bucket, startWindow);
                    if (var <= STABLE_THRESHOLD) {
                        found = true;
                        finalMean = mean;
                        finalVariance = var;
                        selectedBucket = &bucket;
                        break;
                    }
                }
            }
        }

        if (found && selectedBucket) {
            reportedFlows.push_back(flowValue);
            stage3.processSteadySubflow(flowValue, currentWindow + 1 - SUBFLOW_WINDOWS, finalVariance, finalMean);
        }

        return found;
    }
};

#endif
