#ifndef STAGE3_H
#define STAGE3_H
using namespace std;
#include "parm.h"
#include "MurmurHash3.h"
#include <random>
#include <string>
#include <vector>

// Statistics structure: stores mean frequency and frequency variance of stable flows
struct Statistics {
    float mean = 0.0f;
    float variance = 0.0f;
    Statistics() = default;
    Statistics(float m, float v) : mean(m), variance(v) {}
};

// Stage3Cell: stores merged information of a stable flow
struct Stage3Cell {
    char ID[KEY_LEN]{};
    uint16_t window = 0;
    Statistics s;
    uint16_t number = 0;

    bool empty() const { return ID[0] == 0; }
    void clear() {
        fill(begin(ID), end(ID), 0);
        window = 0;
        s = Statistics();
        number = 0;
    }
};

// Stage3: stable subflow merger
class Stage3Merger {
private:
    vector<vector<Stage3Cell>> buckets;
    uint32_t hashSeed;
    mt19937 gen;
    uniform_real_distribution<float> dist;
    size_t l = 0;
    size_t b = 0;

    // Check if new subflow can be merged: incremental variance calculation
    static bool canMergeVariance(const Stage3Cell& cell, float newVar, float newMean) {
        const uint32_t C = cell.number;
        const float mean_star = cell.s.mean;
        const float var_star = cell.s.variance;
        const float mean = newMean;
        const float var = newVar;

        // Calculate merged mean
        const float mu_star = (C * mean_star + mean) / (C + 1);
        // Calculate merged variance using incremental formula
        const float term1 = C * (var_star + (mean_star - mu_star) * (mean_star - mu_star)) / (C + 1);
        const float term2 = (var + (mean - mu_star) * (mean - mu_star)) / (C + 1);
        const float V_star = term1 + term2;

        return V_star <= STABLE_THRESHOLD;
    }

    void clearCell(Stage3Cell& cell) {
        if (!cell.empty() && cell.number >= 1) {
            if (cell.number >= Q) {
                string flowIDStr(cell.ID, strnlen(cell.ID, KEY_LEN));
                float V_star = cell.s.variance;

                if (V_star > STABLE_THRESHOLD) {
                    cell.clear();
                    return;
                }

                uint32_t endWindow = cell.window + cell.number * MIN_SUBFLOWS - 1;
            }
        }
        cell.clear();
    }

    static void initNewCell(Stage3Cell& cell, const char* flowID, uint32_t startW, float var, float mean) {
        memcpy(cell.ID, flowID, KEY_LEN);
        cell.window = startW;
        cell.s = Statistics(mean, var);
        cell.number = 1;
    }

    // Merge new subflow into cell: update mean and variance incrementally
    static void mergeCell(Stage3Cell& cell, float newVar, float newMean) {
        if (cell.number < static_cast<uint32_t>(P)) {
            const uint32_t C = cell.number;
            const float mean_star = cell.s.mean;
            const float var_star = cell.s.variance;
            const float mean = newMean;
            const float var = newVar;
            const float mu_star = (C * mean_star + mean) / (C + 1);

            // Update variance using incremental merging formula
            const float term1 = C * (var_star + (mean_star - mu_star) * (mean_star - mu_star)) / (C + 1);
            const float term2 = (var + (mean - mu_star) * (mean - mu_star)) / (C + 1);
            const float V_star = term1 + term2;

            cell.number++;
            cell.s.mean = mu_star;
            cell.s.variance = V_star;
        }
    }

public:
    explicit Stage3Merger(): hashSeed(0x300),gen(random_device{}()),dist(0.0f, 1.0f)
    {
        l = STAGE3_BUCKETS;
        size_t cellSize = sizeof(Stage3Cell);
        size_t perBucketBytes = (l > 0) ? (STAGE3_MEMORY_BYTES / l) : 0;
        b = (cellSize > 0) ? max<size_t>(1, perBucketBytes / cellSize) : 1;

        buckets.resize(l);
        for (size_t i = 0; i < l; ++i) {
            buckets[i].resize(b);
        }
    }

    ~Stage3Merger() {
        finalize();
    }

    // Process stable subflow: merge or insert based on bucket state
    void processSteadySubflow(const char* flowID, uint32_t startW, float var, float mean) {
        uint32_t h = 0;
        MurmurHash3_x86_32(flowID, KEY_LEN, hashSeed, &h);
        size_t u = h % l;
        auto& bucket = buckets[u];

        Stage3Cell* targetCell = nullptr;
        int emptyIndex = -1;
        vector<int> discontinuousIndices;

        // Scan bucket to find matching cell, empty slot, and discontinuous cells
        for (int a = 0; a < static_cast<int>(b); ++a) {
            if (bucket[a].empty()) {
                if (emptyIndex < 0) emptyIndex = a;
                continue;
            }
            if (memcmp(bucket[a].ID, flowID, KEY_LEN) == 0) {
                targetCell = &bucket[a];
            } else {
                uint32_t lastwin = bucket[a].window + bucket[a].number * MIN_SUBFLOWS;
                if (startW != lastwin) {
                    discontinuousIndices.push_back(a);
                }
            }
        }

        // Case 1: Empty slot available
        if (!targetCell && emptyIndex != -1) {
            initNewCell(bucket[emptyIndex], flowID, startW, var, mean);
        }
        // Case 2: Matching cell found
        else if (targetCell) {
            uint32_t lastwin = targetCell->window + targetCell->number * MIN_SUBFLOWS;

            if (startW != lastwin) {
                // Window discontinuity: report and reset
                clearCell(*targetCell);
                initNewCell(*targetCell, flowID, startW, var, mean);
            } else {
                // Continuous windows: try to merge
                if (canMergeVariance(*targetCell, var, mean)) {
                    mergeCell(*targetCell, var, mean);
                    if (targetCell->number >= static_cast<uint32_t>(P)) {
                        // Max segments reached: report and reset
                        clearCell(*targetCell);
                        initNewCell(*targetCell, flowID, startW, var, mean);
                    }
                } else {
                    // Merge failed: report and reset
                    clearCell(*targetCell);
                    initNewCell(*targetCell, flowID, startW, var, mean);
                }
            }
        }
        // Case 3: No match, no empty slot
        else if (!targetCell) {
            if (!discontinuousIndices.empty()) {
                // Replace discontinuous cell (prefer smallest number)
                int victimIndex = discontinuousIndices[0];
                for (int idx : discontinuousIndices) {
                    if (bucket[idx].number < bucket[victimIndex].number) {
                        victimIndex = idx;
                    }
                }
                clearCell(bucket[victimIndex]);
                initNewCell(bucket[victimIndex], flowID, startW, var, mean);
            } else {
                // All cells continuous: probabilistic replacement
                int victimIndex = -1;
                uint32_t minStableWindows = UINT32_MAX;
                
                for (int a = 0; a < static_cast<int>(b); ++a) {
                    uint32_t stableWindows = bucket[a].number;
                    if (stableWindows < minStableWindows) {
                        minStableWindows = stableWindows;
                        victimIndex = a;
                    }
                }
                
                if (victimIndex != -1) {
                    // Replacement probability: 1/(wmin-s+1)
                    uint32_t totalStableWindows = minStableWindows * MIN_SUBFLOWS;
                    float replaceProb = 1.0f / max(1.0f, static_cast<float>(totalStableWindows - MIN_SUBFLOWS + 1));
                    if (dist(gen) <= replaceProb) {
                        clearCell(bucket[victimIndex]);
                        initNewCell(bucket[victimIndex], flowID, startW, var, mean);
                    }
                }
            }
        }
    }

    void finalize() {
        for (auto& bucket : buckets) {
            for (auto& cell : bucket) {
                clearCell(cell);
            }
        }
    }
};

#endif