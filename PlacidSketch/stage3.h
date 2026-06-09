#ifndef STAGE3_H
#define STAGE3_H

#include "parm.h"
#include "MurmurHash3.h"
#include <cmath>
#include <random>
#include <unordered_set>
#include <vector>

struct Stage3Cell {
    uint32_t flowID;
    uint16_t window;
    uint16_t number;
    uint8_t variance_scaled;
    uint8_t mean_scaled;

    Stage3Cell() { clear(); }

    bool empty() const { return flowID == 0; }
    void clear() {
        flowID = 0;
        window = 0;
        mean_scaled = 0;
        variance_scaled = 0;
        number = 0;
    }
    void setStatistics(float m, float v) {
        mean_scaled = static_cast<uint8_t>(std::roundf(m * 10.0f));
        variance_scaled = static_cast<uint8_t>(std::roundf(v * 10.0f));
    }
    float getMean() const { return mean_scaled / 10.0f; }
    float getVariance() const { return variance_scaled / 10.0f; }
};

class Stage3Merger {
private:
    std::vector<std::vector<Stage3Cell>> buckets;
    uint32_t hashSeed;
    std::mt19937 gen;
    std::uniform_real_distribution<float> dist;
    size_t l = 0;
    size_t b = 0;
    std::vector<uint32_t>* detectedFlowsCallback = nullptr;
    std::unordered_set<uint32_t> reportedFlows;

    static uint32_t getDurationWindows(const Stage3Cell& cell) {
        return static_cast<uint32_t>(cell.number) * SUBFLOW_WINDOWS;
    }

    void tryReportStableFlow(Stage3Cell& cell) {
        if (!detectedFlowsCallback || cell.empty()) {
            return;
        }
        if (reportedFlows.count(cell.flowID)) {
            return;
        }
        if (getDurationWindows(cell) >= static_cast<uint32_t>(Steady_FLOWS)) {
            reportedFlows.insert(cell.flowID);
            detectedFlowsCallback->push_back(cell.flowID);
        }
    }

    static bool canMergeVariance(const Stage3Cell& cell, float newVar, float newMean) {
        const uint32_t C = cell.number;
        const float mean_star = cell.getMean();
        const float var_star = cell.getVariance();
        const float mean = newMean;
        const float var = newVar;

        const float mu_star = (C * mean_star + mean) / (C + 1);
        const float term1 = C * (var_star + (mean_star - mu_star) * (mean_star - mu_star)) / (C + 1);
        const float term2 = (var + (mean - mu_star) * (mean - mu_star)) / (C + 1);
        const float V_star = term1 + term2;

        return V_star <= STABLE_THRESHOLD;
    }

    void clearCell(Stage3Cell& cell) {
        if (!cell.empty() && cell.number >= Q) {
            float V_star = cell.getVariance();

            if (V_star > STABLE_THRESHOLD) {
                cell.clear();
                return;
            }
        }
        cell.clear();
    }

    static void initNewCell(Stage3Cell& cell, uint32_t flowID, uint32_t startW, float var, float mean) {
        cell.flowID = flowID;
        cell.window = static_cast<uint16_t>(startW);
        cell.setStatistics(mean, var);
        cell.number = 1;
    }

    static void mergeCell(Stage3Cell& cell, float newVar, float newMean) {
        if (cell.number < static_cast<uint32_t>(P)) {
            const uint32_t C = cell.number;
            const float mean_star = cell.getMean();
            const float var_star = cell.getVariance();
            const float mean = newMean;
            const float var = newVar;
            const float mu_star = (C * mean_star + mean) / (C + 1);

            const float term1 = C * (var_star + (mean_star - mu_star) * (mean_star - mu_star)) / (C + 1);
            const float term2 = (var + (mean - mu_star) * (mean - mu_star)) / (C + 1);
            const float V_star = term1 + term2;

            cell.number++;
            cell.setStatistics(mu_star, V_star);
        }
    }

public:
    explicit Stage3Merger() : hashSeed(0x300), gen(std::random_device{}()), dist(0.0f, 1.0f) {
        l = STAGE3_BUCKETS;
        size_t cellSize = sizeof(Stage3Cell);
        size_t perBucketBytes = (l > 0) ? (STAGE3_MEMORY_BYTES / l) : 0;
        b = (cellSize > 0) ? std::max<size_t>(1, perBucketBytes / cellSize) : 1;

        buckets.resize(l);
        for (size_t i = 0; i < l; ++i) {
            buckets[i].resize(b);
        }
    }

    ~Stage3Merger() {
        finalize();
    }

    void setDetectedFlowsCallback(std::vector<uint32_t>* callback) {
        detectedFlowsCallback = callback;
    }

    void processSteadySubflow(uint32_t flowID, uint32_t startW, float var, float mean) {
        uint32_t h = 0;
        MurmurHash3_x86_32(&flowID, sizeof(uint32_t), hashSeed, &h);
        size_t u = h % l;
        auto& bucket = buckets[u];

        Stage3Cell* targetCell = nullptr;
        int emptyIndex = -1;
        std::vector<int> discontinuousIndices;

        for (int a = 0; a < static_cast<int>(b); ++a) {
            if (bucket[a].empty()) {
                if (emptyIndex < 0) emptyIndex = a;
                continue;
            }
            if (bucket[a].flowID == flowID) {
                targetCell = &bucket[a];
            } else {
                uint32_t lastwin = bucket[a].window + bucket[a].number * MIN_SUBFLOWS;
                if (startW != lastwin) {
                    discontinuousIndices.push_back(a);
                }
            }
        }

        if (!targetCell && emptyIndex != -1) {
            initNewCell(bucket[emptyIndex], flowID, startW, var, mean);
        } else if (targetCell) {
            uint32_t lastwin = targetCell->window + targetCell->number * MIN_SUBFLOWS;

            if (startW != lastwin) {
                clearCell(*targetCell);
                initNewCell(*targetCell, flowID, startW, var, mean);
            } else {
                if (canMergeVariance(*targetCell, var, mean)) {
                    mergeCell(*targetCell, var, mean);
                    tryReportStableFlow(*targetCell);
                    if (targetCell->number >= static_cast<uint32_t>(P)) {
                        clearCell(*targetCell);
                        initNewCell(*targetCell, flowID, startW, var, mean);
                    }
                } else {
                    clearCell(*targetCell);
                    initNewCell(*targetCell, flowID, startW, var, mean);
                }
            }
        } else if (!targetCell) {
            if (!discontinuousIndices.empty()) {
                int victimIndex = discontinuousIndices[0];
                for (int idx : discontinuousIndices) {
                    if (bucket[idx].number < bucket[victimIndex].number) {
                        victimIndex = idx;
                    }
                }
                clearCell(bucket[victimIndex]);
                initNewCell(bucket[victimIndex], flowID, startW, var, mean);
            } else {
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
                    uint32_t totalStableWindows = minStableWindows * MIN_SUBFLOWS;
                    float replaceProb = 1.0f / std::max(1.0f, static_cast<float>(totalStableWindows - MIN_SUBFLOWS + 1));
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
                tryReportStableFlow(cell);
                clearCell(cell);
            }
        }
    }
};

#endif
