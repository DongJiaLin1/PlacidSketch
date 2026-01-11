#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <cstdint>
#include <cstring>

constexpr int KEY_LEN = 16;
constexpr int KEY_LEN1 = 64;

constexpr size_t STAGE1_2_TOTAL_MEMORY_BYTES = 600ull * 1024;
constexpr double STAGE1_MEMORY_RATIO = 50.0 / 600.0;
constexpr size_t STAGE1_MEMORY_BYTES = static_cast<size_t>(STAGE1_2_TOTAL_MEMORY_BYTES * STAGE1_MEMORY_RATIO);
constexpr size_t STAGE2_MEMORY_BYTES = STAGE1_2_TOTAL_MEMORY_BYTES - STAGE1_MEMORY_BYTES;
constexpr int STAGE1_ROWS = 3;
constexpr int STAGE2_ROWS = 2;

constexpr size_t STAGE3_MEMORY_BYTES = 200ull * 1024;
constexpr int STAGE3_BUCKETS = 4;

constexpr int SUBFLOW_WINDOWS = 5;
constexpr int COUNTER_BITS = 8;
constexpr int ALPHA_THRESHOLD = 10;
constexpr int MIN_SUBFLOWS = 5;
constexpr int Steady_FLOWS = 200;
constexpr int P = 400;
constexpr int Q = 40;
constexpr float STABLE_THRESHOLD = 5.0f;

struct Packet {
    char flowID[KEY_LEN];
    char quintuple[KEY_LEN1];
    uint32_t windowNumber;

    Packet() : windowNumber(0) {
        memset(flowID, 0, KEY_LEN);
        memset(quintuple, 0, KEY_LEN1);
    }

    Packet(const char* fingerprint, const char* quintupleStr, uint32_t window) : windowNumber(window) {
        strncpy(flowID, fingerprint, KEY_LEN);
        flowID[KEY_LEN - 1] = '\0';
        if (quintupleStr) {
            strncpy(quintuple, quintupleStr, KEY_LEN1);
            quintuple[KEY_LEN1 - 1] = '\0';
        } else {
            memset(quintuple, 0, KEY_LEN1);
        }
    }
};

#endif