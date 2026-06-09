#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <cstdint>
#include <cstring>

constexpr int KEY_LEN = 4;

// Memory configuration for three-stage sketch
constexpr size_t STAGE1_2_TOTAL_MEMORY_BYTES = 300ull * 1024;
constexpr double STAGE1_MEMORY_RATIO = 20.0 / 100.0;  // r%
constexpr size_t STAGE1_MEMORY_BYTES = static_cast<size_t>(STAGE1_2_TOTAL_MEMORY_BYTES * STAGE1_MEMORY_RATIO);
constexpr size_t STAGE2_MEMORY_BYTES = STAGE1_2_TOTAL_MEMORY_BYTES - STAGE1_MEMORY_BYTES;
constexpr int STAGE1_ROWS = 3;
constexpr int STAGE2_ROWS = 2;

constexpr size_t STAGE3_MEMORY_BYTES = 100ull * 1024;
constexpr int STAGE3_BUCKETS = 4;

// Stability detection parameters
constexpr int SUBFLOW_WINDOWS = 5;      // Number of windows in a subflow
constexpr int COUNTER_BITS = 8;          // Counter bit width (max value = 256)
constexpr int MIN_SUBFLOWS = 5;          // Window span per subflow in Stage3 continuity tracking
constexpr int Steady_FLOWS = 200;        // Minimum window duration to report a stable flow
constexpr int P = 400;                  // Max subflows per cell before reset
constexpr int Q = 40;                   // Subflow count threshold
constexpr float STABLE_THRESHOLD = 3.0f; // Variance threshold for stability

// Packet structure: represents a network flow packet with 32-bit flow fingerprint
struct Packet {
    char flowID[KEY_LEN];
    uint32_t windowNumber;

    Packet() : windowNumber(0) {
        memset(flowID, 0, KEY_LEN);
    }

    // Constructor for 32-bit flow fingerprint (binary .dat files)
    Packet(uint32_t flowID32, uint32_t window) : windowNumber(window) {
        memcpy(flowID, &flowID32, sizeof(uint32_t));
    }
};

#endif
