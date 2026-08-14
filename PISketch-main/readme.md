# PISketch

A memory-efficient algorithm for detecting **stable (persistent) flows** in high-speed network streams.

## Overview

Given a stream of packet data, PISketch identifies flows that:

1. Appear in **at least 200 time windows** (configurable)
2. Exhibit **low inter-arrival variance** (gap between consecutive appearances ≤ 3.0)

These stable flows typically correspond to long-lived connections, periodic beaconing traffic, or persistent service interactions.

The algorithm uses a compact memory budget (default: **400 KB**) with a combination of a Bloom Filter and a hash-based counter structure to score and rank flows.

## Algorithm

The sketch maintains two components:

- **Bloom Filter**: Tracks which flows have been seen in the current time cycle. A cleared Bloom Filter signals a new cycle.
- **PISketch**: A per-bucket counter structure where each cell holds `<flow_id, w, window_count>`. Insert weight `w = L` on first appearance, `w = -1` on repeat. Full buckets evict the smallest-`w` cell.

After processing the stream, a flow is considered stable if:
- `window_count ≥ 200` (number of cycles where it appeared)
- `variance(gaps) ≤ 3.0` (variance of inter-arrival gaps between appearances)

## Build

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make
```

## Run

Place your binary flow file as `flow.dat` in the working directory, then:

```bash
./pisketch
```

Each 32-bit (4-byte) record in `flow.dat` is treated as a flow identifier. Records are grouped into windows of 200 items each.

## Output

| File | Description |
|------|-------------|
| `output.csv` | Aggregate metrics: ARE, Precision, Recall, F1, MSE, Throughput, stable flow count |
| `stable_flows.csv` | Detected stable flows: hex ID, window count, gap variance, sketch score |

## Configuration

Key constants in `mrun.cpp`:

```cpp
constexpr int    WINDOW_THRESHOLD   = 200;    // min active windows
constexpr double VARIANCE_THRESHOLD = 3.0;    // max gap variance
constexpr int    TOTAL_MEMORY       = 400*1024; // bytes
constexpr int    DEFAULT_L          = 10;     // initial insert weight
constexpr int    DEFAULT_P          = 5;      // number of hash buckets
```

## File Structure

```
.
├── mrun.cpp      # Entry point, data loading, evaluation
├── PISketch.h    # BloomFilter + Bucket + PISketch core
├── data.h        # Data/Packet structures, BOBHash32 implementation
├── CMakeLists.txt
└── README.md
```
