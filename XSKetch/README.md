# XSketch

A two-stage Sketch data structure for network traffic anomaly detection and stable flow identification.

## Overview

XSketch is designed to detect stable network flows in high-speed network traffic using a two-stage architecture:
- **Stage1**: TowerSketch for quickly recording flow frequency information
- **Stage2**: Hash bucket structure for precise anomaly detection

## Build

```bash
mkdir build && cd build
cmake ..
make
```

## Usage

```bash
./XSketch <data_file_path> <memory_size(MB)> <run_length>
```

Parameters:
- `data_file_path`: Path to the binary data file (CAIDA format: 16 bytes per record, 8 bytes timestamp + 8 bytes flow ID)
- `memory_size`: Memory allocation in MB
- `run_length`: Number of packets to process

## Parameters

The algorithm uses the following default parameters (defined in `Param.h`):

| Parameter | Value | Description |
|-----------|-------|-------------|
| window_size | 10000 | Time window size |
| P | 7 | Number of consecutive windows for linear regression |
| K | 0 | Polynomial degree (0 = constant model) |
| S | 4 | Number of windows recorded in Stage1 |
| var_thres | 0.1 | Variance threshold for detection |
| error_thres | 0.5 | Mean squared error threshold |
| bucket_size | 4 | Number of slots per hash bucket |
| stage_ratio | 0.8 | Memory ratio for Stage1 |

## Files

```
.
├── Param.h           # Algorithm parameters and regression functions
├── XSketch.h         # Core XSketch data structure
├── main.cc           # Test program
├── Common/
│   ├── hash.h        # Hash function wrapper
│   ├── BOBHash.h     # BOBHash implementation
│   ├── Matrix.h      # Matrix operations
│   └── Mmap.h        # Memory-mapped file I/O
├── CMakeLists.txt
└── README.md
```

## Algorithm

1. **Flow Insertion**: Each packet's flow ID is hashed and inserted into the appropriate stage
2. **Window Transition**: Every N packets, the algorithm transitions to a new time window
3. **Stability Detection**: Flows with consistent packet counts across windows (low variance) are identified as stable
4. **Linear Regression**: A simplified K=0 model (mean calculation) is used for stability analysis

## License

MIT License
