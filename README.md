# PlacidSketch

PlacidSketch is a three-stage streaming algorithm for detecting stable flows in network traffic data. It uses a multi-stage filtering and monitoring approach to efficiently identify flows with stable frequency patterns.

## Overview

PlacidSketch consists of three main stages:

- **Stage 1 (Filter)**: Filters candidate stable flows using continuity-based detection
- **Stage 2 (Monitor)**: Monitors flows for stability using frequency variance analysis
- **Stage 3 (Merger)**: Merges stable subflows using incremental multi-segment merging technique

## Requirements

- C++17 or later
- CMake 3.26 or later
- A C++ compiler with C++17 support

## Building

```bash
mkdir build
cd build
cmake ..
make
```

Or use CMake directly:

```bash
cmake -B build
cmake --build build
```

## Usage

1. Prepare your input data as CSV files in a directory, where each CSV file represents a time window
2. Modify the `folderPath` variable in `main.cpp` to point to your data directory
3. Run the compiled executable

```bash
./main
```

## Configuration

Parameters can be modified in `parm.h`:

- `STAGE1_MEMORY_BYTES`: Memory allocation for Stage 1
- `STAGE2_MEMORY_BYTES`: Memory allocation for Stage 2
- `STAGE3_MEMORY_BYTES`: Memory allocation for Stage 3
- `SUBFLOW_WINDOWS`: Number of windows for stability detection
- `STABLE_THRESHOLD`: Variance threshold for stability
