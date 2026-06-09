# PlacidSketch

A three-stage sketch-based system for detecting steady (long-lived) network flows. It processes traffic data window-by-window using a multi-stage pipeline, keeping memory usage low through hash-based sketch techniques.

## Overview

PlacidSketch uses a three-stage pipeline:

```
.dat files (one per time window)
        |
        v
   Stage1Filter    -- Quick candidate filtering via per-window arrival tracking
        |
        v
   Stage2Monitor   -- Stability monitoring with ring-buffer counters
        |
        v
   Stage3Merger    -- Merges consecutive stable subflows, reports final results
```

Each flow is represented by a 32-bit identifier. Stage1 uses 3 hash rows; Stage2 uses 2 rows with a ring buffer of 6 counters per row; Stage3 uses a CM-sketch-style merger with variance-based stability merging.

## Building

```bash
mkdir build
cd build
cmake ..
cmake --build .
```

Or with a single command:

```bash
cmake -B build -S . && cmake --build build
```

## Running

```bash
./build/placidsketch <data_folder_path>
```

## Input Format

Each `.dat` file represents one time window. Files are processed in lexicographic order by filename. Each file is treated as a stream of 32-bit binary records, where each record is a flow identifier. File size must be a multiple of 4 bytes.

## Configuration

All key parameters are in `parm.h`:


| Parameter             | Default | Description                            |
| --------------------- | ------- | -------------------------------------- |
| `STAGE1_MEMORY_BYTES` | 80 KB   | Stage1 bucket memory                   |
| `STAGE2_MEMORY_BYTES` | 220 KB  | Stage2 bucket memory                   |
| `STAGE3_MEMORY_BYTES` | 100 KB  | Stage3 merger memory                   |
| `SUBFLOW_WINDOWS`     | 5       | Windows per subflow                    |
| `MIN_SUBFLOWS`        | 40      | Subflows needed for stability          |
| `Steady_FLOWS`        | 200     | Total subflows to report a stable flow |
| `STABLE_THRESHOLD`    | 3.0     | Variance threshold for merging         |
| `ALPHA_THRESHOLD`     | 5       | Counter difference threshold           |


## Project Structure

```
.
├── CMakeLists.txt
├── MurmurHash3.h          # MurmurHash3_x86_32 implementation
├── parm.h                 # All configurable constants and Packet struct
├── stage1.h               # Stage1: candidate filter
├── stage2.h               # Stage2: stability monitor
├── stage3.h               # Stage3: subflow merger
└── main.cpp               # Entry point, data loader, sketch orchestration
```

## Dependencies

A C++17 compiler and CMake 3.15+. No third-party libraries required.