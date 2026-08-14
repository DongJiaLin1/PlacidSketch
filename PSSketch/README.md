# PSSketch

Persistent-item (PI) detection sketches implemented in C++.

The main entry point is `PSSketch.cpp`, which runs the PISketch + BloomFilter
algorithm over a flat binary stream of flow keys and reports precision,
recall, F1, throughput, MSE and ARE.

## Input format

The program reads a single binary file. Each 4 bytes (`uint32_t`, host
endianness) is interpreted as one flow id. No headers, no per-key metadata.

By default the file path is `stream.bin` in the current working directory.
Override it at compile time:

```bash
cmake -DDATA_BIN=/path/to/your/stream.bin -B build
cmake --build build
./build/PSSketch
```

A window consists of `WindowSize` keys (default `1000`). Any trailing tail
(`< WindowSize`) still forms the final window.

## Build

```bash
cmake -B build
cmake --build build -j
./build/PSSketch
```

The build requires a compiler with AVX2 and BMI support (GCC/Clang on a
recent x86_64 CPU). Flags `-mavx2 -mbmi` are added automatically by
`CMakeLists.txt`.

## Tunables

All knobs live in `para.h`. The most commonly adjusted ones:

| Macro          | Meaning                                                  |
| -------------- | -------------------------------------------------------- |
| `DATA_BIN`     | Path to the binary input file.                           |
| `WindowSize`   | Number of keys per window.                               |
| `P_thr`        | Minimum window-presence count to qualify as a PI.        |
| `D_thr`        | Maximum average-per-window frequency to qualify as a PI. |
| `MEM_CONFIG`   | Memory budgets (rows × cols + BloomFilter size).         |

## Project layout

| File           | Purpose                                                   |
| -------------- | --------------------------------------------------------- |
| `PSSketch.cpp` | Experiment driver: loads data, runs PISketch, reports.    |
| `PISketch.h`   | PISketch + BloomFilter used for PI detection.             |
| `class.h`      | Earlier L1/L2 filter-based PSSketch (LPSSketch, legacy).  |
| `class2.h`     | Small/Large-layer PSSketch (LPSSketch, current).          |
| `strawman.h`   | OO_PE + Count-Min baseline for comparison.                |
| `CMSketch.h`   | Count-Min Sketch.                                         |
| `OO_PE.h`      | On-Off persistence estimator.                             |
| `BOBHash32.h`  | 32-bit BOB hash.                                          |
| `bitset.h`     | Compact bitset used by On-Off estimator.                  |
| `data.h`       | Legacy fixed-size-packet container.                       |
| `definition.h` | Global type / size constants.                             |
| `hash.h`       | Hash helper declarations.                                 |
| `para.h`       | All tunable parameters and memory configurations.         |