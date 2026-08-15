/**
 * @file main.cc
 * @brief XSketch test program
 */

#include <bits/stdc++.h>
#include "XSketch.h"
#include "Param.h"
#include "Common/Mmap.h"

using namespace std;

typedef uint64_t ID_TYPE;

int main(int argc, char** argv) {
    init_matrix();

    if (argc < 4) {
        cout << "Usage: " << argv[0] << " <data_file_path> <memory_size(MB)> <run_length>" << endl;
        return 1;
    }

    string PATH = argv[1];
    uint32_t memory = atoi(argv[2]);
    uint32_t run_length = atoi(argv[3]);

    LoadResult load_result = Load(PATH.c_str());
    CAIDA_Tuple* dataset = (CAIDA_Tuple*)load_result.start;
    uint64_t length = load_result.length / sizeof(CAIDA_Tuple);

    cout << "Dataset length: " << length << endl;
    cout << "Memory size: " << memory << " MB" << endl;
    cout << "Run length: " << run_length << endl;

    XSketch<ID_TYPE>* x_sketch = new XSketch<ID_TYPE>(
        memory, 
        var_thres_p,
        error_thres_p,
        stage_ratio_p,
        bucket_size_p,
        S_p,
        potential_thres_p
    );

    auto start = std::chrono::high_resolution_clock::now();

    for (uint32_t i = 0; i < run_length; ++i) {
        x_sketch->insert(dataset[i % length].id, i);

        if (i > 0 && i % 10000 == 0) {
            cout << "Processed " << i << " records..." << endl;
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    chrono::duration<double, milli> tm = end - start;

    cout << "Insertion completed!" << endl;
    cout << "Time elapsed: " << tm.count() << " ms" << endl;
    cout << "Throughput: " << run_length / (1.0 * tm.count() * 1000) << " M/s" << endl;

    vector<pair<ID_TYPE, uint32_t>> result = x_sketch->query();
    vector<Report_Slot<ID_TYPE>> report = x_sketch->report();

    cout << "Number of detected flows: " << result.size() << endl;
    cout << "Number of reported Top-K flows: " << report.size() << endl;

    delete x_sketch;
    UnLoad(load_result);

    return 0;
}
