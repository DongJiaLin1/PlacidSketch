#include <unordered_set>
#include <set>
#include <map>
#include <unordered_map>
#include <ctime>
#include <time.h>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cerrno>
#include <cctype>
#include <float.h>
#include "PISketch.h"
#include "para.h"

using namespace std;

uint32_t insert_data[MAX_INSERT_PACKAGE];
uint32_t query_data[MAX_INSERT_PACKAGE];

int total_pkt_num = 0;
int total_window_count = 0;

int PI_X;
int PI_Y;
int BF_len;
int L;
int W_thr;

// Reads a flat binary stream of flow keys: each 4 bytes is interpreted as
// a uint32_t flow id (host endianness). Windows are sliced purely by key
// count: every WindowSize keys form one window; any tail (< WindowSize)
// still counts as the final window. No headers, no per-key metadata.
int load_data() {
    const char* path = DATA_BIN;
    FILE* pf = fopen(path, "rb");
    if (!pf) {
        fprintf(stderr, "Error: Cannot open binary input '%s': %s\n",
                path, strerror(errno));
        exit(1);
    }

    ground_truth_f.clear();
    ground_truth_p.clear();
    ground_truth_w.clear();

    int ret = 0;
    int now_w = 0;
    int keys_in_window = 0;
    uint32_t key;

    while (ret < MAX_INSERT_PACKAGE &&
           fread(&key, sizeof(uint32_t), 1, pf) == 1) {
        insert_data[ret] = key;
        ground_truth_f[key]++;

        if (ground_truth_w.find(key) == ground_truth_w.end()) {
            ground_truth_p[key]++;
            ground_truth_w[key] = now_w;
        } else if (ground_truth_w[key] < now_w) {
            ground_truth_p[key]++;
            ground_truth_w[key] = now_w;
        }

        ret++;
        keys_in_window++;
        if (keys_in_window >= WindowSize) {
            keys_in_window = 0;
            now_w++;
        }
    }
    if (keys_in_window > 0) now_w++;

    fclose(pf);

    int i = 0;
    for (auto itr : ground_truth_f) query_data[i++] = itr.first;

    printf("Loaded binary input '%s': %d keys, %d distinct flows, %d windows\n\n",
           path, ret, (int)ground_truth_f.size(), now_w);

    total_pkt_num = ret;
    total_window_count = now_w;
    return ret;
}

inline double ABS(double a, double b) {
    if (a > b) return a - b;
    return b - a;
}

struct ExpResult {
    double memoryKB;
    double throughputMpps;
    double PR, CR, F1;
    double MSE, ARE;
};

ExpResult run_experiment(int mem_idx) {
    PI_X   = MEM_CONFIG[mem_idx][0];
    PI_Y   = MEM_CONFIG[mem_idx][1];
    BF_len = MEM_CONFIG[mem_idx][2];

    PISketch PISketch_Model(PI_X, PI_Y);
    BloomFilter BF_SM(3, BF_len);

    double pi_mem = (double)PI_X * PI_Y * PI_size / 1024.0;
    double bf_mem = (double)BF_len * BF_size / 1024.0;
    double total_mem = pi_mem + bf_mem;

    int PI_num = 0;
    for (auto itr : ground_truth_f) {
        uint32_t k = itr.first;
        double den = (double)ground_truth_f[k] / ground_truth_p[k];
        if (ground_truth_p[k] > P_thr && den < D_thr) PI_num++;
    }

    int pkt_num = total_pkt_num;
    int w_cnt = 0;

    clock_t t0, t1;

    t0 = clock();
    for (int i = 0; i < pkt_num; i++) {
        int w = BF_SM.query(insert_data[i]) ? -1 : L;
        BF_SM.insert(insert_data[i]);
        PISketch_Model.insert(insert_data[i], w);
        w_cnt++;
        if (w_cnt % WindowSize == 0) BF_SM.clear();
    }
    t1 = clock();
    double time_sec = (double)(t1 - t0) / CLOCKS_PER_SEC;
    if (time_sec <= 0.0) time_sec = 1e-6;
    double throughput = pkt_num / time_sec / 1e6;

    unordered_set<uint32_t> pi_reported;
    for (int pix = 0; pix < PI_X; pix++) {
        for (int piy = 0; piy < PI_Y; piy++) {
            int tmp_p = PISketch_Model.bucket[pix].p[piy];
            int tmp_f = tmp_p * (L + 1) - PISketch_Model.bucket[pix].w[piy];
            double tmp_d = (double)tmp_f / tmp_p;
            if (tmp_p > P_thr && tmp_d < D_thr) {
                uint32_t findkey = PISketch_Model.getfp(pix, piy);
                pi_reported.insert(findkey);
            }
        }
    }

    int TP = 0, FP = 0;
    for (uint32_t k : pi_reported) {
        double den = (double)ground_truth_f[k] / ground_truth_p[k];
        if (ground_truth_p[k] > P_thr && den < D_thr) TP++;
        else FP++;
    }

    double PR = (TP + FP) == 0 ? 0.0 : (double)TP / (TP + FP);
    double CR = PI_num == 0 ? 0.0 : (double)TP / PI_num;
    double F1 = (PR + CR) == 0.0 ? 0.0 : 2.0 * PR * CR / (PR + CR);

    double errsum = 0.0, mse_sum = 0.0;
    int fk = 0;
    for (int pix = 0; pix < PI_X; pix++) {
        for (int piy = 0; piy < PI_Y; piy++) {
            int pi_p = PISketch_Model.bucket[pix].p[piy];
            int pi_f = (L + 1) * pi_p - PISketch_Model.bucket[pix].w[piy];
            double pi_den = (double)pi_f / pi_p;
            uint32_t findkey = PISketch_Model.getfp(pix, piy);
            double true_den = (double)ground_truth_f[findkey] / ground_truth_p[findkey];
            double err = ABS(pi_den, true_den);
            errsum += err;
            mse_sum += err * err;
            fk++;
        }
    }
    double ARE = fk == 0 ? 0.0 : errsum / fk;
    double MSE = fk == 0 ? 0.0 : mse_sum / fk;

    ExpResult r;
    r.memoryKB = total_mem;
    r.throughputMpps = throughput;
    r.PR = PR;
    r.CR = CR;
    r.F1 = F1;
    r.MSE = MSE;
    r.ARE = ARE;
    return r;
}

void print_sep() {
    printf("======================================================================\n");
}

void print_table_header() {
    printf("%-8s %10s %10s %10s %8s %10s %10s\n",
           "Memory", "Throughput", "PR", "CR", "F1", "MSE", "ARE");
}

void print_result(const char* name, const ExpResult& r) {
    printf("%-8s %10.4f %10.4f %10.4f %8.4f %10.6f %10.6f\n",
           name, r.throughputMpps, r.PR, r.CR, r.F1, r.MSE, r.ARE);
}

int main(int argc, char** argv) {
    L = 10;
    W_thr = 200;

    printf("Loading data from %s...\n", DATA_BIN);
    int pkt_num = load_data();

    int PI_num = 0;
    for (auto itr : ground_truth_f) {
        uint32_t k = itr.first;
        double den = (double)ground_truth_f[k] / ground_truth_p[k];
        if (ground_truth_p[k] > P_thr && den < D_thr) PI_num++;
    }

    print_sep();
    printf("  PSSketch PI Detection Experiment (PISketch Algorithm)\n");
    printf("  P_thr=%d (windows>=%d), D_thr=%.1f, WindowSize=%d, Total_Windows=%d\n",
           P_thr, P_thr, D_thr, WindowSize, total_window_count);
    printf("  Total packets=%d, Total flows=%d, Ground-truth PI flows=%d\n",
           pkt_num, (int)ground_truth_f.size(), PI_num);
    print_sep();

    printf("\n");
    print_table_header();

    ExpResult results[MEM_COUNT];
    for (int i = 0; i < MEM_COUNT; i++) {
        printf("\n--- %s ---\n", MEM_NAMES[i]);
        printf("  PI_X=%d, PI_Y=%d, BF_len=%d, PI_mem=%.1fKB, BF_mem=%.1fKB\n",
               MEM_CONFIG[i][0], MEM_CONFIG[i][1], MEM_CONFIG[i][2],
               (double)MEM_CONFIG[i][0] * MEM_CONFIG[i][1] * PI_size / 1024.0,
               (double)MEM_CONFIG[i][2] * BF_size / 1024.0);
        results[i] = run_experiment(i);
        print_result(MEM_NAMES[i], results[i]);
    }

    printf("\n");
    print_sep();
    printf("  Summary (PISketch)\n");
    print_sep();
    print_table_header();
    for (int i = 0; i < MEM_COUNT; i++) print_result(MEM_NAMES[i], results[i]);
    printf("\nDone.\n");
    return 0;
}