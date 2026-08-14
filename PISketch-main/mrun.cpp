#include "PISketch.h"
#include <set>
#include <map>
#include <vector>
#include <fstream>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <ctime>
#include <iostream>

using namespace std;

constexpr int WINDOW_THRESHOLD   = 200;
constexpr double VARIANCE_THRESHOLD = 3.0;
constexpr int TOTAL_MEMORY      = 400 * 1024;
constexpr int WINDOW_SIZE      = WINDOW_THRESHOLD;
constexpr int DEFAULT_L         = 10;
constexpr int DEFAULT_P         = 5;
constexpr int DEFAULT_AVEV      = 5000;

struct Parameter {
    uint64_t input_max;
    uint64_t max_mem;
    int      threshold;
    int      aveV;
    int      L;
    int      p;
};

struct DetectedFlow {
    Data            flow;
    unsigned int    window_count;
    double          variance;
    int             sketch_w;
};

struct Result {
    double            MEM;
    double            ARE;
    double            PR;
    double            RR;
    double            F1;
    double            CR;
    double            MSE;
    double            THP;
    Parameter         para;
    uint64_t         Npi;
    uint64_t          stable_flow_count;
    double            avg_variance;
    vector<DetectedFlow> detected_flows;

    void print() const {
        printf("%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%ld,%.4f\n",
               MEM, ARE, PR, RR, CR, F1, MSE, THP,
               static_cast<long>(stable_flow_count), avg_variance);
    }
};

Packet* load_dat(const char* path, int& total_packets, int& total_windows) {
    ifstream fin(path, ios::binary);
    if (!fin.is_open()) {
        cerr << "Cannot open file: " << path << endl;
        return nullptr;
    }

    vector<Packet> all;
    unsigned char buf[4];
    int window_id = 1;
    int in_window = 0;

    while (fin.read(reinterpret_cast<char*>(buf), 4)) {
        Packet pkt;
        memset(pkt.dat.str, 0, DATA_LEN);
        memcpy(pkt.dat.str, buf, 4);
        pkt.timestamp = window_id;
        all.push_back(pkt);

        if (++in_window >= WINDOW_SIZE) {
            ++window_id;
            in_window = 0;
        }
    }
    fin.close();

    total_windows = window_id - 1;
    total_packets = static_cast<int>(all.size());
    Packet* P = new Packet[total_packets];
    for (int i = 0; i < total_packets; ++i) P[i] = all[i];

    sort(P, P + total_packets);
    cerr << "Packets loaded: " << total_packets << endl;
    cerr << "Windows: " << total_windows << endl;
    return P;
}

void PI_Test(Parameter para, Packet* P, Result* R, int total_windows) {
    R->para = para;
    R->detected_flows.clear();

    const uint64_t input_max = para.input_max;
    const uint64_t max_mem  = para.max_mem;
    const int      thr      = para.threshold;
    const int      L        = para.L;
    const int      p        = para.p;

    const double alpha = PI_alpha;
    const int BF_L  = static_cast<int>(max_mem * alpha / 3 * 8);
    const int SK_L  = static_cast<int>(max_mem * (1 - alpha) / p /
                                       (sizeof(int) + sizeof(int) + sizeof(Data) + sizeof(bool)));

    BloomFilter* bf = new BloomFilter(3, BF_L);
    PISketch*    sk = new PISketch(SK_L, p);

    const ts_t T_interval = P[input_max - 1].timestamp - P[0].timestamp;
    ts_t preTime = P[0].timestamp;
    const double timeCycle = static_cast<double>(T_interval) / input_max * para.aveV;
    int currentCycle = 1;

    set<Data>            s, s_all;
    map<Data, int>       mp_std, mpwn_std;
    map<Data, vector<int>> flow_windows;

    for (uint64_t i = 0; i < input_max; ++i) {
        const Packet& pkt = P[i];
        const Data& item  = pkt.dat;

        if (pkt.timestamp - preTime > timeCycle) {
            s.clear();
            bf->clear();
            preTime = pkt.timestamp;
            ++currentCycle;
        }

        int w = bf->query(item) ? -1 : L;
        if (w == L) bf->insert(item);
        sk->insert(item, w);

        s_all.insert(item);
        if (s.find(item) != s.end()) {
            mp_std[item] -= 1;
        } else {
            mp_std[item] += L;
            ++mpwn_std[item];
        }
        s.insert(item);

        if (flow_windows[item].empty() || flow_windows[item].back() != currentCycle)
            flow_windows[item].push_back(currentCycle);
    }

    set<Data> s_std;
    for (const auto& kv : mp_std)
        if (kv.second >= thr) s_std.insert(kv.first);

    map<Data, double>   flow_variance;
    map<Data, unsigned> flow_wn;
    double total_var = 0;
    uint64_t stable_count = 0;

    for (const auto& kv : flow_windows) {
        const Data& item    = kv.first;
        const vector<int>& wins = kv.second;
        const int wn = static_cast<int>(wins.size());

        double var = 0;
        if (wn >= 2) {
            double sum = 0;
            for (int k = 1; k < wn; ++k) sum += wins[k] - wins[k - 1];
            const double mean = sum / (wn - 1);
            double ss = 0;
            for (int k = 1; k < wn; ++k) {
                const double d = (wins[k] - wins[k - 1]) - mean;
                ss += d * d;
            }
            var = ss / (wn - 1);
        }

        if (wn >= thr && var <= VARIANCE_THRESHOLD) {
            ++stable_count;
            total_var += var;
            flow_variance[item] = var;
            flow_wn[item] = static_cast<unsigned>(wn);
        }
    }

    set<Data> s_my;
    for (const auto& item : s_all) {
        const auto res = sk->query(item);
        if (res.first >= thr && flow_variance.count(item))
            s_my.insert(item);
    }

    double are = 0;
    int cnt_are = 0;
    for (const auto& item : s_my) {
        if (mpwn_std[item] > 0) {
            const auto res = sk->query(item);
            are += fabs(static_cast<double>(res.second) - mpwn_std[item]) / mpwn_std[item];
            ++cnt_are;
        }
    }
    if (cnt_are > 0) are /= cnt_are;

    int TP = 0, FP = 0, FN = 0;
    for (const auto& item : s_my)
        (s_std.count(item) ? TP : FP++);
    for (const auto& item : s_std)
        if (!s_my.count(item)) ++FN;

    const double PR = (TP + FP > 0) ? static_cast<double>(TP) / (TP + FP) : 0;
    const double RR = (TP + FN > 0) ? static_cast<double>(TP) / (TP + FN) : 0;
    const double F1 = (PR + RR > 0) ? 2 * PR * RR / (PR + RR) : 0;
    const double CR = s_all.empty() ? 0 : static_cast<double>(s_my.size()) / s_all.size();

    double mse = 0;
    int cnt_mse = 0;
    for (const auto& item : s_my) {
        if (mpwn_std[item] > 0) {
            const auto res = sk->query(item);
            const double e = res.second - mpwn_std[item];
            mse += e * e;
            ++cnt_mse;
        }
    }
    if (cnt_mse > 0) mse /= cnt_mse;

    R->MEM             = max_mem;
    R->ARE             = are;
    R->PR              = PR;
    R->RR              = RR;
    R->F1              = F1;
    R->CR              = CR;
    R->MSE             = mse;
    R->Npi             = s_std.size();
    R->stable_flow_count = stable_count;
    R->avg_variance    = stable_count > 0 ? total_var / stable_count : 0;

    for (const auto& kv : flow_variance) {
        DetectedFlow df;
        df.flow        = kv.first;
        df.window_count = flow_wn[kv.first];
        df.variance    = kv.second;
        df.sketch_w    = sk->query(kv.first).first;
        R->detected_flows.push_back(df);
    }

    delete bf;
    delete sk;
}

void PI_Throughput(Parameter para, Packet* P, Result* R) {
    R->para = para;
    const uint64_t input_max = para.input_max;
    const uint64_t max_mem  = para.max_mem;
    const int      L        = para.L;
    const int      p        = para.p;
    const double   alpha    = PI_alpha;

    const int BF_L = static_cast<int>(max_mem * alpha / 3 * 8);
    const int SK_L = static_cast<int>(max_mem * (1 - alpha) / p /
                                       (sizeof(int) + sizeof(int) + sizeof(Data) + sizeof(bool)));

    BloomFilter* bf = new BloomFilter(3, BF_L);
    PISketch*    sk = new PISketch(SK_L, p);

    const ts_t T_interval = P[input_max - 1].timestamp - P[0].timestamp;
    ts_t preTime = P[0].timestamp;
    const double timeCycle = static_cast<double>(T_interval) / input_max * para.aveV;

    const clock_t t0 = clock();
    for (uint64_t i = 0; i < input_max; ++i) {
        const Packet& pkt = P[i];
        const Data& item  = pkt.dat;
        if (pkt.timestamp - preTime > timeCycle) {
            bf->clear();
            preTime = pkt.timestamp;
        }
        const int w = bf->query(item) ? -1 : L;
        if (w == L) bf->insert(item);
        sk->insert(item, w);
    }
    const clock_t t1 = clock();
    const double sec = static_cast<double>(t1 - t0) / CLOCKS_PER_SEC;
    R->THP = sec > 0 ? input_max / sec / 1'000'000.0 : 0;

    delete bf;
    delete sk;
}

int main() {
    const char* input_file = "flow.dat";
    const char* output_csv = "output.csv";
    const char* stable_csv = "stable_flows.csv";

    int total_packets = 0, total_windows = 0;
    Packet* P = load_dat(input_file, total_packets, total_windows);
    if (!P || total_packets == 0) {
        cerr << "Failed to load data." << endl;
        return 1;
    }

    const int max_input = min(total_packets, static_cast<int>(MAX_INPUT_SIZE));
    const Parameter para{static_cast<uint64_t>(max_input),
                         static_cast<uint64_t>(TOTAL_MEMORY),
                         WINDOW_THRESHOLD,
                         DEFAULT_AVEV,
                         DEFAULT_L,
                         DEFAULT_P};

    Result R;
    PI_Test(para, P, &R, total_windows);
    PI_Throughput(para, P, &R);

    FILE* fout = fopen(output_csv, "w");
    if (fout) {
        fprintf(fout, "Total_Packets,%d\n", max_input);
        fprintf(fout, "Total_Windows,%d\n", total_windows);
        fprintf(fout, "Window_Threshold,%d\n", WINDOW_THRESHOLD);
        fprintf(fout, "Variance_Threshold,%.1f\n", VARIANCE_THRESHOLD);
        fprintf(fout, "Total_Memory_Bytes,%d\n", TOTAL_MEMORY);
        fprintf(fout, "MEM,ARE,PR,RR,CR,F1,MSE,THP,Stable_Flow_Count,Avg_Variance\n");
        R.print();
        fclose(fout);
    }

    FILE* sfout = fopen(stable_csv, "w");
    if (sfout) {
        fprintf(sfout, "flow_hex,window_count,variance,sketch_w\n");
        for (const auto& df : R.detected_flows) {
            uint32_t val = 0;
            memcpy(&val, df.flow.str, sizeof(uint32_t));
            fprintf(sfout, "%08x,%u,%.6f,%d\n",
                    val, df.window_count, df.variance, df.sketch_w);
        }
        fclose(sfout);
    }

    cerr << "Stable flows: " << R.stable_flow_count << endl;
    cerr << "Results written to " << output_csv << " and " << stable_csv << endl;

    delete[] P;
    return 0;
}
