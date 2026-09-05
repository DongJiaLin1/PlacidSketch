#ifndef STEADYSKETCH_H
#define STEADYSKETCH_H

#include <cstdint>
#include <cstring>
#include <vector>
#include <map>
#include <utility>

#include "MurmurHash.h"
#include "parm.h"

namespace steady {

// ---------------------------------------------------------------------------
// SteadySketch: detects steady (smooth) flows in a packet stream.
//
//   - Stable    : variance of the per-window count <= kStableThreshold
//                 over the last kPeriod windows.
//   - Steady    : a flow that is stable for >= kSmoothThreshold consecutive
//                 windows.
//
// Construction:
//   SteadySketch ss(kMemoryKB);     // pick sizes so total memory <= kMemoryKB
//
// Stream interface:
//   ss.Insert(item);                // item.id is a 32-bit flow ID,
//                                   // item.window is its time window
//   ss.Flush();                     // finalize any pending steady flows
//   for (auto& f : ss.reports()) {  // iterate detected steady flows
//       f.id; f.variance; f.duration; f.start; f.end;
//   }
// ---------------------------------------------------------------------------
class SteadySketch {
public:
    struct Report {
        std::uint32_t id       = 0;
        double        variance = 0.0;
        int           duration = 0;
        int           start    = 0;
        int           end      = 0;
    };

    // Build a sketch whose total memory (RBF + GS + USS) is at most `memory_kb` KB.
    explicit SteadySketch(int memory_kb = kMemoryKB) {
        // USS is fixed by the bucket parameters; subtract it from the budget.
        const int uss_bytes =
            kNumHeavyHitterBuckets * kElementPerBucket * static_cast<int>(sizeof(StableElement));
        const int remaining = memory_kb * 1024 - uss_bytes;
        const int filter_budget = static_cast<int>(kFilterRatio * remaining);
        const int sketch_budget = remaining - filter_budget;

        d_     = kSketchArray;
        L_     = sketch_budget / (kSketchArray * static_cast<int>(sizeof(unsigned char)) * kWindow);
        rbf_n_ = kFilterArray;
        rbf_l_ = filter_budget / (kFilterArray * static_cast<int>(sizeof(short)));

        rbf_ = new short[rbf_n_ * rbf_l_]();
        gs_  = new unsigned char[d_ * L_ * kWindow]();
        heavy_hitter_ = new StableElement*[kNumHeavyHitterBuckets];
        for (int i = 0; i < kNumHeavyHitterBuckets; ++i) {
            heavy_hitter_[i] = new StableElement[kElementPerBucket]();
        }

        rbf_bytes_ = sizeof(short) * rbf_n_ * rbf_l_;
        gs_bytes_  = sizeof(unsigned char) * d_ * L_ * kWindow;
        uss_bytes_ = sizeof(StableElement) * kNumHeavyHitterBuckets * kElementPerBucket;
    }

    SteadySketch(const SteadySketch&) = delete;
    SteadySketch& operator=(const SteadySketch&) = delete;

    ~SteadySketch() {
        for (int i = 0; i < kNumHeavyHitterBuckets; ++i) delete[] heavy_hitter_[i];
        delete[] heavy_hitter_;
        delete[] gs_;
        delete[] rbf_;
    }

    // Insert one packet into the sketch.
    void Insert(const Item& e) {
        const int last_ts = static_cast<int>(timestamp_);
        timestamp_ = static_cast<unsigned>(e.window % kWindow);
        const int ts  = static_cast<int>(timestamp_);
        int next_ts   = ts + 1;
        if (next_ts == kWindow) next_ts = 0;

        if (ts != last_ts) {
            for (int i = 0; i < rbf_n_ * rbf_l_; ++i) rbf_[i] &= static_cast<short>(~(1 << next_ts));
            for (int p = next_ts; p < d_ * L_ * kWindow; p += kWindow) gs_[p] = 0;
        }

        const std::uint64_t hash = MurmurHash32(&e.id, sizeof(e.id), kHashSeed);
        UpdateRBF(hash, ts, next_ts);
        if (rbf_tag() >= rbf_n_) return;

        const Group grp = UpdateGS(hash, ts, next_ts);

        int ex[2]  = {0, 0};
        int ex2[2] = {0, 0};
        for (int i = 0; i < kWindow; ++i) {
            ex [0] += grp.v[0][i];
            ex2[0] += grp.v[0][i] * grp.v[0][i];
            ex [1] += grp.v[1][i];
            ex2[1] += grp.v[1][i] * grp.v[1][i];
        }

        const int span   = kWindow - 2;
        const int thresh = static_cast<int>(kStableThreshold * span * span);
        const int d0 = ex2[0] * span - ex[0] * ex[0];
        const int d1 = ex2[1] * span - ex[1] * ex[1];

        if (d0 > thresh && d1 > thresh) return;

        const double variance = static_cast<double>(
            (d0 <= thresh ? d0 : d1)) / (span * span);

        UpdateHeavyHitter(hash, e.id, e.window, variance);
    }

    // Drain any still-active steady flows into the report buffer.
    void Flush() {
        for (int i = 0; i < kNumHeavyHitterBuckets; ++i) {
            for (int j = 0; j < kElementPerBucket; ++j) {
                const StableElement& s = heavy_hitter_[i][j];
                if (s.recent_time - s.start_time >= kSmoothThreshold) {
                    PushReport(s.id, s.start_time, s.recent_time);
                }
            }
        }
    }

    // Sorted (descending by duration) reports of detected steady flows.
    std::vector<Report> reports() const {
        std::vector<Report> out;
        out.reserve(report_buffer_.size());
        for (const auto& r : report_buffer_) {
            Report x;
            x.id       = r.id;
            x.duration = r.end - r.start;
            x.start    = r.start;
            x.end      = r.end;
            auto it = variance_map_.find(r.id);
            x.variance = (it != variance_map_.end()) ? it->second : 0.0;
            out.push_back(x);
        }
        std::sort(out.begin(), out.end(),
                  [](const Report& a, const Report& b) { return a.duration > b.duration; });
        return out;
    }

    // Approximate memory usage in bytes.
    int rbf_bytes() const { return rbf_bytes_; }
    int gs_bytes()  const { return gs_bytes_; }
    int uss_bytes() const { return uss_bytes_; }
    int total_bytes() const { return rbf_bytes_ + gs_bytes_ + uss_bytes_; }

private:
    static constexpr unsigned kHashSeed = 0x100;

    struct Group { unsigned char v[2][kWindow]; };

    struct RawReport {
        std::uint32_t id;
        int           start;
        int           end;
    };

    void UpdateRBF(std::uint64_t hash, int ts, int next_ts) {
        unsigned long long h = hash;
        const unsigned tmp   = ((rbf_l_ >> kRBFAlpha) - 1) << kRBFAlpha;
        unsigned index       = (h >> (rbf_n_ * kRBFAlpha)) % (tmp ? tmp : 1);
        rbf_tag_ = 0;
        for (int i = 0; i < rbf_n_; ++i) {
            const unsigned addr = (index + (h & ((1u << kRBFAlpha) - 1))) % ((i + 1) * rbf_l_);
            rbf_tag_ += (rbf_[addr] & (1 << ts)) >> ts;
            rbf_[addr] = static_cast<short>(rbf_[addr] | (1 << ts));
            const int masked = (~(rbf_[addr] | (1 << next_ts))) & ((1 << kWindow) - 1);
            if (masked != 0) rbf_tag_ = rbf_n_;
            index += rbf_l_;
            h >>= kRBFAlpha;
        }
    }

    Group UpdateGS(std::uint64_t hash, int ts, int next_ts) {
        unsigned long long h = hash;
        const unsigned tmp = ((L_ >> kGSAlpha) - 1) << kGSAlpha;
        unsigned index     = (h >> (d_ * kGSAlpha)) % (tmp ? tmp : 1);

        Group g;
        for (int i = 0; i < kWindow; ++i) g.v[0][i] = g.v[1][i] = 0xff;

        for (int i = 0; i < d_; ++i) {
            const unsigned addr = ((index + (h & ((1u << kGSAlpha) - 1))) % ((i + 1) * L_)) * kWindow;
            gs_[addr + ts]      = static_cast<unsigned char>(gs_[addr + ts] + 1);
            gs_[addr + next_ts] = 0;
            index += L_;
            h >>= kGSAlpha;

            if (rbf_tag_ >= rbf_n_) continue;
            for (int j = 0; j < kWindow; ++j) {
                const unsigned char c = gs_[addr + j];
                if (g.v[0][j] > c) g.v[0][j] = c;
                const unsigned char c1 = static_cast<unsigned char>(c + 0x80);
                if (g.v[1][j] > c1) g.v[1][j] = c1;
            }
            g.v[0][ts] = g.v[1][ts] = 0;
        }
        return g;
    }

    int rbf_tag() const { return rbf_tag_; }

    void UpdateHeavyHitter(std::uint64_t hash, std::uint32_t id, int window, double variance) {
        const unsigned bucket = hash % kNumHeavyHitterBuckets;
        StableElement* slot = heavy_hitter_[bucket];

        for (int i = 0; i < kElementPerBucket; ++i) {
            StableElement& s = slot[i];

            if (s.recent_time == 0) {  // empty slot
                s.start_time  = window;
                s.recent_time = window + 1;
                s.id          = id;
                variance_map_[id] = variance;
                return;
            }
            if (s.recent_time < window) {  // interrupted -> replace
                if (s.recent_time - s.start_time >= kSmoothThreshold) {
                    PushReport(s.id, s.start_time, s.recent_time);
                }
                s.start_time  = window;
                s.recent_time = window + 1;
                s.id          = id;
                variance_map_[id] = variance;
                return;
            }
            if (s.id == id) {  // existing flow -> extend
                s.recent_time = window + 1;
                variance_map_[id] = variance;
                return;
            }
        }
    }

    void PushReport(std::uint32_t id, int start, int end) {
        report_buffer_.push_back({id, start, end});
    }

    short*                  rbf_             = nullptr;
    unsigned char*          gs_              = nullptr;
    StableElement**         heavy_hitter_    = nullptr;
    unsigned                d_               = 0;  // #GS arrays
    unsigned                L_               = 0;  // per-GS-array length
    unsigned                rbf_n_           = 0;  // #RBF arrays
    unsigned                rbf_l_           = 0;  // per-RBF-array length
    int                     rbf_tag_         = 0;
    unsigned                timestamp_       = 0;
    int                     rbf_bytes_       = 0;
    int                     gs_bytes_        = 0;
    int                     uss_bytes_       = 0;
    std::vector<RawReport>  report_buffer_;
    std::map<std::uint32_t, double> variance_map_;
};

}  // namespace steady

#endif  // STEADYSKETCH_H