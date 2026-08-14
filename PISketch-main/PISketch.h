#pragma once

#include "data.h"
#include <cstring>
#include <cstdlib>
#include <queue>
#include <map>
#include <iostream>

class BloomFilter {
    int k_;
    int m_;
    BOBHash32* bh_;
    bool* bits_;

public:
    BloomFilter(int k, int m) : k_(k), m_(m) {
        bh_ = new BOBHash32[k];
        for (int i = 0; i < k; ++i)
            bh_[i].initialize(rand() % 10 + 10 * i);
        bits_ = new bool[m];
        memset(bits_, 0, sizeof(bool) * m);
    }

    ~BloomFilter() {
        delete[] bh_;
        delete[] bits_;
    }

    void insert(const Data& x) {
        for (int i = 0; i < k_; ++i) {
            int hi = bh_[i].run(reinterpret_cast<const char*>(x.str), DATA_LEN) % m_;
            bits_[hi] = true;
        }
    }

    bool query(const Data& x) const {
        for (int i = 0; i < k_; ++i) {
            int hi = bh_[i].run(reinterpret_cast<const char*>(x.str), DATA_LEN) % m_;
            if (!bits_[hi]) return false;
        }
        return true;
    }

    void clear() {
        memset(bits_, 0, sizeof(bool) * m_);
    }

    int memUsed() const { return m_; }
};

struct Bucket {
    int* w_;
    int* winNum_;
    Data* dat_;
    bool* valid_;
    int capacity_;

    Bucket() : w_(nullptr), winNum_(nullptr), dat_(nullptr), valid_(nullptr), capacity_(0) {}

    void init(int cap) {
        capacity_ = cap;
        valid_  = new bool[cap];
        w_       = new int[cap];
        winNum_  = new int[cap];
        dat_     = new Data[cap];
        memset(valid_, 0, sizeof(bool) * cap);
    }

    void insert(const Data& d, int wi) {
        int evict = -1;
        for (int i = 0; i < capacity_; ++i) {
            if (!valid_[i]) continue;
            if (dat_[i] == d) {
                w_[i] += wi;
                if (wi > 0) ++winNum_[i];
                return;
            }
            if (evict < 0 || w_[i] < w_[evict]) evict = i;
        }

        for (int i = 0; i < capacity_; ++i) {
            if (!valid_[i]) {
                dat_[i]  = d;
                w_[i]     = wi;
                winNum_[i] = 1;
                valid_[i]  = true;
                return;
            }
        }

        w_[evict]--;
        if (w_[evict] < 0) {
            dat_[evict]  = d;
            w_[evict]    = wi + 1;
            winNum_[evict] = 1;
        }
    }

    std::pair<int, int> query(const Data& d) const {
        for (int i = 0; i < capacity_; ++i) {
            if (valid_[i] && dat_[i] == d)
                return {w_[i], winNum_[i]};
        }
        return {-1, -1};
    }
};

class PISketch {
    Bucket* buckets_;
    int numBuckets_;
    int perBucket_;
    BOBHash32 bh_;

public:
    PISketch(int nb, int ppc) : numBuckets_(nb), perBucket_(ppc) {
        buckets_ = new Bucket[nb];
        for (int i = 0; i < nb; ++i) buckets_[i].init(ppc);
        bh_.initialize(rand() % 50);
    }

    void insert(const Data& d, int wi) {
        int pos = bh_.run(reinterpret_cast<const char*>(d.str), DATA_LEN) % numBuckets_;
        buckets_[pos].insert(d, wi);
    }

    std::pair<int, int> query(const Data& d) const {
        int pos = bh_.run(reinterpret_cast<const char*>(d.str), DATA_LEN) % numBuckets_;
        return buckets_[pos].query(d);
    }

    int memUsed() const {
        return numBuckets_ * perBucket_ *
               (sizeof(int) + sizeof(int) + sizeof(Data) + sizeof(bool));
    }
};
