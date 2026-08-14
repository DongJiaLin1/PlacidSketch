#pragma once
#define DATA_LEN 4

#include <cstdint>
#include <cstring>
#include <functional>

typedef uint32_t uint;
typedef uint64_t ulong;
typedef uint8_t  uchar;
typedef double   ts_t;

constexpr uint MAX_INPUT_SIZE = 8000000;
constexpr double PI_alpha = 0.1;

uint BOBHash32_run(const uchar* str, uint len, uint primeIdx);

class Data {
public:
    uchar str[DATA_LEN];

    Data() { memset(str, 0, DATA_LEN); }

    Data& operator=(const Data& other) {
        memcpy(str, other.str, DATA_LEN);
        return *this;
    }

    uint Hash(uint num = 0) const {
        return BOBHash32_run(str, DATA_LEN, num);
    }
};

bool operator<(const Data& a, const Data& b) {
    for (int i = 0; i < DATA_LEN; ++i) {
        if (a.str[i] < b.str[i]) return true;
        if (a.str[i] > b.str[i]) return false;
    }
    return false;
}

bool operator==(const Data& a, const Data& b) {
    for (int i = 0; i < DATA_LEN; ++i)
        if (a.str[i] != b.str[i]) return false;
    return true;
}

struct Packet {
    Data dat;
    ts_t timestamp;
    bool operator<(const Packet& other) const {
        return timestamp < other.timestamp;
    }
};

struct Element {
    Data item;
    int w;
    int window_number;
    bool operator<(const Element& other) const {
        return w < other.w || (w == other.w && item < other.item);
    }
    bool operator==(const Element& other) const {
        return w == other.w && item == other.item;
    }
    bool operator>(const Element& other) const {
        return w > other.w || (w == other.w && other.item < item);
    }
};

struct HashRet {
    int* ret;
    explicit HashRet(int k) { ret = new int[k]; }
};

struct BOBHash32 {
    uint32_t primeIdx;
    BOBHash32() : primeIdx(0) {}
    explicit BOBHash32(uint32_t p) : primeIdx(p) {}
    void initialize(uint32_t p) { primeIdx = p; }
    uint32_t run(const char* str, uint32_t len) const {
        return BOBHash32_run(reinterpret_cast<const uchar*>(str), len, primeIdx);
    }
};

#define mix(a, b, c)                \
    do {                            \
        a -= b; a -= c; a ^= (c >> 13); \
        b -= c; b -= a; b ^= (a << 8);  \
        c -= a; c -= b; c ^= (b >> 13); \
        a -= b; a -= c; a ^= (c >> 12); \
        b -= c; b -= a; b ^= (a << 16); \
        c -= a; c -= b; c ^= (b >> 5);  \
        a -= b; a -= c; a ^= (c >> 3);  \
        b -= c; b -= a; b ^= (a << 10); \
        c -= a; c -= b; c ^= (b >> 15); \
    } while (0)

static const uint32_t bobhash_prime32[10] = {
    20177, 20201, 20249, 20269, 20297, 20323, 20341, 20353, 20369, 20393
};

uint BOBHash32_run(const uchar* str, uint len, uint primeIdx) {
    uint32_t a = 0x9e3779b9;
    uint32_t b = 0x9e3779b9;
    uint32_t c = bobhash_prime32[primeIdx % 10];

    while (len >= 12) {
        a += (str[0]       ) + ((uint32_t)str[1] <<  8) + ((uint32_t)str[2] << 16) + ((uint32_t)str[3] << 24);
        b += (str[4]       ) + ((uint32_t)str[5] <<  8) + ((uint32_t)str[6] << 16) + ((uint32_t)str[7] << 24);
        c += (str[8]       ) + ((uint32_t)str[9] <<  8) + ((uint32_t)str[10] << 16) + ((uint32_t)str[11] << 24);
        mix(a, b, c);
        str += 12; len -= 12;
    }

    c += len;
    switch (len) {
        case 11: c += ((uint32_t)str[10] << 24);
        case 10: c += ((uint32_t)str[ 9] << 16);
        case  9: c += ((uint32_t)str[ 8] <<  8);
        case  8: b += ((uint32_t)str[ 7] << 24);
        case  7: b += ((uint32_t)str[ 6] << 16);
        case  6: b += ((uint32_t)str[ 5] <<  8);
        case  5: b += str[4];
        case  4: a += ((uint32_t)str[ 3] << 24);
        case  3: a += ((uint32_t)str[ 2] << 16);
        case  2: a += ((uint32_t)str[ 1] <<  8);
        case  1: a += str[0];
        case  0: break;
    }
    mix(a, b, c);
    return c;
}
