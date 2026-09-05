#ifndef STEADYSKETCH_MURMURHASH_H
#define STEADYSKETCH_MURMURHASH_H

#include <cstdint>
#include <cstddef>

namespace steady {

// MurmurHash2 32-bit. Returns a 32-bit hash of `len` bytes at `key`.
// Public-domain reference implementation by Austin Appleby.
inline std::uint32_t MurmurHash32(const void* key, int len, std::uint32_t seed) {
    const std::uint32_t m = 0x5bd1e995u;
    const int           r = 24;

    std::uint32_t h = seed ^ static_cast<std::uint32_t>(len);
    const unsigned char* data = static_cast<const unsigned char*>(key);

    while (len >= 4) {
        std::uint32_t k;
        std::memcpy(&k, data, sizeof(k));
        k *= m;
        k ^= k >> r;
        k *= m;
        h *= m;
        h ^= k;
        data += 4;
        len  -= 4;
    }
    switch (len) {
        case 3: h ^= static_cast<std::uint32_t>(data[2]) << 16; [[fallthrough]];
        case 2: h ^= static_cast<std::uint32_t>(data[1]) << 8;  [[fallthrough]];
        case 1: h ^= data[0];
                h *= m;
    }
    h ^= h >> 13;
    h *= m;
    h ^= h >> 15;
    return h;
}

}  // namespace steady

#endif  // STEADYSKETCH_MURMURHASH_H