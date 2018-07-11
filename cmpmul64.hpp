#pragma once

#include "common.hpp"

#if defined(_MSC_VER) && defined(_M_X64) // Needed for MSVC cmpMul64
#   include <intrin.h>
#   pragma intrinsic(_mul128)
#endif

namespace tinygeom2d {

// Portable implementation of multiplication of 64-bit unsigned integers.
// The result is the pair (high bits, low bits).
inline pair<uint64_t, uint64_t> portableMulU64(uint64_t x, uint64_t y) {
    const uint64_t low = ((uint64_t)1 << 32) - (uint64_t)1;
    
    uint64_t x0 = x & low;
    uint64_t x1 = x >> 32;
    
    uint64_t y0 = y & low;
    uint64_t y1 = y >> 32;
    
    uint64_t p00 = x0 * y0;
    uint64_t p01 = x0 * y1;
    uint64_t p10 = x1 * y0;
    uint64_t p11 = x1 * y1;
    
    uint64_t mid = (p00 >> 32) + (p01 & low) + (p10 & low);
    
    uint64_t res0 = (p00 & low) | (mid << 32);
    uint64_t res1 = p11 + (mid >> 32) + (p01 >> 32) + (p10 >> 32);
    
    return {res1, res0};
}

// Returns 1, 0 or -1, if a * b is greater than, equal to or less than c * d,
// respectively. This is the portable implementation - the function cmpMul64
// may be faster.
inline int portableCmpMul64(int64_t a, int64_t b, int64_t c, int64_t d) {
    pair<uint64_t, uint64_t> x = portableMulU64((uint64_t)a, (uint64_t)b);
    pair<uint64_t, uint64_t> y = portableMulU64((uint64_t)c, (uint64_t)d);
    
    // Transform to signed multiplication
    x.first -= a < 0 ? (uint64_t)b : 0;
    x.first -= b < 0 ? (uint64_t)a : 0;
    y.first -= c < 0 ? (uint64_t)d : 0;
    y.first -= d < 0 ? (uint64_t)c : 0;
    
    // Flip top bits so we can use unsigned comparison for signed comparison
    const uint64_t top = (uint64_t)1 << 63;
    x.first ^= top;
    y.first ^= top;
    
    if(x == y) {
        return 0;
    } else if(x > y) {
        return 1;
    } else {
        return -1;
    }
}

// Equivalent to portableCmpMul64, but may use compiler- or platform-specific
// optimizations.
#if __SIZEOF_INT128__ == 16 // gcc, clang and icc on x86_64
    inline int cmpMul64(int64_t a, int64_t b, int64_t c, int64_t d) {
        __int128 x = (__int128)a * (__int128)b;
        __int128 y = (__int128)c * (__int128)d;
        
        if(x == y) {
            return 0;
        } else {
            return x > y ? 1 : -1;
        }
    }
#elif defined(_MSC_VER) && defined(_M_X64) // MSVC on x86_64
    inline int cmpMul64(int64_t a, int64_t b, int64_t c, int64_t d) {
        pair<int64_t, uint64_t> x, y;
        x.second = _mul128(a, b, &x.first);
        y.second = _mul128(c, d, &y.first);
        
        if(x == y) {
            return 0;
        } else {
            return x > y ? 1 : -1;
        }
    }
#else
    inline int cmpMul64(int64_t a, int64_t b, int64_t c, int64_t d) {
        return portableCmpMul64(a, b, c, d);
    }
#endif

}
