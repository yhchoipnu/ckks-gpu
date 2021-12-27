#include <random>

#include "clmp/clmp.h"

clmp_uint32 clmp_rand_bits_ui(clmp_uint32 _bits) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0x80000000, 0xFFFFFFFF);

    clmp_uint32 res;

    if (_bits >= 32) {
        res = __CLMP_CAST(clmp_uint32, rand()/*dis(gen)*/);
    }
    else {
        res = __CLMP_CAST(clmp_uint32, (rand()/*dis(gen)*/ & (0x7FFFFFFF >> (sizeof(clmp_uint32) * 8 - _bits))) | 0x40000000 >> (sizeof(clmp_uint32) * 8 - _bits - 1));
    }

    return res;
}
