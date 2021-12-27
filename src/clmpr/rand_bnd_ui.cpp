#include <random>

#include "clmp/clmp.h"

clmp_uint32 clmp_rand_bnd_ui(clmp_uint32 _num) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0x0, _num);

    clmp_uint32 res = __CLMP_CAST(clmp_uint32, dis(gen));

    return res;
}
