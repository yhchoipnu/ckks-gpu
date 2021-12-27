#ifndef CKKS_PARAM_H
#define CKKS_PARAM_H

#include "clmp/clmp.h"

namespace CKKS {
    static const clmp_uint64 logN = 4;
    static const clmp_uint64 logQ = 40;

    static const double sigma = 3.2;
    static const clmp_uint64 h = 1;
    static const clmp_uint64 pbnd = 59.0;
    static const clmp_uint64 kbar = pbnd + 1;
    static const clmp_uint64 kbar2 = 2 * kbar;
    static const clmp_uint64 logNh = (logN - 1);
    static const clmp_uint64 logQQ = (2 * logQ);
    static const clmp_uint64 N = (1 << logN);
    static const clmp_uint64 Nh = (1 << logNh);
    static const clmp_uint64 M = (N << 1);
    static const clmp_uint64 nprimes = (2 + logN + 4 * logQ + pbnd - 1) / pbnd;
    static const clmp_uint64 Nnprimes = (nprimes << logN);
    static const clmp_uint64 bignum = 0xfffffff;
}

#endif
