#ifndef CKKS_KEY_H
#define CKKS_KEY_H

#include "clmp/clmp.h"
#include "ckks/param.h"

namespace CKKS {
    class key {
    public:
        clmp_uint64 *rax = new clmp_uint64[Nnprimes];
        clmp_uint64 *rbx = new clmp_uint64[Nnprimes];

        key();
        ~key();
    };
}

#endif
