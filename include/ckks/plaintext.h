#ifndef CKKS_PLAINTEXT_H
#define CKKS_PLAINTEXT_H

#include "clmp/clmp.h"
#include "ckks/param.h"

namespace CKKS {
    class plaintext{
    public:

        clmpz_t *msg;

        clmp_uint64 log_p;
        clmp_uint64 log_q;
        clmp_uint64 len;
        clmp_uint64 num;

        plaintext(clmp_uint64 _log_p = 0, clmp_uint64 _log_q = 0, clmp_uint64 _n = 0, clmp_uint64 _num = 0);

        ~plaintext();

        void free();
        void init();
    };
}

#endif
