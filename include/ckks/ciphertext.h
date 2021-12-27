#ifndef CKKS_CIPHERTEXT_H
#define CKKS_CIPHERTEXT_H

#include "clmp/clmp.h"
#include "ckks/param.h"

namespace CKKS {
    class ciphertext{
    public:

        clmpz_t *ax;
        clmpz_t *bx;

        clmp_uint64 log_p;
        clmp_uint64 log_q;
        clmp_uint64 len;
        clmp_uint64 num;

        ciphertext(clmp_uint64 _log_p = 0, clmp_uint64 _log_q = 0, clmp_uint64 _len = 0, clmp_uint64 _num = 0);
        ciphertext(const ciphertext &_origin);

        ~ciphertext();

        void free();
        void init();
    };
}

#endif
