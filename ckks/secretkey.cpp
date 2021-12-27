#include "ckks/secretkey.h"

namespace CKKS{
    secretkey::secretkey(context &_ctx) {
        for (clmp_uint64 i = 0; i < N; i++) {
            clmpz_set_i(sx[i], 0);
        }

        _ctx.sampleHWT(sx);
    }

    secretkey::~secretkey() {
        delete [] sx;
    }
}
