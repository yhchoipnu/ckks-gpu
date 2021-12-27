#include "clmp/clmp.h"

void clmpz_bit_xor(clmpz_ptr _res, clmpz_srcptr _lhs, clmpz_srcptr _rhs) {
    for (clmp_size_t i = 0; i < CLMP_LIMB_LEN; i++) {
        _res->__limb[i] = _lhs->__limb[i] ^ _rhs->__limb[i];
    }
}
