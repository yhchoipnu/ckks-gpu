#include "clmp/clmp.h"

clmp_limb_t clmpz_bit(clmpz_srcptr _lhs, clmp_uint32 _rhs) {
    clmp_uint32 shift = _rhs / 0x20;
    clmp_uint32 bit = _rhs % 0x20;

    clmp_limb_t target = 0;

    if (shift < CLMP_LIMB_LEN) {
        target = __CLMP_CAST(clmp_limb_t, (__CLMP_CAST(clmp_uint64, _lhs->__limb[shift]) >> (0x20 - bit)));
    }

    if (target) {
        return 1;
    }

    return 0;
}
