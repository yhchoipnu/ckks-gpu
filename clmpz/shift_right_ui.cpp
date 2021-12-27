#include "clmp/clmp.h"

void clmpz_shift_right_ui(clmpz_ptr _res, clmpz_srcptr _lhs, clmp_uint32 _rhs) {
    clmp_uint32 limb_shift = _rhs >> 5;
    clmp_limb_t bit_shift = _rhs & 0x1F;

    if (_rhs >= CLMP_LIMB_BITS) {
        memset(_res->__limb, 0, CLMP_LIMB_BYTES);
        return;
    }

    _res->__limb[0] = _lhs->__limb[limb_shift] >> bit_shift;

    for (clmp_int64 i = 1; i < CLMP_LIMB_LEN - limb_shift; i++) {
        if (bit_shift == 0x0) {
            _res->__limb[i - 1] |= 0x0;
        }
        else {
            _res->__limb[i - 1] |= _lhs->__limb[i + limb_shift] << (0x20 - bit_shift);
        }

        _res->__limb[i] = _lhs->__limb[i + limb_shift] >> bit_shift;
    }

    memset(_res->__limb + CLMP_LIMB_LEN - 1 - limb_shift, 0, (limb_shift + 1) * sizeof(clmp_limb_t));
}
