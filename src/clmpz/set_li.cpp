#include "clmp/clmp.h"

void clmpz_set_li(clmpz_ptr _res, clmp_int64 _rhs) {
    memset(_res->__limb, 0, CLMP_LIMB_BYTES);

    if (_rhs >= 0) {
        _res->__limb[0] = __CLMP_CAST(clmp_limb_t, _rhs);
        _res->__limb[1] = __CLMP_CAST(clmp_limb_t, _rhs >> 32);
    }
    else {
        _res->__limb[0] = __CLMP_CAST(clmp_limb_t, -_rhs);
        _res->__limb[1] = __CLMP_CAST(clmp_limb_t, (-_rhs) >> 32);

        clmpz_neg(_res, _res);
    }
}
