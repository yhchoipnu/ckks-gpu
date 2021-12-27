#include "clmp/clmp.h"

void clmpz_set_i(clmpz_ptr _res, clmp_int32 _rhs) {
    memset(_res->__limb, 0, CLMP_LIMB_BYTES);

    if (_rhs >= 0) {
        _res->__limb[0] = __CLMP_CAST(clmp_limb_t, _rhs);
    }
    else {
        _res->__limb[0] = __CLMP_CAST(clmp_limb_t, -_rhs);

        clmpz_neg(_res, _res);
    }
}
