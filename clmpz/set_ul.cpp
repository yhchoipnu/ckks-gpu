#include "clmp/clmp.h"

void clmpz_set_ul(clmpz_ptr _res, clmp_uint64 _rhs) {
    memset(_res->__limb, 0, CLMP_LIMB_BYTES);

    _res->__limb[0] = __CLMP_CAST(clmp_limb_t, _rhs);
    _res->__limb[1] = __CLMP_CAST(clmp_limb_t, _rhs >> 32);
}
