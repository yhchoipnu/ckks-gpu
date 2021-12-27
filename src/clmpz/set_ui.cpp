#include "clmp/clmp.h"

void clmpz_set_ui(clmpz_ptr _res, clmp_uint32 _rhs) {
    memset(_res->__limb, 0, CLMP_LIMB_BYTES);

    _res->__limb[0] = __CLMP_CAST(clmp_limb_t, _rhs);
}
