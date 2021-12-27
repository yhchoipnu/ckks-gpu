#include "clmp/clmp.h"

void clmpz_neg(clmpz_ptr _res, clmpz_srcptr _rhs) {
    for (clmp_size_t i = 0; i < CLMP_LIMB_LEN; i++) {
        _res->__limb[i] = ~_rhs->__limb[i];
    }

    clmpz_add_ui(_res, _res, 1);
}
