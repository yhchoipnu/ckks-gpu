#include "clmp/clmp.h"
#include <iostream>
clmp_int32 clmpz_is_zero(clmpz_srcptr _rhs) {
    for (clmp_size_t i = 0; i < CLMP_LIMB_LEN; i++) {
        if (_rhs->__limb[i] != 0) {
            return CLMP_BOOL_FALSE;
        }
    }

    return CLMP_BOOL_TRUE;
}
