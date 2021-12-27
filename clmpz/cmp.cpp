#include "clmp/clmp.h"

clmp_int32 clmpz_cmp(clmpz_srcptr _lhs, clmpz_srcptr _rhs) {
    clmp_int32 sign_lhs = clmpz_sign(_lhs);
    clmp_int32 sign_rhs = clmpz_sign(_rhs);

    if (sign_lhs * sign_rhs == -1) {
        if (sign_lhs == CLMP_SIGN_POS) {
            return 1;
        }
        else {
            return -1;
        }
    }
    else {
        for (clmp_int64 i = CLMP_LIMB_LEN - 1; i > -1; i--) {
            if (_lhs->__limb[i] > _rhs->__limb[i]) {
                return 1;
            }
            if (_lhs->__limb[i] < _rhs->__limb[i]) {
                return -1;
            }
        }

        return 0;
    }
}
