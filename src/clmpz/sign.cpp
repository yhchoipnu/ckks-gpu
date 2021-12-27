#include "clmp/clmp.h"

clmp_int32 clmpz_sign(clmpz_srcptr _rhs) {
    if ((_rhs->__limb[CLMP_LIMB_LEN-1] & 0x80000000)) {
        return CLMP_SIGN_NEG;
    }
    else {
        return CLMP_SIGN_POS;
    }
}
