#include "clmp/clmp.h"

void clmpz_mul_ui(clmpz_ptr _res, clmpz_srcptr _lhs, clmp_uint32 _rhs) {
    clmpz_t lhs, rhs;
    clmp_limb_t res[CLMP_LIMB_LEN * 2];

    clmp_int32 sign_lhs = clmpz_sign(_lhs);

    clmpz_set(lhs, _lhs);
    clmpz_set_ui(rhs, _rhs);
    memset(res, 0, CLMP_LIMB_LEN * 2 * sizeof(clmp_limb_t));

    if (sign_lhs == CLMP_SIGN_NEG) {
        clmpz_neg(lhs, lhs);
    }

    for (clmp_size_t i = 0; i < CLMP_LIMB_LEN; i++) {
        clmp_uint64 carry = 0;
        for (clmp_size_t j = 0; j < CLMP_LIMB_LEN; j++) {
            carry += __CLMP_CAST(clmp_uint64, lhs->__limb[i]) * __CLMP_CAST(clmp_uint64, rhs->__limb[j]) + __CLMP_CAST(clmp_uint64, res[i + j]);
            res[i + j] = __CLMP_CAST(clmp_limb_t, carry);
            carry >>= 32;
        }
    }

    memcpy(_res->__limb, res, CLMP_LIMB_BYTES);

    if (sign_lhs == CLMP_SIGN_NEG) {
        clmpz_neg(_res, _res);
    }
}
