#include "clmp/clmp.h"

void clmpz_mul(clmpz_ptr _res, clmpz_srcptr _lhs, clmpz_srcptr _rhs) {
    clmpz_t lhs, rhs;
    clmp_limb_t res[CLMP_LIMB_LEN * 2];

    clmp_int32 sign_lhs = clmpz_sign(_lhs);
    clmp_int32 sign_rhs = clmpz_sign(_rhs);

    clmpz_set(lhs, _lhs);
    clmpz_set(rhs, _rhs);
    memset(res, 0, CLMP_LIMB_LEN * 2 * sizeof(clmp_limb_t));

    if (sign_lhs == CLMP_SIGN_NEG) {
        clmpz_neg(lhs, lhs);
    }

    if (sign_rhs == CLMP_SIGN_NEG) {
        clmpz_neg(rhs, rhs);
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

    if (sign_lhs * sign_rhs == CLMP_SIGN_NEG) {
        clmpz_neg(_res, _res);
    }
}
