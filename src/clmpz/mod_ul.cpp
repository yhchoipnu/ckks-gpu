#include "clmp/clmp.h"


void clmpz_mod_ul(clmpz_ptr _res, clmpz_srcptr _lhs, clmp_uint64 _rhs) {
    clmpz_t lhs, rhs;

    clmp_int32 sign_lhs = clmpz_sign(_lhs);

    clmpz_set(lhs, _lhs);
    clmpz_set_ul(rhs, _rhs);

    if (sign_lhs == CLMP_SIGN_NEG) {
        clmpz_neg(lhs, lhs);
    }

    clmpz_t quot, restore;
    clmpz_t clmp_one;

    clmpz_set_ui(quot, 0);
    clmpz_set_ui(restore, 0);

    clmpz_set_ui(clmp_one, 1);

    for (clmp_bitcnt_t i = clmpz_bits(lhs); i > -1; i--) {
        clmpz_t shift, tmp_r;

        clmpz_shift_left_ui(shift, rhs, i);
        clmpz_add(tmp_r, restore, shift);

        if (clmpz_cmp(lhs, tmp_r) > -1) {
            clmpz_t tmp_q;

            clmpz_add(restore, restore, shift);
            clmpz_shift_left_ui(tmp_q, clmp_one, i);
            clmpz_bit_or(quot, quot, tmp_q);
        }
    }

    clmpz_sub(_res, lhs, restore);

    if (sign_lhs == CLMP_SIGN_NEG) {
        clmpz_neg(_res, _res);

        if (clmpz_is_zero(_res) == CLMP_BOOL_FALSE) {
            clmpz_add(_res, _res, rhs);
        }
    }
}
