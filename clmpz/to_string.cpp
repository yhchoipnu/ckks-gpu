#include "clmp/clmp.h"
#include <iostream>
std::string clmpz_to_string(clmpz_srcptr _clmpz) {
    std::string res = "";

    if (clmpz_is_zero(_clmpz) == CLMP_BOOL_TRUE) {
        res = "0";
    }
    else {
        clmpz_t tmp, mod;
        clmp_int32 sign;

        clmpz_set(tmp, _clmpz);

        sign = clmpz_sign(tmp);
        if (sign == CLMP_SIGN_NEG) {
            clmpz_neg(tmp, tmp);
        }

        do {
            clmpz_mod_ui(mod, tmp, 10);
            
            res = "0123456789"[mod->__limb[0]] + res;

            clmpz_div_ui(tmp, tmp, 10);
        } while (clmpz_is_zero(tmp) == CLMP_BOOL_FALSE);

        if (sign == CLMP_SIGN_NEG) {
            res = "-" + res;
        }
    }

    return  res;
}
