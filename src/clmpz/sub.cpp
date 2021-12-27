#include "clmp/clmp.h"

void clmpz_sub(clmpz_ptr _res, clmpz_srcptr _lhs, clmpz_srcptr _rhs) {
    if (_res != _lhs) {
        clmpz_set(_res, _lhs);
    }

    clmp_limb_t carry = 1;
    clmp_limb_t tmp;

    for (clmp_size_t i = 0; i < CLMP_LIMB_LEN; i++) {
        tmp = ~_rhs->__limb[i];
        _res->__limb[i] += tmp;

        if (carry) {
            if (_res->__limb[i] >= tmp) {
                carry = 0;
            }

            if (_res->__limb[i] == 0xFFFFFFFF) {
                carry = 1;
            }

            _res->__limb[i]++;
        }
        else {
            if (_res->__limb[i] < tmp) {
                carry = 1;
            }
        }
    }
}
