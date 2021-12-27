#include "clmp/clmp.h"

void clmpz_sub_ul(clmpz_ptr _res, clmpz_srcptr _lhs, clmp_uint64 _rhs) {
    if (_res != _lhs) {
        clmpz_set(_res, _lhs);
    }

    clmp_limb_t carry = 1;
    clmp_limb_t tmp;

    tmp = ~(__CLMP_CAST(clmp_limb_t, _rhs));
    _res->__limb[0] += tmp;
    _res->__limb[0]++;

    if (_res->__limb[0] >= tmp) {
        carry = 1;
    }

    tmp = ~(__CLMP_CAST(clmp_limb_t, _rhs >> 32));
    _res->__limb[1] += tmp;

    if (carry) {
        if (_res->__limb[1] >= tmp) {
            carry = 0;
        }

        if (_res->__limb[1] == 0xFFFFFFFF) {
            carry = 1;
        }

        _res->__limb[1]++;
    }
    else {
        if (_res->__limb[1] < tmp) {
            carry = 1;
        }
    }

    for (clmp_size_t i = 2; i < CLMP_LIMB_LEN; i++) {
        tmp = ~0x00000000;
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
