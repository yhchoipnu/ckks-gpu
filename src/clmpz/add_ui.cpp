#include "clmp/clmp.h"

void clmpz_add_ui(clmpz_ptr _res, clmpz_srcptr _lhs, clmp_uint32 _rhs) {
    if (_res != _lhs) {
        clmpz_set(_res, _lhs);
    }

    clmp_limb_t carry = 0;

    _res->__limb[0] = _res->__limb[0] + __CLMP_CAST(clmp_limb_t, _rhs);

    if (_res->__limb[0] < __CLMP_CAST(clmp_limb_t, _rhs)) {
        carry = 1;
    }

    for (clmp_size_t i = 0; i < CLMP_LIMB_LEN; i++) {
        if (carry) {
            carry = 0;

            if (_res->__limb[i] == 0xFFFFFFFF) {
                carry = 1;
            }

            _res->__limb[i]++;
        }
    }
}
