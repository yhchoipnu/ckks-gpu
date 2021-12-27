#include "clmp/clmp.h"

void clmpz_add(clmpz_ptr _res, clmpz_srcptr _lhs, clmpz_srcptr _rhs) {
    clmp_limb_t carry = 0;

    for (clmp_size_t i = 0; i < CLMP_LIMB_LEN; i++) {
        _res->__limb[i] = _lhs->__limb[i] + _rhs->__limb[i];

        if (carry) {
            if ((_res->__limb[i] >= _lhs->__limb[i]) && (_res->__limb[i] >= _rhs->__limb[i])) {
                carry = 0;
            }

            if (_res->__limb[i] == 0xFFFFFFFF) {
                carry = 1;
            }

            _res->__limb[i]++;
        }
        else {
            if ((_res->__limb[i] < _lhs->__limb[i]) || (_res->__limb[i] < _rhs->__limb[i])) {
                carry = 1;
            }
        }
    }
}
