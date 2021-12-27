#include "clmp/clmp.h"

clmp_bitcnt_t clmpz_bits(clmpz_srcptr _rhs) {
    for (clmp_int64 i = CLMP_LIMB_LEN - 1; i > -1; i--) {
        if (_rhs->__limb[i] != 0) {
            if (_rhs->__limb[i] & 0x80000000) {
                return (i + 1) * sizeof(clmp_limb_t) * 8;
            }
            else {
                clmp_limb_t tmp = _rhs->__limb[i];
                clmp_bitcnt_t cnt = 0;

                while (tmp) {
                    tmp >>= 1;
                    cnt++;
                }

                return i * sizeof(clmp_limb_t) * 8 + cnt;
            }
        }
    }

    return 0;
}
