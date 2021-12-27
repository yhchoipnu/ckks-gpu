#include "clmp/clmp.h"

clmp_uint64 clmpz_to_ul(clmpz_srcptr _clmpz) {
    clmp_uint64 res;

    res = __CLMP_CAST(clmp_uint64, _clmpz->__limb[1]) << 32;
    res += __CLMP_CAST(clmp_uint64, _clmpz->__limb[0]);

    return res;
}
