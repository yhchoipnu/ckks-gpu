#include "clmp/clmp.h"

clmp_uint32 clmpz_to_ui(clmpz_srcptr _clmpz) {
    return __CLMP_CAST(clmp_uint32, _clmpz->__limb[0]);
}
