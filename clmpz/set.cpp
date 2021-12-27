#include "clmp/clmp.h"

void clmpz_set(clmpz_ptr _res, clmpz_srcptr _rhs) {
    memcpy(_res->__limb, _rhs->__limb, CLMP_LIMB_BYTES);
}
