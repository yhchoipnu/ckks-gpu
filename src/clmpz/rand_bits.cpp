#include "clmp/clmp.h"

void clmpz_rand_bits(clmpz_ptr _res, clmp_uint64 _bits) {
    clmpz_set_ui(_res, 0);

    for (clmp_size_t i = 0, j = 0; j < CLMP_LIMB_BITS && j < _bits; i++, j += sizeof(clmp_limb_t) * 8) {
        if (j + sizeof(clmp_limb_t) * 8 >= _bits) {
            _res->__limb[i] = __CLMP_CAST(clmp_limb_t, clmp_rand_bits_ui(_bits - j));
        }
        else {
            _res->__limb[i] = __CLMP_CAST(clmp_limb_t, clmp_rand_bits_ui(32));
        }
    }
}
