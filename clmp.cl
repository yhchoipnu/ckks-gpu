#include "clmp/clmp.clh"

void clmp_memcpy(__global void *_dst, __global const void *_src, clmp_size_t _size) {
    __global const clmp_uint8 *src = __CLMP_CAST(__global const clmp_uint8 *, _src);
    __global clmp_uint8 *dst = __CLMP_CAST(__global clmp_uint8 *, _dst);

    for (clmp_size_t i = 0; i < _size; i++) {
        dst[i] = src[i];
    }
}

void clmp_memcpy_global_to_local(void *_dst, __global const void *_src, clmp_size_t _size) {
    __global const clmp_uint8 *src = __CLMP_CAST(__global const clmp_uint8 *, _src);
    clmp_uint8 *dst = __CLMP_CAST(clmp_uint8 *, _dst);

    for (clmp_size_t i = 0; i < _size; i++) {
        dst[i] = src[i];
    }
}

void clmp_memcpy_local_to_global(__global void *_dst, const void *_src, clmp_size_t _size) {
    const clmp_uint8 *src = __CLMP_CAST(const clmp_uint8 *, _src);
    __global clmp_uint8 *dst = __CLMP_CAST(__global clmp_uint8 *, _dst);

    for (clmp_size_t i = 0; i < _size; i++) {
        dst[i] = src[i];
    }
}

void clmp_memcpy_local_to_local(void *_dst, const void *_src, clmp_size_t _size) {
    const clmp_uint8 *src = __CLMP_CAST(const clmp_uint8 *, _src);
    clmp_uint8 *dst = __CLMP_CAST(clmp_uint8 *, _dst);

    for (clmp_size_t i = 0; i < _size; i++) {
        dst[i] = src[i];
    }
}

clmp_uint32 clmp_rand(clmp_uint64 _seed) {
    _seed = (_seed * 0x5DEECE66DL + 0xBL) & ((1L << 48) - 1);
    return __CLMP_CAST(clmp_uint32, _seed >> 16);
}

void clmp_memset(__global void *_dst, clmp_uint8 _byte, clmp_size_t _size) {
    __global clmp_uint8 *dst = __CLMP_CAST(__global clmp_uint8 *, _dst);

    for (clmp_size_t i = 0; i < _size; i++) {
        dst[i] = _byte;
    }
}

void clmp_memset_local(void *_dst, clmp_uint8 _byte, clmp_size_t _size) {
    clmp_uint8 *dst = __CLMP_CAST(clmp_uint8 *, _dst);

    for (clmp_size_t i = 0; i < _size; i++) {
        dst[i] = _byte;
    }
}

void clmpz_set(clmpz_ptr _res, clmpz_srcptr _rhs) {
    if (_res != _rhs) {
        clmp_memcpy(_res->__limb, _rhs->__limb, CLMP_LIMB_BYTES);
    }
}

void clmpz_set_ui(clmpz_ptr _res, clmp_uint32 _rhs) {
    clmp_memset(_res->__limb, 0, CLMP_LIMB_BYTES);

    _res->__limb[0] = __CLMP_CAST(clmp_limb_t, _rhs);
}

void clmpz_set_ul(clmpz_ptr _res, clmp_uint64 _rhs) {
    clmp_memset(_res->__limb, 0, CLMP_LIMB_BYTES);

    _res->__limb[0] = __CLMP_CAST(clmp_limb_t, _rhs);
    _res->__limb[1] = __CLMP_CAST(clmp_limb_t, _rhs >> 32);
}

void clmpz_set_i(clmpz_ptr _res, clmp_int32 _rhs) {
    clmp_memset(_res->__limb, 0, CLMP_LIMB_BYTES);

    if (_rhs >= 0) {
        _res->__limb[0] = __CLMP_CAST(clmp_limb_t, _rhs);
    }
    else {
        _res->__limb[0] = __CLMP_CAST(clmp_limb_t, -_rhs);

        clmpz_neg(_res, _res);
    }
}

void clmpz_set_li(clmpz_ptr _res, clmp_int64 _rhs) {
    clmp_memset(_res->__limb, 0, CLMP_LIMB_BYTES);

    if (_rhs >= 0) {
        _res->__limb[0] = __CLMP_CAST(clmp_limb_t, _rhs);
        _res->__limb[1] = __CLMP_CAST(clmp_limb_t, _rhs >> 32);
    }
    else {
        _res->__limb[0] = __CLMP_CAST(clmp_limb_t, -_rhs);
        _res->__limb[1] = __CLMP_CAST(clmp_limb_t, (-_rhs) >> 32);

        clmpz_neg(_res, _res);
    }
}

clmp_bitcnt_t clmpz_bits_limb(const clmp_limb_t *_rhs) {
    for (clmp_int64 i = CLMP_LIMB_LEN - 1; i > -1; i--) {
        if (_rhs[i] != 0) {
            if (_rhs[i] & 0x80000000) {
                return (i + 1) * sizeof(clmp_limb_t) * 8;
            }
            else {
                clmp_limb_t tmp = _rhs[i];
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

clmp_bitcnt_t clmpz_bits(clmpz_srcptr _rhs) {
    clmp_limb_t rhs[CLMP_LIMB_LEN];

    clmp_memcpy_global_to_local(rhs, _rhs->__limb, CLMP_LIMB_BYTES);

    return clmpz_bits_limb(rhs);
}

void clmpz_bit_and_limb(clmp_limb_t *_res, const clmp_limb_t *_lhs, const clmp_limb_t *_rhs) {
    for (clmp_size_t i = 0; i < CLMP_LIMB_LEN; i++) {
        _res[i] = _lhs[i] & _rhs[i];
    }
}

void clmpz_bit_or_limb(clmp_limb_t *_res, const clmp_limb_t *_lhs, const clmp_limb_t *_rhs) {
    for (clmp_size_t i = 0; i < CLMP_LIMB_LEN; i++) {
        _res[i] = _lhs[i] | _rhs[i];
    }
}

void clmpz_bit_xor_limb(clmp_limb_t *_res, const clmp_limb_t *_lhs, const clmp_limb_t *_rhs) {
    for (clmp_size_t i = 0; i < CLMP_LIMB_LEN; i++) {
        _res[i] = _lhs[i] ^ _rhs[i];
    }
}

void clmpz_bit_and(clmpz_ptr _res, clmpz_srcptr _lhs, clmpz_srcptr _rhs) {
    clmp_limb_t res[CLMP_LIMB_LEN], lhs[CLMP_LIMB_LEN], rhs[CLMP_LIMB_LEN];

    clmp_memcpy_global_to_local(lhs, _lhs->__limb, CLMP_LIMB_BYTES);
    clmp_memcpy_global_to_local(rhs, _rhs->__limb, CLMP_LIMB_BYTES);

    clmpz_bit_and_limb(res, lhs, rhs);

    clmp_memcpy_local_to_global(_res->__limb, res, CLMP_LIMB_BYTES);
}

void clmpz_bit_or(clmpz_ptr _res, clmpz_srcptr _lhs, clmpz_srcptr _rhs) {
    clmp_limb_t res[CLMP_LIMB_LEN], lhs[CLMP_LIMB_LEN], rhs[CLMP_LIMB_LEN];

    clmp_memcpy_global_to_local(lhs, _lhs->__limb, CLMP_LIMB_BYTES);
    clmp_memcpy_global_to_local(rhs, _rhs->__limb, CLMP_LIMB_BYTES);

    clmpz_bit_or_limb(res, lhs, rhs);

    clmp_memcpy_local_to_global(_res->__limb, res, CLMP_LIMB_BYTES);
}

void clmpz_bit_xor(clmpz_ptr _res, clmpz_srcptr _lhs, clmpz_srcptr _rhs) {
    clmp_limb_t res[CLMP_LIMB_LEN], lhs[CLMP_LIMB_LEN], rhs[CLMP_LIMB_LEN];

    clmp_memcpy_global_to_local(lhs, _lhs->__limb, CLMP_LIMB_BYTES);
    clmp_memcpy_global_to_local(rhs, _rhs->__limb, CLMP_LIMB_BYTES);

    clmpz_bit_xor_limb(res, lhs, rhs);

    clmp_memcpy_local_to_global(_res->__limb, res, CLMP_LIMB_BYTES);
}

void clmpz_shift_left_limb_ui(clmp_limb_t *_res, const clmp_limb_t *_lhs, clmp_uint32 _rhs) {
    clmp_uint32 limb_shift = _rhs >> 5;
    clmp_limb_t bit_shift = _rhs & 0x1F;

    if (_rhs >= CLMP_LIMB_BITS) {
        clmp_memset_local(_res, 0, CLMP_LIMB_BYTES);
        return;
    }

    _res[CLMP_LIMB_LEN - 1] = _lhs[CLMP_LIMB_LEN - 1 - limb_shift] << bit_shift;

    for (clmp_int64 i = CLMP_LIMB_LEN - 2; i >= limb_shift; i--) {
        if (bit_shift == 0x0) {
            _res[i + 1] |= 0x0;
        }
        else {
            _res[i + 1] |= _lhs[i - limb_shift] >> (0x20 - bit_shift);
        }

        _res[i] = _lhs[i - limb_shift] << bit_shift;
    }

    clmp_memset_local(_res, 0, limb_shift * sizeof(clmp_limb_t));
}

void clmpz_shift_right_limb_ui(clmp_limb_t *_res, const clmp_limb_t *_lhs, clmp_uint32 _rhs) {
    clmp_uint32 limb_shift = _rhs >> 5;
    clmp_limb_t bit_shift = _rhs & 0x1F;

    if (_rhs >= CLMP_LIMB_BITS) {
        clmp_memset_local(_res, 0, CLMP_LIMB_BYTES);
        return;
    }

    _res[0] = _lhs[limb_shift] >> bit_shift;

    for (clmp_int64 i = 1; i < CLMP_LIMB_LEN - limb_shift; i++) {
        if (bit_shift == 0x0) {
            _res[i - 1] |= 0x0;
        }
        else {
            _res[i - 1] |= _lhs[i + limb_shift] << (0x20 - bit_shift);
        }

        _res[i] = _lhs[i + limb_shift] >> bit_shift;
    }

    clmp_memset_local(_res + CLMP_LIMB_LEN - 1 - limb_shift, 0, (limb_shift + 1) * sizeof(clmp_limb_t));
}

void clmpz_shift_left_ui(clmpz_ptr _res, clmpz_srcptr _lhs, clmp_uint32 _rhs) {
    clmp_limb_t res[CLMP_LIMB_LEN], lhs[CLMP_LIMB_LEN];

    clmp_memcpy_global_to_local(lhs, _lhs->__limb, CLMP_LIMB_BYTES);

    clmpz_shift_left_limb_ui(res, lhs, _rhs);

    clmp_memcpy_local_to_global(_res->__limb, res, CLMP_LIMB_BYTES);
}

void clmpz_shift_right_ui(clmpz_ptr _res, clmpz_srcptr _lhs, clmp_uint32 _rhs) {
    clmp_limb_t res[CLMP_LIMB_LEN], lhs[CLMP_LIMB_LEN];

    clmp_memcpy_global_to_local(lhs, _lhs->__limb, CLMP_LIMB_BYTES);

    clmpz_shift_right_limb_ui(res, lhs, _rhs);

    clmp_memcpy_local_to_global(_res->__limb, res, CLMP_LIMB_BYTES);
}

clmp_int32 clmpz_sign(clmpz_srcptr _rhs) {
    if ((_rhs->__limb[CLMP_LIMB_LEN-1] & 0x80000000)) {
        return CLMP_SIGN_NEG;
    }
    else {
        return CLMP_SIGN_POS;
    }
}

clmp_int32 clmpz_cmp_limb(const clmp_limb_t *_lhs, const clmp_limb_t *_rhs) {
    for (clmp_int64 i = CLMP_LIMB_LEN - 1; i > -1; i--) {
        if (_lhs[i] > _rhs[i]) {
            return 1;
        }
        if (_lhs[i] < _rhs[i]) {
            return -1;
        }
    }

    return 0;
}

clmp_int32 clmpz_cmp(clmpz_srcptr _lhs, clmpz_srcptr _rhs) {
    clmp_int32 sign_lhs = clmpz_sign(_lhs);
    clmp_int32 sign_rhs = clmpz_sign(_rhs);

    if (sign_lhs * sign_rhs == -1) {
        if (sign_lhs == CLMP_SIGN_POS) {
            return 1;
        }
        else {
            return -1;
        }
    }
    else {
        for (clmp_int64 i = CLMP_LIMB_LEN - 1; i > -1; i--) {
            if (_lhs->__limb[i] > _rhs->__limb[i]) {
                return 1;
            }
            if (_lhs->__limb[i] < _rhs->__limb[i]) {
                return -1;
            }
        }

        return 0;
    }
}

clmp_int32 clmpz_is_zero_limb(clmp_limb_t *_rhs) {
    for (clmp_size_t i = 0; i < CLMP_LIMB_LEN; i++) {
        if (_rhs[i] != 0) {
            return CLMP_BOOL_FALSE;
        }
    }

    return CLMP_BOOL_TRUE;
}

clmp_int32 clmpz_is_zero(clmpz_srcptr _rhs) {
    for (clmp_size_t i = 0; i < CLMP_LIMB_LEN; i++) {
        if (_rhs->__limb[i] != 0) {
            return CLMP_BOOL_FALSE;
        }
    }

    return CLMP_BOOL_TRUE;
}

void clmpz_neg_limb(clmp_limb_t *_res , const clmp_limb_t *_rhs) {
    for (clmp_size_t i = 0; i < CLMP_LIMB_LEN; i++) {
        _res[i] = ~_rhs[i];
    }

    clmpz_add_limb_ui(_res, _res, 1);
}

void clmpz_neg(clmpz_ptr _res, clmpz_srcptr _rhs) {
    clmp_limb_t res[CLMP_LIMB_LEN], rhs[CLMP_LIMB_LEN];

    clmp_memcpy_global_to_local(rhs, _rhs->__limb, CLMP_LIMB_BYTES);

    clmpz_neg_limb(res, rhs);

    clmp_memcpy_local_to_global(_res->__limb, res, CLMP_LIMB_BYTES);
}

void clmpz_add_limb(clmp_limb_t *_res, const clmp_limb_t *_lhs, const clmp_limb_t *_rhs) {
    clmp_uint64 carry = 0;

    for (clmp_size_t i = 0; i < CLMP_LIMB_LEN; i++) {
        carry += (__CLMP_CAST(clmp_uint64, _lhs[i]) & 0xFFFFFFFF) + (__CLMP_CAST(clmp_uint64, _rhs[i]) & 0xFFFFFFFF);
        _res[i] = __CLMP_CAST(clmp_limb_t, carry);
        carry >>= 0x20;
    }
}

void clmpz_add_limb_ui(clmp_limb_t *_res, const clmp_limb_t *_lhs, const clmp_uint32 _rhs) {
    clmp_uint64 carry = __CLMP_CAST(clmp_uint64, _rhs) & 0xFFFFFFFF;

    for (clmp_size_t i = 0; i < CLMP_LIMB_LEN; i++) {
        carry += (__CLMP_CAST(clmp_uint64, _lhs[i]) & 0xFFFFFFFF);
        _res[i] = __CLMP_CAST(clmp_limb_t, carry);
        carry >>= 0x20;
    }
}

void clmpz_add_limb_ul(clmp_limb_t *_res, const clmp_limb_t *_lhs, const clmp_uint64 _rhs) {
    clmp_uint64 carry = (__CLMP_CAST(clmp_uint64, _lhs[0]) & 0xFFFFFFFF) + (_rhs & 0xFFFFFFFF);
    _res[0] = __CLMP_CAST(clmp_limb_t, carry);
    carry >>= 0x20;

    carry += (__CLMP_CAST(clmp_uint64, _lhs[1]) & 0xFFFFFFFF) + ((_rhs >> 0x20) & 0xFFFFFFFF);
    _res[1] = __CLMP_CAST(clmp_limb_t, carry);
    carry >>= 0x20;

    for (clmp_size_t i = 2; i < CLMP_LIMB_LEN; i++) {
        carry += (__CLMP_CAST(clmp_uint64, _lhs[i]) & 0xFFFFFFFF);
        _res[i] = __CLMP_CAST(clmp_limb_t, carry);
        carry >>= 0x20;
    }
}

void clmpz_add(clmpz_ptr _res, clmpz_srcptr _lhs, clmpz_srcptr _rhs) {
    clmp_limb_t res[CLMP_LIMB_LEN], lhs[CLMP_LIMB_LEN], rhs[CLMP_LIMB_LEN];

    clmp_memcpy_global_to_local(lhs, _lhs->__limb, CLMP_LIMB_BYTES);
    clmp_memcpy_global_to_local(rhs, _rhs->__limb, CLMP_LIMB_BYTES);

    clmpz_add_limb(res, lhs, rhs);

    clmp_memcpy_local_to_global(_res->__limb, res, CLMP_LIMB_BYTES);
}

void clmpz_add_ui(clmpz_ptr _res, clmpz_srcptr _lhs, const clmp_uint32 _rhs) {
    clmp_limb_t res[CLMP_LIMB_LEN], lhs[CLMP_LIMB_LEN];

    clmp_memcpy_global_to_local(lhs, _lhs->__limb, CLMP_LIMB_BYTES);

    clmpz_add_limb_ui(res, lhs, _rhs);

    clmp_memcpy_local_to_global(_res->__limb, res, CLMP_LIMB_BYTES);
}

void clmpz_add_ul(clmpz_ptr _res, clmpz_srcptr _lhs, const clmp_uint64 _rhs) {
    clmp_limb_t res[CLMP_LIMB_LEN], lhs[CLMP_LIMB_LEN];

    clmp_memcpy_global_to_local(lhs, _lhs->__limb, CLMP_LIMB_BYTES);

    clmpz_add_limb_ul(res, lhs, _rhs);

    clmp_memcpy_local_to_global(_res->__limb, res, CLMP_LIMB_BYTES);
}

void clmpz_sub_limb(clmp_limb_t *_res, const clmp_limb_t *_lhs, const clmp_limb_t *_rhs) {
    if (_res != _lhs) {
        clmp_memcpy_local_to_local(_res, _lhs, CLMP_LIMB_BYTES);
    }

    clmp_limb_t carry = 1;
    clmp_limb_t tmp;

    for (clmp_size_t i = 0; i < CLMP_LIMB_LEN; i++) {
        tmp = ~_rhs[i];
        _res[i] += tmp;

        if (carry) {
            if (_res[i] >= tmp) {
                carry = 0;
            }

            if (_res[i] == 0xFFFFFFFF) {
                carry = 1;
            }

            _res[i]++;
        }
        else {
            if (_res[i] < tmp) {
                carry = 1;
            }
        }
    }
}

void clmpz_sub_limb_ui(clmp_limb_t *_res, const clmp_limb_t *_lhs, const clmp_uint32 _rhs) {
    if (_res != _lhs) {
        clmp_memcpy_local_to_local(_res, _lhs, CLMP_LIMB_BYTES);
    }

    clmp_limb_t carry = 1;
    clmp_limb_t tmp;

    tmp = ~(__CLMP_CAST(clmp_limb_t, _rhs));
    _res[0] += tmp;
    _res[0]++;

    if (_res[0] >= tmp) {
        carry = 0;
    }

    if (_res[0] == 0xFFFFFFFF) {
        carry = 1;
    }

    for (clmp_size_t i = 1; i < CLMP_LIMB_LEN; i++) {
        tmp = ~0x00000000;
        _res[i] += tmp;

        if (carry) {
            if (_res[i] >= tmp) {
                carry = 0;
            }

            if (_res[i] == 0xFFFFFFFF) {
                carry = 1;
            }

            _res[i]++;
        }
        else {
            if (_res[i] < tmp) {
                carry = 1;
            }
        }
    }
}
void clmpz_sub_limb_ul(clmp_limb_t *_res, const clmp_limb_t *_lhs, const clmp_uint64 _rhs) {
    if (_res != _lhs) {
        clmp_memcpy_local_to_local(_res, _lhs, CLMP_LIMB_BYTES);
    }

    clmp_limb_t carry = 1;
    clmp_limb_t tmp;

    tmp = ~(__CLMP_CAST(clmp_limb_t, _rhs));
    _res[0] += tmp;
    _res[0]++;

    if (_res[0] >= tmp) {
        carry = 0;
    }

    if (_res[0] == 0xFFFFFFFF) {
        carry = 1;
    }

    tmp = ~(__CLMP_CAST(clmp_limb_t, _rhs >> 32));
    _res[1] += tmp;

    if (carry) {
        if (_res[1] >= tmp) {
            carry = 0;
        }

        if (_res[1] == 0xFFFFFFFF) {
            carry = 1;
        }

        _res[1]++;
    }
    else {
        if (_res[1] < tmp) {
            carry = 1;
        }
    }

    for (clmp_size_t i = 2; i < CLMP_LIMB_LEN; i++) {
        tmp = ~0x00000000;
        _res[i] += tmp;

        if (carry) {
            if (_res[i] >= tmp) {
                carry = 0;
            }

            if (_res[i] == 0xFFFFFFFF) {
                carry = 1;
            }

            _res[i]++;
        }
        else {
            if (_res[i] < tmp) {
                carry = 1;
            }
        }
    }
}

void clmpz_sub(clmpz_ptr _res, clmpz_srcptr _lhs, clmpz_srcptr _rhs) {
    clmp_limb_t res[CLMP_LIMB_LEN], lhs[CLMP_LIMB_LEN], rhs[CLMP_LIMB_LEN];

    clmp_memcpy_global_to_local(lhs, _lhs->__limb, CLMP_LIMB_BYTES);
    clmp_memcpy_global_to_local(rhs, _rhs->__limb, CLMP_LIMB_BYTES);

    clmpz_sub_limb(res, lhs, rhs);

    clmp_memcpy_local_to_global(_res->__limb, res, CLMP_LIMB_BYTES);
}

void clmpz_sub_ui(clmpz_ptr _res, clmpz_srcptr _lhs, clmp_uint32 _rhs) {
    clmp_limb_t res[CLMP_LIMB_LEN], lhs[CLMP_LIMB_LEN];

    clmp_memcpy_global_to_local(lhs, _lhs->__limb, CLMP_LIMB_BYTES);

    clmpz_sub_limb_ui(res, lhs, _rhs);

    clmp_memcpy_local_to_global(_res->__limb, res, CLMP_LIMB_BYTES);
}

void clmpz_sub_ul(clmpz_ptr _res, clmpz_srcptr _lhs, clmp_uint64 _rhs) {
    clmp_limb_t res[CLMP_LIMB_LEN], lhs[CLMP_LIMB_LEN];

    clmp_memcpy_global_to_local(lhs, _lhs->__limb, CLMP_LIMB_BYTES);

    clmpz_sub_limb_ul(res, lhs, _rhs);

    clmp_memcpy_local_to_global(_res->__limb, res, CLMP_LIMB_BYTES);
}

void clmpz_mul_limb(clmp_limb_t *_res, const clmp_limb_t *_lhs, const clmp_limb_t *_rhs) {
    clmp_limb_t res[CLMP_LIMB_LEN * 2];

    clmp_memset_local(res, 0, CLMP_LIMB_BYTES * 2);

    for (clmp_size_t i = 0; i < CLMP_LIMB_LEN; i++) {
        clmp_uint64 carry = 0;
        for (clmp_size_t j = 0; j < CLMP_LIMB_LEN; j++) {
            carry += __CLMP_CAST(clmp_uint64, _lhs[i]) * __CLMP_CAST(clmp_uint64, _rhs[j]) + __CLMP_CAST(clmp_uint64, res[i + j]);
            res[i + j] = __CLMP_CAST(clmp_limb_t, carry);
            carry >>= 32;
        }
    }

    clmp_memcpy_local_to_local(_res, res, CLMP_LIMB_BYTES);
}

void clmpz_mul(clmpz_ptr _res, clmpz_srcptr _lhs, clmpz_srcptr _rhs) {
    clmp_limb_t lhs[CLMP_LIMB_LEN], rhs[CLMP_LIMB_LEN], res[CLMP_LIMB_LEN];

    clmp_int32 sign_lhs = clmpz_sign(_lhs);
    clmp_int32 sign_rhs = clmpz_sign(_rhs);

    clmp_memcpy_global_to_local(lhs, _lhs->__limb, CLMP_LIMB_BYTES);
    clmp_memcpy_global_to_local(rhs, _rhs->__limb, CLMP_LIMB_BYTES);

    if (sign_lhs == CLMP_SIGN_NEG) {
        clmpz_neg_limb(lhs, lhs);
    }

    if (sign_rhs == CLMP_SIGN_NEG) {
        clmpz_neg_limb(rhs, rhs);
    }

    clmpz_mul_limb(res, lhs, rhs);

    if (sign_lhs * sign_rhs == CLMP_SIGN_NEG) {
        clmpz_neg_limb(res, res);
    }

    clmp_memcpy_local_to_global(_res->__limb, res, CLMP_LIMB_BYTES);
}

void clmpz_mul_ui(clmpz_ptr _res, clmpz_srcptr _lhs, clmp_uint32 _rhs) {
    clmp_limb_t lhs[CLMP_LIMB_LEN], rhs[CLMP_LIMB_LEN], res[CLMP_LIMB_LEN];

    clmp_int32 sign_lhs = clmpz_sign(_lhs);
    clmp_int32 sign_rhs = CLMP_SIGN_POS;

    clmp_memcpy_global_to_local(lhs, _lhs->__limb, CLMP_LIMB_BYTES);
    clmp_memset_local(rhs, 0, CLMP_LIMB_BYTES);
    rhs[0] = __CLMP_CAST(clmp_limb_t, _rhs);

    if (sign_lhs == CLMP_SIGN_NEG) {
        clmpz_neg_limb(lhs, lhs);
    }

    clmpz_mul_limb(res, lhs, rhs);

    if (sign_lhs * sign_rhs == CLMP_SIGN_NEG) {
        clmpz_neg_limb(res, res);
    }

    clmp_memcpy_local_to_global(_res->__limb, res, CLMP_LIMB_BYTES);
}


void clmpz_mul_ul(clmpz_ptr _res, clmpz_srcptr _lhs, clmp_uint64 _rhs) {
    clmp_limb_t lhs[CLMP_LIMB_LEN], rhs[CLMP_LIMB_LEN], res[CLMP_LIMB_LEN];

    clmp_memset_local(res, 0, CLMP_LIMB_BYTES * 2);

    clmp_int32 sign_lhs = clmpz_sign(_lhs);
    clmp_int32 sign_rhs = CLMP_SIGN_POS;

    clmp_memcpy_global_to_local(lhs, _lhs->__limb, CLMP_LIMB_BYTES);
    clmp_memset_local(rhs, 0, CLMP_LIMB_BYTES);
    rhs[0] = __CLMP_CAST(clmp_limb_t, _rhs);
    rhs[1] = __CLMP_CAST(clmp_limb_t, _rhs >> 32);

    if (sign_lhs == CLMP_SIGN_NEG) {
        clmpz_neg_limb(lhs, lhs);
    }

    clmpz_mul_limb(res, lhs, rhs);

    if (sign_lhs * sign_rhs == CLMP_SIGN_NEG) {
        clmpz_neg_limb(res, res);
    }

    clmp_memcpy_local_to_global(_res->__limb, res, CLMP_LIMB_BYTES);
}

void clmpz_div_limb(clmp_limb_t *_res, const clmp_limb_t *_lhs, const clmp_limb_t *_rhs) {
    clmp_limb_t quot[CLMP_LIMB_LEN], restore[CLMP_LIMB_LEN];
    clmp_limb_t clmp_one[CLMP_LIMB_LEN];

    clmp_memset_local(quot, 0, CLMP_LIMB_BYTES);
    clmp_memset_local(restore, 0, CLMP_LIMB_BYTES);

    clmp_memset_local(clmp_one, 0, CLMP_LIMB_BYTES);
    clmp_one[0] = 1;

    for (clmp_bitcnt_t i = clmpz_bits_limb(_lhs); i > -1 ; i--) {
        clmp_limb_t shift[CLMP_LIMB_LEN], tmp_r[CLMP_LIMB_LEN];

        clmpz_shift_left_limb_ui(shift, _rhs, i);
        clmpz_add_limb(tmp_r, restore, shift);

        if (clmpz_cmp_limb(_lhs, tmp_r) > -1) {
            clmp_limb_t tmp_q[CLMP_LIMB_LEN];

            clmpz_add_limb(restore, restore, shift);
            clmpz_shift_left_limb_ui(tmp_q, clmp_one, i);
            clmpz_bit_or_limb(quot, quot, tmp_q);
        }
    }

    clmp_memcpy_local_to_local(_res, quot, CLMP_LIMB_BYTES);
}

void clmpz_div(clmpz_ptr _res, clmpz_srcptr _lhs, clmpz_srcptr _rhs) {
    clmp_limb_t lhs[CLMP_LIMB_LEN], rhs[CLMP_LIMB_LEN], res[CLMP_LIMB_LEN];

    clmp_int32 sign_lhs = clmpz_sign(_lhs);
    clmp_int32 sign_rhs = clmpz_sign(_rhs);

    clmp_memcpy_global_to_local(lhs, _lhs->__limb, CLMP_LIMB_BYTES);
    clmp_memcpy_global_to_local(rhs, _rhs->__limb, CLMP_LIMB_BYTES);

    if (sign_lhs == CLMP_SIGN_NEG) {
        clmpz_neg_limb(lhs, lhs);
    }

    if (sign_rhs == CLMP_SIGN_NEG) {
        clmpz_neg_limb(rhs, rhs);
    }

    clmpz_div_limb(res, lhs, rhs);

    if (sign_lhs * sign_rhs == CLMP_SIGN_NEG) {
        clmpz_neg_limb(res, res);
    }

    clmp_memcpy_local_to_global(_res->__limb, res, CLMP_LIMB_BYTES);
}

void clmpz_div_ui(clmpz_ptr _res, clmpz_srcptr _lhs, clmp_uint32 _rhs) {
    clmp_limb_t lhs[CLMP_LIMB_LEN], rhs[CLMP_LIMB_LEN], res[CLMP_LIMB_LEN];

    clmp_int32 sign_lhs = clmpz_sign(_lhs);

    clmp_memcpy_global_to_local(lhs, _lhs->__limb, CLMP_LIMB_BYTES);
    clmp_memset_local(rhs, 0, CLMP_LIMB_BYTES);
    rhs[0] = __CLMP_CAST(clmp_limb_t, _rhs);

    if (sign_lhs == CLMP_SIGN_NEG) {
        clmpz_neg_limb(lhs, lhs);
    }

    clmpz_div_limb(res, lhs, rhs);

    if (sign_lhs == CLMP_SIGN_NEG) {
        clmpz_neg_limb(res, res);
    }

    clmp_memcpy_local_to_global(_res->__limb, res, CLMP_LIMB_BYTES);
}

void clmpz_div_ul(clmpz_ptr _res, clmpz_srcptr _lhs, clmp_uint64 _rhs) {
    clmp_limb_t lhs[CLMP_LIMB_LEN], rhs[CLMP_LIMB_LEN], res[CLMP_LIMB_LEN];

    clmp_int32 sign_lhs = clmpz_sign(_lhs);

    clmp_memcpy_global_to_local(lhs, _lhs->__limb, CLMP_LIMB_BYTES);
    clmp_memset_local(rhs, 0, CLMP_LIMB_BYTES);
    rhs[0] = __CLMP_CAST(clmp_limb_t, _rhs);
    rhs[1] = __CLMP_CAST(clmp_limb_t, _rhs >> 32);

    if (sign_lhs == CLMP_SIGN_NEG) {
        clmpz_neg_limb(lhs, lhs);
    }

    clmpz_div_limb(res, lhs, rhs);

    if (sign_lhs == CLMP_SIGN_NEG) {
        clmpz_neg_limb(res, res);
    }

    clmp_memcpy_local_to_global(_res->__limb, res, CLMP_LIMB_BYTES);
}

void clmpz_mod_limb(clmp_limb_t *_res, const clmp_limb_t *_lhs, const clmp_limb_t *_rhs) {
    clmp_limb_t quot[CLMP_LIMB_LEN], restore[CLMP_LIMB_LEN];
    clmp_limb_t clmp_one[CLMP_LIMB_LEN];

    clmp_memset_local(quot, 0, CLMP_LIMB_BYTES);
    clmp_memset_local(restore, 0, CLMP_LIMB_BYTES);

    clmp_memset_local(clmp_one, 0, CLMP_LIMB_BYTES);
    clmp_one[0] = 1;

    for (clmp_int64 i = clmpz_bits_limb(_lhs); i > -1 ; i--) {
        clmp_limb_t shift[CLMP_LIMB_LEN], tmp_r[CLMP_LIMB_LEN];

        clmpz_shift_left_limb_ui(shift, _rhs, i);
        clmpz_add_limb(tmp_r, restore, shift);

        if (clmpz_cmp_limb(_lhs, tmp_r) > -1) {
            //printf("%ld\n", i);
            clmp_limb_t tmp_q[CLMP_LIMB_LEN];

            clmpz_add_limb(restore, restore, shift);
            clmpz_shift_left_limb_ui(tmp_q, clmp_one, i);
            clmpz_bit_or_limb(quot, quot, tmp_q);
        }
    }

    //printf("%u\n", restore[0]);

    clmpz_sub_limb(_res, _lhs, restore);
}

void clmpz_mod(clmpz_ptr _res, clmpz_srcptr _lhs, clmpz_srcptr _rhs) {
    clmp_limb_t lhs[CLMP_LIMB_LEN], rhs[CLMP_LIMB_LEN], res[CLMP_LIMB_LEN];

    clmp_int32 sign_lhs = clmpz_sign(_lhs);
    clmp_int32 sign_rhs = clmpz_sign(_rhs);

    clmp_memcpy_global_to_local(lhs, _lhs->__limb, CLMP_LIMB_BYTES);
    clmp_memcpy_global_to_local(rhs, _rhs->__limb, CLMP_LIMB_BYTES);

    if (sign_lhs == CLMP_SIGN_NEG) {
        clmpz_neg_limb(lhs, lhs);
    }

    if (sign_rhs == CLMP_SIGN_NEG) {
        clmpz_neg_limb(rhs, rhs);
    }

    clmpz_mod_limb(res, lhs, rhs);

    if (sign_lhs == CLMP_SIGN_NEG) {
        clmpz_neg_limb(res, res);

        if (clmpz_is_zero_limb(res) == CLMP_BOOL_FALSE) {
            clmpz_add_limb(res, res, rhs);
        }
    }

    clmp_memcpy_local_to_global(_res->__limb, res, CLMP_LIMB_BYTES);
}

void clmpz_mod_ui(clmpz_ptr _res, clmpz_srcptr _lhs, clmp_uint32 _rhs) {
    clmp_limb_t lhs[CLMP_LIMB_LEN], rhs[CLMP_LIMB_LEN], res[CLMP_LIMB_LEN];

    clmp_int32 sign_lhs = clmpz_sign(_lhs);

    clmp_memcpy_global_to_local(lhs, _lhs->__limb, CLMP_LIMB_BYTES);
    clmp_memset_local(rhs, 0, CLMP_LIMB_BYTES);
    rhs[0] = __CLMP_CAST(clmp_limb_t, _rhs);

    if (sign_lhs == CLMP_SIGN_NEG) {
        clmpz_neg_limb(lhs, lhs);
    }

    clmpz_mod_limb(res, lhs, rhs);

    if (sign_lhs == CLMP_SIGN_NEG) {
        clmpz_neg_limb(res, res);

        if (clmpz_is_zero_limb(res) == CLMP_BOOL_FALSE) {
            clmpz_add_limb(res, res, rhs);
        }
    }

    clmp_memcpy_local_to_global(_res->__limb, res, CLMP_LIMB_BYTES);
}

void clmpz_mod_ul(clmpz_ptr _res, clmpz_srcptr _lhs, clmp_uint64 _rhs) {
    clmp_limb_t lhs[CLMP_LIMB_LEN], rhs[CLMP_LIMB_LEN], res[CLMP_LIMB_LEN];

    clmp_int32 sign_lhs = clmpz_sign(_lhs);

    clmp_memcpy_global_to_local(lhs, _lhs->__limb, CLMP_LIMB_BYTES);
    clmp_memset_local(rhs, 0, CLMP_LIMB_BYTES);
    rhs[0] = __CLMP_CAST(clmp_limb_t, _rhs);
    rhs[1] = __CLMP_CAST(clmp_limb_t, _rhs >> 32);

    if (sign_lhs == CLMP_SIGN_NEG) {
        clmpz_neg_limb(lhs, lhs);
    }

    clmpz_mod_limb(res, lhs, rhs);

    if (sign_lhs == CLMP_SIGN_NEG) {
        clmpz_neg_limb(res, res);

        if (clmpz_is_zero_limb(res) == CLMP_BOOL_FALSE) {
            clmpz_add_limb(res, res, rhs);
        }
    }

    clmp_memcpy_local_to_global(_res->__limb, res, CLMP_LIMB_BYTES);
}

void clmpz_accum_mul_add(clmpz_ptr _res, clmpz_srcptr _lhs, clmp_uint64 _rhs) {
    clmp_limb_t lhs[CLMP_LIMB_LEN], rhs[CLMP_LIMB_LEN], res[CLMP_LIMB_LEN];

    clmp_int32 sign_lhs = clmpz_sign(_lhs);

    clmp_memcpy_global_to_local(lhs, _lhs->__limb, CLMP_LIMB_BYTES);

    clmp_memset_local(rhs, 0, CLMP_LIMB_BYTES);
    rhs[0] = __CLMP_CAST(clmp_limb_t, _rhs);
    rhs[1] = __CLMP_CAST(clmp_limb_t, _rhs >> 32);

    clmp_memcpy_global_to_local(res, _res->__limb, CLMP_LIMB_BYTES);

    if (sign_lhs == CLMP_SIGN_NEG) {
        clmpz_neg_limb(lhs, lhs);
    }

    clmpz_mul_limb(lhs, lhs, rhs);

    if (sign_lhs == CLMP_SIGN_NEG) {
        clmpz_neg_limb(lhs, lhs);
    }

    clmpz_add_limb(res, res, lhs);
    /*
    for(size_t i = 0; i < CLMP_LIMB_LEN; i++) {
        printf("%llu ", res[i]);
    }
    printf("\n");*/
    clmp_memcpy_local_to_global(_res->__limb, res, CLMP_LIMB_BYTES);
}

clmp_uint64 clmpz_mod_ul_to_ul(clmpz_srcptr _lhs, clmp_uint64 _rhs) {
    clmp_limb_t lhs[CLMP_LIMB_LEN], rhs[CLMP_LIMB_LEN], res[CLMP_LIMB_LEN];

    clmp_int32 sign_lhs = clmpz_sign(_lhs);

    clmp_memcpy_global_to_local(lhs, _lhs->__limb, CLMP_LIMB_BYTES);
    clmp_memset_local(rhs, 0, CLMP_LIMB_BYTES);
    rhs[0] = __CLMP_CAST(clmp_limb_t, _rhs);
    rhs[1] = __CLMP_CAST(clmp_limb_t, _rhs >> 32);

    if (sign_lhs == CLMP_SIGN_NEG) {
        clmpz_neg_limb(lhs, lhs);
    }

    clmpz_mod_limb(res, lhs, rhs);

    if (sign_lhs == CLMP_SIGN_NEG) {
        clmpz_neg_limb(res, res);

        if (clmpz_is_zero_limb(res) == CLMP_BOOL_FALSE) {
            clmpz_add_limb(res, res, rhs);
        }
    }

    return __CLMP_CAST(clmp_uint64, res[0]) + (__CLMP_CAST(clmp_uint64, res[1]) << 32);
}

clmp_uint64 clmpz_mod_ul_to_ul_barrett(clmpz_srcptr _lhs, clmp_uint64 _rhs, clmp_uint64 _rhs_r, clmp_uint64 _k_bar_d) {
    clmp_limb_t lhs[CLMP_LIMB_LEN], rhs[CLMP_LIMB_LEN], rhs_r[CLMP_LIMB_LEN], res[CLMP_LIMB_LEN], tmp[CLMP_LIMB_LEN];

    clmp_int32 sign_lhs = clmpz_sign(_lhs);

    clmp_memcpy_global_to_local(lhs, _lhs->__limb, CLMP_LIMB_BYTES);

    clmp_memset_local(rhs, 0, CLMP_LIMB_BYTES);
    rhs[0] = __CLMP_CAST(clmp_limb_t, _rhs);
    rhs[1] = __CLMP_CAST(clmp_limb_t, _rhs >> 32);

    clmp_memset_local(rhs_r, 0, CLMP_LIMB_BYTES);
    rhs_r[0] = __CLMP_CAST(clmp_limb_t, _rhs_r);
    rhs_r[1] = __CLMP_CAST(clmp_limb_t, _rhs_r >> 32);

    if (sign_lhs == CLMP_SIGN_NEG) {
        clmpz_neg_limb(lhs, lhs);
    }

    clmpz_mul_limb(tmp, lhs, rhs_r);
    clmpz_shift_right_limb_ui(tmp, tmp, _k_bar_d);
    clmpz_mul_limb(tmp, tmp, rhs);
    clmpz_sub_limb(res, lhs, tmp);

    if (clmpz_cmp_limb(res, rhs) > -1) {
        clmpz_sub_limb(res, res, rhs);
    }

    if (sign_lhs == CLMP_SIGN_NEG) {
        clmpz_neg_limb(res, res);

        if (clmpz_is_zero_limb(res) == CLMP_BOOL_FALSE) {
            clmpz_add_limb(res, res, rhs);
        }
    }

    return __CLMP_CAST(clmp_uint64, res[0]) + (__CLMP_CAST(clmp_uint64, res[1]) << 32);
}
