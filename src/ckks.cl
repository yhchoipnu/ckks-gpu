#include "ckks/ckks.clh"

__constant uint128_t uint128_0 = {0, 0};
__constant uint128_t uint128_1 = {0, 1};

complex_double_t mult_complex_double(complex_double_t _lhs, complex_double_t _rhs) {
    return (complex_double_t)(_lhs.x * _rhs.x - _lhs.y * _rhs.y, _lhs.x * _rhs.y + _lhs.y * _rhs.x);
}

clmp_uint64 to_uint64_uint128(const uint128_t _rhs) {
    return _rhs.lower;
}

uint128_t to_uint128_uint64(const clmp_uint64 _rhs) {
    uint128_t res;

    res.lower = _rhs;
    res.upper = 0;

    return res;
}

uint128_t shift_right_uint128_with_uint64(const uint128_t _lhs, const clmp_uint64 _rhs) {
    uint128_t res;

    if (_rhs >= 128) {
        return uint128_0;
    }
    else if (_rhs == 64) {
        res.upper = 0;
        res.lower = _lhs.upper;

        return res;
    } else if (_rhs == 0) {
        return _lhs;
    } else if (_rhs < 64) {
        res.upper = (_lhs.upper >> _rhs);
        res.lower = (_lhs.upper << (64 - _rhs)) + (_lhs.lower >> _rhs);

        return res;
    } else if (128 > _rhs && _rhs > 64) {
        res.upper = 0;
        res.lower = _lhs.upper >> (_rhs - 64);

        return res;
    } else {
        return uint128_0;
    }
}

uint128_t add_uint128_with_uint128(const uint128_t _lhs, const uint128_t _rhs) {
    uint128_t res;
    res.upper = _lhs.upper + _rhs.upper;

    if ((_lhs.lower + _rhs.lower) < _lhs.lower) {
        res.upper += 1;
    }

    res.lower = _lhs.lower + _rhs.lower;

    return res;
}

uint128_t sub_uint128_with_uint128(const uint128_t _lhs, const uint128_t _rhs) {
    uint128_t res;
    res.upper = _lhs.upper - _rhs.upper;

    if ((_lhs.lower - _rhs.lower) > _lhs.lower) {
        res.upper -= 1;
    }

    res.lower = _lhs.lower - _rhs.lower;

    return res;
}

uint128_t mult_uint128_with_uint128(const uint128_t _lhs, const uint128_t _rhs) {
    uint128_t res;

    clmp_uint64 top[4] = {_lhs.upper >> 32, _lhs.upper & 0xffffffff, _lhs.lower >> 32, _lhs.lower & 0xffffffff};
    clmp_uint64 bottom[4] = {_rhs.upper >> 32, _rhs.upper & 0xffffffff, _rhs.lower >> 32, _rhs.lower & 0xffffffff};
    clmp_uint64 products[4][4];


    for (int i = 3; i > -1; i--) {
        for (int j = 3; j > -1; j--) {
            products[3-j][i] = top[j] * bottom[i];
        }
    }

    clmp_uint64 fourth = products[0][3] & 0xffffffff;
    clmp_uint64 third = (products[0][2] & 0xffffffff) + (products[0][3] >> 32);
    clmp_uint64 second = (products[0][1] & 0xffffffff) + (products[0][2] >> 32);
    clmp_uint64 first = (products[0][0] & 0xffffffff) + (products[0][1] >> 32);

    third += (products[1][3] & 0xffffffff);
    second += (products[1][2] & 0xffffffff) + (products[1][3] >> 32);
    first += (products[1][1] & 0xffffffff) + (products[1][2] >> 32);

    second += (products[2][3] & 0xffffffff);
    first += (products[2][2] & 0xffffffff) + (products[2][3] >> 32);

    first  += (products[3][3] & 0xffffffff);

    third += fourth >> 32;
    second += third >> 32;
    first += second >> 32;

    fourth &= 0xffffffff;
    third &= 0xffffffff;
    second &= 0xffffffff;
    first &= 0xffffffff;

    res.upper = first << 32 | second;
    res.lower = third << 32 | fourth;

    return res;
}

void scale_up_to_clmpz(clmpz_ptr _res, const double _val, const long _log_p) {
    clmp_double tmp = { .d = _val };

    if (_val == 0.0) {
        clmpz_set_i(_res, 0);
        return;
    }

    clmp_int32 d_sign = tmp.l & 0x8000000000000000 ? CLMP_SIGN_NEG : CLMP_SIGN_POS;
    clmp_exp_t d_exp = ((tmp.l >> 52) & 0x7FF) - 0x3ff;
    clmp_uint64 d_mant = 0x10000000000000 + (tmp.l & 0xFFFFFFFFFFFFF);

    clmp_exp_t shift = d_exp - 52 + _log_p;

    clmpz_set_ul(_res, d_mant);
    if (shift > 0) {
        clmpz_shift_left_ui(_res, _res, shift);
    }
    else if (shift < 0) {
        clmpz_shift_right_ui(_res, _res, -shift);
    }


    if (d_sign == CLMP_SIGN_NEG) {
        clmpz_neg(_res, _res);
    }
}

double scale_down_to_double(clmpz_srcptr _val, const long _log_p) {
    clmp_double tmp;
    clmp_limb_t val[CLMP_LIMB_LEN];

    if (clmpz_is_zero(_val)) {
        return tmp.d = 0.0;
    }

    clmp_memcpy_global_to_local(val, _val->__limb, CLMP_LIMB_BYTES);

    clmp_int32 d_sign = clmpz_sign(_val);

    if (d_sign == CLMP_SIGN_NEG) {
        clmpz_neg_limb(val, val);
    }

    clmp_exp_t d_exp = clmpz_bits_limb(val) - 1;
    clmp_exp_t shift = d_exp - 52;
    clmp_uint64 d_mant;

    if (shift > 0) {
        clmpz_shift_right_limb_ui(val, val, shift);
    }
    else if (shift < 0) {
        clmpz_shift_left_limb_ui(val, val, -shift);
    }

    d_mant = (__CLMP_CAST(clmp_uint64, val[1] & 0xfffff) << 32) + __CLMP_CAST(clmp_uint64, val[0]);

    tmp.l = d_sign == CLMP_SIGN_POS ? 0 : 0x8000000000000000;
    tmp.l += (d_exp + 0x3ff - _log_p) << 52;
    tmp.l += d_mant & 0xFFFFFFFFFFFFF;

    return tmp.d;
}

void permute_array(__global complex_double_t *_vals, const long _len) {
    for (long i = 1, j = 0; i < _len; i++) {
        long bit = _len >> 1;

        for (; j >= bit; bit >>= 1) {
            j ^= bit;
        }

        j ^= bit;

        if (i < j) {
            complex_double_t temp = _vals[i];
            _vals[i] = _vals[j];
            _vals[j] = temp;
        }

    }
}


void fft_special(__global complex_double_t *_vals, const long _len, const long _M, __global const long *_rot_group, __global const complex_double_t *_ksi_pows) {
    permute_array(_vals, _len);

    for (long l = 2; l <= _len; l <<= 1) {
        long lh = l >> 1;
        long lq = l << 2;
        long gap = _M / lq;

        for (long i = 0; i < _len; i += l) {
            for (long j = 0; j < lh; j++) {
                long index = (_rot_group[j] % lq) * gap;

                complex_double_t u = _vals[i + j];
                complex_double_t v = mult_complex_double(_vals[i + j + lh], _ksi_pows[index]);

                _vals[i + j] = u + v;
                _vals[i + j + lh] = u - v;
            }
        }
    }
}

void fft_special_inv(__global complex_double_t *_vals, const long _len, const long _M, __global const long *_rot_group, __global const complex_double_t *_ksi_pows) {
    for (long l = _len; l >= 1; l >>= 1) {
        long lh = l >> 1;
        long lq = l << 2;
        long gap = _M / lq;

        for (long i = 0; i < _len; i += l) {
            for (long j = 0; j < lh; j++) {
                long index = (lq - (_rot_group[j] % lq)) * gap;

                complex_double_t u = _vals[i + j] + _vals[i + j + lh];
                complex_double_t v = mult_complex_double(_vals[i + j] - _vals[i + j + lh], _ksi_pows[index]);

                _vals[i + j] = u;
                _vals[i + j + lh] = v;
            }
        }
    }

    permute_array(_vals, _len);

    for (long i = 0; i < _len; i++) {
        _vals[i] /= _len;
    }
}

__kernel void encode(clmpz_t *_res, __global complex_double_t *_vals, const long _len, const long _log_p, const long _Nh, const long _M, __global const long *_rot_group, __global const complex_double_t *_ksi_pows) {
    long g_id = get_global_id(0);

    long i, idx, jdx;
    long gap = _Nh / _len;
    long N = _Nh << 1;

    fft_special_inv(_vals + (g_id * _len), _len, _M, _rot_group, _ksi_pows);

    for (i = 0; i < N; i++) {
        clmpz_set_ui(_res[i + (g_id * N)], 0);
    }

    for (i = 0, idx = 0, jdx = _Nh; i < _len; i++, idx += gap, jdx += gap) {
        scale_up_to_clmpz(_res[idx + (g_id * N)], _vals[i + (g_id * _len)].x, _log_p);
        scale_up_to_clmpz(_res[jdx + (g_id * N)], _vals[i + (g_id * _len)].y, _log_p);
    }

}

__kernel void decode(__global complex_double_t *_res, clmpz_t *_vals, const long _len, const long _log_p, const long _log_q, const long _Nh, const long _M, __global const long *_rot_group, __global const complex_double_t *_ksi_pows, clmpz_t * _q_pows) {
    long g_id = get_global_id(0);

    long gap = _Nh / _len;
    long N = _Nh << 1;

    for (long i = 0, idx = 0; i < _len; i++, idx += gap) {
        clmpz_mod(_vals[idx + (g_id * N)], _vals[idx + (g_id * N)], _q_pows[_log_q]);
        if (clmpz_bits(_vals[idx + (g_id * N)]) == _log_q) {
            clmpz_sub(_vals[idx + (g_id * N)], _vals[idx + (g_id * N)], _q_pows[_log_q]);
        }
        _res[i + (g_id * _len)].x = scale_down_to_double(_vals[idx + (g_id * N)], _log_p);


        clmpz_mod(_vals[idx + _Nh + (g_id * N)], _vals[idx + _Nh + (g_id * N)], _q_pows[_log_q]);
        if (clmpz_bits(_vals[idx + _Nh + (g_id * N)]) == _log_q) {
            clmpz_sub(_vals[idx + _Nh + (g_id * N)], _vals[idx + _Nh + (g_id * N)], _q_pows[_log_q]);
        }
        _res[i + (g_id * _len)].y = scale_down_to_double(_vals[idx + _Nh + (g_id * N)], _log_p);
    }

    fft_special(_res + (g_id * _len), _len, _M, _rot_group, _ksi_pows);
}

clmp_uint64 mul_mod_barrett(const clmp_uint64 _rhs, const clmp_uint64 _lhs, const clmp_uint64 _p, const clmp_uint64 _p_r, const long _k_bar_d) {
    clmp_uint64 res;

    uint128_t mul = mult_uint128_with_uint128(to_uint128_uint64(_lhs), to_uint128_uint64(_rhs));

    clmp_uint64 bot = to_uint64_uint128(mul);
    clmp_uint64 top = to_uint64_uint128(shift_right_uint128_with_uint64(mul, 64));

    uint128_t tmp = mult_uint128_with_uint128(to_uint128_uint64(bot), to_uint128_uint64(_p_r));
    tmp = shift_right_uint128_with_uint64(tmp, 64);
    tmp = add_uint128_with_uint128(tmp, mult_uint128_with_uint128(to_uint128_uint64(top), to_uint128_uint64(_p_r)));
    tmp = shift_right_uint128_with_uint64(tmp, _k_bar_d - 64);
    tmp = mult_uint128_with_uint128(tmp, to_uint128_uint64(_p));
    tmp = sub_uint128_with_uint128(mul, tmp);

    res = to_uint64_uint128(tmp);

    if (res >= _p) {
        res -= _p;
    }

    return res;
}

void butt(__global clmp_uint64 *_lhs, __global clmp_uint64 *_rhs, const clmp_uint64 _W, const clmp_uint64 _p, const clmp_uint64 _p_inv) {
    uint128_t U = mult_uint128_with_uint128(to_uint128_uint64(*_rhs), to_uint128_uint64(_W));

    clmp_uint64 U0 = to_uint64_uint128(U);
    clmp_uint64 U1 = to_uint64_uint128(shift_right_uint128_with_uint64(U, 64));

    clmp_uint64 Q = U0 * _p_inv;

    uint128_t Hx = mult_uint128_with_uint128(to_uint128_uint64(Q), to_uint128_uint64(_p));
    clmp_uint64 H = to_uint64_uint128(shift_right_uint128_with_uint64(Hx, 64));

    clmp_uint64 V = U1 < H ? U1 + _p - H : U1 - H;

    *_rhs = *_lhs < V ? *_lhs + _p - V : *_lhs - V;
    *_lhs += V;
    if (*_lhs > _p) {
        *_lhs -= _p;
    }
}

void ibutt(__global clmp_uint64 *_lhs, __global clmp_uint64 *_rhs, const clmp_uint64 _W, const clmp_uint64 _p, const clmp_uint64 _p_inv) {
    clmp_uint64 T = *_lhs < *_rhs ? *_lhs + _p - *_rhs : *_lhs - *_rhs;

    *_lhs += *_rhs;
    if (*_lhs > _p) {
        *_lhs -= _p;
    }

    uint128_t UU = mult_uint128_with_uint128(to_uint128_uint64(T), to_uint128_uint64(_W));
    clmp_uint64 U0 = to_uint64_uint128(UU);
    clmp_uint64 U1 = to_uint64_uint128(shift_right_uint128_with_uint64(UU, 64));

    clmp_uint64 Q = U0 * _p_inv;

    uint128_t Hx = mult_uint128_with_uint128(to_uint128_uint64(Q), to_uint128_uint64(_p));
    clmp_uint64 H = to_uint64_uint128(shift_right_uint128_with_uint64(Hx, 64));

    *_rhs = (U1 < H) ? U1 + _p - H: U1 - H;
}

void idivN(__global clmp_uint64 *_val, const clmp_uint64 _scaled_N_inv, const clmp_uint64 _p, const clmp_uint64 _p_inv) {
    uint128_t U = mult_uint128_with_uint128(to_uint128_uint64(*_val), to_uint128_uint64(_scaled_N_inv));
    clmp_uint64 U0 = to_uint64_uint128(U);
    clmp_uint64 U1 = to_uint64_uint128(shift_right_uint128_with_uint64(U, 64));

    clmp_uint64 Q = U0 * _p_inv;

    uint128_t Hx = mult_uint128_with_uint128(to_uint128_uint64(Q), to_uint128_uint64(_p));
    clmp_uint64 H = to_uint64_uint128(shift_right_uint128_with_uint64(Hx, 64));

    *_val = (U1 < H) ? U1 + _p - H : U1 - H;
}

void NTT(__global clmp_uint64 *_vals, const clmp_uint64 _p, const clmp_uint64 _p_inv, __global clmp_uint64 *_scaled_root_pows, const long _N, const long _log_N) {
    long t = _N;
    long logt1 = _log_N + 1;

    for (long m = 1; m < _N; m <<= 1) {
        t >>= 1;
        logt1 -= 1;

        for (long i = 0; i < m; i++) {
            long j1 = i << logt1;
            long j2 = j1 + t - 1;
            clmp_uint64 W = _scaled_root_pows[m + i];

            for (long j = j1; j <= j2; j++) {
                butt(&_vals[j], &_vals[j+t], W, _p, _p_inv);
            }
        }
    }
}

void iNTT(__global clmp_uint64 *_vals, const clmp_uint64 _p, const clmp_uint64 _p_inv, __global clmp_uint64 *_scaled_root_inv_pows, const clmp_uint64 _scaled_N_inv, const long _N) {
    long t = 1;

    for (long m = _N; m > 1; m >>= 1) {
        long j1 = 0;
        long h = m >> 1;

        for (long i = 0; i < h; i++) {
            long j2 = j1 + t -1;
            clmp_uint64 W = _scaled_root_inv_pows[h + i];

            for (long j = j1; j <= j2; j++) {
                ibutt(&_vals[j], &_vals[j+t], W, _p, _p_inv);
            }

            j1 += (t << 1);
        }
        t <<= 1;
    }

    for (long i = 0; i < _N; i++) {
        idivN(&_vals[i], _scaled_N_inv, _p, _p_inv);
    }
}

__kernel void add_n_equal(clmpz_t *_lhs, clmpz_t *_rhs, clmpz_t *_mod, clmp_uint64 _N) {
    long n = get_global_id(0);
    long b = get_global_id(1);

    clmpz_add(_lhs[b * _N + n], _lhs[b * _N + n], _rhs[b * _N + n]);
    clmpz_mod(_lhs[b * _N + n], _lhs[b * _N + n], *_mod);
}

__kernel void right_shift_n_equal(clmpz_t *_vals, clmp_uint64 _bits, clmp_uint64 _N, clmpz_t *_tmp) {
    long n = get_global_id(0);
    long b = get_global_id(1);

    clmpz_set_ui(_tmp[b * _N + n], 1);
    clmpz_shift_left_ui(_tmp[b * _N + n], _tmp[b * _N + n], __CLMP_CAST(clmp_uint32, _bits - 1));
    clmpz_add(_vals[b * _N + n], _vals[b * _N + n], _tmp[b * _N + n]);
    clmpz_shift_right_ui(_vals[b * _N + n], _vals[b * _N + n], __CLMP_CAST(clmp_uint32, _bits));
    /*
    clmp_limb_t vals[CLMP_LIMB_LEN], tmp[CLMP_LIMB_LEN];

    clmp_memcpy_global_to_local(vals, _vals[n]->__limb, CLMP_LIMB_BYTES);
    clmp_memset_local(tmp, 0, CLMP_LIMB_BYTES);
    tmp[0] = 1;

    clmpz_shift_left_limb_ui(tmp, tmp, __CLMP_CAST(clmp_uint32, _bits - 1));
    clmpz_add_limb(vals, vals, tmp);
    clmpz_shift_right_limb_ui(vals, vals,  __CLMP_CAST(clmp_uint32, _bits));

    clmp_memcpy_local_to_global(_vals[n]->__limb, vals, CLMP_LIMB_BYTES);
    */
}

__kernel void crt_decomposition(__global clmp_uint64 *_res, clmpz_t *_vals, const long _np, const clmp_uint64 _log_N, const clmp_uint64 _N, __global clmp_uint64 *_p_vec, __global clmp_uint64 *_p_r_vec, clmp_uint64 _k_bar_d) {
    long n = get_global_id(0);
    long b = get_global_id(1);
    //long i = get_global_id(2);

    for (long i = 0; i < _np; i++) {

        __global clmp_uint64 *res_i = _res + b * (_np << _log_N) + (i << _log_N);

        clmp_uint64 p_i = _p_vec[i];
        clmp_uint64 p_r_i = _p_r_vec[i];

        res_i[n] = clmpz_mod_ul_to_ul_barrett(_vals[b * _N + n], p_i, p_r_i, _k_bar_d);
        barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
    }
}

__kernel void crt_ntt(__global clmp_uint64 *_res, const long _log_N, const long _N, __global clmp_uint64 *_p_vec, __global clmp_uint64 *_p_inv_vec, const long _k_bar_d, __global clmp_uint64 *_scaled_root_pows) {
    long i = get_global_id(0);

    __global clmp_uint64 *res_i = _res + (i << _log_N);

    clmp_uint64 p_i = _p_vec[i];
    clmp_uint64 p_inv_i = _p_inv_vec[i];

    NTT(res_i, p_i, p_inv_i, _scaled_root_pows + i * _N, _N, _log_N);
}

__kernel void mult_decomposition(__global clmp_uint64 *_r_lhs, __global clmp_uint64 *_r_rhs, clmpz_t *_lhs, clmpz_t *_rhs, const clmp_uint64 _log_N, __global clmp_uint64 *_p_vec, __global clmp_uint64 *_p_r_vec, clmp_uint64 _k_bar_d) {
    long n = get_global_id(0);
    long i = get_global_id(1);

    __global clmp_uint64 *lhs_i = _r_lhs + (i << _log_N);
    __global clmp_uint64 *rhs_i = _r_rhs + (i << _log_N);

    clmp_uint64 p_i = _p_vec[i];
    clmp_uint64 p_r_i = _p_r_vec[i];

    lhs_i[n] = clmpz_mod_ul_to_ul_barrett(_lhs[n], p_i, p_r_i, _k_bar_d);
    rhs_i[n] = clmpz_mod_ul_to_ul_barrett(_rhs[n], p_i, p_r_i, _k_bar_d);
}

__kernel void reconstruct(clmpz_t *_res, __global clmp_uint64 *_vals, const long _np, clmpz_srcptr _q, const long _nprimes, clmpz_t *_p_hat, __global clmp_uint64 *_p_hat_inv_mod_p, clmpz_t *_p_prod, clmpz_t * _p_prod_h, __global clmp_uint64 *_p_vec, __global clmp_uint64 *_p_r_vec, const long _k_bar_d, const long _N, const long _log_N) {
    clmpz_t *p_hat_np = _p_hat + (_np - 1) * _nprimes;
    __global clmp_uint64 *p_hat_inv_mod_p_np = _p_hat_inv_mod_p + (_np - 1) * _nprimes;
    clmpz_t *p_prod_np = _p_prod + _np - 1;
    clmpz_t *p_prod_h_np = _p_prod_h + _np - 1;

    long n = get_global_id(0);

    clmpz_set_ui(_res[n], 0);

    for (long i = 0; i < _np; i++) {
        //printf("%lu\n", _vals[n + (i << _log_N)]);
        clmp_uint64 s = mul_mod_barrett(_vals[n + (i << _log_N)], *(p_hat_inv_mod_p_np + i), _p_vec[i], _p_r_vec[i], _k_bar_d);
        /*
        if (n == 0) {
            printf("%ld %ld %lu\n", i, n, s);
        }
        */
        clmpz_accum_mul_add(_res[n], *(p_hat_np + i), s);
    }
    /*
    if (n == 0) {
        for (clmp_uint64 i = 0; i < CLMP_LIMB_LEN; i++) {
            printf("%x ", _res[n]->__limb[i]);
        }
        printf("\n");
    }
    if (n == 0) {
        for (clmp_uint64 i = 0; i < CLMP_LIMB_LEN; i++) {
            printf("%x ", (*p_prod_np)->__limb[i]);
        }
        printf("\n");
    }*/
    //printf("%ld %lu \n", n, _res[n]->__limb[0]);
    //printf("%ld %u %u \n", n, _res[n]->__limb[0], (*p_prod_np)->__limb[0]);
    clmpz_mod(_res[n], _res[n], *p_prod_np);
    //printf("%ld %lu \n", n, _res[n]->__limb[0]);


    if (clmpz_cmp(_res[n], *p_prod_h_np) == 1) {
        clmpz_sub(_res[n], _res[n], *p_prod_np);
    }


    clmpz_mod(_res[n], _res[n], _q);

printf("");
}

__kernel void mult_ntt_decomposition(__global clmp_uint64 *_res, clmpz_t *_vals, const long _np, const clmp_uint64 _log_N, const clmp_uint64 _N, __global clmp_uint64 *_p_vec, __global clmp_uint64 *_p_r_vec, clmp_uint64 _k_bar_d) {
    long n = get_global_id(0);
    long b = get_global_id(1);
    //long i = get_global_id(2);

    for (long i = 0; i < _np; i++) {

        __global clmp_uint64 *res_i = _res + (i << _log_N);

        clmp_uint64 p_i = _p_vec[i];
        clmp_uint64 p_r_i = _p_r_vec[i];

        res_i[n] = clmpz_mod_ul_to_ul_barrett(_vals[b * _N + n], p_i, p_r_i, _k_bar_d);
        barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
    }
}

__kernel void reconstruct_batch(clmpz_t *_res, __global clmp_uint64 *_vals, const long _np, clmpz_srcptr _q, const long _nprimes, clmpz_t *_p_hat, __global clmp_uint64 *_p_hat_inv_mod_p, clmpz_t *_p_prod, clmpz_t * _p_prod_h, __global clmp_uint64 *_p_vec, __global clmp_uint64 *_p_r_vec, const long _k_bar_d, const long _N, const long _log_N) {
    clmpz_t *p_hat_np = _p_hat + (_np - 1) * _nprimes;
    __global clmp_uint64 *p_hat_inv_mod_p_np = _p_hat_inv_mod_p + (_np - 1) * _nprimes;
    clmpz_t *p_prod_np = _p_prod + _np - 1;
    clmpz_t *p_prod_h_np = _p_prod_h + _np - 1;

    long n = get_global_id(0);
    long b = get_global_id(1);

    clmpz_set_ui(_res[b * _N + n], 0);

    for (long i = 0; i < _np; i++) {
        clmp_uint64 s = mul_mod_barrett(_vals[b * (_np << _log_N) + n + (i << _log_N)], *(p_hat_inv_mod_p_np + i), _p_vec[i], _p_r_vec[i], _k_bar_d);
        clmpz_accum_mul_add(_res[b * _N + n], *(p_hat_np + i), s);
    }


    clmpz_mod(_res[b * _N + n], _res[b * _N + n], *p_prod_np);

    if (clmpz_cmp(_res[b * _N + n], *p_prod_h_np) == 1) {
        clmpz_sub(_res[b * _N + n], _res[b * _N + n], *p_prod_np);
    }
    clmpz_mod(_res[b * _N + n], _res[b * _N + n], _q);
}

__kernel void ctx_mult(__global clmp_uint64 *_res, __global clmp_uint64 *_lhs, __global clmp_uint64 *_rhs, const long _log_N, const long _N, __global clmp_uint64 *_p_vec, __global clmp_uint64 *_p_inv_vec, __global clmp_uint64 *_p_r_vec, const long _k_bar_d, __global clmp_uint64 *_scaled_root_pows, __global clmp_uint64 *_scaled_root_inv_pows, __global clmp_uint64 *_scaled_N_inv) {
    long i = get_global_id(0);

    __global clmp_uint64 *lhs_i = _lhs + (i << _log_N);
    __global clmp_uint64 *rhs_i = _rhs + (i << _log_N);
    __global clmp_uint64 *res_i = _res + (i << _log_N);

    clmp_uint64 p_i = _p_vec[i];
    clmp_uint64 p_inv_i = _p_inv_vec[i];
    clmp_uint64 p_r_i = _p_r_vec[i];


    NTT(lhs_i, p_i, p_inv_i, _scaled_root_pows + i * _N, _N, _log_N);
    NTT(rhs_i, p_i, p_inv_i, _scaled_root_pows + i * _N, _N, _log_N);

    for (long n = 0; n < _N; n++) {
        res_i[n] = mul_mod_barrett(lhs_i[n], rhs_i[n], p_i, p_r_i, _k_bar_d);
        //printf("%lu %lu %lu\n", lhs_i[n], rhs_i[n],res_i[n]);
        //printf("%lu %lu %lu\n", lhs_i[n], rhs_i[n], res_i[n]);
    }
    /*
    for (clmp_uint64 n = 0; n < _N; n++) {
        printf("%lu\n", res_i[n]);
    }
    */
    iNTT(res_i, p_i, p_inv_i, _scaled_root_inv_pows + i * _N, _scaled_N_inv[i], _N);
    //printf("%lu\n", res_i[0]);
    /*
    for (clmp_uint64 n = 0; n < _N; n++) {
        printf("%lu\n", res_i[n]);
    }
    */

}

__kernel void ctx_mult_decrypt(__global clmp_uint64 *_res, __global clmp_uint64 *_lhs, __global clmp_uint64 *_rhs, const long _np, const long _log_N, const long _N, __global clmp_uint64 *_p_vec, __global clmp_uint64 *_p_inv_vec, __global clmp_uint64 *_p_r_vec, const long _k_bar_d, __global clmp_uint64 *_scaled_root_pows, __global clmp_uint64 *_scaled_root_inv_pows, __global clmp_uint64 *_scaled_N_inv) {
    long i = get_global_id(0);
    long b = get_global_id(1);

    __global clmp_uint64 *lhs_i = _lhs + b * (_np << _log_N) + (i << _log_N);
    __global clmp_uint64 *rhs_i = _rhs + (i << _log_N);
    __global clmp_uint64 *res_i = _res + b * (_np << _log_N) + (i << _log_N);

    clmp_uint64 p_i = _p_vec[i];
    clmp_uint64 p_inv_i = _p_inv_vec[i];
    clmp_uint64 p_r_i = _p_r_vec[i];


    NTT(lhs_i, p_i, p_inv_i, _scaled_root_pows + i * _N, _N, _log_N);
    NTT(rhs_i, p_i, p_inv_i, _scaled_root_pows + i * _N, _N, _log_N);

    for (long n = 0; n < _N; n++) {
        res_i[n] = mul_mod_barrett(lhs_i[n], rhs_i[n], p_i, p_r_i, _k_bar_d);
    }

    iNTT(res_i, p_i, p_inv_i, _scaled_root_inv_pows + i * _N, _scaled_N_inv[i], _N);
}

__kernel void ctx_mult_partial_ntt(__global clmp_uint64 *_res, __global clmp_uint64 *_lhs, __global clmp_uint64 *_rhs, const long _np, const long _log_N, const long _N, __global clmp_uint64 *_p_vec, __global clmp_uint64 *_p_inv_vec, __global clmp_uint64 *_p_r_vec, const long _k_bar_d, __global clmp_uint64 *_scaled_root_pows, __global clmp_uint64 *_scaled_root_inv_pows, __global clmp_uint64 *_scaled_N_inv) {
    long i = get_global_id(0);
    long b = get_global_id(1);

    __global clmp_uint64 *lhs_i = _lhs + (_np << _log_N) * b + (i << _log_N);
    __global clmp_uint64 *rhs_i = _rhs  + (i << _log_N);
    __global clmp_uint64 *res_i = _res + (_np << _log_N) * b + (i << _log_N);

    clmp_uint64 p_i = _p_vec[i];
    clmp_uint64 p_inv_i = _p_inv_vec[i];
    clmp_uint64 p_r_i = _p_r_vec[i];
    if (i == 0) {
        for (clmp_uint64 n = 0; n < _N; n++) {
            printf("%lu %lu %lu\n", b, i, lhs_i[n]);
        }
        printf("\n");
    }
    if (i == 1) {
        for (clmp_uint64 n = 0; n < _N; n++) {
            printf("%lu %lu %lu\n", b, i, lhs_i[n]);
        }
        printf("\n");
    }


    NTT(lhs_i, p_i, p_inv_i, _scaled_root_pows + i * _N, _N, _log_N);

    for (clmp_uint64 n = 0; n < _N; n++) {
        printf("%lu %lu %lu\n", b, i, lhs_i[n]);
    }
    printf("\n");

    for (long n = 0; n < _N; n++) {
        res_i[n] = mul_mod_barrett(lhs_i[n], rhs_i[n], p_i, p_r_i, _k_bar_d);
        printf("%lu\n", res_i[n]);
    }
    printf("\n");


    iNTT(res_i, p_i, p_inv_i, _scaled_root_inv_pows + i * _N, _scaled_N_inv[i], _N);
    for (clmp_uint64 n = 0; n < _N; n++) {
        printf("%lu\n", res_i[n]);
    }
    printf("\n");
}

__kernel void sub_from_gauss_n_equal(clmpz_t *_res, clmpz_t *_mod, clmp_uint64 _seed, double _sigma) {
    long i = get_global_id(0);

    //for (clmp_uint64 i = 0; i < N; i += 2) {
        double r1 = (1 + (clmp_rand(_seed * (i + 1)) & 0xFFFFFFF)) / (__CLMP_CAST(double, 0xFFFFFFF) + 1);
        double r2 = (1 + (clmp_rand(_seed + (i + 1)) & 0xFFFFFFF)) / (__CLMP_CAST(double, 0xFFFFFFF) + 1);
        double theta = 2 * M_PI * r1;
        double rr = sqrt(-2.0 * log(r2)) * _sigma;

        if (i % 2 == 0) {
            clmpz_neg(_res[i], _res[i]);
            //clmpz_add_ul(_res[i], _res[i], __CLMP_CAST(clmp_int64, floor(rr * cos(theta) + 0.5)));
            clmpz_mod(_res[i], _res[i], *_mod);

            clmpz_neg(_res[i + 1], _res[i + 1]);
            //clmpz_add_ul(_res[i + 1], _res[i + 1], __CLMP_CAST(clmp_int64, floor(rr * sin(theta) + 0.5)));
            clmpz_mod(_res[i + 1], _res[i + 1], *_mod);
        }
    //}
}

__kernel void add_gauss_n_equal(clmpz_t *_res, clmpz_t *_mod, clmp_uint64 _seed, double _sigma, clmp_uint64 _N) {
    long i = get_global_id(0);
    long b = get_global_id(1);

    //for (clmp_uint64 i = 0; i < N; i += 2) {
        double r1 = (1 + (clmp_rand(_seed * (i + 1)) & 0xFFFFFFF)) / (__CLMP_CAST(double, 0xFFFFFFF) + 1);
        double r2 = (1 + (clmp_rand(_seed + (i + 1)) & 0xFFFFFFF)) / (__CLMP_CAST(double, 0xFFFFFFF) + 1);
        double theta = 2 * M_PI * r1;
        double rr = sqrt(-2.0 * log(r2)) * _sigma;

        if (i % 2 == 0) {
            clmpz_add_ul(_res[b * _N + i], _res[b * _N + i], __CLMP_CAST(clmp_int64, floor(rr * cos(theta) + 0.5)));
            clmpz_mod(_res[b * _N + i], _res[b * _N + i], *_mod);

            clmpz_add_ul(_res[b * _N + i + 1], _res[b * _N + i + 1], __CLMP_CAST(clmp_int64, floor(rr * sin(theta) + 0.5)));
            clmpz_mod(_res[b * _N + i + 1], _res[b * _N + i + 1], *_mod);
        }
    //}
}
/*
__kernel void reconstruct(clmpz_t *_res, __global clmp_uint64 *_vals, const long _np, clmpz_srcptr _q, const long _nprimes, clmpz_t *_p_hat, __global clmp_uint64 *_p_hat_inv_mod_p, clmpz_t *_p_prod, clmpz_t *_p_prod_h, __global clmp_uint64 *_p_vec, __global clmp_uint64 *_p_r_vec, const long _k_bar_d, const long _N, const long _log_N) {
    clmpz_t *p_hat_np = _p_hat + (_np - 1) * _nprimes;
    __global clmp_uint64 *p_hat_inv_mod_p_np = _p_hat_inv_mod_p + (_np - 1) * _nprimes;
    clmpz_t *p_prod_np = _p_prod + _np - 1;
    clmpz_t *p_prod_h_np = _p_prod_h + _np - 1;

    int n = get_global_id(0);

    clmpz_set_ui(_res[n], 0);

    for (long i = 0; i < _np; i++) {
        clmp_uint64 s = mul_mod_barrett(_vals[n + (i << _log_N)], *(p_hat_inv_mod_p_np + i), _p_vec[i], _p_r_vec[i], _k_bar_d);
        clmpz_accum_mul_add(_res[n], *(p_hat_np + i), s);
    }

    clmpz_mod(_res[n], _res[n], *p_prod_np);

    if (clmpz_cmp(_res[n], *p_prod_h_np) == 1) {
        clmpz_sub(_res[n], _res[n], *p_prod_np);
    }
    clmpz_mod(_res[n], _res[n], _q);
}*/
