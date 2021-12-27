#ifndef CKKS_CONTEXT_H
#define CKKS_CONTEXT_H

#include <iostream>
#include <fstream>
#include <complex>
#include <vector>
#include <cmath>
#include <random>

#ifdef __APPLE__
#include "OpenCL/cl.hpp"
#else
#include <CL/cl.hpp>
#endif
#include "clmp/clmp.h"
#include "ckks/param.h"

namespace CKKS {
    class context {
    public:
        clmpz_t Q;
        clmpz_t QQ;

        cl::Context *cl_ctx;
        cl::Program cl_program;
        cl::CommandQueue cl_queue;

        clmpz_t *q_pows = new clmpz_t[logQQ + 1];
        clmp_uint32 *q_pows_k_d = new clmp_uint32[logQQ + 1];
        clmpz_t *q_pows_r = new clmpz_t[logQQ + 1];
        clmp_uint64 *rot_group = new clmp_uint64[Nh];
        std::complex<double> *ksi_pows = new std::complex<double>[M + 1];

        clmp_uint64 *p_vec = new clmp_uint64[nprimes];
        clmp_uint64 *p_r_vec = new clmp_uint64[nprimes];
        clmp_uint64 *p_inv_vec = new clmp_uint64[nprimes];
        clmp_uint64 *scaled_root_pows = new clmp_uint64[nprimes * N];
        clmp_uint64 *scaled_root_inv_pows = new clmp_uint64[nprimes * N];
        clmp_uint64 *scaled_N_inv = new clmp_uint64[nprimes];

        clmpz_t *p_prod = new clmpz_t[nprimes];
        clmpz_t *p_prod_h = new clmpz_t[nprimes];
        clmpz_t *p_hat = new clmpz_t[nprimes * nprimes];
        clmp_uint64 *p_hat_inv_mod_p = new clmp_uint64[nprimes * nprimes];

        cl::Buffer cl_q_pows;
        cl::Buffer cl_rot_group;
        cl::Buffer cl_ksi_pows;

        cl::Buffer cl_p_vec;
        cl::Buffer cl_p_r_vec;
        cl::Buffer cl_p_inv_vec;
        cl::Buffer cl_scaled_root_pows;
        cl::Buffer cl_scaled_root_inv_pows;
        cl::Buffer cl_scaled_N_inv;

        cl::Buffer cl_p_prod;
        cl::Buffer cl_p_prod_h;
        cl::Buffer cl_p_hat;
        cl::Buffer cl_p_hat_inv_mod_p;

        context(const cl::Context *_cl_ctx = NULL);
        ~context();

        clmp_uint32 reverse_bit(const clmp_uint32 &_rhs);

        clmp_uint64 pow(const clmp_uint64 &_lhs, const clmp_uint64 &_rhs);
        clmp_uint64 inv(const clmp_uint64 &_rhs);

        clmp_uint64 mul_mod(const clmp_uint64 &_lhs, const clmp_uint64 &_rhs, const clmp_uint64 &_mod);
        clmp_uint64 pow_mod(const clmp_uint64 &_lhs, const clmp_uint64 &_rhs, const clmp_uint64 &_mod);
        clmp_uint64 inv_mod(const clmp_uint64 &_rhs, const clmp_uint64 &_mod);

        void find_prime_factors(std::vector<clmp_uint64> &_res, clmp_uint64 _rhs);
        clmp_uint64 find_primitive_root(const clmp_uint64 &_mod);
        clmp_uint64 find_Nth_root_of_unity(const long &_N, const clmp_uint64 &_mod);

        bool is_prime(const clmp_uint64 &_rhs, const clmp_size_t &_n_iter);

        void fft_special(std::complex<double> *_vals, const long _len);
        void fft_special_inv(std::complex<double> *_vals, const long _len);

        void encode(clmpz_t *_res, const std::complex<double> *_vals, const long _len, const long _num, const long _log_p);
        void decode(std::complex<double> *_res, const clmpz_t *_vals, const long _len, const long _num, const long _log_p, const long _log_q);

        void NTT(clmp_uint64 *_vals);
        void iNTT(clmp_uint64 *_vals);

        void add_n_equal(clmpz_t *_lhs, clmpz_t *_rhs, clmp_uint64 _n_batch, clmpz_srcptr _mod);
        void right_shift_n_equal(clmpz_t *_vals, clmp_uint64 _n_batch, clmp_uint64 _bits);

        void CRT(clmp_uint64 *_res, clmpz_t *_vals, const clmp_uint64 _n_batch, const long _np);

        void mult(clmpz_t *_res, const clmpz_t *_lhs, const clmpz_t *_rhs, const long _np, clmpz_srcptr _mod);
        void mult_decrypt(clmpz_t *_res, const clmpz_t *_lhs, const clmpz_t *_rhs, const long _np, clmp_uint64 _n_batch, clmpz_srcptr _mod);
        void mult_ntt(clmpz_t *_res, clmpz_t *_lhs, clmp_uint64 *_rhs, const clmp_uint64 _n_batch, const long _np, clmpz_srcptr _mod);

        void test(clmpz_t * x, clmpz_t * a, clmpz_t * b, long np, clmpz_srcptr mod);
        void reconstruct_cpu(clmpz_t *_res, clmp_uint64* _vals, const long _np, clmpz_srcptr _q);

        void NTT_cpu(clmp_uint64 *_vals, long);
        void iNTT_cpu(clmp_uint64 *_vals, long);
        void butt(clmp_uint64& a, clmp_uint64& b, clmp_uint64 W, clmp_uint64 p, clmp_uint64 pInv);
        void ibutt(clmp_uint64& a, clmp_uint64& b, clmp_uint64 W, clmp_uint64 p, clmp_uint64 pInv);
        void idivN(clmp_uint64& a, clmp_uint64 NScale, clmp_uint64 p, clmp_uint64 pInv);

        void sub_from_gauss_n_equal(clmpz_t *_res, clmpz_srcptr _mod);
        void sub_from_gauss_n_equal_cpu(clmpz_t *_res, clmpz_srcptr _mod);
        void add_gauss_n_equal(clmpz_t *_res, clmp_uint64 _n_batch, clmpz_srcptr _mod);
        void add_gauss_n_equal_cpu(clmpz_t *_res, clmpz_srcptr _mod);
        void sample_ZO(clmpz_t *_res);
        void sampleHWT(clmpz_t *_res);
        void sample_uniform(clmpz_t *_res, long _bits);
    };
}

#endif
