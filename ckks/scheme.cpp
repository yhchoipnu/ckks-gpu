#include "ckks/scheme.h"

#include <complex>
#include <iostream>

namespace CKKS{
    scheme::scheme(context &_ctx, secretkey &_secretkey, bool _is_serialized): ctx(_ctx), is_serialized(_is_serialized) {
        add_enc_key(_secretkey);
    }

    scheme::~scheme() {
        for (auto const & e: key_map) {
            delete e.second;
        }
    }

    void scheme::add_enc_key(secretkey &_secretkey) {
        clmpz_t *ax = new clmpz_t[N];
        clmpz_t *bx = new clmpz_t[N];

        long np = ceil((1 + logQQ + logN + 2) / (double)pbnd);
        ctx.sample_uniform(ax, logQQ);
        ctx.mult(bx, _secretkey.sx, ax, np, ctx.QQ);
        ctx.sub_from_gauss_n_equal(bx, ctx.QQ);

        key *pk = new key();

        ctx.CRT(pk->rax, ax, 1, nprimes);
        ctx.CRT(pk->rbx, bx, 1, nprimes);

        delete [] ax;
        delete [] bx;

        key_map.insert(std::pair<long, key *>(ENCRYPTION, pk));
    }

    void scheme::encode(plaintext &_plain, std::complex<double> *_vals, clmp_uint64 _len, clmp_uint64 _num, clmp_uint64 _log_p, clmp_uint64 _log_q) {
        _plain.log_p = _log_p;
        _plain.log_q = _log_q;
        _plain.len = _len;
        _plain.num = _num;

        _plain.init();

        ctx.encode(_plain.msg, _vals, _len, _num, _log_p + logQ);
    }

    std::complex<double> * scheme::decode(plaintext &_plain) {
        std::complex<double> *res = new std::complex<double>[_plain.len * _plain.num];
        ctx.decode(res, _plain.msg, _plain.len, _plain.num, _plain.log_p, _plain.log_q);
        return res;
    }

    void scheme::encrypt_plaintext(ciphertext &_cipher, plaintext &_plain) {
        _cipher.log_p = _plain.log_p;
        _cipher.log_q = _plain.log_q;
        _cipher.len = _plain.len;
        _cipher.num = _plain.num;

        _cipher.init();

        clmpz_t qQ, qQ_r;
        clmp_uint32 qQ_k_d = ctx.q_pows_k_d[_plain.log_q + logQ];
        clmpz_set(qQ, ctx.q_pows[_plain.log_q + logQ]);
        clmpz_set(qQ_r, ctx.q_pows_r[_plain.log_q + logQ]);

        clmpz_t *vx = new clmpz_t[N];
        ctx.sample_ZO(vx);
        /*
        for (clmp_uint64 i = 0; i < N; i++) {
            std::cout << clmpz_to_string(vx[i]) << std::endl;
        }
        std::cout << std::endl;
        */
        key *pk = key_map.at(ENCRYPTION);

        long np = ceil((1 + logQQ + logN + 2)/(double)pbnd);
        ctx.mult_ntt(_cipher.ax, vx, pk->rax, _plain.num, np, qQ);
        /*
        for (clmp_uint64 i = 0; i < N; i++) {
            std::cout << clmpz_to_string(_cipher.ax[i]) << std::endl;
        }
        std::cout << std::endl;
        */
        ctx.add_gauss_n_equal(_cipher.ax, _cipher.num, qQ);
        /*
        for (clmp_uint64 i = 0; i < N; i++) {
            std::cout << clmpz_to_string(_cipher.ax[i]) << std::endl;
        }
        std::cout << std::endl;
        */
        ctx.mult_ntt(_cipher.bx, vx, pk->rbx, _plain.num, np, qQ);
        ctx.add_gauss_n_equal(_cipher.bx, _cipher.num, qQ);
        delete [] vx;

        ctx.add_n_equal(_cipher.bx, _plain.msg, _cipher.num, qQ);

        ctx.right_shift_n_equal(_cipher.ax, _cipher.num, logQ);
        ctx.right_shift_n_equal(_cipher.bx, _cipher.num, logQ);
    }

    void scheme::decrypt_ciphertext(plaintext &_plain, secretkey &_s_key, ciphertext &_cipher) {
        clmpz_t q;
        clmpz_set(q, ctx.q_pows[_cipher.log_q]);

        _plain.log_p = _cipher.log_p;
        _plain.log_q = _cipher.log_q;
        _plain.len = _cipher.len;
        _plain.num = _cipher.num;

        _plain.init();

        long np = ceil((1 + _cipher.log_q + logN + 2)/(double)pbnd);
        ctx.mult_decrypt(_plain.msg, _cipher.ax, _s_key.sx, np, _plain.num, q);
        ctx.add_n_equal(_plain.msg, _cipher.bx, _plain.num, q);
    }

    void scheme::encrypt(ciphertext &_cipher, std::complex<double> *_vals, clmp_uint64 _len, clmp_uint64 _num, clmp_uint64 _log_p, clmp_uint64 _log_q) {
        plaintext plain;
        encode(plain, _vals, _len, _num, _log_p, _log_q);
        encrypt_plaintext(_cipher, plain);
    }

    std::complex<double> *scheme::decrypt(secretkey &_s_key, ciphertext &_cipher) {
        plaintext plain;
        decrypt_ciphertext(plain, _s_key, _cipher);
        return  decode(plain);
    }
}
