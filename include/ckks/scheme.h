#ifndef CKKS_SCHEME_H
#define CKKS_SCHEME_H

#include <map>
#include <string>
#include <cmath>

#include "clmp/clmp.h"
#include "ckks/context.h"
#include "ckks/key.h"
#include "ckks/secretkey.h"
#include "ckks/param.h"
#include "ckks/plaintext.h"
#include "ckks/ciphertext.h"

namespace CKKS {
    static long ENCRYPTION = 0;
    static long MULTIPLICATION  = 1; 
    static long CONJUGATION = 2;

    class scheme{
    public:

        context &ctx;

        bool is_serialized;

        std::map<long, key *> key_map;

        std::map<long, std::string> ser_key_map;

        scheme(context &_ctx, secretkey &_secretkey, bool _is_serialized = false);
        ~scheme();

        void add_enc_key(secretkey &_secretkey);

        void encode(plaintext &_plain, std::complex<double> *_vals, clmp_uint64 _len, clmp_uint64 _num, clmp_uint64 _log_p, clmp_uint64 _log_q);
        std::complex<double>* decode(plaintext &_plain);

        void encrypt_plaintext(ciphertext &_cipher, plaintext &_plain);
        void decrypt_ciphertext(plaintext &_plain, secretkey &_s_key, ciphertext &_cipher);

        void encrypt(ciphertext &_cipher, std::complex<double> *_vals, clmp_uint64 _len, clmp_uint64 _num, clmp_uint64 _log_p, clmp_uint64 _log_q);
        std::complex<double> *decrypt(secretkey &_s_key, ciphertext &_cipher);
    };
}

#endif
