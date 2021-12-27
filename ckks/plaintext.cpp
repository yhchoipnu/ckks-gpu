#include "ckks/plaintext.h"

namespace CKKS {
    plaintext::plaintext(clmp_uint64 _log_p, clmp_uint64 _log_q, clmp_uint64 _len, clmp_uint64 _num): log_p(_log_p), log_q(_log_q), len(_len), num(_num) {
        msg = new clmpz_t[N * this->num];
    }

    plaintext::~plaintext() {
        delete [] msg;
    }

    void plaintext::free() {
        delete [] msg;
    }

    void plaintext::init() {
        if (msg != NULL) {
            this->free();
        }

        msg = new clmpz_t[N * this->num];
    }
}
