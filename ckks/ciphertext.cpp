#include "ckks/ciphertext.h"

namespace CKKS {
    ciphertext::ciphertext(clmp_uint64 _log_p, clmp_uint64 _log_q, clmp_uint64 _len, clmp_uint64 _num): log_p(_log_p), log_q(_log_q), len(_len), num(_num) {
        ax = new clmpz_t[N * this->num];
        bx = new clmpz_t[N * this->num];
    }

    ciphertext::ciphertext(const ciphertext &_origin): log_p(_origin.log_p), log_q(_origin.log_q), len(_origin.len), num(_origin.num) {
        for (clmp_uint64 i = 0; i < N * num; i++) {
            clmpz_set(ax[i], _origin.ax[i]);
            clmpz_set(bx[i], _origin.bx[i]);
        }
    }

    ciphertext::~ciphertext() {
        delete [] ax;
        delete [] bx;
    }

    void ciphertext::free() {
        delete [] ax;
        delete [] bx;
    }

    void ciphertext::init() {
        if ((ax != NULL) || (bx != NULL)) {
            free();
        }

        ax = new clmpz_t[N * this->num];
        bx = new clmpz_t[N * this->num];
    }
}
