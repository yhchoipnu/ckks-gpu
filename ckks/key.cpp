#include "ckks/key.h"

namespace CKKS{
    key::key() { }

    key::~key() {
        delete [] rax;
        delete [] rbx;
    }
}
