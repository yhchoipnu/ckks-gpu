#ifndef CKKS_SECRETKEY_H
#define CKKS_SECRETKEY_H

#include "clmp/clmp.h"
#include "ckks/context.h"
#include "ckks/param.h"

namespace CKKS {
    class secretkey{
    public:

        clmpz_t *sx = new clmpz_t[N];

        secretkey(context &_ctx);
        ~secretkey();
    };
}

#endif
