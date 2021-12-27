#ifndef CLMP_H_
#define CLMP_H_

#include <stddef.h>
#include <string.h>

#include <string>

#include "clmp_param.h"

#define CLMP_LIMB_BYTES (CLMP_LIMB_LEN * sizeof(clmp_limb_t))
#define CLMP_LIMB_BITS  (CLMP_LIMB_BYTES * 8)
#define CLMP_NAIL_BITS  0
#define CLMP_NUMB_BITS  (CLMP_LIMB_BITS - CLMP_NAIL_BITS)
#define CLMP_NUMB_MASK  ((~ __CLMP_CAST (CLmp_limb_t, 0)) >> CLMP_NAIL_BITS)
#define CLMP_NUMB_MAX   CLMP_NUMB_MASK
#define CLMP_NAIL_MASK  (~ CLMP_NUMB_MASK)

#define CLMP_SIGN_POS   1
#define CLMP_SIGN_NEG   -1

#define CLMP_BOOL_TRUE  1
#define CLMP_BOOL_FALSE 0

#if defined (__cplusplus)
    #define __CLMP_CAST(type, expr) static_cast <type> (expr)
    #define __CLMP_THOTHROW
#else
    #define __CLMP_CAST(type, expr) ((type) (expr))
    #define __CLMP_NOTRHORW
#endif

typedef char                clmp_uint8;
typedef int                 clmp_int32;
typedef long                clmp_int64;
typedef unsigned int        clmp_uint32;
typedef unsigned long       clmp_uint64;

typedef clmp_uint32         clmp_limb_t;
typedef clmp_int32          clmp_limb_signed_t;

typedef clmp_int64          clmp_bitcnt_t;

typedef clmp_uint64         clmp_size_t;
typedef clmp_int64          clmp_exp_t;

typedef clmp_limb_t         *clmp_ptr;
typedef const clmp_limb_t   *clmp_srcptr;

typedef struct __attribute__ ((packed)) {
    clmp_limb_t __limb[CLMP_LIMB_LEN];
} __clmpz_struct;

typedef __clmpz_struct          *clmpz_ptr, clmpz_t[1];
typedef const __clmpz_struct    *clmpz_srcptr;

typedef union {
    double d;
    struct {
        clmp_uint32 mantissa1:32;
        clmp_uint32 mantissa0:20;
        clmp_uint32 exponent:11;
        clmp_uint32 negative:1;
    } parts ;
} __clmp_double_cast_struct;

typedef __clmp_double_cast_struct clmp_double;

clmp_uint32 clmp_rand_bits_ui(clmp_uint32);
clmp_uint32 clmp_rand_bnd_ui(clmp_uint32);

void clmpz_set(clmpz_ptr, clmpz_srcptr);
void clmpz_set_ui(clmpz_ptr, clmp_uint32);
void clmpz_set_ul(clmpz_ptr, clmp_uint64);
void clmpz_set_i(clmpz_ptr, clmp_int32);
void clmpz_set_li(clmpz_ptr, clmp_int64);

clmp_limb_t clmpz_bit(clmpz_srcptr, clmp_uint32);

clmp_bitcnt_t clmpz_bits(clmpz_srcptr);

void clmpz_bit_and(clmpz_ptr, clmpz_srcptr, clmpz_srcptr);
void clmpz_bit_or(clmpz_ptr, clmpz_srcptr, clmpz_srcptr);
void clmpz_bit_xor(clmpz_ptr, clmpz_srcptr, clmpz_srcptr);

void clmpz_shift_left_ui(clmpz_ptr, clmpz_srcptr, clmp_uint32);
void clmpz_shift_right_ui(clmpz_ptr, clmpz_srcptr, clmp_uint32);

clmp_int32 clmpz_sign(clmpz_srcptr);

clmp_int32 clmpz_cmp(clmpz_srcptr, clmpz_srcptr);

clmp_int32 clmpz_is_zero(clmpz_srcptr);

void clmpz_neg(clmpz_ptr, clmpz_srcptr);

void clmpz_add(clmpz_ptr, clmpz_srcptr, clmpz_srcptr);
void clmpz_add_ui(clmpz_ptr, clmpz_srcptr, clmp_uint32);
void clmpz_add_ul(clmpz_ptr, clmpz_srcptr, clmp_uint64);

void clmpz_sub(clmpz_ptr, clmpz_srcptr, clmpz_srcptr);
void clmpz_sub_ui(clmpz_ptr, clmpz_srcptr, clmp_uint32);
void clmpz_sub_ul(clmpz_ptr, clmpz_srcptr, clmp_uint64);

void clmpz_mul(clmpz_ptr, clmpz_srcptr, clmpz_srcptr);
void clmpz_mul_ui(clmpz_ptr, clmpz_srcptr, clmp_uint32);
void clmpz_mul_ul(clmpz_ptr, clmpz_srcptr, clmp_uint64);

void clmpz_div(clmpz_ptr, clmpz_srcptr, clmpz_srcptr);
void clmpz_div_ui(clmpz_ptr, clmpz_srcptr, clmp_uint32);
void clmpz_div_ul(clmpz_ptr, clmpz_srcptr, clmp_uint64);

void clmpz_mod(clmpz_ptr, clmpz_srcptr, clmpz_srcptr);
void clmpz_mod_ui(clmpz_ptr, clmpz_srcptr, clmp_uint32);
void clmpz_mod_ul(clmpz_ptr, clmpz_srcptr, clmp_uint64);

clmp_uint32 clmpz_to_ui(clmpz_srcptr);
clmp_uint64 clmpz_to_ul(clmpz_srcptr);

void clmpz_rand_bits(clmpz_ptr, clmp_uint64);

std::string clmpz_to_string(clmpz_srcptr);

#endif
