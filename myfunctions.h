/*
 * functions.h
 *
 *  Created on: Sep 4, 2016
 *      Author: hongyanwang
 */

#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>

#include <liboptarith/gcd/gcd_binary_l2r.h>
#include <liboptarith/gcd/gcd_brent.h>
#include <liboptarith/gcd/gcd_divrem.h>

#include <liboptarith/math32.h>
#include <liboptarith/math64.h>
#include <liboptarith/math_mpz.h>
#include <liboptarith/mpz_xgcd.h>
#include <liboptarith/primes.h>
#include <liboptarith/sqrtmodp_list.h>
#include <liboptarith/group.h>
#include <liboptarith/group_pow.h>
#include <liboptarith/s128_t.h>
#include <liboptarith/square_free.h>
#include <liboptarith/u128_t.h>


#include <libqform/s64_qform.h>
#include <libqform/s128_qform.h>
#include <libqform/qform_group.h>
#include <libqform/mpz_qform.h>
#include <libqform/gen_qform.h>
#include <libqform/dbreps/s64_pow_reps.h>

#define WITH_INDICES 1
#define FAC_TOTAL 100000
#define BLOCKSIZE 268435456

#ifndef ABS
#define ABS(X) (((X) < 0) ? (-(X)) : (X))
#endif

#ifndef MIN
#define MIN(X,Y) (((X) < (Y)) ? (X) : (Y))
#endif


#define xgcd_s32(s, t, u, v) xgcd_divrem_s32(s, t, u, v)
#define xgcd_left_s32(s, u, v) xgcd_left_divrem_s32(s, u, v)
#define xgcd_partial_s32(R1, R0, C1, C0, bound) xgcd_partial_divrem_s32(R1, R0, C1, C0, bound)


void prime_sieve(const long max_prime, int * primes);
int h_upper_bound(const long D);
void regular_sieve(const int max_prime, const long blocksize, int ** factors, const int * primes, const int flags);
void get_ideal_p(s64_qform_group_t* group, s64_qform_t* pform, const int p);
int legendre_symbol(long a, int p);
int order_p(long D, int p, int h, int *factors_of_h);
void hensel(mpz_t r0,long D,int p, int exp);
void get_unit(mpz_t *unit,long D,int p,int exp);
int isprime(long D, int* primes);

//int *dec_bin(int h,int* bi);
