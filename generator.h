/*
 * generator.h
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

#ifndef ABS
#define ABS(X) (((X) < 0) ? (-(X)) : (X))
#endif

#ifndef MIN
#define MIN(X,Y) (((X) < (Y)) ? (X) : (Y))
#endif

#define xgcd_s32(s, t, u, v) xgcd_divrem_s32(s, t, u, v)
#define xgcd_left_s32(s, u, v) xgcd_left_divrem_s32(s, u, v)
#define xgcd_partial_s32(R1, R0, C1, C0, bound) xgcd_partial_divrem_s32(R1, R0, C1, C0, bound)

typedef struct{
 int64_t a;
 int64_t b;
 int32_t c;
}rela_gen; // the relative generator is equal to a+b*sqrt(D)/2c

static inline void rela_gen_init(rela_gen *gen){
 gen->a=0; 
 gen->b=0;
 gen->c=1;}

static inline void rela_gen_clear(rela_gen *gen){}

static inline void rela_gen_set(rela_gen *g, rela_gen *f){
 g->a=f->a; 
 g->b=f->b;
 g->c=f->c;
}

typedef struct{
 mpz_t a;
 mpz_t b;
 mpz_t c;
}rela_gen_mpz; // the relative generator is equal to a+b*sqrt(D)/2c

static inline void rela_gen_mpz_init(rela_gen_mpz *gen){
 mpz_init(gen->a); 
 mpz_init(gen->b); 
 mpz_init(gen->c); mpz_set_si(gen->c,1);
}

static inline void rela_gen_mpz_clear(rela_gen_mpz *gen){
 mpz_clear(gen->a);
 mpz_clear(gen->b);
 mpz_clear(gen->c);
}

static inline void rela_gen_mpz_set(rela_gen_mpz *g, rela_gen_mpz *f){
 mpz_set(g->a,f->a);
 mpz_set(g->b,f->b);
 mpz_set(g->c,f->c);
}

typedef struct{
 mpz_t a;
 mpz_t b;
}norm_mpz;

static inline void norm_mpz_init(norm_mpz *norm){
 mpz_init(norm->a); 
 mpz_init_set_si(norm->b,1);
}

static inline void norm_mpz_clear(norm_mpz *norm){
 mpz_clear(norm->a);
 mpz_clear(norm->b);
}

//multiply two generators f and g, save the result in r, DD<0
void rela_gen_mpz_mul(rela_gen_mpz *r, rela_gen_mpz *f, rela_gen_mpz *g, mpz_t D);
//find the norm of a generator
void gen_norm(norm_mpz *norm, rela_gen_mpz *gen, mpz_t D);
//multiply two norms and save the result in r
void norm_mul(norm_mpz *r, norm_mpz *m, norm_mpz *n);
//find the generator mod D
void rela_gen_mpz_modD(rela_gen_mpz *gen, mpz_t D);

void mpz_qform_compose_gen(mpz_qform_group_t* group, mpz_qform_t* R, const mpz_qform_t* A, const mpz_qform_t* B, rela_gen_mpz *gen);
void mpz_qform_square_gen(mpz_qform_group_t* group, mpz_qform_t* C, const mpz_qform_t* A, rela_gen_mpz *gen);
void red_gen_mpz(mpz_qform_group_t* group, mpz_qform_t* C, mpz_t s, rela_gen_mpz *gen);
void dec_bin(long ,int* bi);
void AAC_conj(mpz_qform_group_t *group, mpz_qform_t *pform, long p, int order, FILE* wf);
void aac(long D, long p, int order, FILE* wf);
void AAC_unit(mpz_qform_group_t *group, mpz_qform_t *pform, long p, int order, FILE* wf);
void unit(long D, long p, int order, FILE* wf);
long powmodp(long p, long a, long n);
long sqrt_modp(long p, long D);
//int s64_qform_is_primeform(s64_qform_group_t* group,s64_qform_t* form,const int p) ;
int mpz_qform_is_primeform_why(mpz_qform_group_t* group, mpz_qform_t* form, long p);
int order_mpz(long D, long p, int h, int *factors_of_h);
