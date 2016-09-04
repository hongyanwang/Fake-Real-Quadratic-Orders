/*
 * functions.c
 *
 *  Created on: Sep 4, 2016
 *      Author: hongyanwang
 */

#include "myfunctions.h"

// This function is from Anton's C library
// https://github.com/amosunov/class_groups/blob/master/sieve.c
void prime_sieve(const long max_prime, int * primes)
{
	int i, j;

	primes[0] = 0;

	char is_prime[max_prime];
	memset(is_prime, 1, max_prime);

	for (i = 2; i < max_prime; i++)
	{
		if (is_prime[i])
		{
			primes[++primes[0]] = i;

			for (j = (i << 1); j < max_prime; j += i)
			{
				is_prime[j] = 0;
			}
		}
	}

	primes[primes[0] + 1] = 0x7FFFFFFF;
}

// This function is from Anton's C library
// https://github.com/amosunov/class_groups/blob/master/sieve.c
void regular_sieve(int max_prime, long blocksize, int ** factors, const int * primes, const int flags)
{
	long i, j;

	int * f;

	factors[0][0] = 1;
	factors[0][1] = 0;
	factors[1][0] = 1;
	factors[1][1] = 1;

	for (i = 0; i < blocksize; i++)
	{
		factors[i][0] = 0;
	}

	for (int i = 1, p = primes[i]; p < max_prime; p = primes[++i])
	{
		for (j = p; j < blocksize; j += p)
		{
			f = factors[j];
			f[++f[0]] = (flags & WITH_INDICES) ? i : p;
		}
	}
}

// This function is from Anton's C library
// https://github.com/amosunov/class_groups/blob/master/functions.c
int h_upper_bound(const long D)
{
	double E;
	if ((D & 3) == 0)
	{
		E = 0.25 * log(ABS(D)) + 1.25 - 0.5 * log(3);
	}
	else
	{
		E = 0.5 * log(ABS(D)) + 2.5 - log(6);
	}
	
	return (exp(E) * sqrt(ABS(D))) / M_PI;
}

// This function is from Anton's C library
// https://github.com/amosunov/class_groups/blob/master/functions.c
int legendre_symbol(long a, int p)
{
	if (p == 2)
	{
		if (a & 1)
		{
			return ((a & 7) == 1 || (a & 7) == 7) ? 1 : -1;
		}
		else
		{
			return 0;
		}
	}

	int is_negative = 0;
	if (a < 0)
	{
		is_negative = ((p & 3) == 3);
		a = -a;
	}

	long temp;
	int t = 1;

	while (a != 0)
	{
		while ((a & 1) == 0)
		{
			a >>= 1;
			if ((p & 7) == 3 || (p & 7) == 5)
			{
				t = -t;
			}
		}

		if (a < p)
		{
			temp = a;
			a = p;
			p = temp;

			if ((a & 3) == 3 && (p & 3) == 3)
			{
				t = -t;
			}
		}

		a = ((a - p) >> 1);

		if ((p & 7) == 3 || (p & 7) == 5)
		{
			t = -t;
		}
	}

	if (is_negative)
	{
		t = -t;
	}

	return (p == 1) ? t : 0;
}

// Find the reduced ideal that is equivalent to prime ideal P above p 
void get_ideal_p(s64_qform_group_t* group, s64_qform_t* pform, const int p)
{
	if (s64_qform_is_primeform(group, pform, p)==1)
		s64_qform_reduce(group, pform);
}

// Find the order of ideal class [P] in Cl_k
int order_p(long D, int p, int h, int *factors_of_h)
{ 
	s64_qform_group_t group;
	s64_qform_group_init(&group);
	s64_qform_group_set_discriminant_s64(&group,-D);

	group_pow_t gp;
	group_pow_init(&gp,&group.desc.group);
 
	s64_qform_t pform, powform;
	s64_qform_init(&group,&pform);
	s64_qform_init(&group,&powform);
 
	get_ideal_p(&group, &pform, p);
  
	int total_factors=factors_of_h[0];
	uint32_t power=h;  int order_of_p=h;

	for (int j=1;j<=total_factors;j++)
	{ 
		power=order_of_p/factors_of_h[j];
		qform_pow_u32(&gp,&powform,&pform,power);
		while(s64_qform_is_id(&group,&powform)==1)
		{   
			order_of_p=order_of_p/factors_of_h[j];
			if(power%factors_of_h[j]==0)
			{     
				power=power/factors_of_h[j];
				qform_pow_u32(&gp,&powform,&pform,power);   
			}
			else
				break;
		}

	}
	group_pow_clear(&gp);
	s64_qform_group_clear(&group);

	return order_of_p; 
}

// Find r_0 with r_0 = D (mod r^exp) using Hensel's Lifting Lemma
void hensel(mpz_t r0, long D, int p, int exp)
{
	int Dmodp=(-D)%p; 
	if (Dmodp<0)
		Dmodp=p+Dmodp;
	const short* sqrtp;
	sqrtp=sqrtmodp[p];
	//int rr=sqrtp[Dmodp];
	mpz_t r;
	mpz_init(r);
	mpz_set_ui(r,sqrtp[Dmodp]); //r^2=-D (mod p)
 
	for(int j=1;j<exp;j++)
	{
		mpz_t modu,pp,temp; mpz_init(modu);mpz_init(pp);mpz_init(temp);
		mpz_ui_pow_ui(modu,p,j+1);
		mpz_t t,fder; mpz_init(t); mpz_init(fder); 
		mpz_set_ui(pp,p);
		mpz_mul_ui(fder,r,2);     
		mpz_mod(fder,fder,pp);  //fder=f'(r)=2*r%p
		mpz_invert(t,fder,pp);  //t=fder^-1
     
		//r=r-f(r)*f'(r)^-1, f(r)=r^2+D r=r-(r*r+D)*t;   
		mpz_set(temp,r);
		mpz_mul(r,r,r); mpz_add_ui(r,r,D); mpz_mul(r,r,t);mpz_sub(r,temp,r);
		mpz_mod(r,r,modu);

	}

	// we obtained an r s.t. r^2=-D (mod p^n)
	// make sure that r^2= -D (mod 4p^n), we need to check if r%2=-D%2, if not, r=p^exp-r
	mpz_t test; mpz_init(test);
	mpz_mod_ui(test,r,2);  
	if(mpz_cmpabs_ui(test,0)==0)
	{ 
		mpz_t mod;  mpz_init(mod);
		mpz_ui_pow_ui(mod,p,exp);
		mpz_sub(r,mod,r); 
	}

	mpz_set(r0,r);
}

// Find the fundamental unit of O_k,p using Cornacchia's Algorithm
// Save the result in unit
void get_unit(mpz_t *unit, long D, int p, int exp) 
{
	mpz_t m,r0,r1,r2,l; mpz_init(m);  mpz_init(r0);  mpz_init(r1);  mpz_init(r2);  mpz_init(l);
	mpz_ui_pow_ui(m,p,exp); //m=4*p^n  
	mpz_mul_ui(m,m,4);
	mpz_sqrt(l,m);
	if(mpz_cmpabs_ui(m,D)<0) // if D>4*p^n then y=0 and x=sqrt(4*p^n)
	{ 
		mpz_sqrt(unit[0],m); 
		mpz_set_ui(unit[1],0);
		return;
	}
	//find long r0^2 = -D(mod 4p^n)
	hensel(r0,D,p,exp);
 
	if(mpz_cmpabs(r0,l)<=0) 
	{
		mpz_set(unit[0],r0);
		//sol[1]=sqrt((m-r0*r0)/D)
		mpz_mul(r0,r0,r0); mpz_sub(r0,m,r0); mpz_div_ui(r0,r0,D); mpz_sqrt(r0,r0);
		mpz_set(unit[1],r0);
		return;
	} 
 
	mpz_t start; mpz_init(start); mpz_div_ui(start,m,2);
	mpz_mod(r1,start,r0); // start from (2*p^n, r0)
	if(mpz_cmpabs(r1,l)<=0)
	{
		mpz_set(unit[0],r1); 
		mpz_mul(r1,r1,r1); mpz_sub(r1,m,r1); mpz_div_ui(r1,r1,D); mpz_sqrt(r1,r1);
		mpz_set(unit[1],r1);
		return;
	}
 
	mpz_mod(r2,r0,r1);
	while(mpz_cmpabs(r2,l)>0)
	{
		mpz_t temp;  mpz_init_set(temp,r2);
		//mpz_set(temp2,r2);
		mpz_mod(r2,r1,r2);
		mpz_set(r0,r1);
		mpz_set(r1,temp);
	} 
	mpz_set(unit[0],r2); 
	mpz_mul(r2,r2,r2); mpz_sub(r2,m,r2); mpz_div_ui(r2,r2,D); mpz_sqrt(r2,r2);
	mpz_set(unit[1],r2);
}


int isprime(long D, int* primes)
{
	int i=1;
	while(primes[i]<=sqrt(D))
	{ 
		if(D%primes[i]==0) 
			return 0;
	i++;
	}
	return 1;
}

