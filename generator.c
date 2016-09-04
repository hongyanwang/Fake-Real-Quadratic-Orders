/*
 * generator.h
 *
 *  Created on: Sep 4, 2016
 *      Author: hongyanwang
 */

#include "generator.h"
#include "myfunctions.h"

// This is from Max's C library
// https://github.com/maxwellsayles/libqform/blob/master/s64_qform.c
inline int64_t s64_qform_c(const s64_qform_group_t* group, const int32_t a, const int32_t b) 
{
	return (((int64_t)b*(int64_t)b - group->D) >> 2) / a;
}

// This is from Max's C library
// https://github.com/maxwellsayles/libqform/blob/master/s64_qform.c
static inline int32_t avg_s32(const int32_t a, const int32_t b) 
{
	return ((int64_t)a + (int64_t)b) >> 1;
}

// This is from Max's C library
// https://github.com/maxwellsayles/libqform/blob/master/s64_qform.c
static inline uint32_t half_rshift_u32(const uint32_t a) 
{
  // Should not return 0 since this causes weirdness with partial GCD.
	#if defined(__i386) || defined(__x86_64)
	uint32_t r;
	asm("xorl %%ecx, %%ecx\n\t"
	"bsrl %0, %%ecx\n\t"
	"incl %%ecx\n\t"  // round up
	"shrl $1, %%ecx\n\t"
	"shrl %%cl, %0"
	: "=r"(r)
	: "0"(a)
	: "cc", "ecx");
	#else
	uint32_t r = a;
	r = a;
	r >>= ((msb_u32(r) + 1) >> 1;
	#endif
	return r | (!r);  // r == 0 ? 1 : r;
}

// Find the multiplication of two relative generators 'f' and 'g'
// Save the relative generator in 'r'
void rela_gen_mpz_mul(rela_gen_mpz *r, rela_gen_mpz *f, rela_gen_mpz *g, mpz_t D)
{
	mpz_t a1,b1,c1,a2,b2,c2,t;  
	mpz_init(a1); mpz_init(b1); mpz_init(c1); mpz_init(a2); mpz_init(b2); mpz_init(c2); mpz_init(t);
	mpz_set(a1,f->a); mpz_set(b1,f->b); mpz_set(c1,f->c); 
	mpz_set(a2,g->a); mpz_set(b2,g->b); mpz_set(c2,g->c); 
	mpz_mul(r->a,a1,a2); mpz_mul(t,b1,b2); mpz_mul(t,t,D); //r->a=a1*a2+b1*b2*D 
	mpz_add(r->a,r->a,t);
	mpz_mul(r->b,a1,b2); mpz_mul(t,a2,b1); //r->b=a1*b2+a2*b1
	mpz_add(r->b,r->b,t);
	mpz_mul(r->c,c1,c2); //r->c=c1*c2
	mpz_gcd(t,r->a,r->b); mpz_gcd(t,t,r->c);
	mpz_div(r->a,r->a,t); mpz_div(r->b,r->b,t); mpz_div(r->c,r->c,t);
	mpz_clear(a1); mpz_clear(b1); mpz_clear(c1); mpz_clear(a2); mpz_clear(b2); mpz_clear(c2); mpz_clear(t);
} 

// Find the norm of a relative generator 'gen'
// Save the norm in 'norm'
void gen_norm(norm_mpz *norm, rela_gen_mpz *gen, mpz_t D)
{
	mpz_t temp; mpz_init(temp);
	mpz_mul(norm->a,gen->a,gen->a); //norm->a=a^2-b^2*D, norm->b=c^2
	mpz_mul(temp,gen->b,gen->b); mpz_mul(temp,temp,D); 
	mpz_sub(norm->a,norm->a,temp);
	mpz_mul(norm->b,gen->c,gen->c);
	mpz_gcd(temp,norm->a,norm->b); 
	mpz_div(norm->a,norm->a,temp); mpz_div(norm->b,norm->b,temp); 
	mpz_clear(temp);
}

// Find the multiplication of two norms 'm' and 'n'
// Save the result in 'r'
void norm_mul(norm_mpz *r, norm_mpz *m, norm_mpz *n)
{
	mpz_t temp; mpz_init(temp);
	mpz_mul(r->a,m->a,n->a); // r=a1*a2/b1*b2
	mpz_mul(r->b,m->b,n->b);
	mpz_gcd(temp,r->a,r->b);
	mpz_div(r->a,r->a,temp); mpz_div(r->b,r->b,temp);
	mpz_clear(temp);
}

// Generator mod D
void rela_gen_mpz_modD(rela_gen_mpz *gen, mpz_t D)
{
	mpz_mod(gen->a, gen->a, D);
	mpz_mod(gen->b, gen->b, D);
	mpz_mod(gen->c, gen->c, D);
	mpz_t temp; mpz_init(temp); mpz_gcd(temp,gen->a,gen->b); mpz_gcd(temp,temp,gen->c); 
	mpz_div(gen->a, gen->a, temp);  mpz_div(gen->b, gen->b, temp); mpz_div(gen->c, gen->c, temp);
	mpz_clear(temp);
}

// Find the multiplication of two binary quadratic forms R=A*B
// Reduce R and save the relative generator in 'gen'
void mpz_qform_compose_gen(mpz_qform_group_t* group, mpz_qform_t* R, const mpz_qform_t* A, const mpz_qform_t* B, rela_gen_mpz *gen) 
{
	register mpz_qform_compose_t* comp = &group->compose;
  
	if (mpz_cmp(A->a, B->a) > 0) 
	{
		mpz_set(comp->a1, A->a);
	mpz_set(comp->b1, A->b);
	mpz_set(comp->a2, B->a);
	mpz_set(comp->b2, B->b);
	mpz_set(comp->c2, B->c);
	} 
	else 
	{
		mpz_set(comp->a1, B->a);
		mpz_set(comp->b1, B->b);
		mpz_set(comp->a2, A->a);
		mpz_set(comp->b2, A->b);
		mpz_set(comp->c2, A->c);
	}
  
	// s = (b1 + b2)/2, m = (b1 - b2)/2
	mpz_add(comp->ss, comp->b1, comp->b2);
	mpz_fdiv_q_2exp(comp->ss, comp->ss, 1);
  
	mpz_sub(comp->m, comp->b1, comp->b2);
	mpz_fdiv_q_2exp(comp->m, comp->m, 1);
  
	// solve SP = v1 a2 + u1 a1 (only need v1)
	mpz_gcdext(comp->SP, comp->v1, 0, comp->a2, comp->a1);
  
	// K = v1 (b1 - b2) / 2 (mod L)
	mpz_mul(comp->K, comp->m, comp->v1);
	mpz_fdiv_r(comp->K, comp->K, comp->a1);
  
	if (mpz_cmp_ui(comp->SP, 1) != 0) 
	{
		mpz_gcdext(comp->S, comp->u2, comp->v2, comp->SP, comp->ss);
    
		// K = u2 K - v2 c2 (mod L)
		mpz_mul(comp->K, comp->K, comp->u2);
		mpz_mul(comp->temp, comp->v2, comp->c2);
		mpz_sub(comp->K, comp->K, comp->temp);
    
		if (mpz_cmp_ui(comp->S, 1) != 0) 
		{
			mpz_divexact(comp->a1, comp->a1, comp->S);
			mpz_divexact(comp->a2, comp->a2, comp->S);
			mpz_mul(comp->c2, comp->c2, comp->S);
		}
    
		mpz_fdiv_r(comp->K, comp->K, comp->a1);
	}
  
	// T = NK
	mpz_mul(comp->T, comp->a2, comp->K);
    
	// C.a = A.a B.a / d^2 = NL
	mpz_mul(R->a, comp->a2, comp->a1);
    
	// C.b = b2 + 2 a2 K = b2 + 2 T
	mpz_mul_2exp(R->b, comp->T, 1);
	mpz_add(R->b, R->b, comp->b2);
    
	// C.c = (S c2 + K (b2 + T)) / L;
	mpz_add(R->c, comp->b2, comp->T);
	mpz_mul(R->c, R->c, comp->K);
	mpz_add(R->c, R->c, comp->c2);
	mpz_divexact(R->c, R->c, comp->a1);

	red_gen_mpz(group,R,comp->S,gen);
}

// Find the square of a binary quadratic forms R=A*A
// Reduce R and save the relative generator in 'gen'
void mpz_qform_square_gen(mpz_qform_group_t* group, mpz_qform_t* R, const mpz_qform_t* A, rela_gen_mpz *gen) 
{
	register mpz_qform_compose_t* comp = &group->compose;
  
	mpz_set(comp->a1, A->a);
	mpz_set(comp->b1, A->b);
	mpz_set(comp->c1, A->c);
  
	// solve S = v1 b1 + u1 a1 (only need v1)
	mpz_gcdext(comp->S, comp->v1, 0, comp->b1, comp->a1);
	// K = -v1 c1 (mod L)
	mpz_mul(comp->K, comp->v1, comp->c1);
	mpz_neg(comp->K, comp->K);
  
	if (mpz_cmp_ui(comp->S, 1) != 0) 
	{
		mpz_divexact(comp->a1, comp->a1, comp->S);
		mpz_mul(comp->c1, comp->c1, comp->S);
	}
  
	mpz_fdiv_r(comp->K, comp->K, comp->a1);
  
	// T = NK
	mpz_mul(comp->T, comp->a1, comp->K);
    
	// C.a = A.a^2 / S^2 = N^2
	mpz_mul(R->a, comp->a1, comp->a1);
    
	// C.b = b1 + 2 a1 K = b1 + 2 T
	mpz_mul_2exp(R->b, comp->T, 1);
	mpz_add(R->b, R->b, comp->b1);
    
	// C.c = (S c1 + K (b1 + T)) / L;
	mpz_add(R->c, comp->b1, comp->T);
	mpz_mul(R->c, R->c, comp->K);
	mpz_add(R->c, R->c, comp->c1);
	mpz_divexact(R->c, R->c, comp->a1);
    
	red_gen_mpz(group,R,comp->S,gen);
}

// Find the relative generator in the reduced process
// Save the result in 'gen'
void red_gen_mpz(mpz_qform_group_t* group, mpz_qform_t* C, mpz_t s, rela_gen_mpz *gen)
{
	mpz_t na, nb, q, a2, OA, AA, NA, beta_x, temp;  //z=C->a OA=0 AA=1
	mpz_init(na); mpz_init(nb); mpz_init(q); mpz_init(a2); mpz_init(NA); mpz_init(beta_x); mpz_init(temp);
	mpz_init_set_si(AA,1); mpz_init_set_si(OA,0);
	while(mpz_cmp(C->a, C->c)>0)//C->a > C->c
	{ 
		mpz_set(na,C->c); //na=C->c;
		mpz_add(a2,na,na); //a2=na<<1;
		mpz_set(temp,C->b); //temp=-C->b;
		mpz_neg(temp,temp);
		mpz_divmod(q,nb,temp,a2); // divrem_s64(&q,&nb,temp,a2);
		if(mpz_cmp_ui(nb,0)<0) // if(nb<0), that indicates a2<0
		{ 
			mpz_sub(nb,nb,a2); //nb=nb-a2
			mpz_add_ui(q,q,1); //q++
		}
		if(mpz_cmp(nb,na)>0) //(nb>na)
		{
			mpz_sub(nb,nb,a2); //nb=nb-a2
			mpz_add_ui(q,q,1); //q++
		}
		
		mpz_sub(temp,nb,C->b); //temp=(nb-C->b)>>1;
		
		mpz_div_ui(temp,temp,2);
		mpz_mul(temp,temp,q); //temp=temp*q
		mpz_sub(C->c,C->a,temp); //C->c=C->a-temp;
		mpz_set(C->b,nb); //C->b=nb;
		mpz_set(C->a,na); //C->a=na;
		
		mpz_mul(NA,q,AA);//NA=q*AA+OA;
		mpz_add(NA,NA,OA);
		mpz_set(OA,AA); //OA=AA;
		mpz_set(AA,NA); //AA=-NA;
		mpz_neg(AA,AA);
	}
	if(mpz_cmp(C->a,C->c)==0 && mpz_cmp_ui(C->b,0)<0) //(C->a == C->c)&&(C->b < 0))
	{
		mpz_neg(C->b, C->b); //C->b = -C->b;
	}
	mpz_mul(beta_x,AA,C->a);//beta_x=(AA*C->a)<<1;
	mpz_mul_si(beta_x,beta_x,2);
	mpz_mul(temp,OA,C->b); //temp=OA*C->b;
	mpz_sub(beta_x,beta_x,temp);//beta_x=beta_x-temp;
	mpz_mul(gen->a,beta_x,s);//gen->a=s*beta_x; 
	mpz_mul(gen->b,OA,s);//gen->b=-s*OA; 
	mpz_neg(gen->b,gen->b);
	mpz_set(gen->c,C->a); //gen->c=2*C->a;
	mpz_mul_ui(gen->c, gen->c, 2);
	mpz_clear(na); mpz_clear(nb); mpz_clear(q); mpz_clear(a2); mpz_clear(NA); mpz_clear(beta_x); mpz_clear(temp);
	mpz_clear(AA); mpz_clear(OA);
}

// Check if the (D,p) pair is a counterexample of AAC conjecture
void AAC_conj(mpz_qform_group_t *group, mpz_qform_t *pform, long p, int order, FILE *wf)
{
	mpz_t DD; mpz_init_set(DD,group->D);
	mpz_t test;  mpz_init(test); 
	norm_mpz unit_norm; norm_mpz_init(&unit_norm); 
	rela_gen_mpz unitD; rela_gen_mpz_init(&unitD);
	if(order==1)
	{
		mpz_t s; mpz_init(s); mpz_set_si(s,1); red_gen_mpz(group,pform,s,&unitD);
		gen_norm(&unit_norm, &unitD, DD);  rela_gen_mpz_modD(&unitD, DD); 
		mpz_clear(s);
	}
	else
	{
		int order_b[(int)(log(order)/log(2)+2)]; // order_b records the binary sequence of the order of p, order_b[0] is the number of bits of the order
		memset(order_b,0,sizeof(order_b));
		dec_bin(order, order_b);

		int gen_no=order_b[0]-1; // number of generators
		rela_gen_mpz gens[gen_no]; // use gens to record all the generators  
		for(int i=0;i<gen_no;i++)
		{	rela_gen_mpz_init(&gens[i]); 	}
		rela_gen_mpz gen; rela_gen_mpz_init(&gen);
		mpz_qform_t redform; mpz_qform_init(group,&redform); mpz_qform_set(group, &redform, pform); //first, r=p
		for(int i=gen_no;i>=1;i--)
		{	
			mpz_qform_square_gen(group,&redform,&redform, &gen);
			rela_gen_mpz_set(&gens[gen_no-i], &gen);
			if(order_b[i]==1)
			{
				mpz_qform_compose_gen(group,&redform,&redform,pform,&gen);
				rela_gen_mpz_mul(&gens[gen_no-i], &gens[gen_no-i], &gen, DD);
			}
		} // we just found all the generators
		
		// Multiply the norms of all the generators to check the correctness of the unit
		// Multiply the generators modulo D to check if D|b
		norm_mpz temp; norm_mpz_init(&temp); gen_norm(&unit_norm, &gens[0], DD);
		rela_gen_mpz_set(&unitD, &gens[0]); rela_gen_mpz_modD(&unitD, DD);
		for(int k=1;k<gen_no;k++)
		{
			// unit_norm gives us the norm of the unit.
			norm_mul(&unit_norm, &unit_norm,&unit_norm);
			gen_norm(&temp, &gens[k], DD);
			norm_mul(&unit_norm, &unit_norm, &temp);
		
			//here we multiply the generators together
			rela_gen_mpz_mul(&unitD, &unitD, &unitD, DD); 
			rela_gen_mpz_modD(&unitD, DD);
			rela_gen_mpz_mul(&unitD, &unitD, &gens[k], DD);
			rela_gen_mpz_modD(&unitD, DD);
		}
		for(int i=0;i<gen_no;i++)
		{	rela_gen_mpz_clear(&gens[i]); 	}
		rela_gen_mpz_clear(&gen);
		norm_mpz_clear(&temp);	
		mpz_qform_clear(group, &redform);		
	}
	
	mpz_ui_pow_ui(test,p,order); mpz_mul(test,test,unit_norm.b); //To check if a/b=p^n, let test=p^n*b then check if a=test
	if(mpz_cmpabs(unit_norm.a, test)!=0)
	{
		gmp_printf("wrong unit for D=%Zd,p=%ld\n",DD,p);
	}
	else
	{
		if(mpz_cmpabs_ui(DD,4)>0 && mpz_cmpabs_ui(unitD.b,0)==0) //D<-4
		{  
			gmp_printf("Counterexample found for D=%Zd, p=%ld\n",DD,p);
			long D=mpz_get_si(DD); char buf[100]; sprintf(buf,"%ld	%ld\n",D,p);
			fputs(buf,wf);
		}
		if(mpz_cmpabs_ui(DD,3)==0) //(D==-3) // if D=3, there are three different units, we need to check one by one
		{ 
			if(mpz_cmpabs_ui(unitD.b,0)==0)
			{ 
				gmp_printf("Counterexample found for D=-3, p=%ld\n",p);
				long D=-3; char buf[100]; sprintf(buf,"%ld	%ld\n",D,p);
				fputs(buf,wf);
			}
			else    // unit2[0]=(unit[0]-3*unit[1])/2; unit2[1]=(unit[0]+unit[1])/2
			{ 
				rela_gen_mpz unit2; rela_gen_mpz_init(&unit2); 
				mpz_set_si(unit2.a, 1); mpz_set_si(unit2.b, 1); mpz_set_si(unit2.c, 2); 
				rela_gen_mpz_mul(&unitD, &unitD, &unit2, DD); rela_gen_mpz_modD(&unitD,DD);
				if(mpz_cmpabs_ui(unitD.b,0)==0)
				{  
					gmp_printf("Counterexample found for D=-3, p=%ld\n",p);
					long D=-3; char buf[100]; sprintf(buf,"%ld	%ld\n",D,p);
					fputs(buf,wf);
				}
				else
				{ 
					mpz_set_si(unit2.a, 1); mpz_set_si(unit2.b, 1); mpz_set_si(unit2.c, 2); 
					rela_gen_mpz_mul(&unitD, &unitD, &unit2, DD); rela_gen_mpz_modD(&unitD,DD);
					if(mpz_cmpabs_ui(unitD.b,0)==0)
					{  
						gmp_printf("Counterexample found for D=-3, p=%ld\n",p);
						long D=-3; char buf[100]; sprintf(buf,"%ld	%ld\n",D,p);
						fputs(buf,wf);
					}
 				}
			}   
		}
	}
	norm_mpz_clear(&unit_norm);
	rela_gen_mpz_clear(&unitD);
	mpz_clear(DD); mpz_clear(test);
}

// Check the AAC conjecture for D,p
// Call function 'AAC_cnj'
void aac(long D, long p, int order, FILE* wf)
{
	mpz_t DD; mpz_init(DD); mpz_set_si(DD,-D); 
	mpz_qform_group_t group;
	mpz_qform_group_init(&group);
	mpz_qform_group_set_discriminant(&group,DD);

	mpz_qform_t pform;
	mpz_qform_init(&group,&pform); 
	mpz_qform_is_primeform_why(&group, &pform, p);
	AAC_conj(&group, &pform, p, order, wf);
	mpz_qform_clear(&group, &pform);
	mpz_qform_group_clear(&group);
	mpz_clear(DD);
}


// Find the fundamental unit of a fake real quadratic order
void AAC_unit(mpz_qform_group_t *group, mpz_qform_t *pform, long p, int order, FILE *wf)
{
	mpz_t DD; mpz_init_set(DD,group->D);
	rela_gen_mpz unitD; rela_gen_mpz_init(&unitD);
	if(order==1)
	{
		mpz_t s; mpz_init(s); mpz_set_si(s,1); 
		//unitD is  the fundamental unit		
		red_gen_mpz(group,pform,s,&unitD);
		mpz_clear(s);
	}
	else
	{
		int order_b[(int)(log(order)/log(2)+2)]; // order_b records the binary sequence of the order of p, order_b[0] is the number of bits of the order
		memset(order_b,0,sizeof(order_b));
		dec_bin(order, order_b);

		int gen_no=order_b[0]-1; // number of generators
		rela_gen_mpz gens[gen_no]; // use 'gens' to record all the generators  
		for(int i=0;i<gen_no;i++)
			rela_gen_mpz_init(&gens[i]); 	
		rela_gen_mpz gen; rela_gen_mpz_init(&gen);
		mpz_qform_t redform; mpz_qform_init(group,&redform); mpz_qform_set(group, &redform, pform); //first, r=p
		for(int i=gen_no;i>=1;i--)
		{	
			mpz_qform_square_gen(group,&redform,&redform, &gen);
			rela_gen_mpz_set(&gens[gen_no-i], &gen);
			if(order_b[i]==1)
			{
				mpz_qform_compose_gen(group,&redform,&redform,pform,&gen);
				rela_gen_mpz_mul(&gens[gen_no-i], &gens[gen_no-i], &gen, DD);
			}
		} // we just found all the generators
		
		// Multiply the generators to find the fundamental unit
		rela_gen_mpz_set(&unitD, &gens[0]); 
		for(int k=1;k<gen_no;k++)
		{
			rela_gen_mpz_mul(&unitD, &unitD, &unitD, DD); 
			rela_gen_mpz_mul(&unitD, &unitD, &gens[k], DD);
		}
		for(int i=0;i<gen_no;i++)
			rela_gen_mpz_clear(&gens[i]); 	
		rela_gen_mpz_clear(&gen);
		mpz_qform_clear(group, &redform);		
	}

	// output the fundamental unit into the file
	if(mpz_cmpabs_ui(unitD.c,1)==0)
		gmp_fprintf(wf,"%Zd\n\n%Zd",unitD.a,unitD.b);
	else
		if(mpz_cmpabs_ui(unitD.c,2)==0)
			gmp_fprintf(wf,"%Zd/2\n\n%Zd/2",unitD.a,unitD.b);
		else
			gmp_fprintf(wf,"Error! Wrong unit!");

	rela_gen_mpz_clear(&unitD);
	mpz_clear(DD);
}

// Find the fundamental unit of a fake real quadratic order 
// Call function 'AAC_unit'
// Output the fundamental unit into file 'wf'
void unit(long D, long p, int order, FILE* wf)
{
	mpz_t DD; mpz_init(DD); mpz_set_si(DD,-D); 
	mpz_qform_group_t group;
	mpz_qform_group_init(&group);
	mpz_qform_group_set_discriminant(&group,DD);

	mpz_qform_t pform;
	mpz_qform_init(&group,&pform); 
	mpz_qform_is_primeform_why(&group, &pform, p);
	AAC_unit(&group, &pform, p, order, wf);
	mpz_qform_clear(&group, &pform);
	mpz_qform_group_clear(&group);
	mpz_clear(DD);
	
}

// Convert an integer into a binary sequence, bi[0] is the number of bits 11=[4,1,1,0,1] 11=1011
void dec_bin(long h,int* bi)
{
	int i=1;
	while(h!=1)
	{
		if((h&1)==1)
		bi[i]=1;
		h=h>>1;
		i++;
	}
	bi[i]=1; bi[0]=i;
}


// Find a^n (mod p)
long powmodp(long p, long a, long n)
{
	if(a>p)
		a=a%p;
	if(a<0)
		a=a+p;

	//temp=a;
	mpz_t temp; mpz_init(temp);	
	mpz_set_si(temp, a);

	// some special cases
	if(a==0)
		return 0;
	if(n==0 || a==1)
		return 1;

	if(n==1)
		return a;

	int order_n[(int)(log(n)/log(2)+2)];
	memset(order_n,0,sizeof(order_n));
	dec_bin(n, order_n);
	// exponential algorithm
	for(int i=order_n[0]-1;i>=1;i--)
	{	
		//temp=(temp*temp)%p;
		mpz_mul(temp,temp,temp);
		mpz_mod_ui(temp,temp,p);	
		if(order_n[i]==1)
		{
			//temp=(temp*a)%p;
			mpz_mul_ui(temp,temp,a);
			mpz_mod_ui(temp,temp,p);
		}
	}
	long t=mpz_get_ui(temp);
	mpz_clear(temp);
	return t;
}

// Find a with a^2=D (mod p)
long sqrt_modp(long p, long D)
{
	long Dmodp=D%p;

	if(Dmodp<0)
		Dmodp=p+Dmodp;
	long a;
	if(p==2)
	{
		if(Dmodp==1)
			return 1;
		else
			return 0;	
	}
	if(Dmodp==1)
		return 1;
	
	// if p=3(mod 4), a=D^{(p+1)/4} (mod p)
	if(p%4==3)
		//printf("Dmodp=%ld,p=%ld,pow=%ld\n",Dmodp,p,(p+1)/4);	
		return powmodp(p,Dmodp,(p+1)/4);
	
	// if p=5(mod 8), v=(2D)^{(p-5)/8} (mod p), i=2Dv^2 (mod p), a=Dv(i-1) (modp)
	if(p%8==5)
	{	// v=(2D)^{(p-5)/8} (mod p)
		long v=powmodp(p,2*Dmodp,(p-5)/8);
		// i=2Dv^2 (mod p)
		long i=powmodp(p,v,2);
		//i=(2*Dmodp*i)%p;
		mpz_t t1; mpz_init(t1);
		mpz_set_si(t1,Dmodp);
		mpz_mul_ui(t1,t1,2);
		mpz_mul_ui(t1,t1,i);
		mpz_mod_ui(t1,t1,p); // i=t1
		i=mpz_get_ui(t1);
		// a=Dv(i-1) (modp) -> a=(Dmodp*v)%p; a=(a*(i-1))%p;
		mpz_set_si(t1,Dmodp);
		mpz_mul_ui(t1,t1,v);
		mpz_mul_ui(t1,t1,i-1);
		mpz_mod_ui(t1,t1,p);
		a=mpz_get_ui(t1);
		mpz_clear(t1);
		if(a<=p/2)
			return a;
		else
			return p-a;
	}
	
	// if p=1(mod 8), use Shanks algorithm
	if(p%8==1)
	{
		// write p-1=Q*2^s
		long Q=p-1, s=0, z, i;
		while((Q&1)==0) // Q%2=0
		{
			Q=Q>>1;
			s++;
		}
		for(z=2;z<p;z++)
		{
			if(legendre_symbol(z,p)==-1)
				break;
		}
		// c=z^Q(mod p), R=D^{(Q+1)/2} (mod p), t=D^Q (mod p), M=s
		long c=powmodp(p,z,Q); 
		long R=powmodp(p,Dmodp,(Q+1)/2);
		long t=powmodp(p,Dmodp,Q);
		long M=s;
		long b;
		mpz_t t1; mpz_init(t1);
		while(t!=1)
		{	
			// find smallest 0<i<M s.t. t^{2^i}=1 (mod p)
			for(i=1;i<M;i++)
			{
				if(powmodp(p,t,pow(2,i))==1)
					break;
			}
			// b=c^{2^(M-i-1)} (mod p), R=R*b (mod p), t=t*b^2 (mod p), c=b^2 (mod p), M=i
			b=powmodp(p,c,pow(2,M-i-1)); 
			//R=(R*b)%p;			
			mpz_set_si(t1,R);
			mpz_mul_ui(t1,t1,b);
			mpz_mod_ui(t1,t1,p);
			R=mpz_get_ui(t1);
			//t=(t*b*b)%p; 
			mpz_set_si(t1,t);
			mpz_mul_ui(t1,t1,b);
			mpz_mul_ui(t1,t1,b);
			mpz_mod_ui(t1,t1,p);
			t=mpz_get_ui(t1);

			c=powmodp(p,b,2); 
			M=i; 

		}
		mpz_clear(t1);
		if(R<=p/2)
			return R;
		else
			return p-R;
	}
	return 0;
}



// Modificatoin of Max's 's64_qform_is_primeform'
// We do not check table for the square root of D modulo p, instead, we use our function 'sqrt_modp'
//int s64_qform_is_primeform(s64_qform_group_t* group, s64_qform_t* form, const int p) 
//{
	// a = p
//	form->a = p;

//	form->b = sqrt_modp(p,group->D);   
  
	// We know that p | b^2-D for +/- b.
	// Compute c = (b^2 - d) / (4a) if possible.
//	if (p == 2) 
//	{
		// special case if p == 2
//		form->c = (int64_t)form->b * (int64_t)form->b;
//		form->c -= group->D;
//		if ((form->c & 7) == 0) 
//		{
			// 4a | b^2-D, a=2, so we divide by 8
//			form->c >>= 3;
//			return 1;
//		}
//		if (form->b == 1)
//			return 0;
    
		// try -b = 2
//		form->b = 2;
//		form->c = 4 - group->D;
//		if ((form->c & 7) != 0)
//			return 0;
    
		// 4a | b^2-D, a=2, so we divide by 8
//		form->c >>= 3;
//		return 1;
//	}
	// p != 2
  
//	form->c = (int64_t)form->b * (int64_t)form->b;
//	form->c -= group->D;
  
//	if ((form->c & 3) == 0) 
//	{
		// 4a | b^2-D
		// divide by 4
//		form->c >>= 2;
		// divide by a
//		form->c /= p;
//		return 1;
//	}
  
	// try -b
//	form->b = p-form->b;
  
	// make sure that 4 | b^2-D
//	form->c = (int64_t)form->b * (int64_t)form->b;
//	form->c -= group->D;
//	if ((form->c & 3) != 0)
//		return 0;
  
	// divide by 4
//	form->c >>= 2;
  
	// divide by a
//	form->c /= p;
//		return 1;
//}

// Modificatoin of Max's 'mpz_qform_is_primeform'
// We do not check table for the square root of D modulo p, instead, we use our function 'sqrt_modp'
int mpz_qform_is_primeform_why(mpz_qform_group_t* group, mpz_qform_t* form, long p) 
{
	long D=mpz_get_si(group->D); //printf("D=%ld\n",D);
	long sqrtDmodp=sqrt_modp(p,D); //printf("sqrtDmodp=%ld\n",sqrtDmodp);
	// a = p
	mpz_set_si(form->a, p);
	// b = sqrt(D) (mod p)//////////////////////////////////////////////////////////
	mpz_set_si(form->b, sqrtDmodp);
  
	// We know that p | b^2-D for +/- b.
	// Compute c = (b^2 - d) / (4a) if possible.
	if (p == 2) 
	{
		// special case if p == 2
		mpz_mul(form->c, form->b, form->b);
		mpz_sub(form->c, form->c, group->D);

		if ((form->c->_mp_d[0] & 7) == 0) 
	      	{
			// 4a | b^2-D, a=2, so we divide by 8
			mpz_fdiv_q_2exp(form->c, form->c, 3);
        			return 1;
		}

		if (sqrtDmodp == 1) // form->b == 1
			return 0;
      
      		// try -b = 2
		mpz_set_si(form->b, 2);
		mpz_set_si(form->c, 4);
		mpz_sub(form->c, form->c, group->D);
		if ((form->c->_mp_d[0] & 7) != 0)
			return 0;
    
		// 4a | b^2-D, a=2, so we divide by 8
		mpz_fdiv_q_2exp(form->c, form->c, 3);

		return 1;
	} // end p=2

	// p != 2
	mpz_mul(form->c, form->b, form->b);
	mpz_sub(form->c, form->c, group->D);
 
	if ((form->c->_mp_d[0] & 3) == 0) 
	{
		// 4a | b^2-D
		// divide by 4
		mpz_fdiv_q_2exp(form->c, form->c, 2);
    
		// divide by a
		mpz_divexact_ui(form->c, form->c, p);
			return 1;
	}

	// try -b
	mpz_set_si(form->b, p-sqrtDmodp);
  
	// make sure that 4 | b^2-D
	mpz_mul(form->c, form->b, form->b);
	mpz_sub(form->c, form->c, group->D);
	if ((form->c->_mp_d[0] & 3) != 0)
		return 0;
  
  
	// divide by 4
	mpz_fdiv_q_2exp(form->c, form->c, 2);
  
	// divide by a
	mpz_divexact_ui(form->c, form->c, p);
  
	return 1;
}

// Find the order of the ideal class [P] in Cl_k using big integers
int order_mpz(long D, long p, int h, int *factors_of_h)
{ 
	mpz_t DD; mpz_init(DD);
	mpz_init_set_si(DD,-D);
	mpz_qform_group_t group;
	mpz_qform_group_init(&group);
	mpz_qform_group_set_discriminant(&group,DD);

	group_pow_t gp;
	group_pow_init(&gp,&group.desc.group);
 
	//s64_qform_t pform;
	mpz_qform_t pform, powform;
	mpz_qform_init(&group,&pform);
	mpz_qform_init(&group,&powform);

	mpz_qform_is_primeform_why(&group,&pform,p);
	mpz_qform_reduce(&group,&pform);

	//int *factors_of_h=h_factors[h];
	int total_factors=factors_of_h[0];
	uint32_t power=h;  int order_of_p=h;

	for (int j=1;j<=total_factors;j++)
	{ 
		power=order_of_p/factors_of_h[j];
		qform_pow_u32(&gp,&powform,&pform,power);
		while(mpz_qform_is_id(&group,&powform)==1)
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
	mpz_qform_group_clear(&group);
	mpz_clear(DD);

	return order_of_p; 
}
