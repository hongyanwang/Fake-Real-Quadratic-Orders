/*
 * aacD.c
 *
 *  Created on: Sep 4, 2016
 *      Author: hongyanwang
 * Find AAC counterexamples for selected D with (stage-1)*1000000< p <up=stage*1000000
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <sys/time.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <dirent.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>
#include "myfunctions.h"
#include "generator.h"

int main(int argc, char *argv[])
{
	const int stage=atoi(argv[1]);  
	long down=(long)(stage-1)*1000000; long up=(long)stage*1000000; // blocksize = 1000000 (1 million)
	printf("Primes %ld-%ld\n",down,up);
	
	clock_t begin,end;
	double times;
	begin=clock();

	// This is the list of discriminants we want to find counterexamples for
	long D_list[10]={1073741783,4294967291,8589934583,17179869143,34359738319,68719476731,137438953447,274877906899,549755813723,1099511627563};
	int h_list[10]={22991,34805,107849,136023,86757,98601,162855,114227,220527,176451};
	
	// The upper bound of these discriminants is D_max
	long D_max=pow(2,40), D, p;
	// The upper bound of class numbers for these discriminants is h_max
	int h_max = h_upper_bound(-D_max);

	int temp=1, h, order; 
	int *factors_of_h=NULL; factors_of_h = (int *) realloc(factors_of_h, 17* sizeof(int));
	char aac_file[100]; 

	// Use prime sieve to find all primes less than D_root, D_root>sqrt(up)
	int* primes; 
	long D_root=10000000;
	primes = (int *) malloc(((int) (1.25506 * D_root / log(D_root))) * sizeof(int));
	prime_sieve(D_root, primes);
	primes = (int *) realloc(primes, (2 + primes[0]) * sizeof(int));
	
	for(int count=0;count<10;count++)
	{  
		D=D_list[count];
		h=h_list[count];
		// Save the counterexamples in the file ./AAC/D/D_stage
		sprintf(aac_file,"./AAC/%ld/%ld_%d",D,D,stage);
		FILE* wf=fopen(aac_file,"w");

		for(p=down+1;;p++)
		{
			if(p>up)
				break;
			if(isprime(p,primes)==1)
			{	
				if(legendre_symbol(-D,p)==1)
				{  
					if(h==1)
						order=1;
					else
					{ 	
						temp=h;int count=1; factors_of_h[0]=0;  
						for(int s=1;primes[s]<=h_max;s++)
						{    
							if(temp%primes[s]==0)
							{ 
								factors_of_h[count]=primes[s]; count++;  
								factors_of_h[0]=factors_of_h[0]+1;
								while(temp%primes[s]==0)
									temp=temp/primes[s]; 
							}
							if(temp<=1)
								break;
						}
						order=order_mpz(D,p,h,factors_of_h);
					}		
					aac(D,p,order,wf);
				}
			}
		} // run out the primes
		fclose(wf); 
	} // count loop ends
	free(factors_of_h);

	free(primes);

	end=clock();
	times=(double)(end-begin)/CLOCKS_PER_SEC;
	printf("stage %d:%lf seconds\n",stage,times); 

	return 0;
}
