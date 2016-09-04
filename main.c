/*
 * main.c
 *
 *  Created on: Sep 4, 2016
 *      Author: hongyanwang
 * This file is for the combination of Cohen-Lenstra and AAC conjectures
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <time.h>
#include <gmp.h>
#include <dirent.h>

#include "myfunctions.h"
#include "generator.h"


int main(int argc, char *argv[])
{
	const int idx=atoi(argv[1]);  // "index" we need to read

	clock_t begin,end;
	double times;
	begin=clock();
 
	// The following list denotes the primes we want to deal with
	const int ps=7;
	int myprime[7]={2,3,5,7,11,101,1009};

	long D_max=(idx+1)*pow(2,28);
	int h_max = h_upper_bound(-D_max);
	int temp=1, i, h, p;  long D;  

	int aa[4]={3,7,4,8};
	int mm[4]={8,8,16,16};
	int p_1[ps], p_2[ps]; memset(p_1,0,sizeof(p_1)); memset(p_2,0,sizeof(p_2));
	int fold, a, m, D0,order,h_kp;
	char *line;  size_t len = 0; 
	int *factors_of_h=NULL; factors_of_h = (int *) realloc(factors_of_h, 17* sizeof(int));//char *line = NULL; int* factors_of_h=NULL;
	long D_root = sqrt(D_max); 
	
	int* primes; 
	primes = (int *) malloc(((int) (1.25506 * D_root / log(D_root))) * sizeof(int));
	prime_sieve(D_root, primes);
	primes = (int *) realloc(primes, (2 + primes[0]) * sizeof(int));

	//This is the location of data files	
	char name[500]; char folder[500]="../run/data/theta";

	// Create a file to record the counters for Cohen-Lenstra heuristic
	sprintf(name,"./fac_%d_cohen",idx);
	creat(name,0744);
	FILE* wf=fopen(name,"w");
 
	// Create a file to record the counterexamples for AAC conjecture
	sprintf(name,"./fac_%d_aac",idx);
	creat(name,0744);
	FILE* wfaac=fopen(name,"w");

	// Read files from folders cl3mod8,cl7mod8,cl4mod16, cl8mod16
	for (fold=0;fold<4;fold++)
	{
		a=aa[fold]; m=mm[fold];
		D=idx*pow(2,28)+a;

		// unzip the file with index=idx to read data
		sprintf(name,"gunzip %s/cl%dmod%d/cl%dmod%d.%d.gz",folder,a,m,a,m,idx);
		system(name);
		sprintf(name,"%s/cl%dmod%d/cl%dmod%d.%d",folder,a,m,a,m,idx);
		FILE* readf=fopen(name,"r");  //open the file for reading
		printf("reading cl%dmod%d.%d...\n",a,m,idx);

		while(getline(&line, &len, readf) != -1)
		{
			sscanf(line,"%d %d",&D0,&h);
			D=D+D0*m; 
			for(int j=0;j<ps;j++)
			{
				p=myprime[j]; 
				if(legendre_symbol(-D,p)==1)
				{  
					p_1[j]++;
					if(h==1)////////////////////////////////
					{
						order=1;
						p_2[j]++;
					}
					else
					{  
						// factorize the class number h   
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
	
						order=order_p(D,p,h,factors_of_h);
						h_kp=h/order;
						while((h_kp&1)==0) 
 							h_kp=h_kp>>1; 
						if(h_kp==1) 
							p_2[j]++;
					}

						if(fold<2 && isprime(D,primes)==1)
							aac(D,p,order,wfaac);
				}
			}
		} // while getline loop ends

		fclose(readf);
		sprintf(name,"gzip %s/cl%dmod%d/cl%dmod%d.%d",folder,a,m,a,m,idx);
		system(name);

	}// end of fold loop

	fclose(wfaac);/////////////////////////////////////////////////////////////////

	// Write output 
	for (int k=0;k<ps;k++)/////////////////////////////////////////////////////////
		fprintf(wf,"%d	%d\n",p_1[k],p_2[k]);
	fclose(wf);

	free(primes);
	free(factors_of_h);

	end=clock();
	times=(double)(end-begin)/CLOCKS_PER_SEC;
	printf("Index %d: %lf seconds\n",idx,times); 

	return 0;
}
