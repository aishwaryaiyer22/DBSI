#include <xmmintrin.h> //SSE
#include <emmintrin.h> //SSE2
#include <pmmintrin.h> //SSE3
#include <tmmintrin.h> //SSSE3
#include <smmintrin.h> //SSE4.1
#include <nmmintrin.h> //SSE4.2
#include <ammintrin.h> //SSE4A
#include <x86intrin.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "p2randomv2.h"


int levelProbe (/**/int32_t** treeLevels/**/, int32_t currentLevel, int32_t prev_result, int32_t probe, int32_t fanout) {
	
	__m128i key =  _mm_cvtsi32_si128(probe);
	key = _mm_shuffle_epi32(key, 0);
	//-----probe level with fan out of 5---------
	if(fanout==5) {		
		__m128i lvl = _mm_load_si128(/**/ (__m128i*)&treeLevels[currentLevel][prev_result << 2]);
		__m128i cmp_1 = _mm_cmpgt_epi32(lvl, key);
		__m128 new_cmp = _mm_castsi128_ps(cmp_1);
		int mask = _mm_movemask_ps(new_cmp);   // ps: epi32
		int res = _bit_scan_forward(mask ^ 0x1FF);
		res += (res << 2) + res;
		return res;
	}
	//-----probe level with fan out of 9---------
	else if (fanout==9) {
		__m128i lvl_A = _mm_load_si128(( __m128i*)&treeLevels[currentLevel][ prev_result << 3]);
		__m128i lvl_B = _mm_load_si128(( __m128i*)&treeLevels[currentLevel][(prev_result << 3) + 4]);
		__m128i cmp_A = _mm_cmpgt_epi32(lvl_A, key);
		__m128i cmp_B = _mm_cmpgt_epi32(lvl_B, key);
		__m128i cmp = _mm_packs_epi32(cmp_A, cmp_B);
		cmp = _mm_packs_epi16(cmp, _mm_setzero_si128());
		int mask = _mm_movemask_epi8(cmp);
		int res = _bit_scan_forward(mask ^ 0x1FFFF);
		res += (prev_result << 3) + prev_result;
		return res;
	}
	//-----probe level with fan out of 17---------
	else if(fanout==17) {
		
	}
}



int main(int argc, char *argv[]) {
	
	int32_t k,p,nol;

	printf("\n");
	
	if (argc<4){
		printf("You seem to be missing some arguments..\n");
		exit(0);
	}
	
	k = atoi(argv[1]);
	p = atoi(argv[2]);

	nol = argc - 3;
	
	printf("K %d P %d argc %d",k,p,argc);
	int32_t *treeLevels[nol];
	int32_t sizeOfLevel[nol], fanoutOfLevel[nol];
	//Current position to insert at each level
	int32_t cpOfLevel[nol], blockFillLevel[nol], totalTreeKeySpace=0;
	
	int32_t i=0;
	//Pointer for posix mem align
	void *mp;
	while (i<nol) {
		fanoutOfLevel[i] = atoi(argv[i+3]);
		//Check for valid fanout level
		if(fanoutOfLevel[i]<2 || fanoutOfLevel[i]>17) {
			printf("Invalid Fanout for level %d\n", i+1);
			exit(0);
		}

		sizeOfLevel[i] = ( i>0 ? (sizeOfLevel[i-1]+sizeOfLevel[i-1]/(fanoutOfLevel[i-1]-1))*(fanoutOfLevel[i]-1) : fanoutOfLevel[i]-1); 
		totalTreeKeySpace += sizeOfLevel[i];
		//printf("Level - %d -      Size - %d \n",i,sizeOfLevel[i]);
		
		// Allocate each level and point to it
		posix_memalign( &mp, 16, sizeof(int32_t)*sizeOfLevel[i]);
		treeLevels[i] = (int32_t *)mp;
		//Initialize all current insertion point to starting
		cpOfLevel[i] = 0;
		blockFillLevel[i] = 0;
		//Go to next level
		i++;
	}
	
	if(k>totalTreeKeySpace){
		printf("\n Number of keys are more than available space..\n");
		exit(0);
	}
	
	//Insert max int numbers to all space
	for (int32_t j = 0; j < nol; j++) {
		for(int32_t *k = treeLevels[j],l=0; l<sizeOfLevel[j]; k++,l++) {
			//This will automatically pad all the extraspaces with maxint
			*k = INT32_MAX;
		}
	}

	/*
	for (int32_t j = 0; j < nol; j++) {
		for(int32_t *k = treeLevels[j],l=0; l<sizeOfLevel[j]; k++,l++) {
			printf(" %u ",*k);
		}
		printf("\n");
	}
	*/

	//Keys generated and sorted from prandom v2
	rand32_t *gen = rand32_init(time(NULL));
	size_t ii;
	int32_t *a = generate_sorted_unique(k, gen);
	free(gen);

	//Start bulk insertion from bottom most
	int32_t insLevelPosition = nol - 1;
	int32_t *s;
	int t;
	//printf("-------%d-----\n",totalTreeKeySpace);
	//Insert the keys one by one
	for(ii=0;ii<k;ii++) {
		t = nol-1;
		while( blockFillLevel[t]>=(fanoutOfLevel[t]-1) ) {
			blockFillLevel[t] = 0;
			t--;
			if(t<0) {
				printf("\n Woah cant insert more.. Check your code..\n");
				exit(0);
			}
		}
		s = treeLevels[t];
		for(int x=0;x<cpOfLevel[t];x++, s++){}
		
		*s = a[ii];
		//*s = ii;
		//printf("\nInserting at %d value %d\n", *s, a[ii]);
		//Update the positions
		cpOfLevel[t]++;
		blockFillLevel[t]++;
	}
	
	//Print each level
	printf("\nThe Tree: \n");
	for (int32_t j = 0; j < nol; j++) {
		printf("Level %d \n",j);
		for(int32_t *k = treeLevels[j],l=0; l<sizeOfLevel[j]; k++,l++) {
			printf(" %d ",*k);
		}
		printf("\n");
	}
	
	//Check if the root has atleast one element , else return error
	if ( *treeLevels[0] == INT32_MAX) {
		printf("\nToo few keys... :( ..\n");
		exit(0);
	}
	
	//------------------Generate the probes----------
	rand32_t *gen1 = rand32_init(time(NULL));
	int32_t *probes = generate(p, gen1);
	free(gen1);
	int index = levelProbe(treeLevels,1,0,175,5);
	printf("Index found: %d",index);


return (0);

}





