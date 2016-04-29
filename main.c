#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "p2random.h"
#include "tree.h"


#include <xmmintrin.h>
#include <emmintrin.h> 
#include <pmmintrin.h> 
#include <tmmintrin.h> 
#include <smmintrin.h> 
#include <nmmintrin.h> 
#include <ammintrin.h> 
#include <x86intrin.h>

void print128(__m128i x)
{
        int32_t *pt = (int32_t*) &x;
        printf("%d %d %d %d\n", x[0],x[1],x[2],x[3]);

}
//Harcoded probing function
int harcoded_probing(Tree* tree, int32_t* probe) {
        printf("Reached the harcoded_probing function \n");
        
}

//general probing function
int probe_gen_search(Tree* tree, size_t *fanout, int32_t probe) 
{
        // int32_t current_level, int32_t prev_result, int32_t probe, int32_t fanout) {
        printf("Reached the general probing function \n");

        int32_t current_level = 0, prev_result = 0,res;
        int32_t** treeLevels = tree->key_array;
        int32_t num_levels = tree->num_levels;
        __m128i key =  _mm_cvtsi32_si128(probe);
        
        
        key = _mm_shuffle_epi32(key, 0);
        
        printf("\nSearching for key\n");
        print128(key);

        for(current_level=0; current_level < num_levels ; current_level++) {

                
                //-----probe level with fan out of 5---------
                if(fanout[current_level]==5) {         
                        __m128i lvl = _mm_load_si128(/**/ (__m128i*)&treeLevels[current_level][prev_result << 2]);
                        __m128i cmp_1 = _mm_cmpgt_epi32(key,lvl);
                        __m128 new_cmp = _mm_castsi128_ps(cmp_1);
                        int mask = _mm_movemask_ps(new_cmp);   // ps: epi32
                        res = _bit_scan_forward(mask ^ 0x1FF);
                        res += (prev_result << 2) + prev_result;
                }
                //-----probe level with fan out of 9---------
                else if (fanout[current_level]==9) {
                        __m128i lvl_A = _mm_load_si128(( __m128i*)&treeLevels[current_level][ prev_result << 3]);
                        __m128i lvl_B = _mm_load_si128(( __m128i*)&treeLevels[current_level][(prev_result << 3) + 4]);
                        __m128i cmp_A = _mm_cmpgt_epi32(key, lvl_A);
                        __m128i cmp_B = _mm_cmpgt_epi32(key, lvl_B);
                        __m128i cmp = _mm_packs_epi32(cmp_A, cmp_B);
                        cmp = _mm_packs_epi16(cmp, _mm_setzero_si128());
                        int mask = _mm_movemask_epi8(cmp);
                        res = _bit_scan_forward(mask ^ 0x1FFFF);
                        res += (prev_result << 3) + prev_result;
                }
                //-----probe level with fan out of 17---------
                else if(fanout[current_level]==17) {
                        __m128i del_ABCD = _mm_load_si128(( __m128i *) &treeLevels[current_level][prev_result << 4]);
                        __m128i del_EFGH = _mm_load_si128(( __m128i *) &treeLevels[current_level][(prev_result << 4) + 4]);
                        __m128i del_IJKL = _mm_load_si128(( __m128i *) &treeLevels[current_level][(prev_result << 4) + 8]);
                        __m128i del_MNOP = _mm_load_si128(( __m128i *) &treeLevels[current_level][(prev_result << 4) + 12]);
                        
                        // compare with 16 delimiters stored in 4 registers
                        __m128i cmp_ABCD = _mm_cmpgt_epi32(key, del_ABCD);
                        __m128i cmp_EFGH = _mm_cmpgt_epi32(key, del_EFGH);
                        __m128i cmp_IJKL = _mm_cmpgt_epi32(key, del_IJKL);
                        __m128i cmp_MNOP = _mm_cmpgt_epi32(key, del_MNOP);

                        // pack results to 16-bytes in a single SIMD register
                        __m128i cmp_A_to_H = _mm_packs_epi32(cmp_ABCD, cmp_EFGH);
                        __m128i cmp_I_to_P = _mm_packs_epi32(cmp_IJKL, cmp_MNOP);
                        __m128i cmp_A_to_P = _mm_packs_epi16(cmp_A_to_H, cmp_I_to_P);

                        // extract the mask the LSB
                        int mask = _mm_movemask_epi8(cmp_A_to_P);
                        res = _bit_scan_forward(mask ^ 0x1FFFFFFFF); 
                        res += (prev_result << 4) + prev_result;
                }
                prev_result = res;
        }
        return res;
}


//HARD-CODED PROBE IMPLEMENTATION


int main(int argc, char* argv[]) {
        // parsing arguments
        assert(argc > 3);
        size_t num_keys = strtoull(argv[1], NULL, 0);
        size_t num_probes = strtoull(argv[2], NULL, 0);
        size_t num_levels = (size_t) argc - 3;
        size_t* fanout = malloc(sizeof(size_t) * num_levels);
        assert(fanout != NULL);
        for (size_t i = 0; i < num_levels; ++i) {
                fanout[i] = strtoull(argv[i + 3], NULL, 0);
                assert(fanout[i] >= 2 && fanout[i] <= 17);
        }

        // building the tree index
        rand32_t* gen = rand32_init((uint32_t) time(NULL));
        assert(gen != NULL);
        int32_t* delimiter = generate_sorted_unique(num_keys, gen);
        assert(delimiter != NULL);
        Tree* tree = build_index(num_levels, fanout, num_keys, delimiter);
        free(delimiter);
        //free(fanout);
        if (tree == NULL) {
                free(gen);
                exit(EXIT_FAILURE);
        }

        // generate probes
        int32_t* probe = generate(num_probes, gen);
        assert(probe != NULL);
        free(gen);
        uint32_t* result = malloc(sizeof(uint32_t) * num_probes);
        assert(result != NULL);

        if (num_levels == 3 && fanout[0]==9 && fanout[1] == 5 && fanout[2]== 9) {
            // The hardcoded case
            int32_t *r;
            for (size_t i = 0; i < num_probes; i+=4) {
                    r = harcoded_probing(tree, fanout, probe[i]);
                    result[i]   = r[0];
                    result[i+1] = r[1];
                    result[i+2] = r[2];
                    result[i+3] = r[3];
            }
        } else {
            // perform index probing (Phase 2) - The general case
            for (size_t i = 0; i < num_probes; ++i) {
                    result[i] = probe_gen_search(tree, fanout, probe[i]);
            }
        }
        // output results
        for (size_t i = 0; i < num_probes; ++i) {
                fprintf(stdout, "%d %u\n", probe[i], result[i]);
        }

        // cleanup and exit
        free(result);
        free(probe);
        cleanup_index(tree);
        return EXIT_SUCCESS;
}
