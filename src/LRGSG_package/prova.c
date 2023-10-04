#include <stdio.h>
#include <inttypes.h>
#include "prova.h"
#include "SFMT/SFMT.h"
/*////////////////////////////////////////////////////////////////////// RNGS */
//
/** generate and return 64-bit pseudorandom number. init_gen_rand or
 * init_by_array must be called before this function.
 * @return (uint64_t) 64-bit pseudorandom number
 */
uint64_t SFMTrng_u64(void)
{
    return sfmt_genrand_uint64(&sfmt);
}
/** generate and return a random number on [0,1) with 53-bit resolution
 * init_gen_rand or init_by_array must be called before this function.
 * @return (double) number on [0,1) with 53-bit resolution
 */
double SFMTrng_dbl(void)
{
    return sfmt_genrand_res53(&sfmt);
}
/*/////////////////////////////////////////////////////////////////// SEEDERS */
//
/** set the seed of sfmt rng by array.
 */
extern void __set_seed_SFMT(void)
{
    seed_rand = (uint32_t[LENSRND]){SEED, SIID, CEED, CIID};
    sfmt_init_by_array(&sfmt, seed_rand, LENSRND);
}
/** print N_PRIGNG random numbers, integers and floating point types on stdout
 */
extern void __check_RNG(void)
{
    for (int i = 0; i < N_PRIRNG; i++)
        printf(STR_CHECK_RNG, RNG_dbl(), RNG_u64());
}