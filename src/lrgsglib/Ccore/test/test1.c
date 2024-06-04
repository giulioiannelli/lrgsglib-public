#include "sfmtrng.h"
#include <stdlib.h>
sfmt_t sfmt;
uint32_t *seed_rand;
uint32_t seed_array[4] = {3776627104, 436986138, 1717310142, 196696};

int main(int argc, char *argv[])
{
    __set_seed_SFMT();
    // sfmt_init_gen_rand(&sfmt, 4321);
    // sfmt_genrand_uint32(&sfmt);
    // Generate random numbers
    // for (int i = 0; i < 10; i++) {
    //     if (sfmt.idx % 2 != 0) {
    //         fprintf(stderr, "Error: sfmt.idx is not even. Value: %d\n", sfmt.idx);
    //         return 1;
    //     }
    //     uint64_t rand_num = sfmt_genrand_uint64(&sfmt);
    //     printf("Random number: %lu\n", rand_num);
    // }
    uint64_t *random_numbers = __gen_rand_u64_array(atoi(argv[1]));

    // Print the random integers (optional)
    for (size_t i = 0; i < atoi(argv[1]); i++) {
        printf("%lu\n", random_numbers[i]);
    }

    // Free allocated memory
    free(random_numbers);

    return 0;
}
