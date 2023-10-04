#include<stdio.h>
#include<inttypes.h>
#include "LRGSG_utils.h"
#include "LRGSG_customs.h"
#include "sfmtrng.h"


#define BOLTZMANN_FACTOR(DE, T) exp(-DE/T)

sfmt_t sfmt;
uint32_t *seed_rand;

int main(int argc, char *argv[])
{
    __set_seed_SFMT();
    
    FILE *f_init;
    int8_t *s;
    size_t N;


    char ss[10];
    rand_string(ss, 10);;
    // __fopen(f_init, "init_magn", "rb")
    // printf("%d\n", addition(2, 2));

    // s = __chMalloc(100);

}

