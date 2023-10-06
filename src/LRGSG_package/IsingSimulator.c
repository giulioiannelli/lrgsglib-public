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
    
    FILE *f_init, *f_elst;
    int8_t *s;
    size_t N = strtozu(argv[1]);
    double_p *adj;

    s = __chMalloc(N * sizeof(*s));
    __fopen(&f_init, argv[2], "rb"); //argv[2]: init_magn_randomstring
    __fread_check(fread(s, sizeof(*s), N, f_init), N);

    adj = __chMalloc(N * sizeof(*adj));
    __fopen(&f_elst, argv[3], "rb"); //argv[2]: init_magn_randomstring

    for (size_t i = 0; i < N; i++)
    {
        adj[i] = __chMalloc(N * sizeof(**adj));
    }

    __fill_adj__(&f_elst, N, &adj);

    fclose(f_init);
    fclose(f_elst);

}

