#include "LRGSG_customs.h"
#include "sfmtrng.h"
#include <math.h>
#include <sys/types.h>

#ifndef __LRGSGRBIMLIB_H_INC__
#define __LRGSGRBIMLIB_H_INC__

#define BOLTZMANN_FACTOR(DE, T) exp(-DE / T)

void glauber_metropolis_1step(size_t nd, double T, spin_tp s, size_tp nlen,
                              size_tp *neighs, double_p *edgl);
void glauber_metropolis_Nstep(size_t N, double T, spin_tp s, size_tp nlen,
                              size_tp *neighs, double_p *edgl);
double calc_energy_full(size_t N, spin_tp s, size_tp nlen, size_tp *neighs,
                        double_p *edgl);

#endif /* __LRGSGRBIMLIB_H_INC__ */