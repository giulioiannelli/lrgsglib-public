#include "LRGSG_customs.h"
#include "sfmtrng.h"
#include <math.h>
#include <sys/types.h>

#ifndef __LRGSGRBIMLIB_H_INC__
#define __LRGSGRBIMLIB_H_INC__

#define BOLTZMANN_FACTOR(DE, T) exp(-DE / T)

extern sfmt_t sfmt;
extern uint32_t *seed_rand;

void glauber_metropolis_1step(size_t nd, double T, spin_tp s, size_tp nlen,
                              size_tp *neighs, double_p *edgl);
void glauber_metropolis_Nstep(size_t N, double T, spin_tp s, size_tp nlen,
                              size_tp *neighs, double_p *edgl);
void glauber_metropolis_Nstep_mode(size_t N, double T, spin_tp s, size_tp nlen,
                       size_tp *neighs, double_p *edgl, const char *update_mode);
double calc_energy_full(size_t N, spin_tp s, size_tp nlen, size_tp *neighs,
                        double_p *edgl);
void glauberMetropolis_Nstep(size_t N, double T, spin_tp s, size_tp nlen,
                       NodesEdges ne, const char *update_mode);
void glauberMetropolis_1step(size_t nd, double T, spin_tp s, size_t nlen,
                         NodeEdges ne);
#endif /* __LRGSGRBIMLIB_H_INC__ */