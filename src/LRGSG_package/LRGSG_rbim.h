#include <math.h>
#include <sys/types.h>
#include "LRGSG_customs.h"
#include "sfmtrng.h"

#ifndef __LRGSGRBIMLIB_H_INC__
#define __LRGSGRBIMLIB_H_INC__

#define BOLTZMANN_FACTOR(DE, T) exp(-DE / T)

double neigh_weight_magn(size_t nd, size_t n_nn, spin_tp s, size_tp *neighs,
                         double_p *edgl);
double calc_energy_full(size_t N, spin_tp s, size_tp nlen, size_tp *neighs,
                        double_p *edgl);
double calc_ext_magn(size_t N, spin_tp s);
double calc_ext_magn2(size_t N, spin_tp s);
void flip_spin(size_t nd, spin_tp s);
void one_step_metropolis(size_t nd, double T, spin_tp s, size_tp nlen,
                         size_tp *neighs, double_p *edgl);
double calc_clust_magn(size_t cli_l, size_tp cli, spin_tp s);

#endif /* __LRGSGRBIMLIB_H_INC__ */