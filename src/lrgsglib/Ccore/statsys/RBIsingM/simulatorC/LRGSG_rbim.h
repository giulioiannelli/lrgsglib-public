#include "LRGSG_customs.h"
#include "sfmtrng.h"
#include <math.h>
#include <sys/types.h>

#ifndef __LRGSGRBIMLIB_H_INC__
#define __LRGSGRBIMLIB_H_INC__

#define BOLTZMANN_FACTOR(DE, T) exp(-DE / T)

#define ISNG_DIR "%s/ising/%s/"
#define SINI_FNAME ISNG_DIR "s_p=%.3g%s" BINX
#define CLID_FNAME ISNG_DIR "cl%zu_p=%.3g%s" BINX
#define CLOUT_FNAME ISNG_DIR "outcl%zu_p=%.3g_T=%.3g%s%s"
#define ENE_FNAME ISNG_DIR "ene_p=%.3g_T=%.3g_%s" BINX
#define SOUT_FNAME ISNG_DIR "sout_p=%.3g_T=%.3g_%s" BINX
#define MAGN_FNAME ISNG_DIR "m_p=%.3g_T=%.3g_%s" BINX

#define GRPH_DIR "%s/graph/%s/"
#define EDGL_FNAME GRPH_DIR "edgelist_p=%.3g%s" BINX
#define EIGV_FNAME GRPH_DIR "eigV%d_p=%.3g%s" BINX

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