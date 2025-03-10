#include "LRGSG_customs.h"
#include "sfmtrng.h"
#include <math.h>
#include <sys/types.h>

#ifndef __LRGSGRBIMLIB_H_INC__
#define __LRGSGRBIMLIB_H_INC__

#define BOLTZMANN_FACTOR(DE, T) exp(-DE / T)
#define T_THERM_STEP  (size_t)(thrmSTEP * N)
#define T_EQ_STEP (size_t)(eqSTEP * N)
#define T_STEPS (T_THERM_STEP + T_EQ_STEP)

#define ISNG_DIR "%s/ising/%s/"
#define SINI_FNAME ISNG_DIR "s_p=%.3g%s" BINX
#define CLID_FNAME ISNG_DIR "cl%zu_p=%.3g%s" BINX
#define CLOUT_FNAME ISNG_DIR "outcl%zu_p=%.3g_T=%.3g%s%s"
#define ENE_FNAME ISNG_DIR "ene_p=%.3g_T=%.3g_%s" BINX
#define ENEF_FNAME ISNG_DIR "ene_p=%.3g_T=%.3g_%s" TXTX
#define SOUT_FNAME ISNG_DIR "sout_p=%.3g_T=%.3g_%s" BINX
#define MAGN_FNAME ISNG_DIR "m_p=%.3g_T=%.3g_%s" BINX
#define MAGNF_FNAME ISNG_DIR "m_p=%.3g_T=%.3g_%s" TXTX

#define GRPH_DIR "%s/graph/%s/"
#define EDGL_FNAME GRPH_DIR "edgelist_p=%.3g%s" BINX
#define EIGV_FNAME GRPH_DIR "eigV%d_p=%.3g%s" BINX
#define ADJ_FNAME GRPH_DIR "adj_p=%.3g%s" BINX

typedef void (*glauberMetropolisFunc)(size_t, double, spin_tp, size_t, NodeEdges);

extern sfmt_t sfmt;
extern uint32_t *seed_rand;
extern double thrmSTEP;
extern double eqSTEP;
extern glauberMetropolisFunc glauberMetropolis_1step_ptr;


void initialize_glauberMetropolis(double T);
void glauberMetropolis_1step(size_t nd, double T, spin_tp s, size_t nlen,
    NodeEdges ne);
void glauberMetropolis_1step_T0(size_t nd, double T, spin_tp s, size_t nlen,
    NodeEdges ne);
void glauberMetropolis_Nstep(size_t N, double T, spin_tp s, size_tp nlen,
    NodesEdges ne, const char *update_mode);


void glauber_metropolis_1step(size_t nd, double T, spin_tp s, size_tp nlen,
                              size_tp *neighs, double_p *edgl);
void glauber_metropolis_Nstep(size_t N, double T, spin_tp s, size_tp nlen,
                              size_tp *neighs, double_p *edgl);
void glauber_metropolis_Nstep_mode(size_t N, double T, spin_tp s, size_tp nlen,
                       size_tp *neighs, double_p *edgl, const char *update_mode);
double calc_energy_full(size_t N, spin_tp s, size_tp nlen, size_tp *neighs,
                        double_p *edgl);

#endif /* __LRGSGRBIMLIB_H_INC__ */