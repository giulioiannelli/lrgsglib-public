#include "LRGSG_utils.h"
#include "LRGSG_rbim.h"

/** 
 * @brief applies the Glauber-Metropolis dynamics to an Ising system (1 step).
 * 
 * @param nd         (size_t) Index of the spin in the system.
 * @param T          (double) Temperature of the system.
 * @param s          (spin_tp) Array representing the spin configuration.
 * @param nlen       (size_tp) Length of the neighbor list.
 * @param neighs     (size_tp *) Array containing neighbor indices.
 * @param edgl       (double_p *) Array containing edge weights.
 * 
 * @return           None
 */
void glauber_metropolis_1step(size_t nd, double T, spin_tp s, size_tp nlen,
                         size_tp *neighs, double_p *edgl) {
    double nene = neigh_weight_magn(nd, *(nlen + nd), s, neighs, edgl);
    double DE = 2 * s[nd] * nene;
    if (DE < 0) {
        flip_spin(nd, s);
    } else if (RNG_dbl() < BOLTZMANN_FACTOR(DE, T)) {
        flip_spin(nd, s);
    }
}
/** 
 * @brief applies the Glauber-Metropolis dynamics to an Ising system (N steps).
 * 
 * @param N          (size_t) Number of steps to evolve the system.
 * @param T          (double) Temperature of the system.
 * @param s          (spin_tp) Array representing the spin configuration.
 * @param nlen       (size_tp) Length of the neighbor list.
 * @param neighs     (size_tp *) Array containing neighbor indices.
 * @param edgl       (double_p *) Array containing edge weights.
 * 
 * @return           None
 */
void glauber_metropolis_Nstep(size_t N, double T, spin_tp s, size_tp nlen,
                       size_tp *neighs, double_p *edgl) {
    for (size_t i = 0; i < N; i++)
        glauber_metropolis_1step(i, T, s, nlen, neighs, edgl);
}
