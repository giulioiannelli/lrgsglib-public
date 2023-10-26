#include "LRGSG_rbim.h"

double calc_clust_magn(size_t cli_l, size_tp cli, spin_tp s)
{
    double clm = 0.;
    for (size_t i = 0; i < cli_l; i++)
        clm += s[cli[i]];
    return clm / cli_l;
}


double neigh_weight_magn(size_t nd, size_t n_nn, spin_tp s, size_tp *neighs, double_p *edgl)
{
    double sum = 0.;
    for (size_t i = 0; i < n_nn; i++)
        sum += *(*(edgl + nd) + i) * *(s + *(*(neighs + nd) + i));
    return -sum / n_nn;
}
double calc_energy_full(size_t N, spin_tp s, size_tp nlen, size_tp *neighs, double_p *edgl)
{
    double sum = 0.;
    for (size_t i = 0; i < N; i++)
        sum += neigh_weight_magn(i, *(nlen + i), s, neighs, edgl);
    return sum;
}
void flip_spin(size_t nd, spin_tp s)
{
    s[nd] = -s[nd];
}
void one_step_metropolis(size_t nd, double T, spin_tp s, size_tp nlen, size_tp *neighs, double_p *edgl)
{
    double nene, E_old, E_new, DeltaE;
    nene = neigh_weight_magn(nd, *(nlen + nd), s, neighs, edgl);
    E_old = s[nd] * nene;
    E_new = -s[nd] * nene;
    DeltaE = E_new - E_old;
    if (DeltaE < 0)
    {
        flip_spin(nd, s);
    }
    else
    {
        if (RNG_dbl() < BOLTZMANN_FACTOR(DeltaE, T))
        {
            flip_spin(nd, s);
        }
    }
}

double calc_ext_magn(size_t N, spin_tp s)
{
    double m = 0.;
    for (size_t i = 0; i < N; i++)
        m += s[i];
    return m;
}
double calc_ext_magn2(size_t N, spin_tp s)
{
    double m2 = 0.;
    for (size_t i = 0; i < N; i++)
        m2 += s[i] * s[i];
    return m2;
}