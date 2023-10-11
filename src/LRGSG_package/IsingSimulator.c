#include <stdio.h>
#include <math.h>
#include <inttypes.h>
#include "LRGSG_utils.h"
#include "LRGSG_customs.h"
#include "sfmtrng.h"

#define BOLTZMANN_FACTOR(DE, T) exp(-DE / T)
#define T_MAX_STEP 10 * N
#define T_THERM_STEP N
#define DEFAULT_DATA_OUTDIR "data/l2d_sq_ising/"
#define DEFAULT_GRAPH_OUTDIR DEFAULT_DATA_OUTDIR "graphs/"
#define BINX ".bin"
#define TXTX ".txt"

sfmt_t sfmt;
uint32_t *seed_rand;

double neigh_weight_magn(size_t nd, size_t n_nn, spin_tp s, size_tp *neighs,
                         double_p *edgl);
double calc_energy_full(size_t N, spin_tp s, size_tp nlen, size_tp *neighs,
                        double_p *edgl);
double calc_ext_magn(size_t N, spin_tp s);
double calc_ext_magn2(size_t N, spin_tp s);
void flip_spin(size_t nd, spin_tp s);
void one_step_metropolis(size_t nd, double T, spin_tp s, size_tp nlen,
                         size_tp *neighs, double_p *edgl);
void prepend(char *s, const char *t);
double calc_clust_magn(size_t cli_l, size_tp cli, spin_tp s);

int main(int argc, char *argv[])
{
    // char cwd[1024];
    // getcwd(cwd, sizeof(cwd));
    // printf("Current working dir: %s\n", cwd);
    __set_seed_SFMT();
    //
    FILE *f_sini, *f_adj, *f_icl1, *f_icl2, *f_edgel;
    FILE *f_out;
    char buf[STRL512];
    char *ptr;
    double T;
    spin_tp s;
    size_t N, side;
    size_t tmp, cntr, cl1i_l = 0, cl2i_l = 0;
    size_tp neigh_len, cl1i, cl2i;
    size_tp *neighs;
    double_p m1, m2;
    double_p *adj;
    double_p *edgl;
    //
    N = strtozu(argv[1]);
    T = strtod(argv[2], &ptr);
    //
    sprintf(buf, DEFAULT_GRAPH_OUTDIR "adj_%s" BINX, argv[3]);
    __fopen(&f_adj, buf, "rb");
    sprintf(buf, DEFAULT_GRAPH_OUTDIR "s_%s" BINX, argv[3]);
    __fopen(&f_sini, buf, "rb");
    sprintf(buf, DEFAULT_GRAPH_OUTDIR "cl1_%s" BINX, argv[3]);
    __fopen(&f_icl1, buf, "rb");
    sprintf(buf, DEFAULT_GRAPH_OUTDIR "cl2_%s" BINX, argv[3]);
    __fopen(&f_icl2, buf, "rb");
    //
    side = (size_t)sqrt(N);
    //
    adj = __chMalloc(N * sizeof(*adj));
    for (size_t i = 0; i < N; i++)
        adj[i] = __chMalloc(N * sizeof(**adj));
    __fill_adj__(&f_adj, N, &adj);
    //
    s = __chMalloc(N * sizeof(*s));
    __fread_check(fread(s, sizeof(*s), N, f_sini), N);
    //
    cl1i = __chMalloc(1 * sizeof(*cl1i));
    while (fread(&cl1i[cl1i_l++], sizeof(*cl1i), 1, f_icl1) == 1)
        cl1i = realloc(cl1i, (cl1i_l + 1) * sizeof(*cl1i));
    cl1i_l--;
    //
    cl2i = __chMalloc(1 * sizeof(*cl2i));
    while (fread(&cl2i[cl2i_l++], sizeof(*cl2i), 1, f_icl2) == 1)
        cl2i = realloc(cl2i, (cl2i_l + 1) * sizeof(*cl2i));
    cl2i_l--;
    //
    // for (size_t i = 0; i < cl1i_l; i++)
    // {
    //     printf("%zu ", cl1i[i]);
    // }
    // printf("\n\n");
    // for (size_t i = 0; i < cl2i_l; i++)
    // {
    //     printf("%zu ", cl2i[i]);
    // }
    // printf("\n");
    // for (size_t i = 0; i < N; i++)
    // {
    //     printf("%"PRId8" ", s[i]);
    // }
    // printf("\n");

    // for (size_t i = 0; i < N; i++)
    // {
    //     for (size_t j = 0; j < N; j++)
    //     {
    //         printf("%.0lf ", adj[i][j]);
    //     }
    //     printf("\n");
    // }
    // questo forse potremo farlo da python...
    edgl = __chMalloc(N * sizeof(*edgl));
    neighs = __chMalloc(N * sizeof(*neighs));
    neigh_len = __chMalloc(N * sizeof(*neigh_len));
    sprintf(buf, DEFAULT_GRAPH_OUTDIR "edgel_%s" TXTX, argv[3] + 4);
    if (__feexist(buf))
    {
        size_t node_i;
        double w_ij;
        __fopen(&f_edgel, buf, "r+");
        for (size_t i =0; i < N; i++)
        {
            edgl[i] = __chMalloc(1 * sizeof(**edgl));
            neighs[i] = __chMalloc(1 * sizeof(**neighs));
        }
        cntr = 0;
        node_i = 0;
        for(size_t i, j; fscanf(f_edgel, "%zu %zu %lf", &i, &j, &w_ij) != EOF;)
        {
            if (i != node_i)
            {
                neigh_len[i] = cntr;
                node_i++;
                cntr = 0;
            }
            neighs[i][cntr] = j;
            edgl[i][cntr] = w_ij;
            edgl[i] = realloc(edgl[i], (++cntr + 1) * sizeof(**edgl));
            neighs[i] = realloc(neighs[i], (cntr + 1) * sizeof(**edgl));
        }
        fclose(f_edgel);
    }
    else
    {
        __fopen(&f_edgel, buf, "w+");
        for (size_t i = 0; i < N; i++)
        {
            cntr = 0;
            edgl[i] = __chMalloc(1 * sizeof(**edgl));
            neighs[i] = __chMalloc(1 * sizeof(**neighs));
            for (size_t j = 0; j < N; j++)
            {
                if (fabs(adj[i][j]) > 0.)
                {
                    neighs[i][cntr] = j;
                    edgl[i][cntr] = adj[i][j];
                    edgl[i] = realloc(edgl[i], (++cntr + 1) * sizeof(**edgl));
                    neighs[i] = realloc(neighs[i], (cntr + 1) * sizeof(**edgl));
                    // printf("neig of %zu_[%zu] = %zu\n", i, cntr-1, j);
                    fprintf(f_edgel, "%zu %zu %lf\n", i, j, adj[i][j]);
                }
            }
            neigh_len[i] = cntr;
        }
    }

    // printf("energy = %lf\n", calc_energy_full(N, s, neigh_len, neighs, edgl));
    sprintf(buf, DEFAULT_DATA_OUTDIR "out_%s" TXTX, argv[4]);
    __fopen(&f_out, buf, "a+");
    for (size_t t = 0; t < T_MAX_STEP; t++)
    {
        // printf("cycle %zu\r", t);
        for (size_t i = 0; i < N; i++)
            one_step_metropolis(i, T, s, neigh_len, neighs, edgl);
    }
    m1 = __chMalloc(T_MAX_STEP * sizeof(*m1));
    // m2 = __chMalloc(T_MAX_STEP * sizeof(*m2));
    for (size_t t = 0; t < T_THERM_STEP; t++)
    {
        m1[t] = calc_clust_magn(cl1i_l, cl1i, s);
        // m2[t] = calc_clust_magn(cl2i_l, cl2i, s);
        for (size_t i = 0; i < N; i++)
            one_step_metropolis(i, T, s, neigh_len, neighs, edgl);
    }
    fprintf(f_out, "%lf %lf\n", sum_vs(T_THERM_STEP, m1) / T_THERM_STEP,
            sum_vs_2(T_THERM_STEP, m1) / T_THERM_STEP);
    // fprintf(f_out, "%lf %lf\n", sum_vs(T_THERM_STEP, m2) / T_THERM_STEP,
    //         sum_vs_2(T_THERM_STEP, m2) / T_THERM_STEP);
    // printf("energy = %lf\n", calc_energy_full(N, s, neigh_len, neighs, edgl));
    // printf("side = %zu\n", side);

    // for (size_t i = 0; i < N; i++)
    // {
    //     printf("%+" PRId8 " ", s[i]);
    //     if (!((i + 1) % side))
    //     {
    //         printf("\n");
    //     }
    // }
    // printf("\n");

    fclose(f_out);
    fclose(f_sini);
    fclose(f_adj);
    free(neigh_len);
    free(m1);
    // free(m2);
    free(s);

    tmp = N;
    while (tmp)
        free(adj[--tmp]);
    free(adj);

    tmp = N;
    while (tmp)
        free(edgl[--tmp]);
    free(edgl);

    tmp = N;
    while (tmp)
        free(neighs[--tmp]);
    free(neighs);
}
double calc_clust_magn(size_t cli_l, size_tp cli, spin_tp s)
{
    double clm = 0.;
    for (size_t i = 0; i < cli_l; i++)
        clm += s[cli[i]];
    return clm / cli_l;
}

/* Prepends t into s. Assumes s has enough space allocated
** for the combined string.
*/
void prepend(char *s, const char *t)
{
    size_t len = strlen(t);
    memmove(s + len, s, strlen(s) + 1);
    memcpy(s, t, len);
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