#include "LRGSG_utils.h"
#include "LRGSG_customs.h"
#include "LRGSG_rbim.h"
#include "sfmtrng.h"

#define T_MAX_STEP (2 * N)
#define T_THERM_STEP (10 * N)
#define DEFAULT_DATA_OUTDIR "data/l2d_sq_ising/"
#define DEFAULT_GRAPH_OUTDIR DEFAULT_DATA_OUTDIR "graphs/"

sfmt_t sfmt;
uint32_t *seed_rand;

int main(int argc, char *argv[])
{
    /* seed the SFMT RNG */
    __set_seed_SFMT();
    /* variables */
    FILE **f_cl, **f_out, *f_ene;
    FILE *f_sini, *f_adj, *f_edgel;
    char buf[STRL512];
    char *ptr;
    char *code_id, *out_id;
    double T, p;
    spin_tp s;
    size_t N, side, Noclust;
    size_t tmp, *cl_l;
    // sizet cntr;
    size_tp neigh_len, *cl_i;
    size_tp *neighs;
    double *ene;
    double_p *mclus;
    double_p *adj;
    double_p *edgl;
    UNUSED(p);
    UNUSED(argc);
    UNUSED(side);
    /* init variables */
    N = strtozu(argv[1]);
    side = (size_t) sqrt(N);
    T = strtod(argv[2], &ptr);
    p = strtod(argv[3], &ptr);
    code_id = argv[4];
    Noclust = strtozu(argv[5]);
    out_id = argv[6];
    /* open adjacency matrix file */
    sprintf(buf, DEFAULT_GRAPH_OUTDIR "N=%zu/adj_%s" BINX, N, code_id);
    __fopen(&f_adj, buf, "rb");
    /* open magnetization initial condition  */
    sprintf(buf, DEFAULT_GRAPH_OUTDIR "N=%zu/s_%s" BINX, N, code_id);
    __fopen(&f_sini, buf, "rb");
    /* open cluster indices files */
    f_cl = malloc(sizeof(*f_cl) * Noclust);
    for (size_t i = 0; i < Noclust; i++)
    {
        sprintf(buf, DEFAULT_GRAPH_OUTDIR "N=%zu/cl%zu_%s" BINX, N, i, code_id);
        __fopen(&f_cl[i], buf, "rb");
    }
    adj = __chMalloc(N * sizeof(*adj));
    for (size_t i = 0; i < N; i++)
        adj[i] = __chMalloc(N * sizeof(**adj));
    __fill_adj__(&f_adj, N, &adj);
    //
    s = __chMalloc(N * sizeof(*s));
    __fread_check(fread(s, sizeof(*s), N, f_sini), N);
    //
    cl_l = __chCalloc(Noclust, sizeof(*cl_l));
    cl_i = __chMalloc(Noclust * sizeof(*cl_i));
    for (size_t i = 0; i < Noclust; i++)
    {
        cl_i[i] = __chMalloc((cl_l[i] + 1) * sizeof(**cl_i));
        while (fread(&cl_i[i][cl_l[i]++], sizeof(*cl_i[i]), 1, f_cl[i]) == 1)
            cl_i[i] = realloc(cl_i[i], (cl_l[i] + 1) * sizeof(*cl_i[i]));
        cl_l[i]--;
    }
    edgl = __chMalloc(N * sizeof(*edgl));
    neighs = __chMalloc(N * sizeof(*neighs));
    neigh_len = __chMalloc(N * sizeof(*neigh_len));
    sprintf(buf, DEFAULT_GRAPH_OUTDIR "N=%zu/edgel_%s" TXTX, N, code_id);
    if (__feexist(buf))
    {

        __fopen(&f_edgel, buf, "r+");
        for (size_t i =0; i < N; i++)
        {
            edgl[i] = __chMalloc(1 * sizeof(**edgl));
            neighs[i] = __chMalloc(1 * sizeof(**neighs));
        }
        __fill_edgl_read__(&f_edgel, &edgl, &neighs, &neigh_len);
        // size_t node_i;
        // double w_ij;
        // cntr = 0;
        // node_i = 0;
        // for(size_t i, j; fscanf(f_edgel, "%zu %zu %lf", &i, &j, &w_ij) != EOF;)
        // {
        //     if (i != node_i)
        //     {
        //         neigh_len[i] = cntr;
        //         node_i++;
        //         cntr = 0;
        //     }
        //     neighs[i][cntr] = j;
        //     edgl[i][cntr] = w_ij;
        //     edgl[i] = realloc(edgl[i], (++cntr + 1) * sizeof(**edgl));
        //     neighs[i] = realloc(neighs[i], (cntr + 1) * sizeof(**edgl));
        // }
    }
    else
    {
        __fopen(&f_edgel, buf, "w+");
        __fill_edgl_make__(&f_edgel, N, &adj, &edgl, &neighs, &neigh_len);
        //     for (size_t i = 0; i < N; i++)
        //     {
        //         cntr = 0;
        //         edgl[i] = __chMalloc(1 * sizeof(**edgl));
        //         neighs[i] = __chMalloc(1 * sizeof(**neighs));
        //         for (size_t j = 0; j < N; j++)
        //         {
        //             if (fabs(adj[i][j]) > 0.)
        //             {
        //                 neighs[i][cntr] = j;
        //                 edgl[i][cntr] = adj[i][j];
        //                 edgl[i] = realloc(edgl[i], (++cntr + 1) * sizeof(**edgl));
        //                 neighs[i] = realloc(neighs[i], (cntr + 1) * sizeof(**edgl));
        //                 fprintf(f_edgel, "%zu %zu %lf\n", i, j, adj[i][j]);
        //             }
        //         }
        //         neigh_len[i] = cntr;
        //     }
        // }
    }
    fclose(f_edgel);
    // for (size_t i = 0; i < N; i++)
    // {
    //     printf("%zu %zu ", i, neigh_len[i]);
    //     for (size_t j = 0; j < neigh_len[i]; j++)
    //     {
    //         printf("%zu ", neighs[i][j]);
    //     }
    //     for (size_t j = 0; j < neigh_len[i]; j++)
    //     {
    //         printf("%lf ", edgl[i][j]);
    //     }
    //     printf("\n");
    // }
    // printf("\n");
    // printf("energy = %lf\n", calc_energy_full(N, s, neigh_len, neighs, edgl));
    f_out = malloc(sizeof(*f_out) * Noclust);
    for (size_t i = 0; i < Noclust; i++)
    {
        sprintf(buf, DEFAULT_DATA_OUTDIR "N=%zu/outcl%zu_%s_T=%.3g_%s" BINX, N, i, code_id, T, out_id);
        __fopen(&f_out[i], buf, "wb");
    }
    sprintf(buf, DEFAULT_DATA_OUTDIR "N=%zu/ene_%s_T=%.3g_%s" BINX, N, code_id, T, out_id);
    __fopen(&f_ene, buf, "wb");
    ene = malloc(sizeof(*ene) * (T_THERM_STEP + T_MAX_STEP));
    mclus = __chMalloc(Noclust * sizeof(*mclus));
    for (size_t i = 0; i < Noclust; i++)
        mclus[i] = __chMalloc(T_THERM_STEP * sizeof(**mclus));
    
    for (size_t t = 0; t < T_MAX_STEP; t++)
    {
        ene[t] = calc_energy_full(N, s, neigh_len, neighs, edgl);
        for (size_t i = 0; i < N; i++)
            one_step_metropolis(i, T, s, neigh_len, neighs, edgl);
    }
    for (size_t t = 0; t < T_THERM_STEP; t++)
    {
        ene[t] = calc_energy_full(N, s, neigh_len, neighs, edgl);
        for (size_t i = 0; i < Noclust; i++)
            mclus[i][t] = calc_clust_magn(cl_l[i], cl_i[i], s);
        for (size_t i = 0; i < N; i++)
            one_step_metropolis(i, T, s, neigh_len, neighs, edgl);
    }
    for (size_t i = 0; i < Noclust; i++)
        fwrite(mclus[i], sizeof(**mclus), T_THERM_STEP, f_out[i]);
    fwrite(ene, sizeof(*ene), T_MAX_STEP, f_ene);
    // printf("energia positiva = %lf\n", calc_energy_full(N, s, neigh_len, neighs, edgl));
        // fprintf(f_out, "%lf %lf\n", sum_vs(T_THERM_STEP, mclus[i]) / T_THERM_STEP,
        //         sum_vs_2(T_THERM_STEP, mclus[i]) / T_THERM_STEP);

    // printf("energy = %lf\n", calc_energy_full(N, s, neigh_len, neighs, edgl));
    // printf("side = %zu\n", side);

    fclose(f_ene);
    fclose(f_sini);
    fclose(f_adj);
    tmp = Noclust;
    while (tmp)
        fclose(f_cl[--tmp]);
    free(f_cl);
    tmp = Noclust;
    while (tmp)
        fclose(f_out[--tmp]);
    free(f_out);
    tmp = Noclust;
    while (tmp)
        free(cl_i[--tmp]);
    free(neigh_len);
    
    free(s);
    free(cl_l);
    free(ene);
    tmp = Noclust;
    while (tmp)
        free(mclus[--tmp]);
    free(mclus);
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