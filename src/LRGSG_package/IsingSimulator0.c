#include "LRGSG_utils.h"
#include "LRGSG_customs.h"
#include "LRGSG_rbim.h"
#include "sfmtrng.h"

#define DEFAULT_DATA_OUTDIR "data/l2d_sq_ising/"
#define DEFAULT_GRAPH_OUTDIR DEFAULT_DATA_OUTDIR "graphs/"

#define T_MAX_STEP (10 * N)
#define T_THERM_STEP (2 * N)

sfmt_t sfmt;
uint32_t *seed_rand;



int main(int argc, char *argv[])
{
    // char cwd[1024];
    // getcwd(cwd, sizeof(cwd));
    // printf("Current working dir: %s\n", cwd);
    __set_seed_SFMT();
    //
    FILE *f_sini, *f_adj, *f_icl1, *f_icl2, *f_edgel;
    FILE *f_out, *f_out2;
    char buf[STRL512];
    char *ptr;
    double T;
    double p;
    spin_tp s;
    size_t N;//, side;
    size_t tmp, cntr, cl1i_l = 0, cl2i_l = 0;
    size_tp neigh_len, cl1i, cl2i;
    size_tp *neighs;
    double_p m1, m2;
    double_p *adj;
    double_p *edgl;
    //
    if (argc > 5)
    {
        printf("too many arguments!\n");
    }
    N = strtozu(argv[1]);
    T = strtod(argv[2], &ptr);
    p = strtod(argv[3], &ptr);
    //
    sprintf(buf, DEFAULT_GRAPH_OUTDIR "N=%zu/adj_%s" BINX, N, argv[4]);
    __fopen(&f_adj, buf, "rb");
    sprintf(buf, DEFAULT_GRAPH_OUTDIR "N=%zu/s_%s" BINX, N, argv[4]);
    __fopen(&f_sini, buf, "rb");
    sprintf(buf, DEFAULT_GRAPH_OUTDIR "N=%zu/cl1_%s" BINX, N, argv[4]);
    __fopen(&f_icl1, buf, "rb");
    //
    // side = (size_t)sqrt(N);
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
    // qui sostituisci con una condizione su argc se il programma ha trovato piu di un solo cluster
    cl2i = __chMalloc(1 * sizeof(*cl2i));
    if (p > 0.103)
    {
        sprintf(buf, DEFAULT_GRAPH_OUTDIR "N=%zu/cl2_%s" BINX, N, argv[4]);
        __fopen(&f_icl2, buf, "rb");
        while (fread(&cl2i[cl2i_l++], sizeof(*cl2i), 1, f_icl2) == 1)
            cl2i = realloc(cl2i, (cl2i_l + 1) * sizeof(*cl2i));
        cl2i_l--;
    }
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
    sprintf(buf, DEFAULT_GRAPH_OUTDIR "N=%zu/edgel_%s" TXTX, N, argv[4]);
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
                    fprintf(f_edgel, "%zu %zu %lf\n", i, j, adj[i][j]);
                }
            }
            neigh_len[i] = cntr;
        }
    }

    // printf("energy = %lf\n", calc_energy_full(N, s, neigh_len, neighs, edgl));
    sprintf(buf, DEFAULT_DATA_OUTDIR "N=%zu/outcl1_%s_T=%.3g_%s" TXTX, N, argv[4], T, argv[5]);
    __fopen(&f_out, buf, "a+");
    for (size_t t = 0; t < T_MAX_STEP; t++)
    {
        // printf("cycle %zu\r", t);
        for (size_t i = 0; i < N; i++)
            one_step_metropolis(i, T, s, neigh_len, neighs, edgl);
    }
    m1 = __chMalloc(T_THERM_STEP * sizeof(*m1));
    m2 = __chMalloc(T_MAX_STEP * sizeof(*m2));
    if (p < 0.103)
    {
        for (size_t t = 0; t < T_THERM_STEP; t++)
        {
            m1[t] = calc_clust_magn(cl1i_l, cl1i, s);
            for (size_t i = 0; i < N; i++)
                one_step_metropolis(i, T, s, neigh_len, neighs, edgl);
        }
        fprintf(f_out, "%lf %lf\n", sum_vs(T_THERM_STEP, m1) / T_THERM_STEP,
                sum_vs_2(T_THERM_STEP, m1) / T_THERM_STEP);
    }
    else
    {
        sprintf(buf, DEFAULT_DATA_OUTDIR "N=%zu/outcl2_%s_T=%.3g_%s" TXTX, N, argv[4], T, argv[5]);
        __fopen(&f_out2, buf, "a+");
        for (size_t t = 0; t < T_THERM_STEP; t++)
        {
            m1[t] = calc_clust_magn(cl1i_l, cl1i, s);
            m2[t] = calc_clust_magn(cl2i_l, cl2i, s);
            for (size_t i = 0; i < N; i++)
                one_step_metropolis(i, T, s, neigh_len, neighs, edgl);
        }
        fprintf(f_out, "%lf %lf\n", sum_vs(T_THERM_STEP, m1) / T_THERM_STEP,
                sum_vs_2(T_THERM_STEP, m1) / T_THERM_STEP);
        fprintf(f_out2, "%lf %lf\n", sum_vs(T_THERM_STEP, m2) / T_THERM_STEP,
                sum_vs_2(T_THERM_STEP, m2) / T_THERM_STEP);
    }
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

    // for (size_t i = 0; i < N; i++)
    // {
    //     printf("%+"PRId8" ", s[i]);
    //     if (!((i+1) % side))
    //     {
    //         printf("\n");
    //     }
    // }
    // printf("\n");
    fclose(f_out);
    if (p > 0.103)
        fclose(f_out2);
    fclose(f_sini);
    fclose(f_adj);
    free(neigh_len);
    free(m1);
    free(m2);
    free(cl2i);
    free(s);
    free(cl1i);

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
