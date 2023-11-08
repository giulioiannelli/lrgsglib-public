#include "LRGSG_utils.h"
#include "LRGSG_customs.h"
#include "LRGSG_rbim.h"
#include "sfmtrng.h"

#define CLOUT_FNAME "%sN=%zu/outcl%zu_%s_T=%.3g_%s" BINX
#define ENE_FNAME "%sN=%zu/ene_%s_T=%.3g_%s" BINX

#define GRPH_DIR "%sgraphs/"
#define ADJ_FNAME GRPH_DIR "N=%zu/adj_%s" BINX
#define SINI_FNAME GRPH_DIR "N=%zu/s_%s" BINX
#define CLID_FNAME GRPH_DIR "N=%zu/cl%zu_%s" BINX
#define EDGL_FNAME GRPH_DIR "N=%zu/edgel_%s" TXTX

#define T_THERM_STEP (thrmSTEP * N)
#define T_EQ_STEP (eqSTEP * N)

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
    char *ptr, *datdir;
    double T;
    double p;
    spin_tp s;
    size_t N, thrmSTEP, eqSTEP;//, side;
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
    thrmSTEP = strtozu(argv[7]);
    eqSTEP = strtozu(argv[8]);
    datdir = argv[6];
    //
    sprintf(buf, ADJ_FNAME, datdir, N, argv[4]);
    __fopen(&f_adj, buf, "rb");
    sprintf(buf, SINI_FNAME, datdir, N, argv[4]);
    __fopen(&f_sini, buf, "rb");
    sprintf(buf, CLID_FNAME, datdir, N, (size_t)1, argv[4]);
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
        sprintf(buf, CLID_FNAME, datdir, N, (size_t)2, argv[4]);
        __fopen(&f_icl2, buf, "rb");
        while (fread(&cl2i[cl2i_l++], sizeof(*cl2i), 1, f_icl2) == 1)
            cl2i = realloc(cl2i, (cl2i_l + 1) * sizeof(*cl2i));
        cl2i_l--;
    }

    // questo forse potremo farlo da python...
    edgl = __chMalloc(N * sizeof(*edgl));
    neighs = __chMalloc(N * sizeof(*neighs));
    neigh_len = __chMalloc(N * sizeof(*neigh_len));
    sprintf(buf, EDGL_FNAME, datdir, N, argv[4]);
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
    sprintf(buf, CLOUT_FNAME, datdir, N, (size_t)1, argv[4], T, argv[5]);
    __fopen(&f_out, buf, "a+");
    for (size_t t = 0; t < T_THERM_STEP; t++)
    {
        // printf("cycle %zu\r", t);
        for (size_t i = 0; i < N; i++)
            one_step_metropolis(i, T, s, neigh_len, neighs, edgl);
    }
    m1 = __chMalloc(T_EQ_STEP * sizeof(*m1));
    m2 = __chMalloc(T_THERM_STEP * sizeof(*m2));
    if (p < 0.103)
    {
        for (size_t t = 0; t < T_EQ_STEP; t++)
        {
            m1[t] = calc_clust_magn(cl1i_l, cl1i, s);
            for (size_t i = 0; i < N; i++)
                one_step_metropolis(i, T, s, neigh_len, neighs, edgl);
        }
        fprintf(f_out, "%lf %lf\n", sum_vs(T_EQ_STEP, m1) / T_EQ_STEP,
                sum_vs_2(T_EQ_STEP, m1) / T_EQ_STEP);
    }
    else
    {
        sprintf(buf, CLOUT_FNAME, datdir, N, (size_t)2, argv[4], T, argv[5]);
        __fopen(&f_out2, buf, "a+");
        for (size_t t = 0; t < T_EQ_STEP; t++)
        {
            m1[t] = calc_clust_magn(cl1i_l, cl1i, s);
            m2[t] = calc_clust_magn(cl2i_l, cl2i, s);
            for (size_t i = 0; i < N; i++)
                one_step_metropolis(i, T, s, neigh_len, neighs, edgl);
        }
        fprintf(f_out, "%lf %lf\n", sum_vs(T_EQ_STEP, m1) / T_EQ_STEP,
                sum_vs_2(T_EQ_STEP, m1) / T_EQ_STEP);
        fprintf(f_out2, "%lf %lf\n", sum_vs(T_EQ_STEP, m2) / T_EQ_STEP,
                sum_vs_2(T_EQ_STEP, m2) / T_EQ_STEP);
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
