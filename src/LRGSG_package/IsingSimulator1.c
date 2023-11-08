#include "LRGSG_customs.h"
#include "LRGSG_rbim.h"
#include "LRGSG_utils.h"
#include "sfmtrng.h"

#define T_THERM_STEP (thrmSTEP * N)
#define T_EQ_STEP (eqSTEP * N)

#define ISNG_DIR "%sising/"
#define SINI_FNAME ISNG_DIR "N=%zu/s_%s" BINX
#define CLID_FNAME ISNG_DIR "N=%zu/cl%zu_%s" BINX
#define CLOUT_FNAME ISNG_DIR "N=%zu/outcl%zu_%s_T=%.3g_%s" BINX
#define ENE_FNAME ISNG_DIR "N=%zu/ene_%s_T=%.3g_%s" BINX

#define GRPH_DIR "%sgraphs/"
#define ADJ_FNAME GRPH_DIR "N=%zu/adj_%s" BINX
#define EDGL_FNAME GRPH_DIR "N=%zu/edgel_%s" TXTX

sfmt_t sfmt;
uint32_t *seed_rand;

int main(int argc, char *argv[]) {
    /* seed the SFMT RNG */
    __set_seed_SFMT();
    /* variables */
    FILE *f_sini, *f_adj, *f_edgel, *f_ene;
    FILE **f_cl, **f_out;
    char *ptr, *datdir, *code_id, *out_id, buf[STRL512];
    double T, p;
    double *ene;
    spin_tp s;
    size_t N, side, Noclust, tmp, thrmSTEP, eqSTEP;
    size_tp neigh_len, cl_l;
    size_tp *neighs, *cl_i;
    double_p *mclus, *adj, *edgl;
    /* unused variables */
    UNUSED(p);
    UNUSED(argc);
    UNUSED(side);
    /* init variables */
    N = strtozu(argv[1]);
    side = (size_t)sqrt(N);
    T = strtod(argv[2], &ptr);
    p = strtod(argv[3], &ptr);
    code_id = argv[4];
    Noclust = strtozu(argv[5]);
    thrmSTEP = strtozu(argv[6]);
    eqSTEP = strtozu(argv[7]);
    datdir = argv[8];
    out_id = argv[9];
    /* open adjacency matrix file */
    sprintf(buf, ADJ_FNAME, datdir, N, code_id);
    __fopen(&f_adj, buf, "rb");
    /* open magnetization initial condition  */
    sprintf(buf, SINI_FNAME, datdir, N, code_id);
    __fopen(&f_sini, buf, "rb");
    /* open cluster indices files */
    f_cl = malloc(sizeof(*f_cl) * Noclust);
    for (size_t i = 0; i < Noclust; i++) {
        sprintf(buf, CLID_FNAME, datdir, N, i, code_id);
        __fopen((f_cl + i), buf, "rb");
    }
    /* open out files */
    f_out = malloc(sizeof(*f_out) * Noclust);
    for (size_t i = 0; i < Noclust; i++) {
        sprintf(buf, CLOUT_FNAME, datdir, N, i, code_id, T, out_id);
        __fopen((f_out + i), buf, "wb");
    }
    sprintf(buf, ENE_FNAME, datdir, N, code_id, T, out_id);
    __fopen(&f_ene, buf, "wb");
    /* fill adjacency matrix */
    adj = __chMalloc(N * sizeof(*adj));
    for (size_t i = 0; i < N; i++)
        *(adj + i) = __chMalloc(N * sizeof(**adj));
    __fill_adj__(&f_adj, N, &adj);
    /* fill initial condition */
    s = __chMalloc(N * sizeof(*s));
    __fread_check(fread(s, sizeof(*s), N, f_sini), N);
    /* fill clusters indices */
    cl_l = __chCalloc(Noclust, sizeof(*cl_l));
    cl_i = __chMalloc(Noclust * sizeof(*cl_i));
    for (size_t i = 0; i < Noclust; i++) {
        cl_i[i] = __chMalloc((cl_l[i] + 1) * sizeof(**cl_i));
        while (fread(&cl_i[i][cl_l[i]++], sizeof(*cl_i[i]), 1, f_cl[i]) == 1)
            cl_i[i] = realloc(cl_i[i], (cl_l[i] + 1) * sizeof(*cl_i[i]));
        cl_l[i]--;
    }
    /* fill edge list, neighbours list and neighbours lengths */
    edgl = __chMalloc(N * sizeof(*edgl));
    neighs = __chMalloc(N * sizeof(*neighs));
    neigh_len = __chMalloc(N * sizeof(*neigh_len));
    sprintf(buf, EDGL_FNAME, datdir, N, code_id);
    if (__feexist(buf)) {
        __fopen(&f_edgel, buf, "r+");
        __fill_edgl_read__(&f_edgel, N, &edgl, &neighs, &neigh_len);
    } else {
        __fopen(&f_edgel, buf, "w+");
        __fill_edgl_make__(&f_edgel, N, &adj, &edgl, &neighs, &neigh_len);
    }
    /* allocate energy array and mclust */
    ene = __chMalloc(sizeof(*ene) * (T_EQ_STEP + T_THERM_STEP));
    mclus = __chMalloc(sizeof(*mclus) * Noclust);
    for (size_t i = 0; i < Noclust; i++)
        mclus[i] = __chMalloc(sizeof(**mclus) * T_EQ_STEP);

    for (size_t t = 0; t < T_THERM_STEP; t++) {
        ene[t] = calc_energy_full(N, s, neigh_len, neighs, edgl);
        N_step_metropolis(N, T, s, neigh_len, neighs, edgl);
    }
    for (size_t t = 0; t < T_EQ_STEP; t++) {
        ene[t + T_THERM_STEP] = calc_energy_full(N, s, neigh_len, neighs, edgl);
        for (size_t i = 0; i < Noclust; i++)
            mclus[i][t] = calc_clust_magn(cl_l[i], cl_i[i], s);
        N_step_metropolis(N, T, s, neigh_len, neighs, edgl);
    }
    for (size_t i = 0; i < Noclust; i++)
        fwrite(mclus[i], sizeof(**mclus), T_EQ_STEP, f_out[i]);
    fwrite(ene, sizeof(*ene), (T_EQ_STEP + T_THERM_STEP), f_ene);

    fclose(f_edgel);
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

// OLD EDL READ
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

// OLD EGDL MAKE
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