#include "LRGSG_utils.h"
#include "LRGSG_customs.h"
#include "sfmtrng.h"
#include "LRGSG_rbim.h"
//
#define EXPECTED_ARGC 11+1
#define MOD_SAVE 0
//
int main(int argc, char *argv[])
{
    /* check argc */
    if (argc != EXPECTED_ARGC) {
        fprintf(stderr, "Usage: %s N T p Noclust thrmSTEP eqSTEP datdir run_id"\
            " out_id update_mode nSampleLog NeigV\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    /* seed the SFMT RNG */
    __set_seed_SFMT();
    srand(time(NULL));
    /* variables */
    FILE *f_sini;
    FILE *f_sout, *f_ene, *f_m;
    char buf[STRL512];
    char modified_out_id[STRL256];
    char *ptr, *datdir, *syshape, *out_id;
    char *update_mode, *run_id;
    int nSampleLog;
    int* logspc;
    double T, p;
    double *ene, *m;
    size_t N, tmp, freq;
    spin_tp s;
    size_tp neigh_len;
    NodesEdges node_edges;
    Edges edges;
    /* init variables */
    N = strtozu(argv[1]);
    T = strtod(argv[2], &ptr);
    p = strtod(argv[3], &ptr);
    thrmSTEP = strtod(argv[4], &ptr);
    eqSTEP = strtod(argv[5], &ptr);
    datdir = argv[6];
    syshape = argv[7];
    run_id = argv[8];
    out_id = argv[9];
    update_mode = argv[10];
    nSampleLog = atoi(argv[11]);
    freq = (size_t) (T_STEPS / nSampleLog);
    logspc = logspace_int(log10(T_STEPS), &nSampleLog);
    /* init metropolis algorithm */
    initialize_glauberMetropolis(T);
    /* allocate spin, energy and mclust */
    s = __chMalloc(N * sizeof(*s));
    ene = __chMalloc(sizeof(*ene) * T_STEPS);
    m = __chMalloc(sizeof(*m) * T_STEPS);
    /* open spin initial condition  */
    sprintf(buf, SINI_FNAME, datdir, syshape, p, run_id);
    __fopen(&f_sini, buf, "rb");
    __fread_check(fread(s, sizeof(*s), N, f_sini), N);
    sprintf(modified_out_id, "%s_%s", run_id, out_id);
    sprintf(buf, SOUT_FNAME, datdir, syshape, p, T, modified_out_id);
    __fopen(&f_sout, buf, "wb");
    /* fill edge list, neighbours list and neighbours lengths */
    sprintf(buf, EDGL_FNAME, datdir, syshape, p, run_id);
    process_edges(buf, N, &edges, &node_edges, &neigh_len);
    /* simulate Ising (thermalization) */
    for (size_t t = 0; t < T_THERM_STEP; t++) {
        if (t % freq == 0)
            fwrite(s, sizeof(*s), N, f_sout);
        ene[t] = calc_totEnergy(N, s, neigh_len, node_edges);
        m[t] = calc_magn(N, s);
        glauberMetropolis_Nstep(N, T, s, neigh_len, node_edges, update_mode);
    }
    /* simulate Ising (equilibrium) */
    for (size_t t = 0; t < T_EQ_STEP; t++) {
        if ((t + T_THERM_STEP) % freq == 0)
            fwrite(s, sizeof(*s), N, f_sout);
        if (glauber_isStableAtZeroTemp(N, s, neigh_len, node_edges)) {
            fill_array_with_value(ene, T_THERM_STEP-t, T_STEPS, ene[T_THERM_STEP-1]);
            fill_array_with_value(m, T_THERM_STEP-t, T_STEPS, m[T_THERM_STEP-1]);
            break;
        }
        ene[t + T_THERM_STEP] = calc_totEnergy(N, s, neigh_len, node_edges);
        m[t + T_THERM_STEP] = calc_magn(N, s);
        glauberMetropolis_Nstep(N, T, s, neigh_len, node_edges, update_mode);
    }
    switch (MOD_SAVE)
    {
        case 0:
            sprintf(buf, ENE_FNAME, datdir, syshape, p, T, modified_out_id);
            __fopen(&f_ene, buf, "ab");
            fwrite(ene, sizeof(*ene), T_STEPS, f_ene);
            sprintf(buf, MAGN_FNAME, datdir, syshape, p, T, modified_out_id);
            __fopen(&f_m, buf, "ab");
            fwrite(m, sizeof(*m), T_STEPS, f_m);
            break;
        case 1:
            if (out_id[0] != '\0') {
                sprintf(modified_out_id, "_%s", out_id);  // Prepend an underscore if out_id is not empty
            } else {
                modified_out_id[0] = '\0';  // Make sure it's an empty string if out_id is empty
            }
            sprintf(buf, ENEF_FNAME, datdir, syshape, p, T, modified_out_id);
            __fopen(&f_ene, buf, "a+");
            sprintf(buf, MAGNF_FNAME, datdir, syshape, p, T, modified_out_id);
            __fopen(&f_m, buf, "a+");
            fprintf(f_m, "%g %g\n",
                sum_vs(T_EQ_STEP, m + T_THERM_STEP) / T_EQ_STEP,
                sum_vs_2(T_EQ_STEP, m + T_THERM_STEP) / T_EQ_STEP);
            fprintf(f_ene, "%g %g\n",
                sum_vs(T_EQ_STEP, ene + T_THERM_STEP) / T_EQ_STEP,
                sum_vs_2(T_EQ_STEP, ene + T_THERM_STEP) / T_EQ_STEP);
            break;
    }
    fclose(f_ene);
    fclose(f_m);
    fflush(stdout);
    fwrite(s, sizeof(*s), N, stdout);
    /* closing files and freeing arrays*/
    fclose(f_sini);
    free(s);
    free(logspc);
    free(neigh_len);
    free(ene);
    free(edges);
    tmp = N;
    while (tmp)
    {
        free(node_edges[--tmp].neighbors);
        free(node_edges[tmp].weights);
    }
    free(node_edges);
    return 0;
}