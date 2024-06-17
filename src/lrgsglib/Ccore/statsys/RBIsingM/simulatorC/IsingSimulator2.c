#include "LRGSG_utils.h"
#include "LRGSG_customs.h"
#include "sfmtrng.h"
#include "LRGSG_rbim.h"

#define ISNG_DIR "%sising/"
#define S_FNAME ISNG_DIR "N=%zu/sOut_%s_T=%.3g_%s" BINX
#define SINI_FNAME ISNG_DIR "N=%zu/s_%s_%s" BINX
#define ENE_FNAME ISNG_DIR "N=%zu/ene_%s_T=%.3g_%s" BINX
#define MAGN_FNAME ISNG_DIR "N=%zu/magn_%s_T=%.3g_%s" BINX


#define GRPH_DIR "%sgraphs/"
#define ADJ_FNAME GRPH_DIR "N=%zu/adj_%s" BINX
// #define EDGL_FNAME GRPH_DIR "N=%zu/edgel_%s" TXTX
#define EDGL_FNAME GRPH_DIR "N=%zu/edgelist_%s" BINX

#define T_THERM_STEP (thrmSTEP * N)
#define T_EQ_STEP (eqSTEP * N)
#define T_STEPS (T_THERM_STEP + T_EQ_STEP)

#define EXPECTED_ARGC 12
//

sfmt_t sfmt;
uint32_t *seed_rand;

int main(int argc, char *argv[])
{
    if (argc != EXPECTED_ARGC) {
        fprintf(stderr, "Usage: %s N T p code_id <free-arg> thrmSTEP eqSTEP \
            datdir out_id update_mode nSampleLog\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    //
    __set_seed_SFMT();
    srand(time(NULL));
    //
    FILE *f_sini, *f_s, *f_ene, *f_magn;
    char buf[STRL512];
    char *ptr, *datdir, *out_id, *in_id, *update_mode;
    int nSampleLog;
    int* logspc;
    double T, p, thrmSTEP;
    double *m, *ene;
    size_t N, side, eqSTEP, tmp;
    spin_tp s;
    size_tp neigh_len;
    NodesEdges node_edges;
    Edges edges;
    /* unused variables */
    UNUSED(p);
    UNUSED(side);
    /* init variables */
    N = strtozu(argv[1]);
    side = (size_t)sqrt(N);
    T = strtod(argv[2], &ptr);
    p = strtod(argv[3], &ptr);
    in_id = argv[4];
    // we miss argv[5]
    thrmSTEP = strtod(argv[6], &ptr);
    eqSTEP = strtozu(argv[7]);
    datdir = argv[8];
    out_id = argv[9];
    update_mode = argv[10];
    nSampleLog = atoi(argv[11]);
    logspc = logspace_int(log10(T_STEPS), &nSampleLog);
    //
    s = __chMalloc(N * sizeof(*s));
    m = __chMalloc(sizeof(*m) * T_STEPS);
    ene = __chMalloc(sizeof(*ene) * T_STEPS);
    //
    sprintf(buf, ENE_FNAME, datdir, N, in_id, T, out_id);
    __fopen(&f_ene, buf, "ab");
    sprintf(buf, MAGN_FNAME, datdir, N, in_id, T, out_id);
    __fopen(&f_magn, buf, "ab");
    sprintf(buf, SINI_FNAME, datdir, N, in_id, out_id);
    __fopen(&f_sini, buf, "rb");
    sprintf(buf, S_FNAME, datdir, N, in_id, T, out_id);
    __fopen(&f_s, buf, "ab");
    __fread_check(fread(s, sizeof(*s), N, f_sini), N);
    sprintf(buf, EDGL_FNAME, datdir, N, in_id);
    process_edges(buf, N, &edges, &node_edges, &neigh_len);
    //  
    size_t freq = (size_t) (T_STEPS / nSampleLog);
    // printf("%lf %zu %d\n", T_STEPS, freq, nSampleLog);
    for (size_t t = 0; t < T_STEPS; t++)
    {
        if (t % freq == 0)
            fwrite(s, sizeof(*s), N, f_s);
        m[t] = calc_magn(N, s);
        ene[t] = calc_totEnergy(N, s, neigh_len, node_edges);
        glauberMetropolis_Nstep(N, T, s, neigh_len, node_edges, update_mode);
    }
    // printf("\n");
    fwrite(ene, sizeof(*ene), T_STEPS, f_ene);
    fwrite(m, sizeof(*m), T_STEPS, f_magn);
    fflush(stdout);  // Clear any existing buffered output
    fwrite(s, sizeof(*s), N, stdout);
    //
    fclose(f_sini);
    fclose(f_s);
    fclose(f_ene);
    fclose(f_magn);
    //
    free(m);
    free(ene);
    free(s);
    free(logspc);
    free(neigh_len);
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
