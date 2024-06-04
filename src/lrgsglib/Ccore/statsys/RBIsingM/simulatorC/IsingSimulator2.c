#include "LRGSG_utils.h"
#include "LRGSG_customs.h"
#include "sfmtrng.h"
#include "LRGSG_rbim.h"

#define ISNG_DIR "%sising/"
#define S_FNAME ISNG_DIR "N=%zu/sOut_%s_T=%.3g_%s" BINX
#define SINI_FNAME ISNG_DIR "N=%zu/s_%s" BINX
#define ENE_FNAME ISNG_DIR "N=%zu/ene_%s_T=%.3g_%s" BINX
#define MAGN_FNAME ISNG_DIR "N=%zu/magn_%s_T=%.3g_%s" BINX


#define GRPH_DIR "%sgraphs/"
#define ADJ_FNAME GRPH_DIR "N=%zu/adj_%s" BINX
#define EDGL_FNAME GRPH_DIR "N=%zu/edgel_%s" TXTX

#define T_THERM_STEP (thrmSTEP * N)
#define T_EQ_STEP (eqSTEP * N)

sfmt_t sfmt;
uint32_t *seed_rand;

int main(int argc, char *argv[])
{
    __set_seed_SFMT();
    //
    FILE *f_sini, *f_adj, *f_edgel, *f_s;
    char buf[STRL512];
    char *ptr, *datdir, *out_id,  *code_id, *update_mode;
    double T, p;
    double *m, *ene;
    int* logspc;
    spin_tp s;
    size_t N, side, thrmSTEP, eqSTEP;//, side;
    size_t tmp;
    size_tp neigh_len;
    size_tp *neighs;
    double_p *adj;
    double_p *edgl;
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
    // we miss argv[5]
    thrmSTEP = strtozu(argv[6]);
    eqSTEP = strtozu(argv[7]);
    datdir = argv[8];
    out_id = argv[9];
    update_mode = argv[10];
    logspc = logspace_int(0, log10(T_EQ_STEP + T_THERM_STEP), 100);
    //
    sprintf(buf, ADJ_FNAME, datdir, N, argv[4]);
    __fopen(&f_adj, buf, "rb");
    sprintf(buf, SINI_FNAME, datdir, N, argv[4]);
    __fopen(&f_sini, buf, "rb");
    sprintf(buf, S_FNAME, datdir, N, code_id, T, out_id);
    printf("%s\n", buf);
    __fopen(&f_s, buf, "wb");

    /* fill adjacency matrix */
    adj = __chMalloc(N * sizeof(*adj));
    for (size_t i = 0; i < N; i++)
        *(adj + i) = __chMalloc(N * sizeof(**adj));
    __fill_adj__(&f_adj, N, &adj);
    /* fill initial condition */
    s = __chMalloc(N * sizeof(*s));
    __fread_check(fread(s, sizeof(*s), N, f_sini), N);
    // questo forse potremo farlo da python...
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

    m = __chMalloc(sizeof(*m) * (T_THERM_STEP + T_EQ_STEP));
    ene = __chMalloc(sizeof(*ene) * (T_EQ_STEP + T_THERM_STEP));
    
    for (size_t t = 0; t < T_EQ_STEP + T_THERM_STEP; t++)
    {
        if (t+1 == logspc[t])
        {
            fwrite(s, sizeof(*s), (T_EQ_STEP + T_THERM_STEP), f_s);
        }
        glauber_metropolis_Nstep_mode(N, T, s, neigh_len, neighs, edgl, update_mode);
    }

    fclose(f_sini);
    fclose(f_adj);
    fclose(f_s);

    free(neigh_len);
    free(m);
    free(ene);
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
