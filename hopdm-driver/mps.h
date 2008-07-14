typedef struct MPS_struct {

    /* Dimensions. */
    long   m, n, nza, nzq;
    long   max_m, max_n, max_nza, max_nzq;

    /* Control parameters. */
    long   lpqp, presolve, irobj, iolog;
    double objcon, opttol, mult, big;

    /* LP and QP data. */
    long   *clpnts, *qclpts;
    INTS   *rwnmbs, *qrwnbs;
    double *acoeff, *qcoeff, *qdiag;

    /* MPS data. */
    char   *rwname, *clname;
    INTS   *rwstat;
    double *b, *ranges, *lobnd, *upbnd;

} MPS;


MPS*   init_mps (long max_m, long max_n, long max_nza, long max_nzq,
                 long presolve, long iolog, double mult, double big);

void   free_mps (MPS* P);

void   read_mps_data (char* MPS_file, MPS* P);

MPS*   copy_mps (MPS* P);

Lp*    mps2hopdm (MPS* P);

void   FillNames (char C, char *name, long dimension);

int def_MPS_dimensions (MPS* P, long m, long n);

int define_A (MPS* P, long nzA, long irobj,
              long *clpnts, INTS *rwnmbs, double *acoeff);

int define_Q (MPS* P, long nzQ, long lpqp, double *qdiag,
              long *qclpts, INTS *qrwnbs, double *qcoeff);

int def_MPS_data (MPS* P, double objcon, double *b, double *ranges,
                  INTS *rwstat, double *lobnd, double *upbnd);

int define_names (MPS* P, char* rwname, char* clname);
