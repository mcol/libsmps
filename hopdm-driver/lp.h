#ifndef FNAME
#define FNAME 30
#endif

#ifndef MPSNAME
#define MPSNAME 20
#endif


typedef struct Lp_struct {

    /* Dimensions. */
    long   m, n, n_struct, ninmps, nza, m0, mfinal;
    long   max_m, max_n, max_l, max_nza, lauxwk;

    /* Control parameters. */
    long   msglev, presolve, iolog, exit_code, sleep, irobj;
    long   lpqp, mfact, iqp, itslv, kptype, icnvx, freevr;
    long   scaleRows, scaleCols, scaleObj, scaleRhs;
    double mult, objcon, lpobj, qpobj, objective;

    /* Presolve history list. */
    long   lnhist, mxhist;
    INTS   *rwperm, *rwinvp;
    long   *inhsti, *inhstj;
    double *dphist;

    /* LP constraint matrix. */
    long   *clpnts, *rwhead, *rwlink;
    INTS   *rwnmbs, *clnmbs, *lencol;
    double *acoeff;

    /* MPS data. */
    char   *rwname, *clname;
    INTS   *stavar, *vrsign, *rwstat, *starow;
    double *b, *ranges, *p, *q, *rscale;
    double *c, *lobnd, *upbnd, *cscale;
    long   *vused, *vbnded;

    /* A copy of MPS data. */
    long   mcopy, ncopy;
    char   *rwcopy;
    INTS   *cstvar, *cstrow, *crwst, *rwnbcp, *lnclcp;
    double *bcopy, *ccopy, *lobndc, *upbndc, *acffcp;
    long   *clptcp;

    /* A copy of MPS data used by the sparse preconditioner. */
    long   *cl2pts, *rw2hd, *rw2lnk;
    INTS   *rw2nbs, *cl2nbs, *ln2col, *isllt;
    double *a2cff;

    /* QP data. */
    long   nzq, maxnzq;
    INTS   *pvtype, *qrwnbs;
    long   *qclpts;
    double *qcoeff, *qdiag, *objlqp;

    /* Cholesky factorization arrays. */
    long   nzl, maxnzl, iwindw;
    INTS   *perm, *invp, *lrwnbs;
    long   *lclpts, *llinks;
    double *lcoeff, *ldiag, *ldsqrt, *pdReg;

    /* IPM arrays. */
    double *x, *s, *y, *z, *w;
    double *deltax, *deltas, *deltay, *deltaz, *deltaw;
    double *xprox, *yprox, *yorgnl, *theta, *xib, *xic, *xiu;

    /* Work arrays. */
    long   *iltmp1, *iltmp2, *iltmp3, *irow;
    INTS   *intmp2, *intmp3;
    INTS   *header, *linkfd, *linkbk;
    double *rmtmp1, *rmtmp2, *rmtmp3, *rmtmp4, *rmtmp5, *relt;
    double *rntmp1, *rntmp2, *rntmp3, *rntmp4, *rntmp5, *rntmp6;
    double *colnrm, *auxwrk;

} Lp;


typedef struct Lp_optimal_solution {
    long   m, m0, n, n_struct;
    char   *rwsave;          /* a copy of row names           */
    double *xsave;           /* primal       optimal solution */
    double *ssave;           /* primal slack optimal solution */
    double *ysave;           /* dual         optimal solution */
    double *zsave;           /* dual   slack optimal solution */
    double *wsave;           /* dual   slack optimal solution */
    double bsave;            /* barrier in the optimum        */
} Lp_save;


int        save_Lp(Lp* L, Lp_save* L_save);

void       read_MPS(long max_rows, long max_cols, long max_nz, long max_nzq,
             char *rwname, char *clname, INTS *rwstat,
             long *col_ptr, INTS *row_ind, double *coeff,
             double *rhs, double *range, double *low_b, double *upp_b,
             long *m, long *n, long *nz, long *nzq,
             long *lpqp, double *objcon, double *q_diag,
             long *qcol_ptr, INTS *qrow_ind, double *q_coeff,
             long *presolve, double *opt_tol, long *irobj, double *mult,
             long *iolog, double *big, char* file_MPS);

Lp*        load_Lp(char *rwname, char *clname, INTS *rwstat, long *col_ptr,
             INTS *row_ind, double *coeff, double *rhs, double *range,
             double *low_b, double *upp_b, long m, long n, long nz, long nzq,
             long lpqp, double objcon,
             double *q_diag, long *qcol_ptr, INTS *qrow_ind, double *q_coeff,
             long presolve, long irobj, long iolog, double mult);

Lp*        read_Lp(long max_rows, long max_cols, long max_nz, long max_nzq,
             long *presolve, double *opt_tol, long *irobj, double *mult,
             long *iolog, double *big, char* file_MPS);

int        solve_Lp(Lp *L, long restart, Lp_save* L_save, double *opt_tol,
             long *iexit, double *errorb, double *erroru, double *errorc,
             double *asmall, double *alarge, long *innner_iterations);

int        acpdm_Lp(Lp *L, long restart, Lp_save* L_save, double *opt_tol,
             long *iexit, double *errorb, double *erroru, double *errorc,
             double *asmall, double *alarge, long *innner_iterations,
             INTS *go_to_zero);

int        wspdm_Lp(Lp *L, long restart, Lp_save* L_save, double *opt_tol,
             long *iexit, double *errorb, double *erroru, double *errorc,
             double *asmall, double *alarge, long *innner_iterations,
             INTS *go_to_zero);

int        wspdm_Qp(Lp *L, long restart, Lp_save* L_save, double *opt_tol,
             long *iexit, double *errorb, double *erroru, double *errorc,
             double *asmall, double *alarge, long *innner_iterations,
             INTS *go_to_zero);

void       preproc_Lp(Lp* L);

void       get_dual_solution(Lp* L, double *dual);

void       get_opt_solution(Lp* L, long* status, double* obj,
             double* primal, double *dual, double* slack, double* red_costs);

void       print_Lp(Lp* L, char sol_file[FNAME+1], char nammps[MPSNAME+1],
             long iwrite, double *primal);

int        write_Lp(Lp* L, char filename[FNAME+1], char* format);

void       sleep_Lp(Lp* L);

void       wakeup_Lp(Lp* L);

void       free_Lp(Lp* L);

void       free_Lp_save(Lp_save* L_save);


#define    EQUAL              1
#define    GREATER_EQUAL      2
#define    LESS_EQUAL         3
#define    OBJECTIVE          4
#define    MINIMIZE           1
#define    MAXIMIZE          -1

/*
The status of the solution after calling load_Lp routine
(exit_code variable):
  0     Erything OK;
  1     Primal INFEASIBLE (or dual UNBOUNDED);
  2     Primal UNBOUNDED (or dual INFEASIBLE);
  31    Primal INFEASIBLE;
  32    Primal UNBOUNDED;
  33    Dual   INFEASIBLE;
  34    Dual   UNBOUNDED;
  35    Exit due to the lack of storage for presolve history;
  >100  Error.
The status of the solution after calling any solver routine
(solve_Lp, acpdm_Lp or wspdm_Lp):
  0     OPTIMAL solution found;
  1     Primal INFEASIBLE (or dual UNBOUNDED);
  2     Primal UNBOUNDED (or dual INFEASIBLE);
  3     SUBOPTIMAL solution found (accuracy problems);
  4     Excess iterations/time limit.
*/
