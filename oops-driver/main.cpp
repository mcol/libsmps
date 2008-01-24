/*
 *  main.cpp
 *
 *  Driver for the SMPS interface to Oops.
 *
 *  Marco Colombo
 *  School of Mathematics
 *  University of Edinburgh
 *
 */

#include <iostream>
#include "SmpsOops.h"
#include "oops/oopstime.h"
#include "oops/oopsmem.h"
#include "oops/WriteMps.h"


FILE *printout;
double tt_start, tt_end;

static primal_dual_pb* defineProblem(SmpsReturn *Pb);
static void setupOutputFile(OptionsOops &opt);


/** Driver routine for the SMPS interface to OOPS */
int main(const int argc, const char *argv[]) {

  // initialise the parallel stuff
  InitLippPar(argc, argv);

  // parse the command line options
  OptionsOops opt(argc, argv);
  int rv = opt.parse();
  if (rv)
    return 1;

#ifdef WITH_TIME
  initTimes();
#endif

  // decide where the output is going to appear
  setupOutputFile(opt);

  // create an object for the problem
  SmpsOops data(opt.getSmpsFile(), opt.cutoffLevel());

  // read the smps files
  rv = data.read();
  if (rv)
    return 1;

  HopdmOptions *hopdm_options = NewHopdmOptions();

  // solve the problem
  data.solve(opt, hopdm_options);

#ifdef WITH_TIME
  reportTimes(stdout);
#endif

  // close the output file
  if (opt.outputToFile())
    fclose(printout);

  // clean up
  FreeHopdmOptions(hopdm_options);

#ifdef WITH_MPI
  MPI_Finalize();
#endif

  return 0;
}

/**
 *  Call OOPS to solve the problem.
 *
 *  @param opt:
 *         Command line options.
 *  @param hopdm_options:
 *         Pointer to the options for the solver.
 *  @return 1 If something goes wrong; 0 otherwise.
 */
int SmpsOops::solve(const OptionsOops &opt, HopdmOptions *hopdm_options) {

  // generate the problem
  SmpsReturn *prob = generateSmps();

  // check the validity of the arguments
  if (!prob) {
    cout << "No problem defined." << endl;
    return 1;
  }

  // setup the primal-dual problem
  primal_dual_pb *pdProb = defineProblem(prob);

#ifdef REP_MEM
  reportMem();
#endif

  // write the deterministic equivalent in mps format
  if (opt.writeMps()) {
    FILE *fout = fopen("smps.mps", "w");
    Write_MpsFile(fout, pdProb->AlgAug, pdProb->b, pdProb->c, pdProb->u, 0,
		  prob->colnames, prob->rownames);
    fclose(fout);
  }

  printf(" --------------- oopsSolve -----------------\n");

  // options for the complete problem
  hopdm_options->glopt->conv_tol = 1.e-4;

  hopdm_ret *ret;
  PrintOptions *Prt = NewHopdmPrt(PRINT_ITER);

  if (opt.dontSolve()) {
    printf("Problem not solved by request.\n");
    goto TERMINATE;
  }

  // start the clock
  tt_start = oopstime();

  // solve the problem
  ret = hopdm(printout, pdProb, hopdm_options, Prt);

  // stop the clock
  tt_end = oopstime();

  // report statistics
  fprintf(printout, "Elapsed time: %.10g seconds.\n", tt_end - tt_start);

  if (opt.printSolution())
    getSolution(pdProb, prob);

  // clean up
  free(ret);

 TERMINATE:

  // clean up
  FreePDProblem(pdProb);
  freeSmpsReturn(prob);
  delete Prt;

  return 0;
}

/**
 *  Set up the OOPS Algebras and Vectors and build the primal-dual problem.
 *
 *  @param Pb:
 *         Pointer to an SmpsReturn structure.
 *  @return Pointer to a primal_dual_pb structure containing the problem
 *          to be solved by OOPS.
 *
 *  @note
 *  InitAlgebrasNew assumes that a callback function is set up for all
 *  SparseMatrix leaves of A and Q.
 */
primal_dual_pb* defineProblem(SmpsReturn *Pb) {

  Algebra *A = Pb->AlgA;
  Algebra *Q = Pb->AlgQ;

  printf(" --------------- InitAlgebras --------------\n");
  Algebra *AlgAug = InitAlgebrasNew(A, Q);

  Vector *vb = NewVector(A->Trow, "vb");
  Vector *vc = NewVector(A->Tcol, "vc");
  Vector *vu = NewVector(A->Tcol, "vu");

  Vector *vx = NewVector(A->Tcol, "vx");
  Vector *vy = NewVector(A->Trow, "vy");
  Vector *vz = NewVector(A->Tcol, "vz");

  CopyDenseToVector(Pb->b, vb);
  CopyDenseToVector(Pb->c, vc);
  CopyDenseToVector(Pb->u, vu);

  return NewPDProblem(AlgAug, vb, vc, vu, vx, vy, vz);
}

/**
 *  Decide where the output is going to appear.
 *
 *  Sets the global variable @c printout to decide where the output from
 *  OOPS is going to appear.
 *
 *  @param opt:
 *         Command line options.
 */
void setupOutputFile(OptionsOops &opt) {

  char filename[20];

  // redirect the output to a file
  if (opt.outputToFile()) {

#ifdef WITH_MPI
    sprintf(filename, "output%d.dat", MYID_PAR);
#else
    sprintf(filename, "output.dat");
#endif

    // open the file for writing
    printout = fopen(filename, "w");

#ifdef WITH_MPI
    fprintf(printout, "Output from processor %d.\n\n", MYID_PAR);
#endif
  }

  // print the output on the screen
  else
    printout = stdout;
}

/**
 *  Retrieve the solution.
 *
 *  @param Prob:
 *         Pointer to the problem solved by OOPS.
 *  @param Ret:
 *         Pointer to an SmpsReturn structure.
 *
 *  @bug
 *  The value of the slacks is always zero, as the slacks would need to be
 *  computed here.
 */
int SmpsOops::getSolution(primal_dual_pb *Prob, SmpsReturn *Ret) {

  DenseVector *x, *y, *z, *r;
  Tree *Trow = Ret->AlgA->Trow;
  Tree *Tcol = Ret->AlgA->Tcol;
  int nRows  = Trow->end - Trow->begin;
  int nCols  = Tcol->end - Tcol->begin;

  // allocate space for the solution vectors
  x = NewDenseVector(nCols, "dx");
  y = NewDenseVector(nRows, "dy");
  z = NewDenseVector(nCols, "dz");

  // this is just a placeholder, we cannot yet compute the slacks
  r = NewDenseVector(nRows, "slacks");

  // recover initial order on solution vectors
  SmpsVectorToDense(Prob->x, x, Ret, ORDER_COL);
  SmpsVectorToDense(Prob->y, y, Ret, ORDER_ROW);
  SmpsVectorToDense(Prob->z, z, Ret, ORDER_COL);

  // print the solution
  const NodeInfo *info = getNodeInfo();

#ifdef WITH_MPI
  if (IS_ROOT_PAR)
#endif
    printSolution(info, x->elts, y->elts, r->elts, z->elts);

  // clean up
  delete[] info->nRowsNode;
  delete[] info->nColsNode;
  delete info;

  FreeDenseVector(x);
  FreeDenseVector(y);
  FreeDenseVector(z);
  FreeDenseVector(r);

  return 0;
}
