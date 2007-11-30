/*
 *  main.c
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
double tt_end;

static primal_dual_pb* defineProblem(SmpsReturn *Pb);
static void setupOutputFile(opt_st *options);

/** Driver routine for the SMPS interface to OOPS */
int main(const int argc, const char *argv[]) {

  // initialise the parallel stuff
  InitLippPar(argc, argv);

  // parse the command line options
  opt_st *options = parseOptions(argc, argv);
  if (!options)
    return 1;

#ifdef WITH_TIME
  initTimes();
#endif

  // decide where the output is going to appear
  setupOutputFile(options);

  // create an object for the problem
  SmpsOops data(options->smpsFile);

  // read the smps files
  int rv = data.read();
  if (rv) {
    freeOptions(options);
    return 1;
  }

  HopdmOptions *hopdm_options = NewHopdmOptions();

  // solve the problem
  data.solve(options, hopdm_options);

#ifdef WITH_TIME
  reportTimes(stdout);
#endif

  // close the output file
  if (options->OutputToFile == 1)
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
 *  @param options:
 *         Pointer to the options structure.
 *  @param hopdm_options:
 *         Pointer to the options for the solver.
 *  @return 1 If something goes wrong; 0 otherwise.
 */
int SmpsOops::solve(const opt_st *options,
		    HopdmOptions *hopdm_options) {

  // generate the problem
  SmpsReturn *prob = generateSmps();

  // check the validity of the arguments
  if (!prob) {
    cout << "oopsSolve: No problem defined." << endl;
    return 1;
  }

  // setup the primal-dual problem
  primal_dual_pb *pdProb = defineProblem(prob);

#ifdef REP_MEM
  reportMem();
#endif

  // write the deterministic equivalent in mps format
  if (options->WriteMps) {
    FILE *fout = fopen("smps.mps", "w");
    Write_MpsFile(fout, pdProb->AlgAug, pdProb->b, pdProb->c, pdProb->u, 0,
		  prob->colnames, prob->rownames);
    fclose(fout);
  }

  printf(" --------------- oopsSolve -----------------\n");

  // options for the complete problem
  hopdm_options->glopt->conv_tol = 1.e-4;

  if (options->DontSolve) {
    printf("Problem not solved by request.\n");
    goto TERMINATE;
  }

  hopdm_ret *ret;
  hopdm_prt_type *Prt = NewHopdmPrt(1);

  // start the clock
  tt_start = oopstime();

  // solve the problem
  ret = hopdm(printout, pdProb, hopdm_options, Prt);

  // stop the clock
  tt_end = oopstime();

  // report statistics
  fprintf(printout, "Elapsed time: %.10g seconds.\n", tt_end - tt_start);

  if (options->PrintSolution)
    getSolution(pdProb, prob);

  // clean up
  delete Prt;
  free(ret);

 TERMINATE:

  // clean up
  FreePDProblem(pdProb);
  freeSmpsReturn(prob);

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
 *  @param  options: Pointer to an opt_st structure.
 */
void setupOutputFile(opt_st *options) {

  char filename[20];

  // redirect the output to a file
  if (options->OutputToFile == 1) {

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
