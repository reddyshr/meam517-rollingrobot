#include <stdio.h>
#include <string.h>
#include "f2c.h"
#include "snopt.h"
#include "toyfunction.h"

int main()
{
  integer    minrw, miniw, mincw;
  integer    lenrw = 20000, leniw = 10000, lencw = 500;
  doublereal rw[20000];
  integer    iw[10000];
  char       cw[8*500];

  integer    Cold = 0, Basis = 1, Warm = 2;

  doublereal x[2], xlow[2], xupp[2], xmul[2];
  doublereal F[3], Flow[3], Fupp[3], Fmul[3];
  doublereal ObjAdd;
  integer    xstate[2], Fstate[3];

  integer    INFO, ObjRow;
  integer    n, neF;

  integer    lenA = 10, iAfun[10], jAvar[10];
  doublereal A[10];
  integer    lenG = 10, iGfun[10], jGvar[10];

  integer    neA, neG;
  integer    nxname = 1, nFname = 1, npname;
  char       xnames[1*8], Fnames[1*8];
  char       Prob[200];

  integer    iSpecs = 4,  spec_len;
  integer    iSumm  = 6;
  integer    iPrint = 9,  prnt_len;

  char       printname[200];
  char       specname[200];

  integer    nS, nInf;
  doublereal sInf;
  integer    DerOpt, Major, iSum, iPrt, strOpt_len;
  char       strOpt[200];

  printf("\nSolving toy0 without first derivatives ...\n");

  /* open output files using snfilewrappers.[ch] */
  sprintf(specname ,   "%s", "sntoya.spc");   spec_len = strlen(specname);
  sprintf(printname,   "%s", "sntoya.out");   prnt_len = strlen(printname);

  /* Open the print file, fortran style */
  snopenappend_
    ( &iPrint, printname,   &INFO, prnt_len );

  /*     ================================================================== */
  /*     First,  sninit_ MUST be called to initialize optional parameters   */
  /*     to their default values.                                           */
  /*     ================================================================== */

  sninit_
    ( &iPrint, &iSumm, cw, &lencw, iw, &leniw, rw, &lenrw, 8*500 );

  /*     Set up the problem to be solved.                       */
  /*     No derivatives are set in this case.                   */
  /*     NOTE: To mesh with Fortran style coding,               */
  /*           it ObjRow must be treated as if array F          */
  /*           started at 1, not 0.  Hence, if F(0) = objective */
  /*           then ObjRow should be set to 1.                  */

  toy0
    ( &INFO, Prob, &neF, &n, &ObjAdd, &ObjRow, xlow, xupp,
      Flow, Fupp, x, xstate, Fmul );
  npname = strlen(Prob);

  /*     SnoptA will compute the Jacobian by finite-differences.   */
  /*     The user has the option of calling  snJac  to define the  */
  /*     coordinate arrays (iAfun,jAvar,A) and (iGfun, jGvar).     */

  snjac_
    ( &INFO, &neF, &n, toyusrf_,
      iAfun, jAvar, &lenA, &neA, A,
      iGfun, jGvar, &lenG, &neG,
      x, xlow, xupp, &mincw, &miniw, &minrw,
      cw, &lencw, iw, &leniw, rw, &lenrw,
      cw, &lencw, iw, &leniw, rw, &lenrw,
      8*500, 8*500 );

  /*     ------------------------------------------------------------------ */
  /*     Warn SnoptA that userf does not compute derivatives.               */
  /*     The parameters iPrt and iSum may refer to the Print and Summary    */
  /*     file respectively.  Setting them to 0 suppresses printing.         */
  /*     ------------------------------------------------------------------ */

  DerOpt = 0;
  iPrt   = 0;
  iSum   = 0;
  sprintf(strOpt,"%s","Derivative option");
  strOpt_len = strlen(strOpt);
  snseti_
    ( strOpt, &DerOpt, &iPrt, &iSum, &INFO,
      cw, &lencw, iw, &leniw, rw, &lenrw, strOpt_len, 8*500 );

  /*     ------------------------------------------------------------------ */
  /*     Go for it, using a Cold start.                                     */
  /*     ------------------------------------------------------------------ */

  snopta_
    ( &Cold, &neF, &n, &nxname, &nFname,
      &ObjAdd, &ObjRow, Prob, toyusrf_,
      iAfun, jAvar, &lenA, &neA, A,
      iGfun, jGvar, &lenG, &neG,
      xlow, xupp, xnames, Flow, Fupp, Fnames,
      x, xstate, xmul, F, Fstate, Fmul,
      &INFO, &mincw, &miniw, &minrw,
      &nS, &nInf, &sInf,
      cw, &lencw, iw, &leniw, rw, &lenrw,
      cw, &lencw, iw, &leniw, rw, &lenrw,
      npname, 8*nxname, 8*nFname,
      8*500, 8*500);


  printf("\nSolving toy1 using first derivatives ...\n");

  toy1
    ( &INFO, Prob, &neF, &n,
      iAfun, jAvar, &lenA, &neA, A,
      iGfun, jGvar, &lenG, &neG,
      &ObjAdd, &ObjRow, xlow, xupp,
      Flow, Fupp, x, xstate, Fmul );

  /* Read in specs file (optional) */
  /* snfilewrapper_ will open the specs file, fortran style, */
  /* then call snspec_ to read in specs.                        */

  snfilewrapper_
    ( specname, &iSpecs, &INFO, cw, &lencw,
      iw, &leniw, rw, &lenrw, spec_len, 8*lencw);

  if( INFO != 101 )
    {
      printf("Warning: trouble reading specs file %s \n", specname);
    }

  /* Specify any user options not set in the Specs file. */
  DerOpt = 1;
  snseti_
    ( strOpt, &DerOpt, &iPrint, &iSumm, &INFO,
      cw, &lencw, iw, &leniw, rw, &lenrw, strOpt_len, 8*500 );

  Major = 250;
  strcpy( strOpt,"Major Iteration limit");
  strOpt_len = strlen(strOpt);
  snseti_
    ( strOpt, &Major, &iPrint, &iSumm, &INFO,
      cw, &lencw, iw, &leniw, rw, &lenrw, strOpt_len, 8*500 );

  /*     ------------------------------------------------------------------ */
  /*     Solve the problem again, this time with derivatives specified.     */
  /*     ------------------------------------------------------------------ */
  snopta_
    ( &Cold, &neF, &n, &nxname, &nFname,
      &ObjAdd, &ObjRow, Prob, toyusrfg_,
      iAfun, jAvar, &lenA, &neA, A,
      iGfun, jGvar, &lenG, &neG,
      xlow, xupp, xnames, Flow, Fupp, Fnames,
      x, xstate, xmul, F, Fstate, Fmul,
      &INFO, &mincw, &miniw, &minrw,
      &nS, &nInf, &sInf,
      cw, &lencw, iw, &leniw, rw, &lenrw,
      cw, &lencw, iw, &leniw, rw, &lenrw,
      npname, 8*nxname, 8*nFname,
      8*500, 8*500);

  snclose_( &iPrint );
  snclose_( &iSpecs );
}
