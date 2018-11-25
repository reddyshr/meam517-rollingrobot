#include <stdio.h>
#include <string.h>
#include "f2c.h"


void toy0
( integer *inform, char *Prob, integer *neF, integer *n, doublereal *ObjAdd,
  integer *ObjRow, doublereal *xlow, doublereal *xupp,
  doublereal *Flow, doublereal *Fupp, doublereal *x,
  integer *xstate, doublereal *Fmul )
{
  /*     ================================================================== */
  /*     Toy0   defines input data for the toy problem discussed in the     */
  /*     SnoptA Users Guide.                                                */
  /*                                                                        */
  /*        Minimize                      x(2)                              */
  /*                                                                        */
  /*        subject to   x(1)**2      + 4 x(2)**2  <= 4,                    */
  /*                    (x(1) - 2)**2 +   x(2)**2  <= 5,                    */
  /*                     x(1) >= 0.                                         */
  /*                                                                        */
  /*                                                                        */
  /*     On exit:                                                           */
  /*     inform      is 0 if there is enough storage, 1 otherwise.          */
  /*     neF         is the number of problem functions                     */
  /*                 (objective and constraints, linear and nonlinear).     */
  /*     n           is the number of variables.                            */
  /*     xlow        holds the lower bounds on x.                           */
  /*     xupp        holds the upper bounds on x.                           */
  /*     Flow        holds the lower bounds on F.                           */
  /*     Fupp        holds the upper bounds on F.                           */

  /*     xstate(1:n) is a set of initial states for each x  (0,1,2,3,4,5).  */
  /*     x (1:n)     is a set of initial values for x.                      */
  /*     Fmul(1:neF) is a set of initial values for the dual variables.     */
  /*                                                                        */
  /*     ================================================================== */

  /* Give the problem a name.  */
  /*   An 8-character name is expected so fill out with spaces */
  sprintf(Prob,"%s","Toy0    ");

  /* Assign the dimensions of the constraint Jacobian */

  *neF    = 3;
  *n      = 2;

  *ObjRow = 1; /* NOTE: Me must add one to mesh with fortran */
  *ObjAdd = 0;

  /* Set the upper and lower bounds. */
  xlow[0]   =   0.0;  xlow[1] = -1e6;
  xupp[0]   =   1e6;  xupp[1] =  1e6;
  xstate[0] =   0;    xstate[1] = 0;

  Flow[0] = -1e6; Flow[1] = -1e6; Flow[2] = -1e6;
  Fupp[0] =  1e6; Fupp[1] =  4.0; Fupp[2] =  5.0;

  x[0]    = 1.0;
  x[1]    = 1.0;

}

int toyusrf_
( integer    *Status, integer *n,    doublereal x[],
  integer    *needF,  integer *neF,  doublereal F[],
  integer    *needG,  integer *neG,  doublereal G[],
  char       *cu,     integer *lencu,
  integer    iu[],    integer *leniu,
  doublereal ru[],    integer *lenru )
{
  /*     ================================================================== */
  /*     Computes the nonlinear objective and constraint terms for the toy  */
  /*     problem featured in the SnoptA users guide.                        */
  /*     neF = 3, n = 2.                                                    */
  /*                                                                        */
  /*        Minimize                      x(2)                              */
  /*                                                                        */
  /*        subject to   x(1)**2      + 4 x(2)**2  <= 4,                    */
  /*                    (x(1) - 2)**2 +   x(2)**2  <= 5,                    */
  /*                     x(1) >= 0.                                         */
  /*                                                                        */
  /*     ================================================================== */

  F[0] = x[1];
  F[1] = x[0]*x[0] + 4*x[1]*x[1];
  F[2] = (x[0] - 2)*(x[0] - 2) + x[1]*x[1];
  return 0;
}

void toy1
( integer *inform, char *Prob, integer *neF, integer *n,
  integer *iAfun, integer *jAvar, integer *lenA, integer *neA, doublereal *A,
  integer *iGfun, integer *jGvar, integer *lenG, integer *neG, doublereal *ObjAdd,
  integer *ObjRow, doublereal *xlow, doublereal *xupp,
  doublereal *Flow, doublereal *Fupp, doublereal *x,
  integer *xstate, doublereal *Fmul)
{
  /*     ================================================================== */
  /*     Toy1   defines input data for the toy problem discussed in the     */
  /*     SnoptA Users Guide.                                                */
  /*                                                                        */
  /*        Minimize                      x(2)                              */
  /*                                                                        */
  /*        subject to   x(1)**2      + 4 x(2)**2  <= 4,                    */
  /*                    (x(1) - 2)**2 +   x(2)**2  <= 5,                    */
  /*                     x(1) >= 0.                                         */
  /*                                                                        */
  /*                                                                        */
  /*     On exit:                                                           */
  /*        neF  is the number of objective and constraint functions        */
  /*               (including linear and nonlinear)                         */
  /*        n    is the number of variables.                                */
  /*                                                                        */
  /*                                                                        */
  /*        (iGfun(k),jGvar(k)), k = 1,2,...,neG, define the coordinates    */
  /*             of the nonzero problem derivatives.                        */
  /*             If (iGfun(k),jGvar(k)) = (i,j), G(k) is the ijth element   */
  /*             of the problem vector F(i), i = 0,1,2,...,neF,  with       */
  /*             objective function in position 0 and constraint functions  */
  /*             in positions  1  through  m.                               */
  /*                                                                        */
  /*        (iAfun(k),jAvar(k),a(k)), k = 1,2,...,neA, are the coordinates  */
  /*             of the nonzero constant problem derivatives.               */
  /*                                                                        */
  /*             To keep things simple, no constant elements are set here.  */
  /*                                                                        */
  /*     ================================================================== */
  /* Give the problem a name.  */

  strcpy(Prob,"Toy1");

  /* Assign the dimensions of the constraint Jacobian */

  *neF    = 3;
  *n      = 2;

  *ObjRow = 1; /* NOTE: Me must add one to mesh with fortran */
  *ObjAdd = 0;

  /* Set the upper and lower bounds. */
  xlow[0]   =   0.0;  xlow[1] = -1e6;
  xupp[0]   =   1e6;  xupp[1] =  1e6;
  xstate[0] =   0;    xstate[1] = 0;

  Flow[0] = -1e6; Flow[1] = -1e6; Flow[2] = -1e6;
  Fupp[0] =  1e6; Fupp[1] =  4.0; Fupp[2] =  5.0;
  Fmul[0] =    0; Fmul[1] =    0; Fmul[2] =    0;

  x[0]    = 1.0;
  x[1]    = 1.0;

  *inform = 0;
  *neG    = 0;
  iGfun[*neG] = 1;
  jGvar[*neG] = 1;

  *neG    = *neG + 1;
  iGfun[*neG] = 1;
  jGvar[*neG] = 2;

  *neG    = *neG + 1;
  iGfun[*neG] = 2;
  jGvar[*neG] = 1;

  *neG    = *neG + 1;
  iGfun[*neG] = 2;
  jGvar[*neG] = 2;

  *neG    = *neG + 1;
  iGfun[*neG] = 3;
  jGvar[*neG] = 1;

  *neG    = *neG + 1;
  iGfun[*neG] = 3;
  jGvar[*neG] = 2;

  *neG = *neG + 1;
  /* *neG = 6 */

  *neA = 0;
}


int toyusrfg_
( integer    *Status, integer *n,    doublereal x[],
  integer    *needF,  integer *neF,  doublereal F[],
  integer    *needG,  integer *neG,  doublereal G[],
  char       *cu,     integer *lencu,
  integer    iu[],    integer *leniu,
  doublereal ru[],    integer *lenru )
{
  /*     ==================================================================  */
  /*     Computes the nonlinear objective and constraint terms for the toy   */
  /*     problem featured in the SnoptA users guide.                         */
  /*     neF = 3, n = 2.                                                     */
  /*                                                                         */
  /*        Minimize                      x(2)                               */
  /*                                                                         */
  /*        subject to   x(1)**2      + 4 x(2)**2  <= 4,                     */
  /*                    (x(1) - 2)**2 +   x(2)**2  <= 5,                     */
  /*                     x(1) >= 0.                                          */
  /*                                                                         */
  /*     The triples (g(k),iGfun(k),jGvar(k)), k = 1,2,...,neG, define       */
  /*     the sparsity pattern and values of the nonlinear elements           */
  /*     of the Jacobian.                                                    */
  /*     ==================================================================  */

  if( *needF > 0 ) {
    F[0] = x[1]; /* ! The objective row */
    F[1] = x[0]*x[0] + 4*x[1]*x[1];
    F[2] = (x[0] - 2)*(x[0] - 2) + x[1]*x[1];
  }

  if( *needG > 0 ){
    *neG = 0;
    /* iGfun[*neG] = 1 */
    /* jGvar[*neG] = 1 */
    G[*neG] = 0;

    /* iGfun[*neG] = 1 */
    /* jGvar[*neG] = 2 */
    *neG = *neG + 1;
    G[*neG] = 1.0;

    /* iGfun[*neG] = 2 */
    /* jGvar[*neG] = 1 */
    *neG = *neG + 1;
    G[*neG] = 2*x[0];

    /* iGfun[*neG] = 2 */
    /* jGvar[*neG] = 2 */
    *neG = *neG + 1;
    G[*neG] = 8*x[1];

    /* iGfun[*neG] = 3 */
    /* jGvar[*neG] = 1 */
    *neG = *neG + 1;
    G[*neG] = 2*(x[0] - 2);

    /* iGfun[*neG] = 3 */
    /* jGvar[*neG] = 2 */
    *neG = *neG + 1;
    G[*neG] = 2*x[1];
    *neG = *neG + 1;
  }
  return 0;
}
