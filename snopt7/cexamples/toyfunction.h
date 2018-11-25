#include "f2c.h"

void toy0
( integer *inform, char *Prob, integer *neF, integer *n, doublereal *ObjAdd,
  integer *ObjRow, doublereal *xlow, doublereal *xupp,
  doublereal *Flow, doublereal *Fupp, doublereal *x,
  integer *xstate, doublereal *Fmul );

int toyusrf_
( integer    *Status, integer *n,    doublereal x[],
  integer    *needF,  integer *neF,  doublereal F[],
  integer    *needG,  integer *neG,  doublereal G[],
  char       *cu,     integer *lencu,
  integer    iu[],    integer *leniu,
  doublereal ru[],    integer *lenru );

void toy1
( integer *inform, char *Prob, integer *neF, integer *n,
  integer *iAfun, integer *jAvar, integer *lenA, integer *neA, doublereal *A,
  integer *iGfun, integer *jGvar, integer *lenG, integer *neG, doublereal *ObjAdd,
  integer *ObjRow, doublereal *xlow, doublereal *xupp,
  doublereal *Flow, doublereal *Fupp, doublereal *x,
  integer *xstate, doublereal *Fmul );

int toyusrfg_
( integer    *Status, integer *n,    doublereal x[],
  integer    *needF,  integer *neF,  doublereal F[],
  integer    *needG,  integer *neG,  doublereal G[],
  char       *cu,     integer *lencu,
  integer    iu[],    integer *leniu,
  doublereal ru[],    integer *lenru);

