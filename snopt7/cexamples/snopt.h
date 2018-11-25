/* Josh Griffin ... modeled after npsol.h written by:          */
/* Mike Gertz - 2-Aug-98                                       */
/* Function prototypes for functions in the snopt distribution */

#ifndef SNOPT
#define SNOPT

#pragma once

#ifndef F2C_INCLUDE
#include "f2c.h"
#endif

extern int snopta_
( integer *start, integer *nef, integer *n,
  integer *nxname, integer *nfname, doublereal *objadd, integer *objrow,
  char *prob, U_fp usrfun, integer *iafun, integer *javar,
  integer *lena, integer *nea, doublereal *a, integer *igfun,
  integer *jgvar, integer *leng, integer *neg,
  doublereal *xlow, doublereal *xupp,
  char *xnames, doublereal *flow, doublereal *fupp, char *fnames,
  doublereal *x, integer *xstate, doublereal *xmul, doublereal *f,
  integer *fstate, doublereal *fmul, integer *inform, integer *mincw,
  integer *miniw, integer *minrw, integer *ns, integer *ninf,
  doublereal *sinf, char *cu, integer *lencu, integer *iu, integer *leniu,
  doublereal *ru, integer *lenru, char *cw, integer *lencw,
  integer *iw, integer *leniw, doublereal *rw, integer *lenrw, ftnlen
  prob_len, ftnlen xnames_len, ftnlen fnames_len, ftnlen cu_len, ftnlen cw_len );

extern int sninit_
( integer *iPrint, integer *iSumm, char *cw,
  integer *lencw, integer *iw, integer *leniw,
  doublereal *rw, integer *lenrw, ftnlen cw_len );

extern int sngeti_
( char *buffer, integer *ivalue, integer *inform,
  char *cw, integer *lencw, integer *iw,
  integer *leniw, doublereal *rw, integer *lenrw,
  ftnlen buffer_len, ftnlen cw_len);

extern int snset_
( char *buffer, integer *iprint, integer *isumm,
  integer *inform, char *cw, integer *lencw,
  integer *iw, integer *leniw,
  doublereal *rw, integer *lenrw,
  ftnlen buffer_len, ftnlen cw_len );

extern int snseti_
( char *buffer, integer *ivalue, integer *iprint,
  integer *isumm, integer *inform, char *cw,
  integer *lencw, integer *iw, integer *leniw,
  doublereal *rw, integer *lenrw, ftnlen buffer_len,
  ftnlen cw_len );

extern int snsetr_
( char *buffer, doublereal *rvalue, integer * iprint,
  integer *isumm, integer *inform, char *cw,
  integer *lencw, integer *iw, integer *leniw,
  doublereal *rw, integer *lenrw, ftnlen buffer_len,
  ftnlen cw_len );

extern int snspec_
( integer *ispecs, integer *inform, char *cw,
  integer *lencw, integer *iw, integer *leniw,
  doublereal *rw, integer *lenrw, ftnlen cw_len);

extern int snmema_
( integer *iExit, integer *nef, integer *n, integer *nxname,
  integer *nfname, integer *nea, integer *neg,
  integer *mincw, integer *miniw, integer *minrw,
  char *cw, integer *lencw, integer *iw,
  integer *leniw, doublereal *rw, integer *lenrw,
  ftnlen cw_len);

extern int snjac_
( integer *iExit, integer *nef, integer *n, U_fp userfg,
  integer *iafun, integer *javar, integer *lena, integer *nea, doublereal *a,
  integer *igfun, integer *jgvar, integer *leng, integer *neg,
  doublereal *x, doublereal *xlow, doublereal *xupp,
  integer *mincw, integer *miniw, integer *minrw,
  char *cu, integer *lencu, integer *iu, integer *leniu, doublereal *ru, integer *lenru,
  char *cw, integer *lencw, integer *iw, integer *leniw, doublereal *rw, integer *lenrw,
  ftnlen cu_len, ftnlen cw_len);

extern int snopenappend_
( integer *iunit, char *name, integer *inform, ftnlen name_len);

extern int snfilewrapper_
( char *name__, integer *ispec, integer *
  inform__, char *cw, integer *lencw, integer *iw,
  integer *leniw, doublereal *rw, integer *lenrw,
  ftnlen name_len, ftnlen cw_len);

extern int snclose_
( integer *iunit);

extern int snopen_
( integer *iunit);

#endif
