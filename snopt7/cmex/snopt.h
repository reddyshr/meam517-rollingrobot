/* Josh Griffin ... modeled after npsol.h written by:          */
/* Mike Gertz - 2-Aug-98                                       */
/* Function prototypes for functions in the snOpt distribution */

#ifndef SNOPT
#define SNOPT

#pragma once

#ifndef F2C_INCLUDE
#include "f2c.h"
#endif

extern int
snkera_ ( integer *start, integer *nef, integer *n,
          integer *nxname, integer *nfname, doublereal *objadd,
          integer *objrow, char *prob, U_fp usrfun,
          U_fp snLog, U_fp snLog2, U_fp sqLog, U_fp snAbort,
          integer *iafun, integer *javar, integer *lena, integer *nea, doublereal *a,
          integer *igfun, integer *jgvar, integer *leng,
          integer *neg, doublereal *xlow, doublereal *xupp,
          char *xnames, doublereal *flow, doublereal *fupp,
          char *fnames, doublereal *x, integer *xstate, doublereal *xmul,
          doublereal *f, integer *fstate, doublereal *fmul, integer *inform,
          integer *mincw, integer *miniw, integer *minrw, integer *ns, integer *ninf, doublereal *sinf,
          char *cu, integer *lencu, integer *iu, integer *leniu, doublereal *ru, integer *lenru,
          char *cw, integer *lencw, integer *iw, integer *leniw, doublereal *rw, integer *lenrw,
          ftnlen prob_len, ftnlen xnames_len, ftnlen fnames_len, ftnlen cu_len, ftnlen cw_len );

extern int
sninit_ ( integer *iPrint, integer *iSumm,
          char *cw, integer *lencw, integer *iw, integer *leniw, doublereal *rw, integer *lenrw,
          ftnlen cw_len );

extern int
snset_  ( char *buffer,                     integer *iprint, integer *isumm,
          integer *errors,
          char *cw, integer *lencw, integer *iw, integer *leniw, doublereal *rw, integer *lenrw,
          ftnlen buffer_len, ftnlen cw_len);

extern int
snseti_(  char *buffer, integer *ivalue,    integer *iprint, integer *isumm,
          integer *errors,
          char *cw, integer *lencw, integer *iw, integer *leniw, doublereal *rw, integer *lenrw,
          ftnlen buffer_len, ftnlen cw_len);

extern int
snsetr_(  char *buffer, doublereal *rvalue, integer *iprint, integer *isumm,
          integer *errors,
          char *cw, integer *lencw, integer *iw, integer *leniw, doublereal *rw, integer *lenrw,
          ftnlen buffer_len, ftnlen cw_len);

extern integer
snget_ ( char *buffer,
         integer *errors,
         char *cw, integer *lencw, integer *iw, integer *leniw, doublereal *rw, integer *lenrw,
         ftnlen buffer_len, ftnlen cw_len);

extern int
sngeti_( char *buffer, integer *ivalue,
         integer *errors,
         char *cw, integer *lencw, integer *iw, integer *leniw, doublereal *rw, integer *lenrw,
         ftnlen buffer_len, ftnlen cw_len);

extern int
sngetr_( char *buffer, doublereal *rvalue,
         integer *errors,
         char *cw, integer *lencw, integer *iw, integer *leniw, doublereal *rw, integer *lenrw,
         ftnlen buffer_len, ftnlen cw_len);

extern int
snspec_( integer *ispecs, integer *inform,
         char *cw, integer *lencw, integer *iw, integer *leniw, doublereal *rw, integer *lenrw,
         ftnlen cw_len);

extern int
snmema0_(integer *inform, integer *iPrint, integer *iSumm,
         integer *nef, integer *n, integer *nxname, integer *nfname, integer *nea, integer *neg,
         integer *mincw, integer *miniw, integer *minrw,
         char *cw, integer *lencw, integer *iw, integer *leniw, doublereal *rw, integer *lenrw,
         ftnlen cw_len);

extern int
snjac_(  integer *inform,
         integer *nef, integer *n, U_fp userfg,
         integer *iafun, integer *javar, integer *lena,
         integer *nea, doublereal *a, integer *igfun,
         integer *jgvar, integer *leng, integer *neg,
         doublereal *x, doublereal *xlow, doublereal *xupp,
         integer *mincw, integer *miniw, integer *minrw,
         char *cu, integer *lencu, integer *iu, integer *leniu, doublereal *ru, integer *lenru,
         char *cw, integer *lencw, integer *iw, integer *leniw, doublereal *rw, integer *lenrw,
         ftnlen cu_len, ftnlen cw_len);
#endif
