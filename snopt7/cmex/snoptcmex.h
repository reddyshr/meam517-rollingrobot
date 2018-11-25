#ifndef SNOPTCMEX
#define SNOPTCMEX

#pragma once

#ifndef mex_h
#include "mex.h"
#endif

#ifndef F2C_INCLUDE
#include "f2c.h"
#endif

/* sncmxInit() is called on first entry to the mex gateway.   */
/* It sets all optional parameters as "undefined" and         */
/* hence must be called prior to setting SNOPT options.       */
/* Must NEVER be called more than once.                       */
void
sncmxInit ();

/* sncmxAlloc() allocates workspace once the problem dimensions are known. */
int sncmxAlloc  ( integer neF,    integer n,     integer nxname,
		  integer nfname, integer neA,   integer neG );
int sncmxReAlloc( integer mincw,  integer miniw, integer minrw);

/* Functions for setting and getting options */
integer  sncmxGet  ( char option[] );
integer  sncmxGetc ( char option[], char       *cvalue[] );
integer  sncmxGeti ( char option[], integer    *ivalue   );
integer  sncmxGetr ( char option[], doublereal *rvalue   );
integer  sncmxSet  ( char option[] );
integer  sncmxSeti ( char option[], integer     ivalue   );
integer  sncmxSetr ( char option[], doublereal  rvalue   );
integer  sncmxSpecs( char specsfile[] );

void     sncmxJac  ( int nlhs, mxArray *plhs[], int nrhs,
		     const mxArray *prhs[] );
void     sncmxSnopt( int nlhs, mxArray *plhs[], int nrhs,
		     const mxArray *prhs[], integer advancedCallbacks );
void
sncmxFreeWorkSpace ( );

integer
sncmxGetksnoptValue ( char *ksnoptname );

void
sncmxSingleStringArg( int nlhs, mxArray *plhs[], int nrhs,
		      const mxArray *prhs[], int what );

extern FILE   *snPrintFile, *snSummaryFile, *snSpecsFile;

extern integer snSpecsUnit, snPrintUnit, snSummaryUnit;

extern int     callType;
extern int     printFileIsOpen, summaryFileIsOpen, screenIsOn;

extern int     snlog_  ( );
extern int     snlog2_ ( );
extern int     sqlog_  ( );
extern int     snabort_( );
#endif
