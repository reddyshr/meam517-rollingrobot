#include "mex.h"
#include "mexUtility.h"
#include "string.h"
#include "snopt.h"
#include "snoptcmex.h"
#include "snoptprint.h"
#include "f2c.h"
#include "snoptprobdef.h"
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

/* Initialize some global constants. */
enum { ksnoptValue = 0,        /*  0 */
       ksnopt,                 /*  1 */
       ksnoptSetOption,        /*  2 */
       ksnoptSetOptionI,       /*  3 */
       ksnoptSetOptionR,       /*  4 */
       ksnoptGetOption,        /*  5 */
       ksnoptGetOptionC,       /*  6 */
       ksnoptGetOptionI,       /*  7 */
       ksnoptGetOptionR,       /*  8 */
       ksnoptSpecsFile,        /*  9 */
       ksnoptSetPrintFile,     /* 10 */
       ksnoptSetSummaryFile,   /* 11 */
       ksnoptClosePrintFile,   /* 12 */
       ksnoptCloseSummaryFile, /* 13 */
       ksnoptAdvanced,         /* 14 */
       ksnoptScreenOn,         /* 15 */
       ksnoptScreenOff,        /* 16 */
       ksnoptJac,              /* 17 */
       ksnoptSetStatus,        /* 18 */
       ksnoptGetStatus };      /* 19 */

int     printFileIsOpen = 0, summaryFileIsOpen = 0, screenIsOn  = 0;
int     snIsInitialized = 0;
int     memoryIsSet     = 0;
int     systemCall      = 0, userCall          = 1, callType    = 0;

FILE   *snSummaryFile   = NULL;
FILE   *snPrintFile     = NULL;
FILE   *snSpecsFile     = NULL;

integer snSpecsUnit     =  4;
integer snPrintUnit     =  9;
integer snSummaryUnit   = 55;

integer lenrw = 500; /* Default values for options. */
integer leniw = 500;
integer lencw = 500;

doublereal *rusnopt;
doublereal *rwsnopt;
integer    *iwsnopt;
char       *cwsnopt;

integer lrusnopt = 0; /* For snopt/fmincon interface. */

integer lrwsnopt = 0; /* Default values for options. */
integer liwsnopt = 0;
integer lcwsnopt = 0;

doublereal  rw[500];
integer     iw[500];
char        cw[8*500];

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs,
		  const mxArray *prhs[] )
     /* mexFunction - MATLAB gateway routine.                           */
{
  enum { ksimpleCallbacks = 0, kadvancedCallbacks = 1 };

  integer    what, gstatus; /* input argument 1 */
  doublereal *dummyDoublePtr;
  static int firstCall = 1;

  if ( firstCall ) {
    /* Initialize memory for storing the options and */
    /* call snInit_().                               */
    sncmxInit();
    /* Local defaults: 1. Output to screen      */
    /*                 2. No timing information */
    sncmxSeti( "Timing level ",     -1 );
    sncmxSet ( "Sticky parameters yes" );
    firstCall  = 0;
  }

  if ( nrhs < 1 ) {
    mexErrMsgTxt( "Must have at least 1 input argument" );
  }
  /* what :input argument 1 */
  assertScalar( prhs[0], "what" );
  what = *mxGetPr( prhs[0] );

  switch(what) {
  case ksnopt:
    sncmxSnopt( nlhs, plhs, nrhs, prhs, ksimpleCallbacks   );
    break;

  case ksnoptAdvanced:
    sncmxSnopt( nlhs, plhs, nrhs, prhs, kadvancedCallbacks );
    break;

  case ksnoptJac:
    sncmxJac  ( nlhs, plhs, nrhs, prhs );
    break;

  case ksnoptValue:
  case ksnoptSetOption:
  case ksnoptSetOptionI:
  case ksnoptSetOptionR:
  case ksnoptGetOption:
  case ksnoptGetOptionC:
  case ksnoptGetOptionI:
  case ksnoptGetOptionR:
  case ksnoptSpecsFile:
  case ksnoptSetPrintFile:
  case ksnoptSetSummaryFile:
    /* All these require a string argument for input, */
    /* so we pass them all to the same routine.       */
    callType = userCall;
    sncmxSingleStringArg( nlhs, plhs, nrhs, prhs, what );
    callType = systemCall;
    break;
  case ksnoptScreenOn:
    screenIsOn = 1;
    break;
  case ksnoptScreenOff:
    screenIsOn = 0;
    break;
  case ksnoptClosePrintFile:
    /* the print file. */
    if ( 1 != nrhs )
      mexErrMsgTxt( "Wrong number of input arguments" );
    if ( printFileIsOpen ) { /* Only close it if it was open */
      printFileIsOpen = sncmxCloseFile( snPrintFile );
    }
    break;
  case ksnoptCloseSummaryFile:
    /* There are no additional input argument, so we simply close */
    /* the summary file.                                          */
    if ( 1 != nrhs )
      mexErrMsgTxt( "Wrong number of input arguments" );
    if( summaryFileIsOpen ) { /* Only close it if it was open */
      summaryFileIsOpen = sncmxCloseFile( snSummaryFile );
    }
    break;
  case ksnoptSetStatus:
    if ( 2 != nrhs )
      mexErrMsgTxt( "Wrong number of input arguments" );
     /* gstatus :input argument 2 */
    assertScalar( prhs[1], "status" );
    gstatus = *mxGetPr( prhs[1] );
    setStatus( &gstatus );
    break;
  case ksnoptGetStatus:
    if ( 1 != nrhs )
      mexErrMsgTxt( "Wrong number of input arguments" );
    /* gstatus :input argument 2 */
    getStatus( &gstatus );
    plhs[0] = mxCreateDoubleMatrix( 1, 1, mxREAL );
    dummyDoublePtr  = mxGetPr( plhs[0] );
    *dummyDoublePtr = (doublereal) gstatus;
    break;

  default:
    mexErrMsgTxt( "Unrecognized request for the SNOPT mexFunction" );
    break;
  }
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
void sncmxSingleStringArg
( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[], int what )
{
  char *stringArg; /* input argument 2 */

  int         callSave, info;
  doublereal *dummyDoublePtr, rvalue;
  integer     lenStringArg, inform = 0, ivalue, OptionIsSet;
  char        cvalue[8];

  /* Check for correct number of arguments */
  switch( what ) {
  case ksnoptValue:
  case ksnoptSpecsFile:
  case ksnoptSetOption:
  case ksnoptGetOption:
  case ksnoptGetOptionI:
  case ksnoptGetOptionR:
  case ksnoptGetOptionC:
  case ksnoptSetPrintFile:
  case ksnoptSetSummaryFile:
    if ( 2 != nrhs ) {
      mexErrMsgTxt( "Wrong number of input arguments" );
    }
    break;
  case ksnoptSetOptionI:
  case ksnoptSetOptionR:
    if ( 3 != nrhs ) {
      mexErrMsgTxt( "Wrong number of input arguments" );
    }
  }

  /* Read the string argument */
  assertString( prhs[1], "stringArg" );
  lenStringArg = mxGetN( prhs[1] );
  stringArg    = mxCalloc( (lenStringArg + 1), sizeof( char ) );
  mxGetString( prhs[1], stringArg, lenStringArg + 1);

  switch( what ) {
  case ksnoptValue:
    inform = sncmxGetksnoptValue( stringArg );
    break;
  case ksnoptSpecsFile:
    inform = sncmxSpecs( stringArg );
    break;
  case ksnoptSetOption:
    inform = sncmxSet( stringArg );
    /* Set any non-standard Matlab defaults here */
    if ( inform == 0 && (   strcmp(stringArg, "Defaults" ) == 0
                         || strcmp(stringArg, "defaults" ) == 0)) {
      callSave = callType;
      callType = systemCall;
      inform   = sncmxSeti( "Timing level ",     -1 );
      inform   = sncmxSet ( "Sticky parameters yes" );
      /* ...   */
      callType = callSave;
    }
    break;
  case ksnoptSetOptionI:
    assertScalar( prhs[2], "Ivalue" );
    ivalue =    (integer)*mxGetPr( prhs[2] );
    inform =  sncmxSeti( stringArg, ivalue );
    break;
  case ksnoptSetOptionR:
    assertScalar( prhs[2], "Rvalue" );
    rvalue = (doublereal)*mxGetPr( prhs[2] );
    inform =  sncmxSetr( stringArg, rvalue );
    break;
  case ksnoptGetOption:
    OptionIsSet = sncmxGet( stringArg );
    break;
  case ksnoptGetOptionI:
    inform =  sncmxGeti( stringArg, &ivalue );
    break;
  case ksnoptGetOptionR:
    inform =  sncmxGetr( stringArg, &rvalue );
    break;
  case ksnoptSetPrintFile:
    /* Close the current print file */
    if ( printFileIsOpen ) info = sncmxCloseFile( snPrintFile );
    snPrintFile = sncmxOpenReplace( stringArg );
    if ( snPrintFile != NULL ) {
      /* We were able to open the file. */
      printFileIsOpen = 1;
      /* REMEMBER set some options here later */
    }
    break;
  case ksnoptSetSummaryFile:
    /* Close the current summary file */
    if ( summaryFileIsOpen ) info = sncmxCloseFile( snSummaryFile );
    snSummaryFile = sncmxOpenAppend( stringArg );
    if ( snSummaryFile != NULL ) {
      /* We were able to open the file */
      summaryFileIsOpen = 1;
      /* REMEMBER set some options here later */
    }
    break;
  }

  mxFree( stringArg );
  switch (what) {
  case ksnoptValue:
  case ksnoptSetOption:
  case ksnoptSetOptionI:
  case ksnoptSetOptionR:
  case ksnoptSpecsFile:
  case ksnoptSetPrintFile:
  case ksnoptSetSummaryFile:
    /* These all return 0 if successful */
    plhs[0]         = mxCreateDoubleMatrix( 1, 1, mxREAL );
    dummyDoublePtr  = mxGetPr( plhs[0] );
    *dummyDoublePtr = inform;
    break;
  case ksnoptGetOption:
    plhs[0]         = mxCreateDoubleMatrix( 1, 1, mxREAL );
    dummyDoublePtr  = mxGetPr( plhs[0] );
    *dummyDoublePtr = OptionIsSet;
    break;
  case ksnoptGetOptionC:
    plhs[0]         = mxCreateString(cvalue);
    break;
  case ksnoptGetOptionI:
    plhs[0]         = mxCreateDoubleMatrix( 1, 1, mxREAL );
    dummyDoublePtr  = mxGetPr( plhs[0] );
    *dummyDoublePtr = ivalue;
    break;
  case ksnoptGetOptionR:
    plhs[0]         = mxCreateDoubleMatrix( 1, 1, mxREAL );
    dummyDoublePtr  = mxGetPr( plhs[0] );
    *dummyDoublePtr = rvalue;
    break;
  }
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
int sncmxAlloc( integer neF, integer n, integer nxname,
		integer nFname, integer neA, integer neG )
{
  int        k;
  integer    mincw, miniw, minrw;
  integer    inform;

  if (!memoryIsSet){

    /*   Estimate mincw, miniw, minrw  */
    snmema0_( &inform, &snPrintUnit, &snSummaryUnit, &neF, &n,
	      &nxname, &nFname, &neA, &neG,
	      &mincw, &miniw, &minrw,
	      cw, &lencw, iw, &leniw, rw, &lenrw, 8*lencw);

    lrwsnopt = minrw;
    liwsnopt = miniw;
    lcwsnopt = mincw;

    rwsnopt = mxCalloc( lrwsnopt, sizeof(doublereal) );
    iwsnopt = mxCalloc( liwsnopt, sizeof(integer)    );
    cwsnopt = mxCalloc( lcwsnopt, 8*sizeof(char)     );

    if ((rwsnopt != NULL) && (iwsnopt != NULL) && (cwsnopt != NULL)) {

      /* Update cw, iw and rw and copy them into the allocated arrays */

      sncmxSeti( "Total real workspace     ", lrwsnopt );
      sncmxSeti( "Total integer workspace  ", liwsnopt );
      sncmxSeti( "Total character workspace", lcwsnopt );

      for ( k = 0; k < 500; k++ ) {
	rwsnopt[k] = rw[k];
	iwsnopt[k] = iw[k];
      }
      for ( k = 0; k < 8*500; k++ ) {
	cwsnopt[k] = cw[k];
      }
      memoryIsSet = 1;
    }
  }
  return memoryIsSet == 0;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
int sncmxReAlloc( integer mincw, integer miniw, integer minrw)
/* Kluge to avoid using Matlab's buggy mxRealloc */
{
  int        k;

  if ( memoryIsSet ){
    for ( k = 0; k < 500; k++ ) {
      rw[k] = rwsnopt[k];
      iw[k] = iwsnopt[k];
    }
    for ( k = 0; k < 8*500; k++ ) {
      cw[k] = cwsnopt[k];
    }
    mxFree( rwsnopt );
    mxFree( iwsnopt );
    mxFree( cwsnopt );
    sncmxSeti( "Total real workspace     ", minrw );
    sncmxSeti( "Total integer workspace  ", miniw );
    sncmxSeti( "Total character workspace", mincw );
    memoryIsSet = 0;
  }

  lrwsnopt = minrw;
  liwsnopt = miniw;
  lcwsnopt = mincw;

  rwsnopt  = mxCalloc( lrwsnopt, sizeof(doublereal) );
  iwsnopt  = mxCalloc( liwsnopt, sizeof(integer)    );
  cwsnopt  = mxCalloc( lcwsnopt, 8*sizeof(char)     );

  if ((rwsnopt != NULL) && (iwsnopt != NULL) && (cwsnopt != NULL)) {
    for ( k = 0; k < 500; k++ ) {
      rwsnopt[k] = rw[k];
      iwsnopt[k] = iw[k];
    }
    for ( k = 0; k < 8*500; k++ ) {
      cwsnopt[k] = cw[k];
    }

    memoryIsSet = 1;
  }
  return memoryIsSet == 0;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
void sncmxInit()
     /* Should be called at most once per Matlab session                */
     /* Define all options as "unset".  Fortran file unit numbers are   */
     /* set to dummy values as printing is controlled by the global     */
     /* variables printFileisOpen and summaryFileisOpen.                */
{
  if ( snIsInitialized ) {
    mexErrMsgTxt( "SNOPT has already been intialized" );
  } else {
    sninit_( &snPrintUnit, &snSummaryUnit, cw, &lencw,
	     iw, &leniw, rw, &lenrw, 8*lencw );

    snIsInitialized = 1;
  }
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
integer sncmxGet ( char option[] )
{
  integer inform = 0, ivalue;
  integer lenopt = strlen(option);

  ivalue = snget_( option, &inform, cw, &lencw, iw, &leniw, rw, &lenrw,
		   lenopt, 8*lencw);

  return ivalue;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
integer sncmxGetc( char option[], char *cvalue[] )
{
  integer inform = 0;
  integer lenopt = strlen(option);

  sngetc_( option, cvalue, &inform, cw, &lencw, iw, &leniw, rw, &lenrw,
	   lenopt, 8, 8*lencw);

  return inform;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
integer sncmxGeti( char option[], integer *ivalue )
{
  integer inform = 0;
  integer lenopt = strlen(option);

  sngeti_( option, ivalue, &inform, cw, &lencw, iw, &leniw, rw, &lenrw,
           lenopt, 8*lencw);

  return inform;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
integer sncmxGetr( char option[], doublereal *rvalue )
{
  integer inform = 0;
  integer lenopt = strlen(option);

  sngetr_( option, rvalue, &inform, cw, &lencw, iw, &leniw, rw, &lenrw,
	   lenopt, 8*lencw);

  return inform;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
integer sncmxSet ( char option[] )
{
  integer inform = 0;
  integer lenopt = strlen(option);

  snset_( option, &snPrintUnit, &snSummaryUnit, &inform, cw, &lencw,
	 iw, &leniw, rw, &lenrw, lenopt, 8*lencw);

  return inform;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
integer sncmxSeti( char option[], integer iopt )
{
  integer inform = 0;
  integer lenopt = strlen(option);

  snseti_( option, &iopt, &snPrintUnit, &snSummaryUnit, &inform, cw,
	   &lencw, iw, &leniw, rw, &lenrw, lenopt, 8*lencw);

  return inform;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
integer sncmxSetr( char option[], doublereal ropt )
{
  integer inform = 0;
  integer lenopt = strlen(option);

  snsetr_( option, &ropt, &snPrintUnit, &snSummaryUnit, &inform, cw,
	   &lencw, iw, &leniw, rw, &lenrw, lenopt, 8*lencw);

  return inform;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
integer sncmxSpecs( char specsFile[] )
{
  integer inform = 0, lenfile = strlen(specsFile);

  if (specsFile == NULL) specsFile = "\0";
  snSpecsFile = fopen( specsFile, "r");
  if (snSpecsFile == NULL)
    fprintf( stderr,
             "sncmxSpecs(\"%s\") failed: %s\n",
             specsFile, strerror(errno));

  snspec_( &snSpecsUnit, &inform, cw, &lencw,
	   iw, &leniw, rw, &lenrw, 8*lencw );

  sncmxCloseFile( snSpecsFile );

  return inform;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
void sncmxFreeWorkSpace()
{
  int        k;

  if ( memoryIsSet ){
    for ( k = 0; k < 500; k++ ) {
      rw[k] = rwsnopt[k];
      iw[k] = iwsnopt[k];
    }
    for ( k = 0; k < 8*500; k++ ) {
      cw[k] = cwsnopt[k];
    }
    mxFree( rwsnopt );
    mxFree( iwsnopt );
    mxFree( cwsnopt );

    sncmxSeti( "Total real workspace     ", lenrw );
    sncmxSeti( "Total integer workspace  ", leniw );
    sncmxSeti( "Total character workspace", lencw );
    memoryIsSet = 0;
  }
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
void sncmxSnopt( int nlhs, mxArray *plhs[], int nrhs,
		 const mxArray *prhs[], integer advancedCallbacks )
{
  /* ***************Input Arguments*********************** */
  char     *Prob;         /* Input Argument 16 (optional)  */
  double   *x0;           /* Input Argument 1              */
  double   *xlow;         /* Input Argument 2              */
  double   *xupp;         /* Input Argument 3              */
  double   *xmul0;        /* Input Argument 4              */
  /*integer  *xstate0;       Input Argument 5              */
  char     *xnames;       /* Inpute Argument 18 (optional) */

  double   *Flow;         /* Input Argument 6              */
  double   *Fupp;         /* Input Argument 7              */
  double   *Fmul0;        /* Input Argument 8              */
  /*integer  *Fstate0;       Input Argument 9              */
  char     *Fnames;       /* Input Argument  19 (optional) */

  double   ObjAdd;        /* Input Argument 10             */
  integer  ObjRow;        /* Input Argument 11             */

  double   *A;            /* Input Argument 12             */
  integer  *iAfun;        /* Input Argument 13             */
  integer  *jAvar;        /* Input Argument 14             */
  integer  *iGfun;        /* Input Argument 15             */
  integer  *jGvar;        /* Input Argument 16             */
  /* char userfun[32]        Input Argument 17 (declared globally) */
  /********************End Input Arguments******************/

  double   *dummyDoublePtr;
  double   *x, *xmul;
  integer  *xstate;
  mxArray  *mexx, *mexxmul, *mexxstate;
  double   *F, *Fmul;
  integer  *Fstate;
  mxArray  *mexF, *mexFmul, *mexFstate;
  integer  n;              /* length of x0        */
  integer  neF;            /* number of rows of F */

  integer  lenA, neA;      /* number of nonzero elements in A */
  integer  lenG, neG;

  integer  nxname = 1, nFname = 1; /* No names for now. */
  integer  mincw, miniw, minrw, ns, ninf;
  double   sinf;

  integer  Start;
  integer  Cold    = 0;
  integer  Warm    = 2;
  integer  info    = 0;
  integer  Freebie = 0;
  integer  maxColSize = 300, maxRowSize = 300; /* if Freebie == 1 */

  integer  k, memError, m, zzz;
  double   *AA;


  if ( nrhs != 18 && nrhs != 21 ){
    mexErrMsgTxt( "Wrong number of input arguments" );
  }

  /* x0 :input argument 2, read this first to get a value for n */
  n  = mxGetM( prhs[1] );
  if (Freebie > 0){
    if ( n > maxColSize ) {
      mexErrMsgTxt( "Maximum column size of 300 exceeded" );
    }
  }
  assertColDim( prhs[1], 1, "x0" );
  x0 = mxGetPr( prhs[1] );
  /* xlow :input argument 3 */
  assertColDim( prhs[2], 1, "xlow" );
  assertRowDim( prhs[2], n, "xlow" );
  xlow = mxGetPr( prhs[2] );
  /* xupp :input argument 4 */
  assertColDim( prhs[3], 1, "xupp" );
  assertRowDim( prhs[3], n, "xupp" );
  xupp = mxGetPr( prhs[3] );
  /* xmul0 :input argument 5 */
  assertColDim( prhs[4], 1, "xmul0" );
  assertRowDim( prhs[4], n, "xmul0" );
  xmul0 = mxGetPr( prhs[4] );
  /* xstate0 :inpute argument 6 */
  assertColDim( prhs[5], 1, "xstate0" );
  assertRowDim( prhs[5], n, "xstate0" );
  /* xstate0 = mxGetPr( prhs[5] ); */
  /* xstate0 will be copied later into output argument xstate. */

  /* Flow :input argument 7 */
  neF  = mxGetM( prhs[6] );
  if (Freebie > 0){
    if ( neF > maxRowSize ) {
      mexErrMsgTxt( "Maximum row size of 300 exceeded" );
    }
  }
  assertColDim( prhs[6], 1, "Flow" );
  Flow = mxGetPr( prhs[6] );
  /* Fupp :input argument 8 */
  assertColDim( prhs[7],   1, "Fupp" );
  assertRowDim( prhs[7], neF, "Fupp" );
  Fupp = mxGetPr( prhs[7] );
  /* Fmul0 :input argument 9 */
  assertColDim( prhs[8],   1, "Fmul0" );
  assertRowDim( prhs[8], neF, "Fmul0" );
  Fmul0 = mxGetPr( prhs[8] );
  /* Fstate0 :input argument 10 */
  assertColDim( prhs[9],   1, "Fstate0" );
  assertRowDim( prhs[9], neF, "Fstate0" );
  /* Fstate0 = mxGetPr( prhs[9] ); */
  /* Fstate0 will be copied later into output argument Fstate. */

  /* ObjAdd :input argument 11 */
  assertScalar( prhs[10], "ObjAdd" );
  ObjAdd = *mxGetPr( prhs[10] );
  /* ObjRow :input argument 12 */
  assertScalar( prhs[11], "ObjRow" );
  ObjRow = (integer)*mxGetPr( prhs[11] );

  /* A :input argument 13 */
  lenA  = mxGetM( prhs[12] );
  neA   = lenA;
  if ( lenA > 0 ) {
    assertColDim( prhs[12], 1, "A" );
    A     = mxGetPr( prhs[12] );
    /* iAfun :input argument 14 */
    assertColDim( prhs[13],    1, "iAfun" );
    assertRowDim( prhs[13], lenA, "iAfun" );
    iAfun    = mxCalloc( lenA, sizeof( integer ) );
    dble2int( lenA*1, mxGetPr( prhs[13] ), iAfun );
    /* jAvar :input argument 15 */
    assertColDim( prhs[14],    1, "jAvar" );
    assertRowDim( prhs[14], lenA, "jAvar" );
    jAvar    = mxCalloc( lenA, sizeof( integer ) );
    dble2int( lenA*1, mxGetPr( prhs[14] ), jAvar );
  } else {
    A     = 0;
    iAfun = 0;
    jAvar = 0;
  }
  lenG  = mxGetM( prhs[15] );
  neG   = lenG;
  if ( lenG > 0 ){
    /* iGfun :input argument 16  */
    assertColDim( prhs[15], 1, "iGfun" );
    iGfun    = mxCalloc( lenG, sizeof( integer ) );
    dble2int( lenG*1, mxGetPr( prhs[15] ), iGfun );
    /* jGvar :input argument 17 */
    assertColDim( prhs[16],    1, "jGvar" );
    assertRowDim( prhs[16], lenG, "jGvar" );
    jGvar    = mxCalloc( lenG, sizeof( integer ) );
    dble2int( lenG*1, mxGetPr( prhs[16] ), jGvar );
  } else {
    iGfun = 0;
    jGvar = 0;
  }

  /* userfg:input argument 18 */
  assertString( prhs[17], "userfun" );
  mxGetString( prhs[17], userfun, 32 );

  /* userfg:input argument 19-20 */
  /*        optional arguments for "fmincon" version of snopt interface */
  if ( nrhs == 21 ) {
    zzz = 1;

    assertString( prhs[18], "myobj" );
    mxGetString( prhs[18], myobj, 32 );

    assertString( prhs[19], "nonlcon" );
    mxGetString( prhs[19], nonlcon, 32 );

  } else {
    zzz        = 0;

    myobj[0]   = 0;
    nonlcon[0] = 0;
  }

  /* x     : output argument 1 */
  mexx = mxCreateDoubleMatrix( n, 1, mxREAL );
  x    = mxGetPr( mexx );
  memmove( x, x0, n*sizeof(double) );

  /* F     : output argument 2 */
  mexF = mxCreateDoubleMatrix( neF, 1, mxREAL );
  F    = mxGetPr( mexF );

  /* xmul  : output argument 3 */
  mexxmul = mxCreateDoubleMatrix( n, 1, mxREAL);
  xmul    = mxGetPr( mexxmul );
  memmove( xmul, xmul0, n*sizeof(double) );

  /* Fmul  : output argument 4 */
  mexFmul = mxCreateDoubleMatrix( neF, 1, mxREAL );
  Fmul    = mxGetPr( mexFmul );
  memmove( Fmul, Fmul0, neF*sizeof(double) );

  /* xstate: output argument 5 */
  mexxstate = mxCreateDoubleMatrix( n, 1, mxREAL );
  xstate    = (integer *) mxGetPr( mexxstate );
  dble2int( n*1, mxGetPr(prhs[5]), xstate );

  /* Fstate: output argument 6 */
  mexFstate = mxCreateDoubleMatrix( neF, 1, mxREAL );
  Fstate    = (integer *) mxGetPr( mexFstate );
  dble2int( neF*1, mxGetPr(prhs[9]), Fstate );

  /* Make Prob an option later */
  Prob   = mxCalloc(      1, 8*sizeof(char) );
  xnames = mxCalloc( nxname, 8*sizeof(char) );
  Fnames = mxCalloc( nFname, 8*sizeof(char) );

  /* If needed, create memory for userfg_ */
  if ( neG > 0 ) {
    sncmxAllocFuns( &n, &neG, iGfun, jGvar );
  }

  memError = sncmxAlloc( neF, n, nxname, nFname, neA, neG );

  if ( zzz == 1 ) {
    assertColDim( prhs[20], n, "AA" );
    m  = mxGetM( prhs[20] );
    AA = mxGetPr( prhs[20] );

    lrusnopt = m*n + 1;
    rusnopt  = mxCalloc( lrusnopt, sizeof(doublereal) );

    for ( k = 0; k < m*n; k++ ) {
      rusnopt[k] = AA[k];
    }
    rusnopt[lrusnopt-1] = m;
  }
  else {
    lrusnopt = lrwsnopt;
    rusnopt  = rwsnopt;
  }


  if (memError == 0) {
    /* Assign dummy problem name  */
    /*  sncmxSetProbName( Prob ); */
    snprob_( Prob, (ftnlen)8 );

    Start = Cold;

    for ( k = 0; k < 5; k++ ) {

      callType = userCall;

      snkera_( &Start, &neF, &n, &nxname, &nFname,
	       &ObjAdd, &ObjRow, Prob,
	       userfg_, snlog_, snlog2_, sqlog_, snabort_,
	       iAfun, jAvar, &lenA, &neA, A,
	       iGfun, jGvar, &lenG, &neG,
	       xlow, xupp, xnames, Flow, Fupp, Fnames,
	       x, xstate, xmul, F, Fstate, Fmul,
	       &info, &mincw, &miniw, &minrw,
	       &ns, &ninf, &sinf,
	       cwsnopt, &lcwsnopt, iwsnopt, &liwsnopt, rusnopt, &lrusnopt,
	       cwsnopt, &lcwsnopt, iwsnopt, &liwsnopt, rwsnopt, &lrwsnopt,
	       (ftnlen)8, (ftnlen)8,  (ftnlen)8, (ftnlen)8,  (ftnlen)8);

      callType = systemCall;

      if ((info != 82) && (info != 83) && (info != 84)) break;
      /* The LU ran out of memory.  Allocate more space and continue. */
      Start       = Warm;
      iwsnopt[68] = Start;              /* Options override arguments */
      memError    = sncmxReAlloc( mincw, miniw, minrw );
      if (memError > 0) break;
    }
  }

  sncmxFreeFuns( );

  plhs[0] = mexx;
  if ( 1 < nlhs ) {
    plhs[1] = mexF;
  }
  if ( 2 < nlhs ) {
    plhs[2] = mexxmul;
  }
  if ( 3 < nlhs ) {
    plhs[3] = mexFmul;
  }
  if ( 4 < nlhs ) {
    plhs[4] = mxCreateDoubleMatrix( 1, 1, mxREAL );
    dummyDoublePtr  = mxGetPr( plhs[4] );
    *dummyDoublePtr = (double) info;
  }

  if ( 5 < nlhs ) {
    plhs[5] = mexxstate;
    int2dble( n*1, xstate, (double *) xstate);
  }

  if ( 6 < nlhs ) {
    plhs[6] = mexFstate;
    int2dble( neF*1, Fstate, (double *) Fstate);
  }

  if ( 7 < nlhs ) {
    plhs[7] = mxCreateDoubleMatrix( 1, 1, mxREAL );
    dummyDoublePtr  = mxGetPr( plhs[7] );
    *dummyDoublePtr = (double) ns;
  }

 if ( 8 < nlhs ) {
   plhs[8] = mxCreateDoubleMatrix( 1, 1, mxREAL );
   dummyDoublePtr  = mxGetPr( plhs[8] );
   *dummyDoublePtr = (double) ninf;
  }

 if ( 9 < nlhs ) {
   plhs[9] = mxCreateDoubleMatrix( 1, 1, mxREAL );
   dummyDoublePtr  = mxGetPr( plhs[9] );
   *dummyDoublePtr = (double) sinf;
  }

 if ( 10 < nlhs ) {
   plhs[10] = mxCreateDoubleMatrix( 1, 1, mxREAL );
   dummyDoublePtr  = mxGetPr( plhs[10] );
   *dummyDoublePtr = (double) mincw;
  }

 if ( 11 < nlhs ) {
   plhs[11] = mxCreateDoubleMatrix( 1, 1, mxREAL );
   dummyDoublePtr  = mxGetPr( plhs[11] );
   *dummyDoublePtr = (double) miniw;
  }

 if ( 12 < nlhs ) {
   plhs[12] = mxCreateDoubleMatrix( 1, 1, mxREAL );
   dummyDoublePtr  = mxGetPr( plhs[12] );
   *dummyDoublePtr = (double) minrw;
  }

 if ( zzz == 1 ) {
   mxFree( rusnopt );
 }

 sncmxFreeWorkSpace();
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
integer sncmxGetksnoptValue( char *ksnoptname )
{
  if ( strcmp(ksnoptname, "Solve")            == 0 )
    return ksnopt;
  if ( strcmp(ksnoptname, "SetOption")        == 0 )
    return ksnoptSetOption;
  if ( strcmp(ksnoptname, "SetOptionI")       == 0 )
    return ksnoptSetOptionI;
  if ( strcmp(ksnoptname, "SetOptionR")       == 0 )
    return ksnoptSetOptionR;
  if ( strcmp(ksnoptname, "SpecsFile")        == 0 )
    return ksnoptSpecsFile;
  if ( strcmp(ksnoptname, "SetPrintFile")     == 0 )
    return ksnoptSetPrintFile;
  if ( strcmp(ksnoptname, "SetSummaryFile")   == 0 )
    return ksnoptSetSummaryFile;
  if ( strcmp(ksnoptname, "ClosePrintFile")   == 0 )
    return ksnoptClosePrintFile;
  if ( strcmp(ksnoptname, "CloseSummaryFile") == 0 )
    return ksnoptCloseSummaryFile;
  if ( strcmp(ksnoptname, "Advanced")         == 0 )
    return ksnoptAdvanced;

  mexErrMsgTxt("Unknown string name");
  return 0;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
void sncmxJac( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  /* ***************Input Variables**************** */
  /* char userfun[32]              Input Argument 1 */
  doublereal *x0;               /* Input Argument 2 */
  doublereal *xlow;             /* Input Argument 3 */
  doublereal *xupp;             /* Input Argument 4 */
  integer    neF;               /* Input Argument 3 */
  /* ***************Input Variables**************** */
  integer     n;
  integer     lenA, lenG, neA, neG;
  integer     mincw, miniw, minrw;

  doublereal *A;
  integer    *iAfun;
  integer    *jAvar;
  integer    *iGfun;
  integer    *jGvar;

  mxArray    *mexA, *mexiAfun, *mexjAvar;
  mxArray    *mexiGfun, *mexjGvar;

  integer     nxname = 1, nFname = 1, ObjRow = 1;
  integer     inform = 0;
  int         memError;

  if ( nrhs != 6 ){
    mexErrMsgTxt( "Wrong number of input arguments" );
  }

  /* userfg :input argument 1 */
  assertString( prhs[1], "userfun" );
  mxGetString ( prhs[1], userfun, 32 );

  /* x0     :input argument 2 */
  n  = mxGetM( prhs[2] );
  assertColDim( prhs[2], 1, "x0" );
  x0 = mxGetPr( prhs[2] );

  /* xlow   :input argument 3 */
  assertColDim( prhs[3], 1, "xlow" );
  assertRowDim( prhs[3], n, "xlow" );
  xlow = mxGetPr( prhs[3] );

  /* xupp   :input argument 4 */
  assertColDim( prhs[4], 1, "xupp" );
  assertRowDim( prhs[4], n, "xupp" );
  xupp = mxGetPr( prhs[4] );

  /* neF   :input argument 5 */
  assertScalar( prhs[5], "neF" );
  neF = (integer)*mxGetPr( prhs[5] );

  lenA = n*neF; neA = lenA;
  lenG = lenA;  neG = lenA;

  iAfun = mxCalloc( lenA, sizeof(integer)    );
  jAvar = mxCalloc( lenA, sizeof(integer)    );
  A     = mxCalloc( lenA, sizeof(doublereal) );
  iGfun = mxCalloc( lenG, sizeof(integer)    );
  jGvar = mxCalloc( lenG, sizeof(integer)    );

  sncmxAllocJac( &n );
  memError = sncmxAlloc   ( neF, n, nxname, nFname, neA, neG);

  if (memError == 0) {
    snjac_( &inform, &neF, &n, userfg_,
	    iAfun, jAvar, &lenA, &neA, A,
	    iGfun, jGvar, &lenG, &neG,
	    x0, xlow, xupp,
	    &mincw, &miniw, &minrw,
	    cwsnopt, &lcwsnopt, iwsnopt, &liwsnopt, rwsnopt, &lrwsnopt,
	    cwsnopt, &lcwsnopt, iwsnopt, &liwsnopt, rwsnopt, &lrwsnopt,
	    8*lcwsnopt, 8*lcwsnopt);
  }
  sncmxFreeJac();

  /* mexA     : output argument 1 */
  mexA = mxCreateDoubleMatrix( neA, 1, mxREAL );
  memmove( mxGetPr(mexA), A, neA*sizeof(double) );
  plhs[0] = mexA;

  /* mexiAfun : output argument 2 */
  if ( 1 < nlhs ) {
    mexiAfun = mxCreateDoubleMatrix( neA, 1, mxREAL );
    plhs[1] = mexiAfun;
    int2dble( neA*1, iAfun, mxGetPr(mexiAfun));
  }

  /* mexjAvar : output argument 3 */
  if ( 2 < nlhs ) {
    mexjAvar = mxCreateDoubleMatrix( neA, 1, mxREAL );
    plhs[2] = mexjAvar;
    int2dble( neA*1, jAvar, mxGetPr(mexjAvar));
  }

  /* mexiGfun : output argument 4 */
  if ( 3 < nlhs ) {
    mexiGfun = mxCreateDoubleMatrix( neG, 1, mxREAL );
    plhs[3] = mexiGfun;
    int2dble( neG*1, iGfun, mxGetPr(mexiGfun));
  }

  /* mexjGvar : output argument 5 */
  if ( 4 < nlhs ) {
    mexjGvar = mxCreateDoubleMatrix( neG, 1, mxREAL );
    plhs[4] = mexjGvar;
    int2dble( neG*1, jGvar, mxGetPr(mexjGvar));
  }

  mxFree(iAfun);
  mxFree(jAvar);
  mxFree(A);
  mxFree(iGfun);
  mxFree(jGvar);

  sncmxFreeWorkSpace();
}
