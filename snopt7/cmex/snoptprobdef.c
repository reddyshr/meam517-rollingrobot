/* Created  20 July 2002 Josh Griffin  */
/* Modified 18 June 2005 Philip Gill   */

#include <memory.h>
#include "snoptprobdef.h"
#include "mexUtility.h"
#include "snoptcmex.h"
#include "matrix.h"

static mxArray *mexx     = 0;
static mxArray *mexiGfun = 0;
static mxArray *mexjGvar = 0;

static integer derivativeOption; /* Should equal 0 or 1 */
static integer globalStatus = 1;

char userfun[32];
char   myobj[32];
char nonlcon[32];

static void copyIfNotNAN( long rows, long cols, double from[], double to[] );

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int sncmxAllocJac(integer *n)
{
  mexx      = mxCreateDoubleMatrix( *n,   1, mxREAL );
  return 0;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int sncmxFreeJac ( )
{
  if ( mexx ) {
    mxDestroyArray( mexx );
    mexx = 0;
  }
  return 0;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int sncmxAllocFuns( integer *n, integer *neG, integer iGfun[], integer jGvar[] )
{
  integer inform;

  inform = sncmxGeti("Derivative Option", &derivativeOption );
  /*Note: If derivative option is not set, it will be -11111 not 1.*/
  if ( derivativeOption != 0 && derivativeOption != 1 ) {
    derivativeOption = 1;
  }

  mexx      = mxCreateDoubleMatrix( *n,   1, mxREAL );
  mexiGfun  = mxCreateDoubleMatrix( *neG, 1, mxREAL );
  mexjGvar  = mxCreateDoubleMatrix( *neG, 1, mxREAL );
  if( iGfun == 0 || jGvar == 0 )
    mexErrMsgTxt("iGfun and iGvar must be set prior to call of inifns_.");

  int2dble( *neG*1, iGfun, mxGetPr(mexiGfun) );
  int2dble( *neG*1, jGvar, mxGetPr(mexjGvar) );

  return 0;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int sncmxFreeFuns( void )
{
  if ( mexx ) {
    mxDestroyArray( mexx );
    mexx = 0;
  }
  if ( mexiGfun ) {
    mxDestroyArray( mexiGfun );
    mexiGfun = 0;
  }
  if ( mexjGvar ) {
    mxDestroyArray( mexjGvar );
    mexjGvar = 0;
  }

  return 0;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int setStatus( integer *status_value )
{
  globalStatus = *status_value;
  return 0;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

int getStatus( integer *status_value )
{
  *status_value = globalStatus;
  return 0;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

int userfg_
( integer    *Status, integer *n,    doublereal x[],
  integer    *needF,  integer *neF,  doublereal F[],
  integer    *needG,  integer *neG,  doublereal G[],
  char       *cu,     integer *lencu,
  integer    iu[],    integer *leniu,
  doublereal ru[],    integer *lenru )
     /* =========================================================== */
     /* Make callback to the user-supplied m-file:                  */
     /*    function [F,G] = userfun(x)                              */
     /* =========================================================== */
{
  int nlhs = 2;
  int nrhs = 5;
  mxArray *plhs[2];
  /*  mxArray *prhs[2]; */
  mxArray *prhs[5];
  mxArray *mexF, *mexG;

  /* Define the following variables in case a full G is returned */
  int Gnlhs = 1;
  int Gnrhs = 3;
  mxArray *Gplhs[1];
  mxArray *Gprhs[3];

  /* For fmincon interface */
  int     i, j, ij, m, off;
  double  *dummyDoublePtr;


  globalStatus = *Status;
  memmove( mxGetPr( mexx ), x, *n*sizeof(double) );
  prhs[0] = mexx;
  prhs[1] = mxCreateString(userfun);
  prhs[2] = mxCreateDoubleMatrix( 1, 1, mxREAL );
  dummyDoublePtr  = mxGetPr( prhs[2] );
  *dummyDoublePtr = (double) *neF;


  if ( strlen(myobj) > 0 ) {
    nrhs    = 4;
    prhs[3] = mxCreateString(myobj);

    if ( strlen(nonlcon) > 0 ) {
      nrhs    = 5;
      prhs[4] = mxCreateString(nonlcon);
    }
  }
  else {
    nrhs = 3;
  }

  mexCallMATLAB( nlhs, plhs, nrhs, prhs, "snwrapper" );

  if( globalStatus >= 0 ) { /* Functions are defined */

    if( *needF ) {
      mexF = plhs[0];
      assertRowDim( mexF, *neF, "F" );
      assertColDim( mexF,    1, "F" );
      memmove( F, mxGetPr( mexF ), *neF *sizeof(double));
      mxDestroyArray( mexF );

      /* The linear constraints (the matrix A) are held in ru.
         We need to compute A*x  */
      if ( *lenru > 0 ) {
	m   = ru[*lenru - 1];
	off = *neF - m;
	for ( i = 0; i < m; i++ ) {
	  F[i+off] = 0;
	}

	for ( j = 0; j < *n; j++ ) {
	  for ( i = 0; i < m; i++) {
	    ij  = i+j*m;
	    F[i+off] = F[i+off] + ru[ij]*x[j];
	  }
	}
      }

    }

    if( *needG ) {
      mexG = plhs[1];
      if ( mxGetN(mexG) == 0 ) {
	/* No derivatives are defined: do nothing.
	   The input G is left untouched.              */

      } else {
	/* Some derivatives have been set.
	   Copy derivatives that are not NaNs to G.    */

	if( mxGetN(mexG) > 1 ) {
	  /* mexG is a full  neF by n  Jacobian.       */

	  assertRowDim( mexG, *neF, "G" );
	  assertColDim( mexG, *n,   "G" );
	  Gprhs[0] = mexiGfun;
	  Gprhs[1] = mexjGvar;
	  Gprhs[2] = mexG;

	  mexCallMATLAB( Gnlhs, Gplhs, Gnrhs, Gprhs, "snfindG");
	  mxDestroyArray( mexG );

	  mexG = Gplhs[0];
	}

	assertRowDim( mexG, *neG, "G" );
	assertColDim( mexG,    1, "G" );

	if (derivativeOption != 1) {
	  copyIfNotNAN ( 1, *neG, mxGetPr( mexG ), G );
	} else {
	  memmove( G, mxGetPr( mexG ), *neG *sizeof(double));
	}
      } /* End if ( mxGetN(mexG) == 0 ) */
      mxDestroyArray( mexG );
    } /* End if ( *needG ) */
  }

  return 0;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

static void copyIfNotNAN( long rows, long cols, double from[], double to[] )
{
  double *fromi;
  double *toi;

  long i, j;
  for( i = 0; i < rows; i++ ) {
    fromi = &from[i*cols];
    toi   = &to[i*cols];
    for( j = 0; j < cols; j++ ) {
      if (!mxIsNaN(fromi[j]) ) toi[j] = fromi[j];
    }
  }
}
