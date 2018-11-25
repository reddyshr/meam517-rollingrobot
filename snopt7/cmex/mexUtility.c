/* Mike Gertz 2-Aug-98 */
/* Correction: 3-Oct-98 */

#include "mexUtility.h"

/* **************************************************************** *
 * void dble2int ( int n, double * da, integer * ia )
 *     Copy n items from the array double precision array da
 *     to the integer array ia.
 * **************************************************************** */
void dble2int ( int n, double * da, integer * ia )
{
  int i;
  if ( sizeof(integer) <= sizeof(double) ) {
    for ( i = 0; i < n; i++ ) {
      ia[i] = da[i];
    }
  } else {
    for ( i = n - 1; i >= 0; i-- ) {
      ia[i] = da[i];
    }
  }
}

/* **************************************************************** *
 * void int2dble( int n, integer * ia, double * da )
 *     Copy n items from the integer array ia to the double precision
 *     array da.
 * **************************************************************** */
void int2dble( int n, integer * ia, double * da )
{
  int i;
  if ( sizeof(integer) <= sizeof(double) ) {
    for ( i = n - 1; i >= 0; i-- ) {
      da[i] = ia[i];
    }
  } else {
    for ( i = 0; i < n; i++ ) {
      da[i] = ia[i];
    }
  }
}

/* ****************************************************************
 * static void errNotString( char name[8] )
 *     The argument with name "name" should have contained a string,
 *     but did not. Print an error message and abort the mex
 *     function.
 *
 *     This function is only used by assertString(), and exists
 *     to keep assertString as lightweight as possible. assertString
 *     is also a good candidate for conversion to an inline function
 *     or macro.
 * **************************************************************** */
static void errNotString( char name[8] )
{
  char errStr[24];

  sprintf( errStr, "'%.8s' must be a string.", name );
  mexErrMsgTxt( errStr );
}

/* ****************************************************************
 * void assertString( const mxArray * mex, char name[8] )
 *     If "mex" does not contain the MATLAB representation of a string,
 *     prints an error message and aborts the mex function. "name" is
 *     a symbolic name for the argument that should have contained a
 *     string. "name" will be used in displaying the error message.
 * **************************************************************** */
void assertString( const mxArray * mex, char name[8] )
{
  if( !mxIsChar( mex ) ) {
    errNotString( name );
  }
}

/* **************************************************************** *
 * static void errNotScalar( char name[8] )
 *     The argument with name "name" should have contained a scalar
 *     value, but did not. Print an error message and abort the mex
 *     function.
 * **************************************************************** */
static void errNotScalar( char name[8] )
{
  char errStr[40];

  sprintf ( errStr, "Parameter '%.8s' must be scalar", name );
  mexErrMsgTxt( errStr );
}

/* **************************************************************** *
 * void assertScalar( const mxArray * mex, char name[8] )
 *     If "mex" does not contain the MATLAB representation of a scalar
 *     value, prints an error message and aborts the mex
 *     function. "name" is a symbolic name for the argument that
 *     should have contained a scalar value. "name" will be used in
 *     displaying the error message.
 * **************************************************************** */
void assertScalar( const mxArray * mex, char name[8] )
{
  if ( 1 != mxGetM( mex ) | 1 != mxGetN( mex ) ) {
    errNotScalar( name );
  }
}

/* **************************************************************** *
 * static void errRowDim( char name[8] )
 *     The argument with name "name" should have had dim rows, but
 *     instead had m rows. Display an error message and abort the mex
 *     function.
 * **************************************************************** */

static void errRowDim( int dim, int m, char name[8] )
{
  char errStr[64];

  sprintf( errStr, "Expected %4d rows for parameter '%.8s'. Got %4d",
	   dim, name, m );
  mexErrMsgTxt( errStr );
}

/* **************************************************************** *
 * void assertRowDim( const mxArray * mex, int dim, char name[8] )
 *     If "mex" does not contain a matrix with "dim" rows, prints an
 *     error message and aborts the mex function. "name" is a symbolic
 *     name for the argument that should have contained a scalar
 *     value. "name" will be used in displaying the error message.
 * **************************************************************** */
void assertRowDim( const mxArray * mex, int dim, char name[8] )
{
  int m;

  m = mxGetM( mex );
  if ( dim != m ) {
    errRowDim( dim, m, name );
  }
}

/* **************************************************************** *
 * static void errColDim( int dim, int m, char name[8] )
 *     See errRowDim()
 * **************************************************************** */
static void errColDim( int dim, int m, char name[8] )
{
  char errStr[64];

  sprintf( errStr, "Expected %4d columns for parameter '%.8s'. Got %4d",
	   dim, name, m );
  mexErrMsgTxt( errStr );
}

/* **************************************************************** *
 * void assertColDim( const mxArray * mex, int dim, char name[8] )
 *     See assertRowDim()
 * **************************************************************** */
void assertColDim( const mxArray * mex, int dim, char name[8] )
{
  int n;
  n = mxGetN( mex );
  if ( dim != n ) {
    errColDim( dim, n, name );
  }
}

/* **************************************************************** *
 * static void errSparse( char name[8] )
 * See errRowDim()
 * **************************************************************** */
static void errSparse( char name[8] )
{
  char errStr[64];
  sprintf( errStr, "Expected sparse matrix for parameter '%.8s'.",
	   name );
  mexErrMsgTxt( errStr );
}

/* ************************************************************ *
 * void assertSparse( const mxArray * mex, char name[8] )
 * See assertRowDim()
 * ************************************************************ */
void assertSparse( const mxArray * mex, char name[8] )
{
  if( !mxIsSparse(mex) ){
    errSparse( name );
  }
}
