/* Mike Gertz 2-Aug-98 */
#ifndef MEXUTILITY
#define MEXUTILITY

#ifndef mex_h
#include "mex.h"
#endif

#ifndef F2C_INCLUDE
#include "f2c.h"
#endif

#pragma once

void dble2int( int n,  double * da, integer * ia );
void int2dble( int n, integer * ia,  double * da );

void assertScalar( const mxArray * mex, char name[8] );
void assertRowDim( const mxArray * mex, int dim, char name[8] );
void assertColDim( const mxArray * mex, int dim, char name[8] );
void assertString( const mxArray * mex, char name[8] );
void assertSparse( const mxArray * mex, char name[8] );

#endif
