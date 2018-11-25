/* Josh Griffin 9-Jul-02 */
/* Mimics file written by Mike Gertz */

#ifndef SNOPTPROBDEF
#define SNOPTPROBDEF

#pragma once

#ifndef mex_h
#include "mex.h"
#endif

#ifndef F2C_INCLUDE
#include "f2c.h"
#endif

extern char userfun[32];
extern char   myobj[32];
extern char nonlcon[32];

int sncmxAllocFuns( integer *n, integer *neG, integer iGfun[], integer jGvar[] );
int sncmxFreeFuns ( );
int sncmxAllocJac ( integer *n );
int sncmxFreeJac  ( );

int setStatus     ( integer *status_value );
int getStatus     ( integer *status_value );
int userfg_       ( integer *status,   integer *n,    doublereal x[],
		    integer    *needF, integer *neF,  doublereal F[],
		    integer    *needG, integer *neG,  doublereal G[],
		    char       *cu,    integer *lencu,
		    integer    iu[],   integer *leniu,
		    doublereal ru[],   integer *lenru );

#endif
