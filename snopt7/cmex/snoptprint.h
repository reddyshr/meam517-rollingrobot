
#ifndef SNOPTCMEXPRINT
#define SNOPTCMEXPRINT

#pragma once

#ifndef F2C_INCLUDE
#include "f2c.h"
#endif

extern int   sncmxSetProbName  ( char *buff );
extern int   sncmxRead         ( char *buff, int nchar, FILE *stream );
extern int   sncmxReadLine     ( char *buff, int nchar, FILE *stream );
extern int   sncmxCloseFile    ( FILE *f );
extern FILE *sncmxOpenAppend   ( const char fileName[] );
extern FILE *sncmxOpenReplace  ( const char fileName[] );

extern int   sncmxwritescreen_ ( char buff[], integer *lenbuf );
extern int   sncmxwritefile_   ( integer *mode, integer *UnitNumber, char buff[], integer *lenbuf );
extern int   sncmxabort_       ( integer *iabort );
extern int   sncmxfilestatus_  ( integer *scrnOK, integer *summOK, integer *prntOK );
extern int   sncmxprob_        ( char *prob, ftnlen prob_len );
extern int   sncmxread_        ( integer *Unit, char buff[], integer *nchar, integer *endfile );

#endif
