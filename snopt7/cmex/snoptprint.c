/* Philip E. Gill and Josh Griffin  -25-Apr-06 */

#include <stdio.h>
#include <errno.h>
#include "mexUtility.h"
#include "snoptprint.h"
#include "snoptcmex.h"

#ifndef mex_h
#include "mex.h"
#endif

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
FILE *sncmxOpenAppend( const char fileName[] )
     /* Open fileName for appending; return NULL if there is trouble    */
{
  FILE *f;

  if (fileName == NULL) fileName = "\0";
  f = fopen( fileName, "a");
  if (f == NULL)
    fprintf( stderr,
             "openAppend(\"%s\") failed: errno = %d\n",
             fileName, errno );
  return f;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
FILE *sncmxOpenReplace( const char fileName[] )
     /* Open fileName, replacing existing file; return NULL if problem  */
{
  FILE *f;

  if (fileName == NULL) fileName = "\0";
  f = fopen( fileName, "w");
  if (f == NULL)
    fprintf( stderr,
             "openReplace(\"%s\") failed: %d\n",
             fileName, errno );
  return f;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
int sncmxCloseFile( FILE *file )
     /* Close file  */
{
  int   s = 1;
  if (file == NULL) return 0;  /* Ignore this case */
  errno = 0;
  s = fclose( file );
  if (s == EOF) perror("Close failed");
  return s;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
int sncmxabort_( integer *iabort )
     /* Called by fortran routine snABORT              */
     /* On entry,    iAbort =  0                       */
     /* On exit,  if iAbort ne 0, SNOPT is terminated. */
{
  FILE *file;
  char  filename[32];
  sprintf(filename,"snabort.txt");

  file = fopen( filename, "r");
  if (file) {
    if (fscanf(file, "%d", (int*)iabort) != 1) {
      mexPrintf("\n==> Improper data inside %s:\n",filename);
      mexPrintf(  "==> First entry should be an integer.\n");
      mexPrintf(  "==> Continuing with run ...\n");
    }
    if (*iabort != 0) {
      mexPrintf("\n==> Detected nonzero value in %s", filename);
      mexPrintf(  "==> terminating run\n");
    }
  }
  return 0;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
int sncmxwritescreen_( char buffP[], integer *lenbuf )
     /* Called by fortran routine snPRNT  */
{
  if (*lenbuf > 140) mexErrMsgTxt("Print buffer too long for snPRNT");
  buffP[*lenbuf]='\0';
  mexPrintf("%s\n",buffP);
  return 0;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
int sncmxwritefile_( integer *newline, integer *UnitNumber, char buff[],
                     integer *lenbuf )
     /* Called by fortran routine snPRNT  */
{
  if (*lenbuf > 140) mexErrMsgTxt("Line is too long to be printed");
  buff[*lenbuf]='\0';

  if (*UnitNumber == snPrintUnit) {
    if (*newline)
      fprintf( snPrintFile, "\n%s\n", buff );
    else
      fprintf( snPrintFile,   "%s\n", buff );
  }
  else if (*UnitNumber == snSummaryUnit) {
    if (*newline)
      fprintf( snSummaryFile, "\n%s\n", buff );
    else
      fprintf( snSummaryFile,   "%s\n", buff );
  }

  return 0;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
int sncmxReadLine( char buff[], int nchar, FILE *stream )
     /* Reads nchar characters terminated by a '\n' from the file *stream. */
     /* On successful exit, buff[0:nchar-1] is a fortran character string, */
     /* i.e.,  buff includes trailing blanks and has no terminating '\0'.  */
{
  int   c, i = 0;
  char *cs;

  cs = buff;
  while (i < nchar && (c = fgetc(stream)) != EOF)
    if (c != '\n') {
      *cs++ = c;  i++;
    } else
      break;

  /* Flush data to the start of the next newline */
  while (c != '\n' && c != EOF)
    c = fgetc(stream);

  while (i < nchar) {
    *cs++ = ' '; i++;
  }

  return feof(stream);
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
int sncmxread_( integer *UnitNumber, char buff[], integer *nchar, integer *endfile )
     /* Called by fortran routine snPRNT  */
{
  int n, result;
  n        = *nchar;
  result   = sncmxReadLine( buff, n, snSpecsFile );
  *endfile = (integer)result;
  return 0;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
int sncmxfilestatus_( integer *scrnOK, integer *summOK, integer *prntOK )
     /* Called by fortran routine snPRNT  */
{
  *scrnOK = (integer)(screenIsOn        && callType);
  *summOK = (integer)(summaryFileIsOpen && callType);
  *prntOK = (integer)(printFileIsOpen   && callType);

  return 0;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
int sncmxSetProbName( char *Prob )
{
  Prob = "        ";

  return 0;
}
