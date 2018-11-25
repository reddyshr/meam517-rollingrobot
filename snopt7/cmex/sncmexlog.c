/* ./cmex/sncmexlog.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*     File  sncmexlog.f */

/*     snPRNT  snAbort    snPROB */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Subroutine */ int snread_(integer *unitno, char *string, integer *nchar, 
	integer *endfile, ftnlen string_len)
{
    extern /* Subroutine */ int sncmxread_(integer *, char *, integer *, 
	    integer *, ftnlen);

/*     ================================================================== */
/*     snREAD reads a string of length nchar from file  unitno. */

/*     30 Apr 2006: First version of snREAD. */
/*     30 Apr 2006: Matlab version. */
/*     ================================================================== */
    sncmxread_(unitno, string, nchar, endfile, string_len);
    return 0;
} /* snread_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine snREAD */
/* Subroutine */ int snprnt_(integer *mode, char *string, integer *iw, 
	integer *leniw, ftnlen string_len)
{
    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    extern /* Subroutine */ int sncmxfilestatus_(integer *, integer *, 
	    integer *);
    static integer screenok, m;
    extern /* Subroutine */ int sncmxwritescreen_(char *, integer *, ftnlen);
    static integer summaryok;
    static char buff[140];
    static integer isumm;
    extern /* Subroutine */ int s1trim_(char *, integer *, ftnlen);
    static integer length, iprint, lvlsys;
    extern /* Subroutine */ int sncmxwritefile_(integer *, integer *, char *, 
	    integer *, ftnlen);
    static integer newline, printok;

/*     ================================================================== */
/*     snPRNT  prints a trimmed form of "string" on various files. */
/*     If mode = 0,      nothing is output. */
/*     If mode = 1,      string is output to iPrint. */
/*     If mode = 2,      string is output to iSumm. */
/*     If mode = 3 or 4, string is output to iPrint and iSumm. */
/*     If mode = 4,      string is output to the screen. */
/*                       This mode is intended for error messages. */
/*     If mode = 5,      string is output to iStdo (standard output) */
/*                       This mode is to be used when the elements of */
/*                       the integer work array iw cannot be trusted. */

/*     mode 11-15 are the same as mode 1-5 with blank line before output. */

/*     If mode > 15 then nothing is printed unless  lvlSys > 0. */
/*     mode 21-25 are the same as mode 1-5 */
/*     mode 31-35 are the same as mode 11-15 */

/*     25 Sep 2002: First version of snPRNT. */
/*     31 Jul 2003: mode 11-14 added.  form introduced. */
/*     27 Dec 2003: mode 5 added to allow printing before iw is set. */
/*     12 Mar 2004: s1trim called to trim the string. */
/*     22 Jun 2004: System printing option added. */
/*     14 Oct 2004: Matlab version of snPRNT. */
/*     30 Apr 2006: Files opened and closed in C. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --iw;

    /* Function Body */
    lvlsys = iw[71];
/* > 0   => print system info */
    newline = 0;
    m = 0;
    if (*mode <= 0) {
/*        Relax */
    } else if (*mode < 10) {
	m = *mode;
    } else if (*mode < 20) {
/* Blank line first */
	m = *mode - 10;
	newline = 1;
    } else if (lvlsys > 0) {
/* Print system Info */
	if (*mode < 30) {
	    m = *mode - 20;
	} else {
	    m = *mode - 30;
	    newline = 1;
	}
    }
    if (m > 0) {
	sncmxfilestatus_(&screenok, &summaryok, &printok);
/*        length = len_trim(string)     ! An F90 intrinsic */
	s1trim_(string, &length, string_len);
/* The F77 equivalent */
	s_copy(buff, string, (ftnlen)140, string_len);
	if (m == 5) {
	    sncmxwritescreen_(buff, &length, (ftnlen)140);
	} else {
	    iprint = iw[12];
/* Print file */
	    isumm = iw[13];
/* Summary file */
	    if (m == 1 || m >= 3) {
		if (printok > 0) {
		    sncmxwritefile_(&newline, &iprint, buff, &length, (ftnlen)
			    140);
		}
	    }
	    if (m == 2 || m >= 3) {
		if (screenok > 0) {
		    sncmxwritescreen_(buff, &length, (ftnlen)140);
		}
		if (summaryok > 0) {
		    sncmxwritefile_(&newline, &isumm, buff, &length, (ftnlen)
			    140);
		}
	    }
	    if (m == 4) {
		sncmxwritescreen_(buff, &length, (ftnlen)140);
	    }
	}
    }
    return 0;
} /* snprnt_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine snPRNT */
/* Subroutine */ int snabort_(integer *iabort, integer *info, integer *htype, 
	logical *ktcond, integer *mjrprt, integer *minimz, integer *m, 
	integer *maxs, integer *n, integer *nb, integer *nncon0, integer *
	nncon, integer *nnobj0, integer *nnobj, integer *ns, integer *itn, 
	integer *nmajor, integer *nminor, integer *nswap, doublereal *condhz, 
	integer *iobj, doublereal *sclobj, doublereal *objadd, doublereal *
	fmrt, doublereal *pennrm, doublereal *step, doublereal *prinf, 
	doublereal *duinf, doublereal *vimax, doublereal *virel, integer *hs, 
	integer *ne, integer *nlocj, integer *locj, integer *indj, doublereal 
	*jcol, integer *negcon, doublereal *ascale, doublereal *bl, 
	doublereal *bu, doublereal *fcon, doublereal *gcon, doublereal *gobj, 
	doublereal *ycon, doublereal *pi, doublereal *rc, doublereal *rg, 
	doublereal *x, char *cu, integer *lencu, integer *iu, integer *leniu, 
	doublereal *ru, integer *lenru, char *cw, integer *lencw, integer *iw,
	 integer *leniw, doublereal *rw, integer *lenrw, ftnlen cu_len, 
	ftnlen cw_len)
{
    extern /* Subroutine */ int sncmxabort_(integer *);

/*     ================================================================== */
/*     snAbort  is called every major iteration. */
/*     If iAbort > 0 on exit, the run is terminated. */
/*     By specifying a custom version of snSTOP, the user can arrange for */
/*     snopt to be terminated at any given major iteration. */

/*     14 Oct 2004: First version of   snAbort. */
/*     01 Sep 2007: Parameter list expanded. */
/*     ================================================================== */
    /* Parameter adjustments */
    --info;
    --ktcond;
    --pi;
    --rg;
    --x;
    --rc;
    --bu;
    --bl;
    --ascale;
    --hs;
    --ycon;
    --fcon;
    --gobj;
    --jcol;
    --indj;
    --locj;
    --gcon;
    cu -= 8;
    --iu;
    --ru;
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    *iabort = 0;
    sncmxabort_(iabort);
    return 0;
} /* snabort_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine snAbort */
/* Subroutine */ int snprob_(char *prob, ftnlen prob_len)
{
    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

/*     ================================================================== */
/*     Assigns an empty problem name. */

/*     31 Dec 2002: First version of snPROB */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    s_copy(prob, "        ", (ftnlen)8, (ftnlen)8);
    return 0;
} /* snprob_ */

