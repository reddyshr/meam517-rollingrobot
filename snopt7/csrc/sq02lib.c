/* ../snopt7/src/sq02lib.f -- translated by f2c (version 20100827).
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

/* Table of constant values */

static integer c__2 = 2;
static integer c__11 = 11;
static integer c__1 = 1;
static integer c__12 = 12;
static integer c__130 = 130;
static integer c__13 = 13;
static logical c_true = TRUE_;
static logical c_false = FALSE_;
static integer c__0 = 0;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*     File  sq02lib.f */

/*     sqTitle   sqInit   sqSpec   sqHx     sqMem    sqlog */

/*     sqSet     sqSeti       sqSetr */
/*     sqGet     sqGetc       sqGeti   sqGetr */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Subroutine */ int sqtitle_(char *title, ftnlen title_len)
{
    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

/*     ================================================================== */
/*     sqtitl sets the title. */
/*     ================================================================== */
    s_copy(title, "S Q O P T  7.4-1.2  (Feb 2015)", (ftnlen)30, (ftnlen)30);
/* ---------------123456789|123456789|123456789|-------------------------- */
    return 0;
} /* sqtitle_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine sqTitle */
/* Subroutine */ int sqinit_(integer *iprint, integer *isumm, char *cw, 
	integer *lencw, integer *iw, integer *leniw, doublereal *rw, integer *
	lenrw, ftnlen cw_len)
{
    /* Initialized data */

    static char dashes[30] = "==============================";

    /* System generated locals */
    address a__1[2];
    integer i__1[2];
    char ch__1[39], ch__2[31];

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen), s_cat(char *,
	     char **, integer *, integer *, ftnlen);

    /* Local variables */
    extern /* Subroutine */ int s3unsetall_(char *, integer *, integer *, 
	    integer *, doublereal *, integer *, ftnlen);
    static char str[80], str2[80], title[30];
    static integer istdo;
    extern /* Subroutine */ int s1init_(char *, integer *, integer *, 
	    doublereal *, integer *, ftnlen);
    static integer ispecs, inform__;
    static char solver[6];
    extern /* Subroutine */ int snwrap_(integer *, char *, char *, char *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen), snprnt_(integer *, 
	    char *, integer *, integer *, ftnlen);
    extern integer s1outpt_(void);
    extern /* Subroutine */ int sqtitle_(char *, ftnlen);

/*     ================================================================== */
/*     sqInit  is called by the user to do the following: */
/*     1. Open default files (Print, Summary). */
/*     2. Initialize title. */
/*     3. Set options to default values. */

/*     15 Nov 1991: First version. */
/*     14 Jul 1997: Thread-safe version. */
/*     21 Mar 1997: First version based on snopt routine snInit */
/*     14 Jul 1997: Thread-safe version. */
/*     02 Oct 1997: Character workspace added. */
/*     15 Oct 2003: snEXIT and snPRNT added. */
/*     18 Jun 2007: First 500 elements of cw, iw and rw are initialized */
/*     18 Jun 2008: Global call-status values added. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/* start of SNOPT part of rw */
/* end   of SNOPT part of rw */
/* start of SNOPT part of iw */
/* end   of SNOPT part of iw */
/* start of SNOPT part of cw */
/* end   of SNOPT part of cw */
/* # nonlinear Jac, variables */
/* # variables in gObj */
/* nonlinear constraints */
/* nonlinear vars */
/* Timing level */
/* QP user-routine call-status */
/* NP user-routine call-status */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
/*     ------------------------------------------------------------------ */
    s_copy(solver, "SQINIT", (ftnlen)6, (ftnlen)6);
    if (*lencw < 500 || *leniw < 500 || *lenrw < 500) {
/* -------------------------------------------------------------- */
/* Not enough workspace to do ANYTHING! */
/* Print and exit without accessing the work arrays. */
/* -------------------------------------------------------------- */
	inform__ = 81;
/* Work arrays must have at least 500 elements */
	snwrap_(&inform__, solver, str, str2, &iw[1], leniw, (ftnlen)6, (
		ftnlen)80, (ftnlen)80);
	goto L999;
    }
/* ----------------------------------------------------------------- */
/* Initialize cw, iw, rw so that they may be copied safely. */

/* This also sets the options to a specific "undefined" state. */
/* snopt  will check the options later and maybe print them. */
/* ----------------------------------------------------------------- */
    s3unsetall_(cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)8);
/* Writing concatenation */
    i__1[0] = 6, a__1[0] = solver;
    i__1[1] = 2, a__1[1] = "  ";
    s_cat(cw + 8, a__1, i__1, &c__2, (ftnlen)8);
/* ----------------------------------------------------------------- */
/* Initialize some default values. */
/* ----------------------------------------------------------------- */
    ispecs = 0;
    istdo = s1outpt_();
    iw[10] = istdo;
/* Standard Output */
    iw[11] = ispecs;
    iw[12] = *iprint;
/* Print file */
    iw[13] = *isumm;
/* Summary file */
    iw[6] = 500;
    iw[4] = 500;
    iw[2] = 500;
    iw[7] = *lencw;
    iw[5] = *leniw;
    iw[3] = *lenrw;
/* ----------------------------------------------------------------- */
/* These dimensions need to be initialized for an MPS run. */
/* ----------------------------------------------------------------- */
    iw[23] = 0;
    iw[21] = 0;
    iw[22] = 0;
    iw[24] = 0;
    sqtitle_(title, (ftnlen)30);
    s1init_(title, &iw[1], leniw, &rw[1], lenrw, (ftnlen)30);
/* Writing concatenation */
    i__1[0] = 9, a__1[0] = "         ";
    i__1[1] = 30, a__1[1] = dashes;
    s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)39);
    snprnt_(&c__11, ch__1, &iw[1], leniw, (ftnlen)39);
/* Writing concatenation */
    i__1[0] = 9, a__1[0] = "         ";
    i__1[1] = 30, a__1[1] = title;
    s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)39);
    snprnt_(&c__1, ch__1, &iw[1], leniw, (ftnlen)39);
/* Writing concatenation */
    i__1[0] = 9, a__1[0] = "         ";
    i__1[1] = 30, a__1[1] = dashes;
    s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)39);
    snprnt_(&c__1, ch__1, &iw[1], leniw, (ftnlen)39);
/* Writing concatenation */
    i__1[0] = 1, a__1[0] = " ";
    i__1[1] = 30, a__1[1] = dashes;
    s_cat(ch__2, a__1, i__1, &c__2, (ftnlen)31);
    snprnt_(&c__12, ch__2, &iw[1], leniw, (ftnlen)31);
/* Writing concatenation */
    i__1[0] = 1, a__1[0] = " ";
    i__1[1] = 30, a__1[1] = title;
    s_cat(ch__2, a__1, i__1, &c__2, (ftnlen)31);
    snprnt_(&c__2, ch__2, &iw[1], leniw, (ftnlen)31);
/* Writing concatenation */
    i__1[0] = 1, a__1[0] = " ";
    i__1[1] = 30, a__1[1] = dashes;
    s_cat(ch__2, a__1, i__1, &c__2, (ftnlen)31);
    snprnt_(&c__2, ch__2, &iw[1], leniw, (ftnlen)31);
/* ----------------------------------------------------------------- */
/* Initialize some global values. */
/* ----------------------------------------------------------------- */
    iw[235] = -11111;
    iw[236] = -11111;
    iw[182] = 3;
L999:
    return 0;
} /* sqinit_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine sqInit */
/* Subroutine */ int sqspec_(integer *ispecs, integer *iexit, char *cw, 
	integer *lencw, integer *iw, integer *leniw, doublereal *rw, integer *
	lenrw, ftnlen cw_len)
{
    /* System generated locals */
    address a__1[2];
    integer i__1[2];

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen), s_cat(char *,
	     char **, integer *, integer *, ftnlen);

    /* Local variables */
    static char str[80], str2[80];
    extern /* Subroutine */ int s3opt_();
    static integer calls, isumm;
    extern /* Subroutine */ int s3file_(integer *, integer *, integer *, U_fp,
	     char *, integer *, integer *, integer *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, ftnlen, ftnlen);
    static integer iprint;
    static char solver[6];
    extern /* Subroutine */ int snwrap_(integer *, char *, char *, char *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen);
    static integer errors;

/*     ================================================================== */
/*     sqSpec  is called by the user to read the Specs file. */

/*     07 Feb 1998: First version. */
/*     01 Aug 2003: s3file now has a "title" parameter.  Use ' '. */
/*     27 Oct 2003: Current version of sqSpec. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    s_copy(solver, "SQSPEC", (ftnlen)6, (ftnlen)6);
    if (*lencw < 500 || *leniw < 500 || *lenrw < 500) {
/*        --------------------------------------------------------------- */
/*        Not enough workspace to do ANYTHING! */
/*        Print and exit without accessing the work arrays. */
/*        --------------------------------------------------------------- */
	*iexit = 81;
/* Work arrays must have at least 500 elements */
	snwrap_(iexit, solver, str, str2, &iw[1], leniw, (ftnlen)6, (ftnlen)
		80, (ftnlen)80);
	goto L999;
    }
/* Writing concatenation */
    i__1[0] = 6, a__1[0] = solver;
    i__1[1] = 2, a__1[1] = "  ";
    s_cat(cw + 8, a__1, i__1, &c__2, (ftnlen)8);
    if (*ispecs <= 0) {
	*iexit = 131;
	goto L800;
    }
    iw[11] = *ispecs;
    iprint = iw[12];
    isumm = iw[13];
    *iexit = 0;
    calls = 1;
/*     ------------------------------------------------------------------ */
/*     Read the Specs file. */
/*     sqopt  will check the options later and maybe print them. */
/*     ------------------------------------------------------------------ */
    s3file_(iexit, &calls, ispecs, (U_fp)s3opt_, " ", &iprint, &isumm, &
	    errors, cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)1, (
	    ftnlen)8);
L800:
    if (*iexit == 0) {
	*iexit = 101;
/* SPECS file read successfully */
    }
L999:
    return 0;
} /* sqspec_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine sqSpec */
/* Subroutine */ int sqhx_(S_fp usrhx, integer *nnh, integer *neh, integer *
	nloch, integer *loch, integer *indh, doublereal *hcol, doublereal *x, 
	doublereal *hx, integer *sqstat, char *cu, integer *lencu, integer *
	iu, integer *leniu, doublereal *ru, integer *lenru, char *cw, integer 
	*lencw, integer *iw, integer *leniw, doublereal *rw, integer *lenrw, 
	ftnlen cu_len, ftnlen cw_len)
{
    static integer lxscaled, lvlscale;
    extern /* Subroutine */ int s5callstatus_(integer *, integer *, integer *)
	    , ddscl_(integer *, doublereal *, integer *, doublereal *, 
	    integer *), dcopy_(integer *, doublereal *, integer *, doublereal 
	    *, integer *);
    static logical scaled;
    static integer status, lscales;

/*     ================================================================== */
/*     sqHx  computes the user-defined product  Hx  and scales it. */

/*     15 Mar 1999: First   version of sqHx */
/*     16 Jun 2008: Call-status implemented correctly. */
/*     05 Nov 2014: neH, indH, locH, Hcol added as arguments. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --hx;
    --x;
    --hcol;
    --indh;
    --loch;
    cu -= 8;
    --iu;
    --ru;
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    lvlscale = iw[75];
/* scale option */
    lscales = iw[296];
/* scales(nb)  = row and column scales */
    lxscaled = iw[302];
/* xScaled(n)  = copy of scaled x(nnL) */
    scaled = lvlscale > 0;
/* Determine the user-function call-status. */
    s5callstatus_(&status, &iw[1], leniw);
    if (scaled) {
	dcopy_(nnh, &x[1], &c__1, &rw[lxscaled], &c__1);
	ddscl_(nnh, &rw[lscales], &c__1, &x[1], &c__1);
    }
    (*usrhx)(nnh, &x[1], &hx[1], &status, cu + 8, lencu, &iu[1], leniu, &ru[1]
	    , lenru, (ftnlen)8);
    if (scaled) {
	dcopy_(nnh, &rw[lxscaled], &c__1, &x[1], &c__1);
	ddscl_(nnh, &rw[lscales], &c__1, &hx[1], &c__1);
    }
    return 0;
} /* sqhx_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine sqHx */
/* Subroutine */ int sqmem_(integer *iexit, integer *m, integer *n, integer *
	nea, integer *lencobj, integer *ncolh, integer *mincw, integer *miniw,
	 integer *minrw, char *cw, integer *lencw, integer *iw, integer *
	leniw, doublereal *rw, integer *lenrw, ftnlen cw_len)
{
    /* System generated locals */
    address a__1[2];
    integer i__1[2];

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen), s_cat(char *,
	     char **, integer *, integer *, ftnlen);

    /* Local variables */
    static logical printmem;
    extern /* Subroutine */ int s5defaults_(integer *, integer *, integer *, 
	    integer *, char *, integer *, integer *, integer *, doublereal *, 
	    integer *, ftnlen);
    static integer nkx;
    static char str[80], str2[80];
    static integer lenr, maxr, maxs;
    extern /* Subroutine */ int s2mem_(integer *, logical *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *), s5map_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *);
    static integer nnjac, nnobj, nncon, maxcw;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer maxiw;
    extern /* Subroutine */ int icopy_(integer *, integer *, integer *, 
	    integer *, integer *);
    static integer maxrw;
    extern /* Subroutine */ int s2bmap_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *);
    static integer llencw;
    extern /* Subroutine */ int chcopy_(integer *, char *, integer *, char *, 
	    integer *, ftnlen, ftnlen);
    static integer inform__, lleniw, llenrw, liwest;
    static char usercw[8*130];
    static integer nextcw;
    static char solver[6];
    extern /* Subroutine */ int snwrap_(integer *, char *, char *, char *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen);
    static integer nextiw, useriw[130], lrwest, nextrw;
    static doublereal userrw[130];

/*     ================================================================== */
/*     sqMem   estimates the memory requirements for sqopt, */
/*     using the values: */
/*        m       , n    , ne */
/*        lencObj , ncolH */

/*     These values are used to compute the minimum required storage: */
/*     mincw, miniw, minrw. */

/*     Note: */
/*     1. All default parameters must be set before calling sqMem, */
/*        since some values affect the amount of memory required. */

/*     2. The arrays rw and iw hold  constants and work-space addresses. */
/*        They must have dimension at least 500. */

/*     3. This version of sqMem does not allow user accessible */
/*        partitions of cw, iw and rw. */

/*     01 May 1998: First version. */
/*     31 Jul 2003: snPRNT adopted. */
/*     25 Oct 2003: Current version of sqMem. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    s_copy(solver, "SQMEM ", (ftnlen)6, (ftnlen)6);
    *iexit = 0;
    if (*lencw < 500 || *leniw < 500 || *lenrw < 500) {
/*        --------------------------------------------------------------- */
/*        Not enough workspace to do ANYTHING! */
/*        Print and exit without accessing the work arrays. */
/*        --------------------------------------------------------------- */
	*iexit = 81;
/* Work arrays must have at least 500 elements */
	snwrap_(iexit, solver, str, str2, &iw[1], leniw, (ftnlen)6, (ftnlen)
		80, (ftnlen)80);
	goto L999;
    }
/* Writing concatenation */
    i__1[0] = 6, a__1[0] = solver;
    i__1[1] = 2, a__1[1] = "  ";
    s_cat(cw + 8, a__1, i__1, &c__2, (ftnlen)8);
/*     Save the user's option choices  (weird choices get overwritten). */
    chcopy_(&c__130, cw + 408, &c__1, usercw, &c__1, (ftnlen)8, (ftnlen)8);
    icopy_(&c__130, &iw[51], &c__1, useriw, &c__1);
    dcopy_(&c__130, &rw[51], &c__1, userrw, &c__1);
/*     Assign fake values for lencw, leniw, lenrw. */
/*     This will force s5Mem to estimate the memory requirements. */
    llenrw = 500;
    lleniw = 500;
    llencw = 500;
/*     An obligatory call to sqInit has `undefined' all options. */
/*     Check the user-defined values and assign undefined values. */
/*     s5Defaults needs various problem dimensions in iw. */
    nncon = 0;
/* Not used in sqopt */
    nnobj = 0;
/* ditto */
    nnjac = 0;
/* ditto */
    iw[15] = *n;
/* copy of the number of columns */
    iw[16] = *m;
/* copy of the number of rows */
    iw[17] = *nea;
/* copy of the number of nonzeros in Acol */
    iw[21] = nnjac;
/* # nonlinear Jacobian variables */
    iw[22] = nnobj;
/* # variables in gObj */
    iw[23] = nncon;
/* # of nonlinear constraints */
    iw[26] = *lencobj;
/* length of QP constant vector */
    iw[27] = *ncolh;
/* # QP Hessian columns */
    s5defaults_(m, n, lencobj, ncolh, cw + 8, &llencw, &iw[1], &lleniw, &rw[1]
	    , &llenrw, (ftnlen)8);
    nextcw = 501;
    nextiw = 501;
    nextrw = 501;
    maxcw = *lencw;
    maxiw = *leniw;
    maxrw = *lenrw;
    maxr = iw[52];
/* max columns of R. */
    maxs = iw[53];
/* max # of superbasics */
    lenr = maxr * (maxr + 1) / 2 + (maxs - maxr);
    nkx = *n + *m;
    s5map_(m, n, &nkx, lencobj, ncolh, &lenr, &maxr, &maxs, &nextcw, &nextiw, 
	    &nextrw, &iw[1], leniw);
    s2bmap_(m, n, nea, &maxs, &nextiw, &nextrw, &maxiw, &maxrw, &liwest, &
	    lrwest, &iw[1], leniw);
    printmem = FALSE_;
/* Print all messages in s2Mem */
    s2mem_(&inform__, &printmem, &liwest, &lrwest, &nextcw, &nextiw, &nextrw, 
	    &maxcw, &maxiw, &maxrw, &llencw, &lleniw, &llenrw, mincw, miniw, 
	    minrw, &iw[1]);
/*     mincw = mincw */
    *miniw = liwest;
    *minrw = lrwest;
/*     Restore the user's choices of options. */
    chcopy_(&c__130, usercw, &c__1, cw + 408, &c__1, (ftnlen)8, (ftnlen)8);
    icopy_(&c__130, useriw, &c__1, &iw[51], &c__1);
    dcopy_(&c__130, userrw, &c__1, &rw[51], &c__1);
/*     Print the exit conditions. */
    if (*iexit == 0) {
	*iexit = 104;
/* memory requirements estimated */
    }
    snwrap_(iexit, solver, str, str2, &iw[1], leniw, (ftnlen)6, (ftnlen)80, (
	    ftnlen)80);
L999:
    return 0;
} /* sqmem_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine sqMem */
/* Subroutine */ int sqlog_(integer *probtype, char *probtag, logical *
	elastic, logical *gotr, logical *firstfeas, logical *feasible, 
	integer *m, integer *mbs, integer *nnh, integer *ns, integer *jsq, 
	integer *jbr, integer *jsr, integer *linesp, integer *liness, integer 
	*itn, integer *itqp, integer *kprc, integer *lvlobje, doublereal *
	pivot, doublereal *step, integer *ninf, doublereal *sinf, integer *
	ninfe, doublereal *sinfe, doublereal *wtinf, integer *nonopt, 
	doublereal *objprt, doublereal *condzhz, doublereal *djqprt, 
	doublereal *rgnorm, integer *kbs, doublereal *xbs, integer *iw, 
	integer *leniw, ftnlen probtag_len)
{
    /* Format strings */
    static char fmt_8010[] = "(\002 Itn\002,i7,\002: Feasible \002,a)";
    static char fmt_8030[] = "(\002 Itn\002,i7,\002: Elastic Phase 2 -- mini"
	    "mizing\002,\002 elastic variables\002)";
    static char fmt_8040[] = "(\002 Itn\002,i7,\002: Elastic Phase 2 -- mini"
	    "mizing\002,\002 obj + weighted elastics\002)";
    static char fmt_3000[] = "(1p,i7,i3,e9.1,3i7,2e9.1,i7,e16.8,i8,i4,e8.1,i"
	    "6,e8.1,i7)";
    static char fmt_5000[] = "(1p,i7,e9.1,1x,i7,e15.7,i6,e10.1,i7)";
    static char fmt_5010[] = "(1p,i7,e9.1,1x,i7,e15.7)";
    static char fmt_6000[] = "(i7,g17.8)";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static logical printlog, printsum;
    static integer k;
    static logical printhead, printheadp, printheads;
    static integer ncp;
    static char str[80];
    static integer jbrn, lenl, lenu, itns, jsqn, jsrn;
    static char buffp[138], buffs[80];
    static integer cgitn, width;
    extern integer s2varn_(integer *, integer *, integer *);
    static integer lprdbg, qpmode, numinf;
    static doublereal suminf;
    static logical newset;
    static integer minors;
    static doublereal mxwdth;
    static integer printp, prints;
    static doublereal rmxint;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);
    extern integer s1intmx_(void);

    /* Fortran I/O blocks */
    static icilist io___67 = { 0, str, 0, fmt_8010, 80, 1 };
    static icilist io___68 = { 0, str, 0, fmt_8030, 80, 1 };
    static icilist io___69 = { 0, str, 0, fmt_8040, 80, 1 };
    static icilist io___77 = { 0, buffp, 0, fmt_3000, 138, 1 };
    static icilist io___78 = { 0, buffp, 0, fmt_3000, 138, 1 };
    static icilist io___79 = { 0, buffp, 0, fmt_3000, 138, 1 };
    static icilist io___80 = { 0, buffp, 0, fmt_3000, 138, 1 };
    static icilist io___81 = { 0, buffs, 0, fmt_5000, 80, 1 };
    static icilist io___82 = { 0, buffs, 0, fmt_5000, 80, 1 };
    static icilist io___83 = { 0, buffs, 0, fmt_5010, 80, 1 };
    static icilist io___84 = { 0, buffs, 0, fmt_5010, 80, 1 };
    static icilist io___86 = { 0, buffp, 0, fmt_6000, 138, 1 };
    static icilist io___87 = { 0, buffp, 0, fmt_6000, 138, 1 };


/*     ================================================================== */
/*     sqLog  prints the LP/QP iteration log. */

/*     mBS = m + maxS  .ge.  m + nS. */

/*     mnrHdP is  1  if a new heading is required for some reason other */
/*     than frequency (e.g., after a basis factorization). */

/*     The output consists of a number of ``sections'' of one line */
/*     summaries, with each section preceded by a header message. */
/*     linesP and linesS count the number of lines remaining to be */
/*     printed in each section of the print and summary files */
/*     respectively.   They too may force a new heading. */

/*     01 Dec 1991: First version based on Minos routine m5log. */
/*     17 Nov 2000: First version for SQP minor itns. */
/*     28 Dec 2000: Row and column permutations added. */
/*     01 Aug 2003: cgItn added to Print log. */
/*     26 Mar 2005: Reordered and sparsified. */
/*     20 Dec 2005: LP/QP Headings made consistent. */
/*     06 Dec 2014: nonOpt added as an argument. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/* >0 => Minor heading for Log */
/*     ------------------------------------------------------------------ */
/* >0 => Minor heading for Summary */
    /* Parameter adjustments */
    --xbs;
    --kbs;
    --iw;

    /* Function Body */
    lprdbg = iw[85];
/* > 0    => private debug print */
    ncp = iw[176];
/* no. of LU compressions */
    lenl = iw[173];
/* size of current  L */
    lenu = iw[174];
/* size of current  U */
    qpmode = iw[208];
/* Current QP solver   (based on QPslvr) */
    printp = iw[218];
/* (on/off) current log     file output status */
    prints = iw[219];
/* (on/off) current summary file output status */
    cgitn = iw[387];
/* symmlq itns for the last QP minor itn */
    printlog = printp == 1;
    printsum = prints == 1;
    mxwdth = 1e7;
/* Integers printed i7 */
    rmxint = (doublereal) s1intmx_();
/* Largest integer without overflow */
    if (mxwdth > rmxint) {
	mxwdth = rmxint;
    }
    width = (integer) mxwdth;
    itns = *itn % width;
    minors = *itqp % width;
    s_copy(buffp, " ", (ftnlen)138, (ftnlen)1);
    s_copy(buffs, " ", (ftnlen)80, (ftnlen)1);
    if (*feasible) {
	numinf = *ninfe;
	suminf = *sinfe;
    } else {
	numinf = *ninf;
	suminf = *sinf;
    }
/* If  newly feasible, print something. */
    if (*firstfeas) {
	if (! (*elastic)) {
/* Constraints feasible in Normal mode. */
/* Print a message. */
/* probTag is one of the following: */
/* probTag = 'QP problem' */
/* probTag = 'LP problem' */
/* probTag = 'QP subproblem' */
/* probTag = 'Linear constraints' */
	    if (*probtype != 4 && *probtype != 0 && *probtype != 3) {
		s_wsfi(&io___67);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		do_fio(&c__1, probtag, (ftnlen)20);
		e_wsfi();
		snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)80);
	    }
	} else {
/* Elastic mode */
/* Elastic Phase 1 has completed. */
	    if (*lvlobje == 2) {
/* Infinite weight on sumInf. */
/* Minimize the infeasible elastics. */
		s_wsfi(&io___68);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)80);
	    } else if (*lvlobje == 1) {
/* Finite nonzero weight on sumInf */
/* Minimize a weighted objective. */
		s_wsfi(&io___69);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)80);
	    }
	}
	iw[223] = 1;
/* Print the header to the print   file */
	iw[225] = 1;
/* Print the header to the summary file */
    }
    printheadp = iw[223] > 0;
    printheads = iw[225] > 0;
/*     ================================================================== */
/*     Print (log) file */
/*     ================================================================== */
    if (printlog) {
/* Print one line to the print file */
	newset = *linesp == 0;
	printhead = printheadp || newset;
	if (printhead) {
	    iw[223] = 0;
	    *linesp = 40;
	}
	--(*linesp);
	jsqn = s2varn_(jsq, leniw, &iw[1]);
	jsrn = s2varn_(jsr, leniw, &iw[1]);
	jbrn = s2varn_(jbr, leniw, &iw[1]);
	if (*nnh > 0) {
	    if (printhead) {
		s_copy(buffp, "    Itn pp  QP mult   +SBS   -SBS    -BS     "
			"Step    Pivot NonOpt    QP objective     L+U ncp  rg"
			"Norm    nS condZHZ", (ftnlen)138, (ftnlen)115);
		if (*feasible && *elastic) {
		    if (*lvlobje == 1) {
			s_copy(buffp + 68, "Elastic QPobj", (ftnlen)13, (
				ftnlen)13);
		    }
		    if (*lvlobje == 2) {
			s_copy(buffp + 68, "      SumInfE", (ftnlen)13, (
				ftnlen)13);
		    }
		    if (qpmode == 1) {
			s_copy(buffp + 116, "cgItns", (ftnlen)6, (ftnlen)6);
		    }
		} else if (! (*feasible)) {
		    s_copy(buffp + 12, "FP", (ftnlen)2, (ftnlen)2);
		    s_copy(buffp + 59, "NumInf", (ftnlen)6, (ftnlen)6);
		    s_copy(buffp + 68, " ", (ftnlen)7, (ftnlen)1);
		    s_copy(buffp + 75, "SumInf", (ftnlen)6, (ftnlen)6);
		}
		snprnt_(&c__11, buffp, &iw[1], leniw, (ftnlen)138);
	    }
	    if (*feasible) {
		s_wsfi(&io___77);
		do_fio(&c__1, (char *)&itns, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&(*kprc), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&(*djqprt), (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&jsqn, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&jsrn, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&jbrn, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&(*step), (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&(*pivot), (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&(*nonopt), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&(*objprt), (ftnlen)sizeof(doublereal));
		i__1 = lenl + lenu;
		do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ncp, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&(*rgnorm), (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&(*ns), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&(*condzhz), (ftnlen)sizeof(doublereal))
			;
		do_fio(&c__1, (char *)&cgitn, (ftnlen)sizeof(integer));
		e_wsfi();
	    } else {
		s_wsfi(&io___78);
		do_fio(&c__1, (char *)&itns, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&(*kprc), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&(*djqprt), (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&jsqn, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&jsrn, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&jbrn, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&(*step), (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&(*pivot), (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&numinf, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&suminf, (ftnlen)sizeof(doublereal));
		i__1 = lenl + lenu;
		do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ncp, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&(*rgnorm), (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&(*ns), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&(*condzhz), (ftnlen)sizeof(doublereal))
			;
		do_fio(&c__1, (char *)&cgitn, (ftnlen)sizeof(integer));
		e_wsfi();
	    }
	} else {
/* nnH == 0 */
	    if (printhead) {
		s_copy(buffp, "    Itn pp  LP mult   +SBS   -SBS    -BS     "
			"Step    Pivot NonOpt    LP objective     L+U ncp", (
			ftnlen)138, (ftnlen)93);
		if (*feasible && *elastic) {
		    if (*lvlobje == 1) {
			s_copy(buffp + 68, "Elastic LPobj", (ftnlen)13, (
				ftnlen)13);
		    }
		    if (*lvlobje == 2) {
			s_copy(buffp + 68, "      SumInfE", (ftnlen)13, (
				ftnlen)13);
		    }
		} else if (! (*feasible)) {
		    s_copy(buffp + 12, "FP", (ftnlen)2, (ftnlen)2);
		    s_copy(buffp + 59, "NumInf", (ftnlen)6, (ftnlen)6);
		    s_copy(buffp + 68, " ", (ftnlen)7, (ftnlen)1);
		    s_copy(buffp + 75, "SumInf", (ftnlen)6, (ftnlen)6);
		}
		if (*ns > 0) {
		    s_copy(buffp + 93, "  rgNorm    nS", (ftnlen)14, (ftnlen)
			    14);
		}
		snprnt_(&c__11, buffp, &iw[1], leniw, (ftnlen)138);
	    }
	    if (*feasible) {
		s_wsfi(&io___79);
		do_fio(&c__1, (char *)&itns, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&(*kprc), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&(*djqprt), (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&jsqn, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&jsrn, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&jbrn, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&(*step), (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&(*pivot), (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&(*nonopt), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&(*objprt), (ftnlen)sizeof(doublereal));
		i__1 = lenl + lenu;
		do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ncp, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&(*rgnorm), (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&(*ns), (ftnlen)sizeof(integer));
		e_wsfi();
	    } else {
		s_wsfi(&io___80);
		do_fio(&c__1, (char *)&itns, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&(*kprc), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&(*djqprt), (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&jsqn, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&jsrn, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&jbrn, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&(*step), (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&(*pivot), (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&numinf, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&suminf, (ftnlen)sizeof(doublereal));
		i__1 = lenl + lenu;
		do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ncp, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&(*rgnorm), (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&(*ns), (ftnlen)sizeof(integer));
		e_wsfi();
	    }
	}
	if ((doublereal) (*kprc) == 0.) {
	    s_copy(buffp + 7, " ", (ftnlen)3, (ftnlen)1);
	}
	if (*djqprt == 0.) {
	    s_copy(buffp + 10, " ", (ftnlen)9, (ftnlen)1);
	}
	if (*jsq == 0) {
	    s_copy(buffp + 19, " ", (ftnlen)7, (ftnlen)1);
	}
	if (*jsr == 0) {
	    s_copy(buffp + 26, " ", (ftnlen)7, (ftnlen)1);
	}
	if (*jbr == 0) {
	    s_copy(buffp + 33, " ", (ftnlen)7, (ftnlen)1);
	}
	if (*step == 0.) {
	    s_copy(buffp + 40, " ", (ftnlen)9, (ftnlen)1);
	}
	if (*pivot == 0.) {
	    s_copy(buffp + 49, " ", (ftnlen)9, (ftnlen)1);
	}
	if (ncp == 0) {
	    s_copy(buffp + 89, " ", (ftnlen)4, (ftnlen)1);
	}
	if (*rgnorm == 0.) {
	    s_copy(buffp + 94, " ", (ftnlen)7, (ftnlen)1);
	}
	if (*ns == 0) {
	    s_copy(buffp + 101, " ", (ftnlen)6, (ftnlen)1);
	}
	if (*condzhz == 0.) {
	    s_copy(buffp + 108, " ", (ftnlen)7, (ftnlen)1);
	}
/* condZHZ */
	if (cgitn == 0) {
	    s_copy(buffp + 116, " ", (ftnlen)6, (ftnlen)1);
	}
	snprnt_(&c__1, buffp, &iw[1], leniw, (ftnlen)138);
    }
/*     ================================================================== */
/*     Summary file */
/*     ================================================================== */
    if (printsum) {
/* Print one line to the summary file */
	newset = *liness == 0;
	printhead = printheads || newset;
	if (printhead) {
	    iw[225] = 0;
	    *liness = 10;
	}
	--(*liness);
	if (*nnh > 0) {
	    if (printhead) {
		s_copy(buffs, "    Itn  QP mult  NonOpt   QP objective    nS"
			"    rgNorm", (ftnlen)80, (ftnlen)55);
		if (*feasible && *elastic) {
		    if (*lvlobje == 1) {
			s_copy(buffs + 26, "Elastic QPobj", (ftnlen)13, (
				ftnlen)13);
		    }
		    if (*lvlobje == 2) {
			s_copy(buffs + 26, "      SumInfE", (ftnlen)13, (
				ftnlen)13);
		    }
		    if (qpmode == 1) {
			s_copy(buffs + 57, "cgItns", (ftnlen)6, (ftnlen)6);
		    }
		} else if (! (*feasible)) {
		    s_copy(buffs + 9, "FP", (ftnlen)2, (ftnlen)2);
		    s_copy(buffs + 18, "NumInf", (ftnlen)6, (ftnlen)6);
		    s_copy(buffs + 26, " ", (ftnlen)7, (ftnlen)1);
		    s_copy(buffs + 33, "SumInf", (ftnlen)6, (ftnlen)6);
		}
		snprnt_(&c__12, buffs, &iw[1], leniw, (ftnlen)80);
	    }
	    if (*feasible) {
		s_wsfi(&io___81);
		do_fio(&c__1, (char *)&minors, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&(*djqprt), (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&(*nonopt), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&(*objprt), (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&(*ns), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&(*rgnorm), (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&cgitn, (ftnlen)sizeof(integer));
		e_wsfi();
	    } else {
		s_wsfi(&io___82);
		do_fio(&c__1, (char *)&minors, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&(*djqprt), (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&numinf, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&suminf, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&(*ns), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&(*rgnorm), (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&cgitn, (ftnlen)sizeof(integer));
		e_wsfi();
	    }
	} else {
/* nnH == 0 */
	    if (printhead) {
		s_copy(buffs, "    Itn  LP mult  NonOpt   LP objective", (
			ftnlen)80, (ftnlen)39);
		if (*feasible && *elastic) {
		    if (*lvlobje == 1) {
			s_copy(buffs + 26, "Elastic LPobj", (ftnlen)13, (
				ftnlen)13);
		    }
		    if (*lvlobje == 2) {
			s_copy(buffs + 26, "      SumInfE", (ftnlen)13, (
				ftnlen)13);
		    }
		} else if (! (*feasible)) {
		    s_copy(buffs + 9, "FP", (ftnlen)2, (ftnlen)2);
		    s_copy(buffs + 18, "NumInf", (ftnlen)6, (ftnlen)6);
		    s_copy(buffs + 26, " ", (ftnlen)7, (ftnlen)1);
		    s_copy(buffs + 33, "SumInf", (ftnlen)6, (ftnlen)6);
		}
		snprnt_(&c__12, buffs, &iw[1], leniw, (ftnlen)80);
	    }
	    if (*feasible) {
		s_wsfi(&io___83);
		do_fio(&c__1, (char *)&minors, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&(*djqprt), (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&(*nonopt), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&(*objprt), (ftnlen)sizeof(doublereal));
		e_wsfi();
	    } else {
		s_wsfi(&io___84);
		do_fio(&c__1, (char *)&minors, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&(*djqprt), (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&numinf, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&suminf, (ftnlen)sizeof(doublereal));
		e_wsfi();
	    }
	}
	if (*djqprt == 0.) {
	    s_copy(buffs + 8, " ", (ftnlen)8, (ftnlen)1);
	}
	if (*ns == 0) {
	    s_copy(buffs + 39, " ", (ftnlen)6, (ftnlen)1);
	}
	if (*rgnorm == 0.) {
	    s_copy(buffs + 45, " ", (ftnlen)10, (ftnlen)1);
	}
	if (cgitn == 0) {
	    s_copy(buffs + 55, " ", (ftnlen)8, (ftnlen)1);
	}
	snprnt_(&c__2, buffs, &iw[1], leniw, (ftnlen)80);
    }
/*     ------------------------------------------------------------------ */
/*     Debug output. */
/*     ------------------------------------------------------------------ */
    if (lprdbg == 100) {
	snprnt_(&c__11, " BS values...", &iw[1], leniw, (ftnlen)13);
	i__1 = *m;
	for (k = 1; k <= i__1; ++k) {
	    s_wsfi(&io___86);
	    i__2 = s2varn_(&kbs[k], leniw, &iw[1]);
	    do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&xbs[k], (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    snprnt_(&c__1, buffp, &iw[1], leniw, (ftnlen)138);
	}
	snprnt_(&c__11, " SB values...", &iw[1], leniw, (ftnlen)13);
	i__1 = *m + *ns;
	for (k = *m + 1; k <= i__1; ++k) {
	    s_wsfi(&io___87);
	    i__2 = s2varn_(&kbs[k], leniw, &iw[1]);
	    do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&xbs[k], (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    snprnt_(&c__1, buffp, &iw[1], leniw, (ftnlen)138);
	}
    }
    return 0;
/*     Minor log,  Print file. */
/*     Minor log,  Summary file. */
/*                itn  mul    ninf    obj  ns  rg     cg */
} /* sqlog_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine sqLog */
/* Subroutine */ int sqset_(char *buffer, integer *iprint, integer *isumm, 
	integer *iexit, char *cw, integer *lencw, integer *iw, integer *leniw,
	 doublereal *rw, integer *lenrw, ftnlen buffer_len, ftnlen cw_len)
{
    static char key[16], cval[8];
    static integer ival;
    static doublereal rval;
    extern /* Subroutine */ int s3opt_(logical *, char *, char *, char *, 
	    integer *, doublereal *, integer *, integer *, integer *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, ftnlen, 
	    ftnlen, ftnlen, ftnlen);

/*     ================================================================== */
/*     sqSet  decodes the option contained in  buffer. */

/*     The buffer is output to file iPrint, minus trailing blanks. */
/*     Error messages are output to files iPrint and iSumm. */
/*     Buffer is echoed to iPrint but normally not to iSumm. */
/*     It is echoed to iSumm before any error msg. */

/*     On entry, */
/*     iPrint is the print   file.  no output occurs if iPrint .le 0. */
/*     iSumm  is the Summary file.  no output occurs if iSumm  .le 0. */
/*     iExit  is the number of errors so far. */

/*     On exit, */
/*     iExit  is the number of errors so far. */

/*     27 Nov 1991: first version of sqSet. */
/*     20 Sep 1998: current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    s3opt_(&c_true, buffer, key, cval, &ival, &rval, iprint, isumm, iexit, cw 
	    + 8, lencw, &iw[1], leniw, &rw[1], lenrw, buffer_len, (ftnlen)16, 
	    (ftnlen)8, (ftnlen)8);
    return 0;
} /* sqset_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine sqSet */
/* Subroutine */ int sqseti_(char *buffer, integer *ivalue, integer *iprint, 
	integer *isumm, integer *iexit, char *cw, integer *lencw, integer *iw,
	 integer *leniw, doublereal *rw, integer *lenrw, ftnlen buffer_len, 
	ftnlen cw_len)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    , i_len(char *, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static char key[16], cval[8];
    static integer ival;
    static doublereal rval;
    extern /* Subroutine */ int s3opt_(logical *, char *, char *, char *, 
	    integer *, doublereal *, integer *, integer *, integer *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, ftnlen, 
	    ftnlen, ftnlen, ftnlen);
    static char buff72[72];
    static integer lenbuf;

    /* Fortran I/O blocks */
    static icilist io___93 = { 0, key, 0, "(i16)", 16, 1 };


/*     ================================================================== */
/*     sqSeti decodes the option contained in  buffer // ivalue. */
/*     The parameters other than ivalue are as in sqSet. */

/*     27 Nov 1991: first version of sqSeti. */
/*     20 Sep 1998: current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    s_wsfi(&io___93);
    do_fio(&c__1, (char *)&(*ivalue), (ftnlen)sizeof(integer));
    e_wsfi();
    lenbuf = i_len(buffer, buffer_len);
    s_copy(buff72, buffer, (ftnlen)72, buffer_len);
    i__1 = lenbuf;
    s_copy(buff72 + i__1, key, lenbuf + 16 - i__1, (ftnlen)16);
    s3opt_(&c_true, buff72, key, cval, &ival, &rval, iprint, isumm, iexit, cw 
	    + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)72, (ftnlen)16, 
	    (ftnlen)8, (ftnlen)8);
    return 0;
} /* sqseti_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine sqSeti */
/* Subroutine */ int sqsetr_(char *buffer, doublereal *rvalue, integer *
	iprint, integer *isumm, integer *iexit, char *cw, integer *lencw, 
	integer *iw, integer *leniw, doublereal *rw, integer *lenrw, ftnlen 
	buffer_len, ftnlen cw_len)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    , i_len(char *, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static char key[16], cval[8];
    static integer ival;
    static doublereal rval;
    extern /* Subroutine */ int s3opt_(logical *, char *, char *, char *, 
	    integer *, doublereal *, integer *, integer *, integer *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, ftnlen, 
	    ftnlen, ftnlen, ftnlen);
    static char buff72[72];
    static integer lenbuf;

    /* Fortran I/O blocks */
    static icilist io___100 = { 0, key, 0, "(1p, e16.8)", 16, 1 };


/*     ================================================================== */
/*     sqSetr decodes the option contained in  buffer // rvalue. */
/*     The parameters other than rvalue are as in sqSet. */

/*     27 Nov 1991: first version of sqSetr. */
/*     20 Sep 1998: current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    s_wsfi(&io___100);
    do_fio(&c__1, (char *)&(*rvalue), (ftnlen)sizeof(doublereal));
    e_wsfi();
    lenbuf = i_len(buffer, buffer_len);
    s_copy(buff72, buffer, (ftnlen)72, buffer_len);
    i__1 = lenbuf;
    s_copy(buff72 + i__1, key, lenbuf + 16 - i__1, (ftnlen)16);
    s3opt_(&c_true, buff72, key, cval, &ival, &rval, iprint, isumm, iexit, cw 
	    + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)72, (ftnlen)16, 
	    (ftnlen)8, (ftnlen)8);
    return 0;
} /* sqsetr_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine sqSetr */
integer sqget_(char *buffer, integer *iexit, char *cw, integer *lencw, 
	integer *iw, integer *leniw, doublereal *rw, integer *lenrw, ftnlen 
	buffer_len, ftnlen cw_len)
{
    /* System generated locals */
    integer ret_val;

    /* Local variables */
    static char key[16];
    extern /* Subroutine */ int s3opt_(logical *, char *, char *, char *, 
	    integer *, doublereal *, integer *, integer *, integer *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, ftnlen, 
	    ftnlen, ftnlen, ftnlen);
    static char cvalue[8];
    static integer ivalue;
    static doublereal rvalue;

/*     ================================================================== */
/*     sqGet  decodes the option contained in  buffer */
/*     and returns 1 if the option has previously been set, else 0. */
/*     For example, */
/*     i = sqGet ( 'Maximize', iExit, cw, lencw, iw, leniw, rw, lenrw ) */

/*     01 Aug 2003: First version of sqGet.  Needed because */
/*                  sqGetc, sqGeti, sqGetr were not well defined */
/*                  for strings that had no numerical value. */
/*     01 Aug 2003: Current version of sqGet. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    s3opt_(&c_false, buffer, key, cvalue, &ivalue, &rvalue, &c__0, &c__0, 
	    iexit, cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, buffer_len, (
	    ftnlen)16, (ftnlen)8, (ftnlen)8);
    ret_val = ivalue;
    return ret_val;
} /* sqget_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* integer function sqGet */
/* Subroutine */ int sqgetc_(char *buffer, char *cvalue, integer *iexit, char 
	*cw, integer *lencw, integer *iw, integer *leniw, doublereal *rw, 
	integer *lenrw, ftnlen buffer_len, ftnlen cvalue_len, ftnlen cw_len)
{
    static char key[16];
    static integer ival;
    static doublereal rval;
    extern /* Subroutine */ int s3opt_(logical *, char *, char *, char *, 
	    integer *, doublereal *, integer *, integer *, integer *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, ftnlen, 
	    ftnlen, ftnlen, ftnlen);

/*     ================================================================== */
/*     sqGetc gets the value of the option contained in  buffer. */
/*     The parameters other than cvalue are as in sqSet. */

/*     17 May 1998: first version of sqGetc. */
/*     20 Sep 1998: current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    s3opt_(&c_false, buffer, key, cvalue, &ival, &rval, &c__0, &c__0, iexit, 
	    cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, buffer_len, (ftnlen)
	    16, (ftnlen)8, (ftnlen)8);
    return 0;
} /* sqgetc_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine sqGetc */
/* Subroutine */ int sqgeti_(char *buffer, integer *ivalue, integer *iexit, 
	char *cw, integer *lencw, integer *iw, integer *leniw, doublereal *rw,
	 integer *lenrw, ftnlen buffer_len, ftnlen cw_len)
{
    static char key[16], cval[8];
    static doublereal rval;
    extern /* Subroutine */ int s3opt_(logical *, char *, char *, char *, 
	    integer *, doublereal *, integer *, integer *, integer *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, ftnlen, 
	    ftnlen, ftnlen, ftnlen);

/*     ================================================================== */
/*     sqGeti gets the value of the option contained in  buffer. */
/*     The parameters other than ivalue are as in sqSet. */

/*     17 May 1998: first version of sqGeti. */
/*     20 Sep 1998: current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    s3opt_(&c_false, buffer, key, cval, ivalue, &rval, &c__0, &c__0, iexit, 
	    cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, buffer_len, (ftnlen)
	    16, (ftnlen)8, (ftnlen)8);
    return 0;
} /* sqgeti_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine sqGeti */
/* Subroutine */ int sqgetr_(char *buffer, doublereal *rvalue, integer *iexit,
	 char *cw, integer *lencw, integer *iw, integer *leniw, doublereal *
	rw, integer *lenrw, ftnlen buffer_len, ftnlen cw_len)
{
    static char key[16], cval[8];
    static integer ival;
    extern /* Subroutine */ int s3opt_(logical *, char *, char *, char *, 
	    integer *, doublereal *, integer *, integer *, integer *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, ftnlen, 
	    ftnlen, ftnlen, ftnlen);

/*     ================================================================== */
/*     sqGetr gets the value of the option contained in  buffer. */
/*     The parameters other than rvalue are as in sqSet. */

/*     17 May 1998: first version of sqGetr. */
/*     20 Sep 1998: current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    s3opt_(&c_false, buffer, key, cval, &ival, rvalue, &c__0, &c__0, iexit, 
	    cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, buffer_len, (ftnlen)
	    16, (ftnlen)8, (ftnlen)8);
    return 0;
} /* sqgetr_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine sqGetr */
/* Subroutine */ int nullhx_(integer *nnh, doublereal *x, doublereal *hx, 
	integer *status, char *cu, integer *lencu, integer *iu, integer *
	leniu, doublereal *ru, integer *lenru, ftnlen cu_len)
{
    /* Format strings */
    static char fmt_1000[] = "(//\002 XXX  The default (dummy) version of su"
	    "broutine Hx\002,\002     has been called. \002/\002 XXX  A user-"
	    "defined version is required when solving a QP\002)";

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void);

    /* Local variables */
    static integer nout;

    /* Fortran I/O blocks */
    static cilist io___120 = { 0, 0, 0, fmt_1000, 0 };


/*     ================================================================== */
/*     This is the dummy (empty) version of the routine qpHx. */
/*     It should never be called by SQOPT. */

/*     Warn the user (on the standard output) that it has been called. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --hx;
    --x;
    cu -= 8;
    --iu;
    --ru;

    /* Function Body */
    nout = 6;
    if (*status == 1) {
	if (nout > 0) {
	    io___120.ciunit = nout;
	    s_wsfe(&io___120);
	    e_wsfe();
	}
    }
    return 0;
} /* nullhx_ */
