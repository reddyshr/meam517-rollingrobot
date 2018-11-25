/* ../snopt7/src/sn87sopt.f -- translated by f2c (version 20100827).
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

static integer c__1 = 1;
static doublereal c_b4 = 0.;
static integer c__4 = 4;
static integer c__0 = 0;
static integer c__2 = 2;
static integer c__3 = 3;
static integer c_n2 = -2;
static integer c__13 = 13;
static integer c_n3 = -3;
static doublereal c_b198 = .33333;
static integer c__7 = 7;
static doublereal c_b245 = 1.;
static doublereal c_b246 = -1.;
static integer c__23 = 23;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*     File  sn87sopt.f */

/*     s8solve */
/*     s8defaults  s8firstCall */
/*     s8Map       s8SQP        s8callStatus */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Subroutine */ int s8solve_(integer *iexit, char *solver, integer *
	starttype, S_fp funwrapper, U_fp funcon, U_fp funobj, U_fp userhv, 
	U_fp mjrlog, U_fp mnrlog, U_fp snstop, logical *gotr, integer *m, 
	integer *n, integer *nb, integer *nnh, integer *nncon, integer *nnjac,
	 integer *nnobj, integer *nnames, integer *iobj, doublereal *objadd, 
	doublereal *fobj, doublereal *objtrue, integer *ninf, doublereal *
	sinf, integer *ninfe, doublereal *sinfe, integer *nej, integer *nlocj,
	 integer *locj, integer *indj, doublereal *jcol, integer *neh, 
	integer *nloch, integer *loch, integer *indh, doublereal *hcol, 
	doublereal *bl, doublereal *bu, char *names, integer *hs, doublereal *
	x, doublereal *pi, doublereal *rc, integer *majors, integer *ns, char 
	*cu, integer *lencu, integer *iu, integer *leniu, doublereal *ru, 
	integer *lenru, char *cw, integer *lencw, integer *iw, integer *leniw,
	 doublereal *rw, integer *lenrw, ftnlen solver_len, ftnlen names_len, 
	ftnlen cu_len, ftnlen cw_len)
{
    /* Format strings */
    static char fmt_1900[] = "(\002 No. of iterations\002,i20,2x,\002 Object"
	    "ive\002,6x,1p,e22.10)";
    static char fmt_1920[] = "(\002 No. of major iterations\002,i14,2x,\002 "
	    "Linear    obj. term\002,1p,e18.10)";
    static char fmt_1935[] = "(\002 Penalty parameter\002,1p,e20.3,2x,\002 N"
	    "onlinear obj. term\002,1p,e18.10)";
    static char fmt_1915[] = "(\002 Elastic weight\002,12x,1p,e11.1,2x,\002 "
	    "Elastic   objective\002,1p,e18.10)";
    static char fmt_2050[] = "(\002                  \002,22x,\002 Elastic i"
	    "nfeas   \002,1p,e20.10)";
    static char fmt_1930[] = "(\002 Penalty parameter\002,1p,e20.3,2x,\002 N"
	    "orm (x - x0)**2   \002,1p,e18.10)";
    static char fmt_1940[] = "(\002                  \002,22x,\002 Nonlinear"
	    " obj. term\002,1p,e18.10)";
    static char fmt_1952[] = "(\002 User function calls (total)\002,i10,2x"
	    ",\002 Calls with modes 1,2 (known g)\002,i7)";
    static char fmt_1951[] = "(\002 User function calls (total)\002,i10)";
    static char fmt_1953[] = "(\002 Calls for forward differencing\002,i7,"
	    "2x,\002 Calls for central differencing\002,i7)";
    static char fmt_1950[] = "(\002 No. of calls to funobj\002,i15,2x,\002 N"
	    "o. of calls to funcon\002,i15)";
    static char fmt_1955[] = "(\002 Calls with modes 1,2 (known g)\002,i7,"
	    "2x,\002 Calls with modes 1,2 (known g)\002,i7)";
    static char fmt_1960[] = "(\002 Calls for forward differencing\002,i7,"
	    "2x,\002 Calls for forward differencing\002,i7)";
    static char fmt_1962[] = "(\002 Calls for central differencing\002,i7,"
	    "2x,\002 Calls for central differencing\002,i7)";
    static char fmt_1963[] = "(\002                  \002,22x,\002 Calls for"
	    " forward differencing\002,i7)";
    static char fmt_1964[] = "(\002                  \002,22x,\002 Calls for"
	    " central differencing\002,i7)";
    static char fmt_1965[] = "(\002 Calls for forward differencing\002,i7)";
    static char fmt_1966[] = "(\002 Calls for central differencing\002,i7)";
    static char fmt_1970[] = "(\002 No. of superbasics\002,i19,2x,\002 No. o"
	    "f basic nonlinears\002,i14)";
    static char fmt_1973[] = "(\002 No. of CG iterations\002,i17)";
    static char fmt_1975[] = "(\002 No. of degenerate steps\002,i14,2x,\002 "
	    "Percentage\002,f27.2)";
    static char fmt_1910[] = "(\002 No. of infeasibilities\002,i15,2x,\002 S"
	    "um of infeas\002,1p,e24.10)";
    static char fmt_1905[] = "(\002 No. of iterations\002,i20)";
    static char fmt_1980[] = "(\002 No. of infeas elastics\002,i15,2x,\002 E"
	    "lastic infeas   \002,1p,e20.10)";

    /* System generated locals */
    address a__1[2];
    integer i__1, i__2, i__3[2];
    char ch__1[38];

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen), s_cat(char *,
	     char **, integer *, integer *, ftnlen);
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    , s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static doublereal scaleobj;
    static integer minimize;
    static doublereal maxvirel;
    static integer mjrprint;
    static doublereal toloptfp;
    static integer mnrprint;
    static doublereal toloptqp;
    static integer j, k;
    static logical nonlinear;
    static integer lfeastype, lr, ly;
    static logical feasiblelc;
    static integer lx1, ly1, ly2, nx0, ldg, lhd;
    extern /* Subroutine */ int s8getfeaslc_(integer *, integer *, U_fp, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, doublereal *, doublereal *, doublereal *
	    , integer *, integer *, integer *, integer *, doublereal *, 
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *, integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, doublereal *, doublereal *, doublereal *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, ftnlen);
    static integer nnb;
    extern /* Subroutine */ int s2getscales_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, doublereal *, integer *);
    static integer mbs, lrg, ldx, lfv, lfx, itn;
    extern /* Subroutine */ int s8firstcall_(integer *, logical *, logical *, 
	    logical *, logical *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, S_fp, U_fp, U_fp, 
	    U_fp, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, doublereal *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, char *, integer *, integer *, integer 
	    *, doublereal *, integer *, char *, integer *, integer *, integer 
	    *, doublereal *, integer *, ftnlen, ftnlen);
    static integer liy, lkx, nkx;
    static char str[133];
    static integer nnh0, lrg2;
    static doublereal eps0;
    static integer liy1, liy2;
    extern /* Subroutine */ int s8fx_(integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    static char str2[133];
    static integer lgbs, lkbs, lenr, lhdx, lgqp, lpbs, lxbs, maxs, nrhs, ludx,
	     lxqp;
    static doublereal tolx;
    static logical nonlinearobj, nonlinearcon;
    static integer lenx0, nrhs0, lxqp0;
    static logical needb;
    extern /* Subroutine */ int s8sqp_(integer *, S_fp, U_fp, U_fp, U_fp, 
	    U_fp, U_fp, U_fp, logical *, logical *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, integer *, doublereal *, doublereal *, doublereal *
	    , doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, integer *, integer 
	    *, integer *, doublereal *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, char *, integer *, integer *, integer *, doublereal 
	    *, integer *, char *, integer *, integer *, integer *, doublereal 
	    *, integer *, ftnlen, ftnlen);
    static doublereal degen;
    extern /* Subroutine */ int dload_(integer *, doublereal *, doublereal *, 
	    integer *), iload_(integer *, integer *, integer *, integer *), 
	    dddiv_(integer *, doublereal *, integer *, doublereal *, integer *
	    ), ddscl_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer lgobj, lblbs, llocg, lfcon, lgcon, inewb;
    static char mprob[8];
    static integer lblqp, lbubs, lbuqp, lycon, lxpen, nlocg, numlc;
    static doublereal duinf, vilim, maxvi, wtinf;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), s2applyscales_(integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    static doublereal xnorm, supvi;
    static integer lgobj1, lgobj2, lfcon1, lfcon2, lgcon1, lgcon2, nnobj0, 
	    nnobj1, nncon0, nncon1;
    extern /* Subroutine */ int s1page_(integer *, integer *, integer *), 
	    s2amat_(integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, integer *), s1time_(integer *, integer *, integer *,
	     integer *, doublereal *, integer *), s8rand_(integer *, integer *
	    , doublereal *), s5getb_(integer *, integer *, U_fp, logical *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, integer *, integer *, doublereal *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, char *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, char *, integer *, integer *, integer *, doublereal 
	    *, integer *, ftnlen, ftnlen);
    static integer lycon1, lycon2;
    static doublereal wtinf0;
    extern /* Subroutine */ int s4newb_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, char *, char *, integer *, integer *, integer *, 
	    ftnlen, ftnlen);
    static doublereal pnorm1, pnorm2;
    extern /* Subroutine */ int s2vmax_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), s4stat_(
	    integer *, char *, ftnlen), s5fixx_(integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *);
    static integer ndegen, modefg;
    static doublereal infbnd, objinf;
    static integer icrash, negcon, lcrash;
    static doublereal objlin, objmin, tcrash;
    static char istate[4*3];
    extern doublereal dnormi_(integer *, doublereal *, integer *);
    static integer inform__, itnlim, lssave, ldycon, lgconu, letype, imaxvi;
    static logical fponly;
    static integer lqprhs, minors, numliq;
    static doublereal fmerit, pinorm, rgnorm;
    extern /* Subroutine */ int snwrap_(integer *, char *, char *, char *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen), snprnt_(integer *, 
	    char *, integer *, integer *, ftnlen), s2crash_(integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *), s4saveb_(integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, doublereal *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, char *, 
	    doublereal *, doublereal *, doublereal *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, ftnlen, ftnlen), 
	    s8gcopy_(integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *, integer *, integer *, doublereal *);
    static logical elastic;
    static integer lscales, lblsave;
    static doublereal signobj;
    static integer lbusave, lestate;
    static doublereal penparm[4];
    static integer lvlsrch;
    static doublereal pennorm;
    static logical gotfuns;
    extern /* Subroutine */ int s8scaleg_(integer *, doublereal *, doublereal 
	    *, doublereal *, integer *), s8scalej_(integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, integer *, integer 
	    *, integer *, doublereal *, doublereal *, integer *);

    /* Fortran I/O blocks */
    static icilist io___121 = { 0, str, 0, fmt_1900, 133, 1 };
    static icilist io___122 = { 0, str, 0, fmt_1920, 133, 1 };
    static icilist io___123 = { 0, str, 0, fmt_1935, 133, 1 };
    static icilist io___124 = { 0, str, 0, fmt_1915, 133, 1 };
    static icilist io___125 = { 0, str, 0, fmt_2050, 133, 1 };
    static icilist io___126 = { 0, str, 0, fmt_1900, 133, 1 };
    static icilist io___127 = { 0, str, 0, fmt_1920, 133, 1 };
    static icilist io___128 = { 0, str, 0, fmt_1930, 133, 1 };
    static icilist io___129 = { 0, str, 0, fmt_1935, 133, 1 };
    static icilist io___130 = { 0, str, 0, fmt_1940, 133, 1 };
    static icilist io___131 = { 0, str, 0, fmt_1952, 133, 1 };
    static icilist io___132 = { 0, str, 0, fmt_1951, 133, 1 };
    static icilist io___133 = { 0, str, 0, fmt_1953, 133, 1 };
    static icilist io___134 = { 0, str, 0, fmt_1950, 133, 1 };
    static icilist io___135 = { 0, str, 0, fmt_1955, 133, 1 };
    static icilist io___136 = { 0, str, 0, fmt_1960, 133, 1 };
    static icilist io___137 = { 0, str, 0, fmt_1962, 133, 1 };
    static icilist io___138 = { 0, str, 0, fmt_1963, 133, 1 };
    static icilist io___139 = { 0, str, 0, fmt_1964, 133, 1 };
    static icilist io___140 = { 0, str, 0, fmt_1965, 133, 1 };
    static icilist io___141 = { 0, str, 0, fmt_1966, 133, 1 };
    static icilist io___142 = { 0, str, 0, fmt_1970, 133, 1 };
    static icilist io___143 = { 0, str, 0, fmt_1973, 133, 1 };
    static icilist io___144 = { 0, str, 0, fmt_1975, 133, 1 };
    static icilist io___145 = { 0, str, 0, fmt_1910, 133, 1 };
    static icilist io___146 = { 0, str, 0, fmt_1905, 133, 1 };
    static icilist io___147 = { 0, str, 0, fmt_1980, 133, 1 };


/*     ================================================================== */
/*     s8solve solves the current problem. */

/*     On entry, */
/*     the specs file has been read, */
/*     all data items have been loaded (including locJ, indJ, Jcol, ...), */
/*     and workspace has been allocated. */

/*     On exit, */
/*     iExit  =  0 if an optimal solution was found, */
/*            =  1 if the problem was infeasible, */
/*            =  2 if the problem was unbounded, */
/*            =  3 if the Iteration limit was exceeded, */
/*           ge  4 if iterations were terminated by some other */
/*                 error condition (see the SNOPT user's guide). */

/*     15 Nov 1991: First version based on Minos 5.4 routine misolv. */
/*     13 Feb 1994: Eliminated "Cycle" options. */
/*                  Simplified s4getB. */
/*     12 Nov 1994: Integer workspace added. */
/*     25 Jul 1996: Sign of the slacks changed. */
/*     28 Sep 1997: Character workspace added. */
/*     11 Nov 1997: Backtracking for undefined functions. */
/*     26 Dec 1997: Dummy Jacobian scaled in feasibility phase. */
/*     27 Aug 1998: Constant Jacobian elements handled correctly. */
/*     10 Oct 1998: Objective and constraint gradient checking merged. */
/*     11 Oct 1998: Facility to combine funobj and funcon added. */
/*     23 Dec 1999: Suboptimize option added. */
/*     30 Dec 2000: Housekeeping for first function call now in snwrapup. */
/*     03 Aug 2003: snEXIT and snPRNT adopted. */
/*     19 Mar 2006: Fx initialized for s8savB. */
/*     16 Jun 2008: Call-status implemented correctly. */
/*     18 Jun 2008: iy2, pBS and rg2 added to workspace. */
/*     23 Oct 2010: pinorm initialized at 1.0. */
/*     01 Dec 2012: eState initialized in s8solve. */
/*     05 Apr 2014: nnH added to argument list. */
/*     19 Oct 2014: Added user-defined Hessian. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/* = 0,1,2,3 or 4, deriv level */
/* 0,1,2  => LM, FM, Newton */
/* scale option */
/* 1, 0, -1  => MIN, FP, MAX */
/* =1(2) for forwd(cntrl) diffs */
/* > 0 => some exact derivs */
/* number of calls of fCon */
/* number of calls of fCon */
/* number of calls of fCon */
/* number of calls of fCon */
/* number of calls of fObj */
/* number of calls of fObj */
/* number of calls of fObj */
/* number of calls of fObj */
/* =1(0) for pd  QP Hessian */
/* # of LU factorizations */
/* # lines in log     file */
/* # lines in summary file */
/* NP user-routine call-status */
/* Number of symmlq iterations */
/* symmlq itns for last minor */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --pi;
    --rc;
    --x;
    --hs;
    --bu;
    --bl;
    names -= 8;
    --jcol;
    --indj;
    --locj;
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
    inewb = iw[124];
/* new basis file */
    negcon = iw[20];
/* # of nonzero elems in J */
    lenr = iw[28];
/* R(lenR) is the reduced Hessian factor */
    maxs = iw[53];
/* max # of superbasics */
    lvlsrch = iw[76];
/* >0     => use derivatives in the line search */
    icrash = iw[88];
/* Crash option */
    itnlim = iw[89];
/* limit on total iterations */
    mjrprint = iw[92];
/* Major print level */
    mnrprint = iw[93];
/* Minor print level */
    minimize = iw[199];
/* 1 (-1)    => minimize (maximize) */
    nkx = iw[247];
/* dimension of kx and its inverse, kxN */
    eps0 = rw[2];
    toloptfp = rw[51];
/* Minor Phase 1 Opt tol */
    toloptqp = rw[52];
/* Minor Phase 2 Opt tol */
    tolx = rw[56];
/* Minor feasibility tolerance. */
    tcrash = rw[62];
/* crash tolerance. */
    infbnd = rw[70];
/* definition of an infinite bound. */
    vilim = rw[81];
/* violation limit */
    wtinf0 = rw[88];
/* infeasibility weight */
    s_copy(mprob, cw + 408, (ftnlen)8, (ftnlen)8);
/*     Addresses */
/* Problem name */
    lkx = iw[251];
/* j  = kx (jN) => col j of Jcol is variable jN */
    llocg = iw[260];
/* locG(nnJac+1) = column pointers for indG */
    lblqp = iw[271];
/* blQP(nb)    = working lower bounds */
    lbuqp = iw[272];
/* buQP(nb)    = working upper bounds */
    lblbs = iw[273];
/* blBS(mBS)   = lower bounds for xBS */
    lbubs = iw[274];
/* buBS(mBS)   = upper bounds for xBS */
    lblsave = iw[275];
/* blSave(nb)  = bl for LC feasibility */
    lbusave = iw[276];
/* buSave(nb)  = bu for LC feasibility */
    lpbs = iw[277];
/* pBS(nb)     = search direction */
    lqprhs = iw[278];
/* QPrhs(nnCon)=  QP constraint rhs */
    letype = iw[283];
/* eType(nb) list of elastic vars */
    lfeastype = iw[284];
/* feasType(mBS), feasibility types */
    lestate = iw[285];
/* eState(nb), status of elastics */
    ldx = iw[287];
/* dx(nb)      = x1 - x */
    lhdx = iw[288];
/* Hdx(nnH)    = product of H with  x1 - x */
    ldg = iw[289];
/* dg(nnH)     = gradient difference */
    lgqp = iw[290];
/* gQP(ngQP)   = QP gradient */
    lgbs = iw[291];
/* gBS(mBS)    = BS components of g */
    lkbs = iw[292];
/* kBS(mBS), ( B  S ) list */
    lrg = iw[293];
/* rg (maxS)   = reduced gradient */
    lrg2 = iw[294];
/* rg2(maxS)   = reduced gradient */
    lr = iw[295];
/* R(lenR)     = factor of Z'HZ */
    lscales = iw[296];
/* scales(nb)  = row and column scales */
    lgobj = iw[297];
/* gObj(nnObj) = Objective gradient */
    lx1 = iw[300];
/* x1(nb)      = new x, used to store x0 */
    lxbs = iw[301];
/* xBS(mBS)    = basics, superbasics */
    lxpen = iw[304];
/* xPen(nnCon) = penalty params */
    lxqp = iw[305];
/* xQP(nb)     = QP solution */
    lxqp0 = iw[306];
/* xQP0(nb)    = QP feasible pt. */
    liy = iw[308];
/* iy(nb)      =  integer work vector */
    liy1 = iw[309];
/* iy1(nb)     =  integer work vector */
    liy2 = iw[310];
/* iy2(nb)     =  integer work vector */
    ly = iw[311];
/* y(nb)       =  real work vector */
    ly1 = iw[312];
/* y1(nb)      =  real work vector */
    ly2 = iw[313];
/* y2(nb)      =  real work vector */
    lfcon = iw[316];
/* fCon (nnCon) constraints at x */
    lfcon1 = iw[317];
/* fCon1(nnCon) constraints at x1 */
    lfcon2 = iw[318];
/* fCon2(nnCon) work vector */
    lgconu = iw[319];
/* record of unknown derivatives and constants */
    lgcon = iw[320];
/* gCon (negCon)   constraint gradients at x */
    lgcon1 = iw[321];
/* gCon1(negCon)   constraint gradients at x1 */
    lgcon2 = iw[322];
/* gCon2(negCon)   work vector */
    lgobj1 = iw[324];
/* gObj1(nnObj) objective gradients at x1 */
    lgobj2 = iw[325];
/* gObj2(nnObj) work gObj */
    lfx = iw[336];
/* Fx (nnCon)  = F(x) + A(linear)x */
    lfv = iw[337];
/* Fv          = F(x) + A(linear)x - sN */
    ludx = iw[345];
/* Udx(nnH)      = product of U with dx */
    lhd = iw[347];
/* Diagonal of BFGS Hessian */
    lycon = iw[348];
/* yCon (nnCon)  = multipliers for fCon */
    lycon1 = iw[349];
/* yCon1(nnCon)  = yCon at x1 */
    lycon2 = iw[350];
/* yCon2(nnCon)  = work copy of yCon */
    ldycon = iw[351];
/* dyCon(nnCon)  = yCon1 - yCon */
    *iexit = 0;
    feasiblelc = FALSE_;
    gotfuns = FALSE_;
    fponly = iw[87] == 0;
    signobj = (doublereal) minimize;
    nonlinearcon = *nncon > 0;
    nonlinearobj = *nnobj > 0;
    nonlinear = *nnh > 0;
    nnobj0 = max(*nnobj,1);
    nncon0 = max(*nncon,1);
    nnh0 = max(*nnh,1);
    mbs = *m + maxs;
    nlocg = *nnjac + 1;
    numlc = *m - *nncon;
/*     Initialize yCon from pi. */
/*     Zap the pi(i) to prevent them being printed without being set. */
    if (nonlinearcon) {
	dcopy_(nncon, &pi[1], &c__1, &rw[lycon], &c__1);
    }
    if (numlc > 0) {
	dload_(&numlc, &c_b4, &pi[*nncon + 1], &c__1);
    }
/*     Initialize a few things. */
/*     Define the Hessian type for the QP subproblem. */
    if (iw[72] == 0 || iw[72] == 1) {
	if (*nnh < *n) {
	    iw[200] = 0;
	} else {
	    iw[200] = 1;
	}
    }
    iw[181] = 1;
    iw[210] = 0;
    iw[386] = 0;
    iw[387] = 0;
    *ninf = 0;
    *ninfe = 0;
    wtinf = 1.;
    pinorm = 1.;
    elastic = FALSE_;
    duinf = 0.;
    fmerit = 0.;
    *fobj = 0.;
    *sinf = 0.;
    *sinfe = 0.;
    maxvi = 0.;
    maxvirel = 0.;
    dload_(&c__4, &c_b4, penparm, &c__1);
    iw[220] = 0;
/* Line count for the print   file */
    iw[221] = 0;
/* Line count for the summary file */
    iload_(&c__4, &c__0, &iw[189], &c__1);
    iload_(&c__4, &c__0, &iw[194], &c__1);
    itn = 0;
    ndegen = 0;
    *majors = 0;
    minors = 0;
/*     Start recording the solve time. */
    s1page_(&c__1, &iw[1], leniw);
    s1time_(&c__2, &c__0, &iw[1], leniw, &rw[1], lenrw);
/*     Make copies of the upper and lower bounds. */
    dcopy_(nb, &bl[1], &c__1, &rw[lblqp], &c__1);
    dcopy_(nb, &bu[1], &c__1, &rw[lbuqp], &c__1);
/*     ------------------------------------------------------------------ */
/*     Print the matrix statistics before the nonlinear part of Jcol is */
/*     loaded with random elements. */
/*     Find the rowtypes for use in s5getB (they are held in iy2). */
/*     ------------------------------------------------------------------ */
    s2amat_(&c__1, &mjrprint, m, n, nb, nncon, nnjac, nnobj, iobj, &numlc, &
	    numliq, nej, nlocj, &locj[1], &indj[1], &jcol[1], &rw[lblqp], &rw[
	    lbuqp], &iw[liy2], &iw[1], leniw, &rw[1], lenrw);
/*     ------------------------------------------------------------------ */
/*     Make a permanent copy in gConU of the constant Jacobian elements */
/*     stored in J.  Load the nonlinear part of J with random elements. */
/*     ------------------------------------------------------------------ */
    if (nonlinearcon) {
	s8gcopy_(nncon, nnjac, nej, nlocj, &locj[1], &indj[1], nej, nlocj, &
		locj[1], &jcol[1], &negcon, &nlocg, &iw[llocg], &rw[lgconu]);
	s8rand_(&negcon, &negcon, &rw[lgcon]);
	s8gcopy_(nncon, nnjac, nej, nlocj, &locj[1], &indj[1], &negcon, &
		nlocg, &iw[llocg], &rw[lgcon], nej, nlocj, &locj[1], &jcol[1])
		;
    }
/* ================================================================= */
/* Find a basis kBS(1:m) for the linear constraints and bounds. */
/* ================================================================= */
/* s5getB does the following. */
/*  1. The linear constraints are (optionally) scaled. */
/*  2. The bounds bl and bu on the nonlinear rows are relaxed. */
/*  3. Elements x(n+1:n+m) of the initial x are assigned. */
/*  4. An LP is used to find a feasible x for the bounds and */
/*     linear equality constraints. */
/*  5. x(nb) is (optionally) scaled and saved in x1. */

/*  The linear constraints are unscaled in s8getFeasLC after a */
/*  feasible point for the linear constraints has been found. */
    nrhs = 0;
/* No QP rhs when finding the first basis */
    nrhs0 = 1;
    nx0 = *nb;
/* elements in x1(nb), the base point */
    lenx0 = *nb;
    iload_(nb, &c__0, &iw[letype], &c__1);
/* placeholder for s5getB */
    iload_(nb, &c__0, &iw[lestate], &c__1);
    s5getb_(&inform__, starttype, (U_fp)mnrlog, &needb, m, &maxs, &mbs, n, nb,
	     nncon, nnjac, nnobj, nnames, ns, &minors, &itnlim, &itn, &ndegen,
	     &numlc, &numliq, &toloptfp, &toloptqp, &tolx, ninf, sinf, &wtinf,
	     iobj, &scaleobj, &pinorm, &rgnorm, nej, nlocj, &locj[1], &indj[1]
	    , &jcol[1], &iw[letype], &iw[lestate], &iw[liy2], &iw[lfeastype], 
	    &hs[1], &iw[lkbs], names + 8, &bl[1], &bu[1], &rw[lblqp], &rw[
	    lbuqp], &rw[lblbs], &rw[lbubs], &rw[lblsave], &rw[lbusave], &rw[
	    lgbs], &pi[1], &rc[1], &nrhs0, &nrhs, &rw[lqprhs], &rw[lscales], &
	    lenx0, &nx0, &rw[lx1], &x[1], &rw[lxbs], &iw[liy], &iw[liy1], &rw[
	    ly], &rw[ly1], &rw[ly2], cw + 8, lencw, &iw[1], leniw, &rw[1], 
	    lenrw, (ftnlen)8, (ftnlen)8);
/*     Potential inform values are: */
/*        0   basis found. */
/*            nInf = 0 => linear equalities are     satisfied. */
/*            nInf > 0 +> linear equalities are not satisfied and the */
/*                        FP problem is unbounded. */
/*       >0   fatal error. No feasible point */
    if (inform__ > 0) {
	*iexit = inform__;
/* fatal error */
    }
/* ================================================================= */
/* An initial basis has been assigned to hs. The array kBS is */
/* defined if s5LP was used to find a feasible point for the linear */
/* equality constraints. Otherwise kBS is set after the first basis */
/* factorization. */

/* Find a feasible point for all the linear constraints. */
/* The norm of x is minimized via a proximal-point QP. */
/* If there is no feasible point, the linear rows can be elastic. */
/* ================================================================= */
    if (numlc > 0) {
	if (*iexit == 0) {
/* the E rows are feasible, check LG rows */
	    s8getfeaslc_(iexit, starttype, (U_fp)mnrlog, &lenr, m, &maxs, &
		    mbs, n, nb, &nncon0, nncon, &nnh0, nnh, &ndegen, ns, &
		    numlc, &numliq, &itn, &itnlim, &minors, &mnrprint, &
		    scaleobj, &toloptqp, &tolx, ninf, sinf, ninfe, sinfe, &
		    wtinf, &pinorm, &rgnorm, nej, nlocj, &locj[1], &indj[1], &
		    jcol[1], neh, nloch, &loch[1], &indh[1], &hcol[1], &iw[
		    letype], &iw[lestate], &iw[lfeastype], &hs[1], &iw[lkbs], 
		    &bl[1], &bu[1], &rw[lblqp], &rw[lbuqp], &rw[lblbs], &rw[
		    lbubs], &rw[lblsave], &rw[lbusave], &rw[lgbs], &rw[lgqp], 
		    &rw[lhdx], &rw[lpbs], &pi[1], &rw[lr], &rc[1], &rw[lrg], &
		    rw[lqprhs], &rw[lscales], &rw[lx1], &x[1], &rw[lxbs], &iw[
		    liy], &iw[liy1], &rw[ly], &rw[ly1], &rw[ly2], cw + 8, 
		    lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)8);
	}
/* -------------------------------------------------------------- */
/* Do some housekeeping in case we have to exit */
/* -------------------------------------------------------------- */
/* Reinstate the scaled bounds on the nonlinear rows. */
	if (nonlinearcon) {
	    dcopy_(nncon, &rw[lblsave + *n], &c__1, &bl[*n + 1], &c__1);
	    dcopy_(nncon, &rw[lbusave + *n], &c__1, &bu[*n + 1], &c__1);
	}
/* Unscale any linear constraints that were scaled in s5getB. */
/* NonlinearCon */
	if (iw[75] > 0) {
	    s2applyscales_(&c__1, m, n, nb, iobj, &infbnd, &scaleobj, nej, 
		    nlocj, &locj[1], &indj[1], &jcol[1], &rw[lscales], &bl[1],
		     &bu[1], &pi[1], &x[1]);
	}
    }
/*     Exit if the linear constraints are infeasible. */
/* numLC > 0 */
    feasiblelc = *iexit == 0;
    if (! feasiblelc) {
	goto L900;
    }
/*     ------------------------------------------------------------------ */
/*     Copy the constant Jacobian elements into gCon, gCon1 and gCon2. */
/*     Reset eType so that only nonlinear rows are elastic. */
/*     Make sure variables are not outside their bounds */
/*     (in particular, check the nonlinear slacks). */
/*     ------------------------------------------------------------------ */
    if (nonlinearcon) {
	dcopy_(&negcon, &rw[lgconu], &c__1, &rw[lgcon], &c__1);
	dcopy_(&negcon, &rw[lgconu], &c__1, &rw[lgcon1], &c__1);
	dcopy_(&negcon, &rw[lgconu], &c__1, &rw[lgcon2], &c__1);
	iload_(nncon, &c__3, &iw[letype + *n], &c__1);
    }
/* NonlinearCon */
    s5fixx_(&c__0, &c__1, nb, &tolx, &hs[1], &bl[1], &bu[1], &x[1]);
/*     ================================================================== */
/*     ================================================================== */
/*     The linear constraints have been satisfied! */
/*     Compute the problem functions at this all-important point. */
/*     No scaling yet. */
/*     Compute any missing derivatives. */
/*     ================================================================== */
/*     ================================================================== */
    if (nonlinear) {
	s8firstcall_(iexit, &feasiblelc, &gotfuns, &nonlinearcon, &
		nonlinearobj, n, nb, &nncon0, nncon, nnjac, nnh, &nnobj0, 
		nnobj, (S_fp)funwrapper, (U_fp)funcon, (U_fp)funobj, (U_fp)
		userhv, &bl[1], &bu[1], &x[1], &rw[lx1], &rw[lycon], nej, 
		nlocj, &locj[1], &indj[1], neh, nloch, &loch[1], &indh[1], &
		hcol[1], &negcon, &nlocg, &iw[llocg], fobj, &rw[lfcon], &rw[
		lgcon], &rw[lgobj], &rw[lfcon2], &rw[lgcon2], &rw[lgobj2], &
		rw[ly], &rw[ly1], cu + 8, lencu, &iu[1], leniu, &ru[1], lenru,
		 cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)8, (
		ftnlen)8);
/*        --------------------------------------------------------------- */
/*        Load the Jacobian gCon in  J. */
/*        --------------------------------------------------------------- */
	if (gotfuns && nonlinearcon) {
	    s8gcopy_(nncon, nnjac, nej, nlocj, &locj[1], &indj[1], &negcon, &
		    nlocg, &iw[llocg], &rw[lgcon], nej, nlocj, &locj[1], &
		    jcol[1]);
	    dcopy_(nncon, &rw[lfcon], &c__1, &rw[lfx], &c__1);
	}
	if (*iexit != 0) {
	    goto L900;
	}
    }
/*     ================================================================== */
/*     Scale the problem. */
/*     ================================================================== */
    if (iw[75] > 0) {
/*        --------------------------------------------------------------- */
/*        Recompute the vector of row types. */
/*        --------------------------------------------------------------- */
	s2amat_(&c__0, &mjrprint, m, n, nb, nncon, nnjac, nnobj, iobj, &numlc,
		 &numliq, nej, nlocj, &locj[1], &indj[1], &jcol[1], &bl[1], &
		bu[1], &iw[liy2], &iw[1], leniw, &rw[1], lenrw);
	s2getscales_(&mjrprint, m, n, nb, nnh, nncon, nnjac, &iw[liy2], nej, 
		nlocj, &locj[1], &indj[1], &jcol[1], &rw[lscales], &bl[1], &
		bu[1], &rw[ly], &rw[ly2], &iw[1], leniw, &rw[1], lenrw);
	s2applyscales_(&c__0, m, n, nb, iobj, &infbnd, &scaleobj, nej, nlocj, 
		&locj[1], &indj[1], &jcol[1], &rw[lscales], &bl[1], &bu[1], &
		pi[1], &x[1]);
/*        --------------------------------------------------------------- */
/*        The objective and constraint functions haven't been scaled yet. */
/*        Scale the constant elements in gCon1 and gCon2. */
/*        Don't forget the initial pi. */
/*        --------------------------------------------------------------- */
	if (nonlinearcon) {
	    dddiv_(nncon, &rw[lscales + *n], &c__1, &rw[lfcon], &c__1);
	    if (iw[184] > 0) {
		s8scalej_(nncon, nnjac, &negcon, n, &rw[lscales], nej, nlocj, 
			&locj[1], &indj[1], &rw[lgcon], &rw[1], lenrw);
		dcopy_(&negcon, &rw[lgcon], &c__1, &rw[lgcon1], &c__1);
		dcopy_(&negcon, &rw[lgcon], &c__1, &rw[lgcon2], &c__1);
	    }
	    ddscl_(nncon, &rw[lscales + *n], &c__1, &rw[lycon], &c__1);
	}
	if (nonlinearobj && iw[184] > 0) {
	    s8scaleg_(nnobj, &rw[lscales], &rw[lgobj], &rw[1], lenrw);
	}
    }
/*     ================================================================== */
/*     s8Fx computes the nonlinear constraint values Fx. */
/*     Copy these into the slacks x(n+i) and make sure they are feasible. */
/*     Crash uses them to decide which slacks to grab for the basis */
/*     If any nonbasic nonlinear slacks are close to a bound, */
/*     move them exactly onto the bound to avoid very small steps. */
/*     ================================================================== */
/* iw(lvlScale) > 0 */
    if (nonlinearcon) {
	s8fx_(n, nncon, nnjac, &eps0, nej, nlocj, &locj[1], &indj[1], &jcol[1]
		, &rw[lfcon], &x[1], &rw[lfx]);
	s2vmax_(n, nncon, &imaxvi, &maxvi, &bl[1], &bu[1], &rw[lfx]);
	supvi = vilim * max(10.,maxvi);
	dcopy_(nncon, &rw[lfx], &c__1, &x[*n + 1], &c__1);
	i__1 = *n + 1;
	i__2 = *n + *nncon;
	s5fixx_(&c__1, &i__1, &i__2, &tolx, &hs[1], &bl[1], &bu[1], &x[1]);
/*        =============================================================== */
/*        Crash on the nonlinear rows. */
/*        hs(*) already defines a basis for the full problem,  but we */
/*        want to do better by not including all of the slacks. */
/*        =============================================================== */
	if (needb) {
/*           Load  iy2  with the row types. */
/*           s2crash uses kBS as workspace.  It may alter x(n+i) for */
/*           nonlinear slacks. */
	    s2amat_(&c__0, &mjrprint, m, n, nb, nncon, nnjac, nnobj, iobj, &
		    numlc, &numliq, nej, nlocj, &locj[1], &indj[1], &jcol[1], 
		    &bl[1], &bu[1], &iw[liy2], &iw[1], leniw, &rw[1], lenrw);
	    lcrash = 5;
	    s2crash_(&lcrash, &mjrprint, m, n, nb, nncon, &icrash, &tcrash, 
		    nej, nlocj, &locj[1], &indj[1], &jcol[1], &iw[lkbs], &hs[
		    1], &iw[liy2], &bl[1], &bu[1], &x[1], &iw[1], leniw, &rw[
		    1], lenrw);
	    needb = FALSE_;
	}
/* NeedB */
    }
/*     ------------------------------------------------------------------ */
/*     Solve the problem. */
/*     ------------------------------------------------------------------ */
/* NonlinearCon */
    s1page_(&c__1, &iw[1], leniw);
    s8sqp_(iexit, (S_fp)funwrapper, (U_fp)funcon, (U_fp)funobj, (U_fp)userhv, 
	    (U_fp)mjrlog, (U_fp)mnrlog, (U_fp)snstop, &elastic, gotr, 
	    starttype, &itn, &lenr, m, &maxs, &mbs, n, nb, ns, &nncon0, nncon,
	     &nnobj0, nnobj, &nnh0, nnh, majors, &minors, &ndegen, &duinf, &
	    minimize, iobj, &scaleobj, objadd, fobj, &fmerit, &maxvi, &
	    maxvirel, &supvi, ninf, sinf, ninfe, sinfe, &wtinf0, &wtinf, 
	    penparm, &pinorm, &xnorm, nej, nlocj, &locj[1], &indj[1], &jcol[1]
	    , neh, nloch, &loch[1], &indh[1], &hcol[1], &negcon, &nlocg, &iw[
	    llocg], &iw[letype], &iw[lestate], &iw[lfeastype], &hs[1], &iw[
	    lkbs], &bl[1], &bu[1], &rw[lblqp], &rw[lbuqp], &rw[lblbs], &rw[
	    lbubs], &rw[lfv], &rw[lfx], &rw[lfcon], &rw[lgcon], &rw[lgobj], &
	    rw[lfcon1], &rw[lgcon1], &rw[lgobj1], &rw[lfcon2], &rw[lgcon2], &
	    rw[lgobj2], &rw[lgbs], &rw[lgqp], &rw[ldycon], &rw[ldx], &rw[ldg],
	     &rw[ludx], &rw[lhd], &rw[lhdx], &rw[lpbs], &rw[lycon], &rw[
	    lycon1], &rw[lycon2], &pi[1], &rw[lqprhs], &rw[lr], &rc[1], &rw[
	    lrg], &rw[lrg2], &rw[lscales], &x[1], &rw[lx1], &rw[lxbs], &rw[
	    lxqp0], &rw[lxqp], &rw[lxpen], &iw[liy], &iw[liy1], &rw[ly], &rw[
	    ly1], &rw[ly2], cu + 8, lencu, &iu[1], leniu, &ru[1], lenru, cw + 
	    8, lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)8, (ftnlen)8);
/*     ================================================================== */
/*     Exit. */
/*     Set output variables and print a summary of the final solution. */
/*     ================================================================== */
L900:
    snwrap_(iexit, solver, str, str2, &iw[1], leniw, (ftnlen)6, (ftnlen)133, (
	    ftnlen)133);
    s1time_(&c_n2, &c__0, &iw[1], leniw, &rw[1], lenrw);
/*     Set the objective function  (minimized or maximized) */
/*     objTrue = objAdd + fObj + x(n+iObj)*scaleObj */
    if (*iobj == 0) {
	objlin = *objadd;
    } else {
	objlin = *objadd + x[*n + *iobj] * scaleobj;
    }
    *objtrue = objlin;
    objmin = 0.;
    objinf = 0.;
    if (gotfuns) {
	if (elastic) {
	    if (fponly) {
		objmin = wtinf * *sinfe;
		objinf = *sinfe;
	    } else {
		if (*nnobj > 0) {
		    *objtrue += *fobj;
		}
		objmin = signobj * *objtrue + wtinf * *sinfe;
		objinf = *sinfe;
	    }
	} else {
/* Normal mode */
	    if (fponly) {
/*              Relax */
	    } else {
		if (*nnobj > 0) {
		    *objtrue += *fobj;
		}
		objmin = signobj * *objtrue;
	    }
	}
    }
/*     ------------------------------------------------------------------ */
/*     Print statistics. */
/*     ------------------------------------------------------------------ */
    degen = ndegen * 100. / max(itn,1);
    xnorm = dnormi_(n, &x[1], &c__1);
    pennorm = penparm[2];
/*     Count basic nonlinear variables (used only for printing). */
    nnb = 0;
    i__1 = *nnh;
    for (j = 1; j <= i__1; ++j) {
	if (hs[j] == 3) {
	    ++nnb;
	}
    }
    if (inewb > 0 && *iexit / 10 < 8) {
	k = *iexit / 10 + 1;
	s4stat_(&k, istate, (ftnlen)4);
	s4newb_(&c__1, &inewb, &minimize, m, n, nb, ns, &mbs, &itn, ninf, 
		sinf, fobj, &iw[lkbs], &hs[1], &rw[lscales], &bl[1], &bu[1], &
		x[1], &rw[lxbs], istate, cw + 8, lencw, &iw[1], leniw, (
		ftnlen)4, (ftnlen)8);
    }
/* Writing concatenation */
    i__3[0] = 30, a__1[0] = " Problem name                 ";
    i__3[1] = 8, a__1[1] = mprob;
    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)38);
    snprnt_(&c__13, ch__1, &iw[1], leniw, (ftnlen)38);
    if (gotfuns) {
	if (elastic) {
	    s_wsfi(&io___121);
	    do_fio(&c__1, (char *)&itn, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*objtrue), (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)133);
	    s_wsfi(&io___122);
	    do_fio(&c__1, (char *)&(*majors), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&objlin, (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)133);
	    s_wsfi(&io___123);
	    do_fio(&c__1, (char *)&pennorm, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*fobj), (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)133);
	    if (*ninfe > 0) {
/* objective is known */
		s_wsfi(&io___124);
		do_fio(&c__1, (char *)&wtinf, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&objmin, (ftnlen)sizeof(doublereal));
		e_wsfi();
		snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)133);
		s_wsfi(&io___125);
		do_fio(&c__1, (char *)&objinf, (ftnlen)sizeof(doublereal));
		e_wsfi();
		snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)133);
	    }
	} else {
/* Normal mode */
	    s_wsfi(&io___126);
	    do_fio(&c__1, (char *)&itn, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*objtrue), (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)133);
	    s_wsfi(&io___127);
	    do_fio(&c__1, (char *)&(*majors), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&objlin, (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)133);
	    if (*nncon > 0) {
/* Nonlinear constraints */
		if (fponly) {
		    s_wsfi(&io___128);
		    do_fio(&c__1, (char *)&pennorm, (ftnlen)sizeof(doublereal)
			    );
		    do_fio(&c__1, (char *)&(*fobj), (ftnlen)sizeof(doublereal)
			    );
		    e_wsfi();
		    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)133);
		} else {
		    s_wsfi(&io___129);
		    do_fio(&c__1, (char *)&pennorm, (ftnlen)sizeof(doublereal)
			    );
		    do_fio(&c__1, (char *)&(*fobj), (ftnlen)sizeof(doublereal)
			    );
		    e_wsfi();
		    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)133);
		}
	    } else {
		if (fponly) {
/*                 Relax. Nothing is printed */
		} else {
		    s_wsfi(&io___130);
		    do_fio(&c__1, (char *)&(*fobj), (ftnlen)sizeof(doublereal)
			    );
		    e_wsfi();
		    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)133);
		}
	    }
	}
	if (s_cmp(solver, "SNOPTA", (ftnlen)6, (ftnlen)6) == 0 || s_cmp(
		solver, "SNOPTC", (ftnlen)6, (ftnlen)6) == 0) {
/*           SNOPTA and SNOPTC call one user-supplied function */
	    if (iw[70] < 3 || lvlsrch == 0 && nonlinear) {
		s_wsfi(&io___131);
		do_fio(&c__1, (char *)&iw[194], (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&iw[195], (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)133);
	    } else {
		s_wsfi(&io___132);
		do_fio(&c__1, (char *)&iw[194], (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)133);
	    }
	    if (iw[70] < 3) {
		s_wsfi(&io___133);
		do_fio(&c__1, (char *)&iw[196], (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&iw[197], (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)133);
	    }
	} else {
/*           SNOPTB and NPOPT call separate obj and constraint functions. */
	    s_wsfi(&io___134);
	    do_fio(&c__1, (char *)&iw[194], (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&iw[189], (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)133);
	    if (iw[70] < 3 || lvlsrch == 0 && nonlinear) {
		s_wsfi(&io___135);
		do_fio(&c__1, (char *)&iw[195], (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&iw[190], (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)133);
	    }
	    if (iw[70] < 3) {
		if (iw[70] == 0) {
		    s_wsfi(&io___136);
		    do_fio(&c__1, (char *)&iw[196], (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&iw[191], (ftnlen)sizeof(integer));
		    e_wsfi();
		    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)133);
		    s_wsfi(&io___137);
		    do_fio(&c__1, (char *)&iw[197], (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&iw[192], (ftnlen)sizeof(integer));
		    e_wsfi();
		    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)133);
		} else if (iw[70] == 1) {
		    s_wsfi(&io___138);
		    do_fio(&c__1, (char *)&iw[191], (ftnlen)sizeof(integer));
		    e_wsfi();
		    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)133);
		    s_wsfi(&io___139);
		    do_fio(&c__1, (char *)&iw[192], (ftnlen)sizeof(integer));
		    e_wsfi();
		    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)133);
		} else if (iw[70] == 2) {
		    s_wsfi(&io___140);
		    do_fio(&c__1, (char *)&iw[196], (ftnlen)sizeof(integer));
		    e_wsfi();
		    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)133);
		    s_wsfi(&io___141);
		    do_fio(&c__1, (char *)&iw[197], (ftnlen)sizeof(integer));
		    e_wsfi();
		    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)133);
		}
	    }
	}
	if (*ns > 0) {
	    s_wsfi(&io___142);
	    do_fio(&c__1, (char *)&(*ns), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&nnb, (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)133);
	}
	if (iw[386] > 0) {
	    s_wsfi(&io___143);
	    do_fio(&c__1, (char *)&iw[386], (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)133);
	}
	s_wsfi(&io___144);
	do_fio(&c__1, (char *)&ndegen, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&degen, (ftnlen)sizeof(doublereal));
	e_wsfi();
	snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)133);
    } else {
/* No functions computed */
	if (*ninf > 0) {
	    s_wsfi(&io___145);
	    do_fio(&c__1, (char *)&(*ninf), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*sinf), (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)133);
	}
	s_wsfi(&io___146);
	do_fio(&c__1, (char *)&itn, (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)133);
	if (elastic) {
	    if (*ninfe > 0) {
		s_wsfi(&io___147);
		do_fio(&c__1, (char *)&(*ninfe), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&(*sinfe), (ftnlen)sizeof(doublereal));
		e_wsfi();
		snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)133);
	    }
	}
    }
/*     ------------------------------------------------------------------ */
/*     Unscale, compute nonlinear constraint violations, */
/*     save basis files and prepare to print the solution. */
/*     Clock 3 is "Output time". */
/*     ------------------------------------------------------------------ */
    s1time_(&c__3, &c__0, &iw[1], leniw, &rw[1], lenrw);
/*     Skip the functions if we don't have them. */
/*     Skip unscaling everything for infeasible linear constraints, */
/*     they have already been unscaled. */
    lssave = iw[75];
    if (! gotfuns) {
	nncon1 = 0;
	nnobj1 = 0;
	if (! feasiblelc) {
	    iw[75] = 0;
	}
    } else {
	nncon1 = *nncon;
	nnobj1 = *nnobj;
    }
    s4saveb_(iexit, &c__0, &minimize, m, n, nb, &nkx, &nncon0, &nncon1, &nnh0,
	     &nnobj1, nnames, ns, &itn, ninf, sinf, &wtinf, &maxvi, iobj, &
	    scaleobj, objtrue, &pnorm1, &pnorm2, &pinorm, &xnorm, nej, nlocj, 
	    &locj[1], &indj[1], &jcol[1], &iw[lkx], &iw[lestate], &hs[1], &rw[
	    lscales], &bl[1], &bu[1], &rw[lfx], &rw[lgobj], names + 8, &pi[1],
	     &rc[1], &x[1], cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (
	    ftnlen)8, (ftnlen)8);
/*     If task = 'Print', s4saveB prints the solution under the control */
/*     of lprSol (set by the  Solution  keyword in the SPECS file). */
/*     The printed solution may or may not be wanted, as follows: */

/*     lprSol = 0   means      No */
/*            = 2   means      Yes */
    s4saveb_(iexit, &c__1, &minimize, m, n, nb, &nkx, &nncon0, &nncon1, &nnh0,
	     nnobj, nnames, ns, &itn, ninf, sinf, &wtinf, &maxvi, iobj, &
	    scaleobj, objtrue, &pnorm1, &pnorm2, &pinorm, &xnorm, nej, nlocj, 
	    &locj[1], &indj[1], &jcol[1], &iw[lkx], &iw[lestate], &hs[1], &rw[
	    lscales], &bl[1], &bu[1], &rw[lfx], &rw[lgobj], names + 8, &pi[1],
	     &rc[1], &x[1], cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (
	    ftnlen)8, (ftnlen)8);
    iw[75] = lssave;
    s1time_(&c_n3, &c__0, &iw[1], leniw, &rw[1], lenrw);
/*     ------------------------------------------------------------------ */
/*     If the user hasn't already pulled the plug, */
/*     call the functions one last time with  Status .ge. 2. */
/*     Everything has been  unscaled, so we have to disable scaling. */
/*     modefg = 0  requests that no gradients are computed. */
/*     ------------------------------------------------------------------ */
    *fobj = objlin;
    if (gotfuns && *iexit / 10 != 7 && *iexit / 10 != 6) {
/* Computing MIN */
	i__1 = *iexit / 10;
	iw[236] = min(i__1,4) + 2;
	modefg = 0;
	lssave = iw[75];
	iw[75] = 0;
	(*funwrapper)(&inform__, &modefg, &nonlinearcon, &nonlinearobj, n, &
		negcon, &nncon0, nncon, nnjac, nnh, &nnobj0, nnobj, (U_fp)
		funcon, (U_fp)funobj, (U_fp)userhv, &x[1], &rw[lycon], nej, 
		nlocj, &locj[1], &indj[1], neh, nloch, &loch[1], &indh[1], &
		hcol[1], &rw[lfcon], fobj, &rw[lgcon2], &rw[lgobj2], cu + 8, 
		lencu, &iu[1], leniu, &ru[1], lenru, cw + 8, lencw, &iw[1], 
		leniw, &rw[1], lenrw, (ftnlen)8, (ftnlen)8);
	iw[75] = lssave;
    }
/* Save some things needed by solvers calling SNOPT */
    rw[421] = *objtrue;
/* The objective (minimized or maximized) */
    rw[422] = pinorm;
/* Lagrange multiplier norm */
    rw[423] = xnorm;
/* Norm of the variables (for GAMS) */
    rw[424] = wtinf;
/* Infeasibility weight */
    rw[433] = *sinf + *sinfe;
/* Sum of infeasibilities */
    rw[434] = objlin;
/* Linear    objective term */
    iw[421] = itn;
/* Total iteration count */
    iw[423] = maxs;
/*     Nonlinear constraint information. */
/* max # of superbasics */
    rw[432] = maxvi;
/* Inf norm of the constraint violation */
    rw[435] = *fobj;
/* Nonlinear objective term */
    rw[436] = penparm[2];
/* Norm of penalty parameters */
    iw[422] = *majors;
/* Major iterations */
    return 0;
} /* s8solve_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8solve */
/* Subroutine */ int s8defaults_(integer *m, integer *n, integer *nncon, 
	integer *nnjac, integer *nnobj, integer *iobj, char *cw, integer *
	lencw, integer *iw, integer *leniw, doublereal *rw, integer *lenrw, 
	ftnlen cw_len)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    double sqrt(doublereal), pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static doublereal hcondbnd;
    static integer lvlscale;
    static doublereal condumax;
    static integer minimize, nparprlp, nparprqp, mjrprint, stickyop;
    static doublereal toloptfp;
    static integer mnrprint, qpsolver;
    static doublereal toloptnp, toloptqp, c4, c6;
    static logical nonlinear;
    static doublereal proxweight;
    static integer nnh;
    static doublereal eps, eps0, eps1, eps2, eps3, eps4;
    static integer npr1, npr2, kfac;
    static char mbnd[8];
    static integer kchk;
    static char mobj[8];
    static integer klog;
    static char mrng[8];
    static integer ksav, maxr, maxs;
    static char mrhs[8];
    static integer nout;
    static doublereal tolx, dens1, dens2, lmax1;
    static logical nonlinearcon;
    static doublereal lmax2, xpen0, utol1, utol2;
    static integer iback, emode, ioldb;
    static doublereal bigdx, bigfx;
    static integer ipnch;
    static doublereal etarg;
    static integer inewb;
    static doublereal small, tolcg, xdlim;
    static char mprob[8];
    static integer idump, maxmn, never, mskip, isoln;
    static doublereal epsrf, vilim;
    static integer ksumm;
    static doublereal wtmax, fdint1, fdint2;
    static integer jverf1, jverf2, jverf3, jverf4, jverf5, jverf6;
    static doublereal toldj3, wtinf0, utol1m, utol2m;
    static integer iloadb, kdegen;
    static doublereal infbnd, zcndbd, chzbnd;
    static integer icrash;
    static logical linear;
    static integer lprdbg;
    static doublereal tolfac, uspace;
    static logical lincon;
    static integer maxcol;
    static doublereal tcrash;
    static integer mmajor;
    static doublereal toldcp, tolddp;
    static integer lvlder, minmax, minprc, cgitmx, itnlim, deropt, kreset, 
	    lprsch, lprscl, mflush, mminor, lvlpre, iprint, ireprt, mqnmod, 
	    iinsrt, mnewsb;
    static char solver[8];
    static integer lprsol, lprprm, lvlpiv, lvlppm, lvlver, objrow;
    static doublereal scltol, tolcon, toldpp, toldrp;
    static integer tpivot;
    static doublereal toldup, tolpiv, tolrow, tolswp, tolupd, wolfeg;
    static integer lvlsys, lvlobje;
    static doublereal maxtime;
    static integer lvlsrch, lvlhess, nparpru;

/*     ================================================================== */
/*     s8defaults checks the optional parameter values and possibly */
/*     changes them to reasonable values. */

/*     Note that checking occurs before the amount of working storage has */
/*     been defined. */

/*     See  snworkspace.info  for full documentation of cw, iw and rw. */

/*     15 Nov 1991: first version. */
/*     27 Apr 2001: wtMax introduced. */
/*     10 Dec 2002: Added defaults for LU Rook and Diagonal Pivoting. */
/*     31 Dec 2002: Added default MPS character names. */
/*     30 Jul 2003: Added default CG tolerance. */
/*     22 Jun 2004: Added default LU mod singularity tol */
/*     20 Dec 2004: Default LU tols reduced. */
/*     21 Dec 2004: Default LU tols fixed up. */
/*     09 Jun 2005: Default tolpiv same as in s5dflt. */
/*     02 Jul 2005: Default Utol's set back to eps1. */
/*     02 May 2006: lvlTim removed. */
/*     01 Sep 2007: stickyOp added. */
/*     25 Nov 2007: Hessian options added. */
/*     05 Apr 2014: Objective derivatives not verified in FP mode. */
/*     12 Sep 2014: FPonly set if nnObj .eq. 0 and iObj .eq. 0 */
/*     25 Oct 2014: Bogus assignment of LUprnt removed (its set in s2BLU). */
/*     25 Oct 2014: nout set independently of lvlsys. */
/*     15 Nov 2014: No scaling is now the default. */
/*     21 Nov 2014: Default wtInf now 10^5. */
/*     18 Feb 2015: SNOPTA parameters set separately. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
/*     Set some local machine-dependent constants. */
    /* Parameter adjustments */
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    eps = rw[1];
/* machine precision.  IEEE DP  2.22e-16 */
    eps0 = rw[2];
/* eps**(4/5)          IEEE DP  3.00e-13 */
    eps1 = rw[3];
/* eps**(2/3)          IEEE DP  3.67e-11 */
    eps2 = rw[4];
/* eps**(1/2)          IEEE DP  1.49e-08 */
    eps3 = rw[5];
/* eps**(1/3)          IEEE DP  6.05e-06 */
    eps4 = rw[6];
/* eps**(1/4)          IEEE DP  1.22e-04 */
    s_copy(solver, cw + 8, (ftnlen)8, (ftnlen)8);
/*     ------------------------------------------------------------------ */
/*     rw(51)--rw(150): optional parameters set via the specs file. */
/*     ------------------------------------------------------------------ */
    toloptfp = rw[51];
/* Minor Phase 1 Opt tol */
    toloptqp = rw[52];
/* Minor Phase 2 Opt tol */
    toloptnp = rw[53];
/* Major Optimality tolerance */
    tolcg = rw[54];
/* cg tolerance */
    tolx = rw[56];
/* Minor feasibility tolerance. */
    tolcon = rw[57];
/* Major feasibility tolerance. */
    tolpiv = rw[60];
/* excludes small elements of y */
    tolrow = rw[61];
/* tolerance for the row error */
    tcrash = rw[62];
/* crash tolerance */
    utol1m = rw[63];
/* abs tol for small diag of U in LU mod */
    utol2m = rw[64];
/* rel tol for small diag of U in LU mod */
    tolswp = rw[65];
/* LU swap tolerance */
    tolfac = rw[66];
/* User-defined LU factor tolerance */
    tolupd = rw[67];
/* User-defined LU update tolerance */
    infbnd = rw[70];
/* definition of an infinite bound */
    bigfx = rw[71];
/* unbounded objective */
    bigdx = rw[72];
/* unbounded step */
    epsrf = rw[73];
/* relative function precision. */
    fdint1 = rw[76];
/* (1) forwrd diff. interval */
    fdint2 = rw[77];
/* (2) cntrl  diff. interval */
    maxtime = rw[79];
/* Time limit */
    xdlim = rw[80];
/* Step limit */
    vilim = rw[81];
/* violation limit */
    etarg = rw[83];
/* Quasi-Newton QP rg tolerance */
    wolfeg = rw[84];
/* line search tolerance. */
    hcondbnd = rw[85];
/* bound on the condition of Hz */
    zcndbd = rw[86];
/* bound on the condition of Z */
    condumax = rw[87];
/* max cond estimator for U with H = U'U */
    wtinf0 = rw[88];
/* infeasibility weight */
    xpen0 = rw[89];
/* initial penalty parameter. */
    wtmax = rw[90];
/* max     infeasibility weight */
    proxweight = rw[91];
/* Proximal-point weight */
    scltol = rw[92];
/*     ------------------------------------------------------------------ */
/*     rw(151)--rw(180) are parmLU parameters for LUSOL (some optional). */
/*     ------------------------------------------------------------------ */
/* scale tolerance. */
    small = rw[153];
/* defn of small real. */
    utol1 = rw[154];
/* abs tol for small diag of U. */
    utol2 = rw[155];
/* rel tol for small diag of U. */
    uspace = rw[156];
/* limit on waste space in U. */
    dens1 = rw[157];
/* switch to search maxcol columns and no rows */
    dens2 = rw[158];
/*     ------------------------------------------------------------------ */
/*     rw(181)--rw(199) pass parameters into various routines. */
/*     ------------------------------------------------------------------ */
/*     toldj3     = rw(186) ! current optimality tol */
/*     ------------------------------------------------------------------ */
/*     iw(1)--iw(50): I/O file numbers and dimensions. */
/*     ------------------------------------------------------------------ */
/* switch to dense LU. */
    iprint = iw[12];
/*     ------------------------------------------------------------------ */
/*     iw(51)--iw(150): optional parameters set via the specs file. */
/*     ------------------------------------------------------------------ */
/* Print file */
    maxr = iw[52];
/* max columns of R. */
    maxs = iw[53];
/* max # of superbasics */
    mqnmod = iw[54];
/* (ge 0) max # of BFGS updates */
    qpsolver = iw[55];
/* 0(1) => QP(QN) QP solver */
    emode = iw[56];
/* >0    => use elastic mode */
    kchk = iw[58];
/* check (row) frequency */
    kfac = iw[59];
/* factorization frequency */
    ksav = iw[60];
/* save basis map */
    klog = iw[61];
/* log/print frequency */
    ksumm = iw[62];
/* Summary print frequency */
    kdegen = iw[63];
/* max. expansions of featol */
    kreset = iw[64];
/* Hessian frequency */
    mflush = iw[66];
/* Hessian flush */
    mskip = iw[67];
/*     lvlStart   = iw( 69) ! = 0:1:2:3 => cold:basis:warm:hot start */
/* # largest value of nSkip */
    lvlder = iw[70];
/* = 0, 1 or 2, the derivative level */
    lvlsys = iw[71];
/* > 0   => print system info */
    lvlhess = iw[72];
/* 0,1,2  => LM, FM, Exact Hessian */
    lvlobje = iw[73];
/* Elastic option */
    lvlscale = iw[75];
/* scale option */
    lvlsrch = iw[76];
/* >0     => use derivatives in the line searc */
    lvlpre = iw[77];
/* >0    => QN preconditioned CG */
    lvlver = iw[78];
/* Verify level */
    lvlppm = iw[79];
/* 1(2)-norm proximal point method for x0 */
    lvlpiv = iw[80];
/* 0(1 2 3) LU partial(rook complete diagonal) */
    lprprm = iw[81];
/* > 0    => parms are printed */
    lprsch = iw[82];
/* line search debug starting itn */
    lprscl = iw[83];
/* > 0    => print the scales */
    lprsol = iw[84];
/* > 0    => print the solution */
    lprdbg = iw[85];
/* > 0    => private debug print */
    minmax = iw[87];
/* 1, 0, -1  => MIN, FP, MAX */
    icrash = iw[88];
/* Crash option */
    itnlim = iw[89];
/* limit on total iterations */
    mmajor = iw[90];
/* limit on major iterations */
    mminor = iw[91];
/* limit on minor iterations */
    mjrprint = iw[92];
/* Major print level */
    mnrprint = iw[93];
/* Minor print level */
    mnewsb = iw[95];
/* maximum # of new superbasics per major */
    cgitmx = iw[97];
/* CG iteration limit */
    nparpru = iw[101];
/* # of partial pricing sections */
    objrow = iw[103];
/* Objective row of user-defined F */
    deropt = iw[104];
/* 0, 1, 2 => derivative option */
    jverf1 = iw[110];
/* col # to start derivative checking */
    jverf2 = iw[111];
/* col # to stop  derivative checking */
    jverf3 = iw[112];
/* col # to start derivative checking */
    jverf4 = iw[113];
/* col # to stop  derivative checking */
    jverf5 = iw[114];
/* start col for Hessian checking */
    jverf6 = iw[115];
/* stop  col for Hessian checking */
    stickyop = iw[116];
/* > 0 => optional parameters are sticky */
    iback = iw[120];
/* backup file */
    idump = iw[121];
/* dump file */
    iloadb = iw[122];
/* load file */
    inewb = iw[124];
/* new basis file */
    iinsrt = iw[125];
/* insert file */
    ioldb = iw[126];
/* old basis file */
    ipnch = iw[127];
/* punch file */
    ireprt = iw[130];
/* report file */
    isoln = iw[131];
/*     ------------------------------------------------------------------ */
/*     iw(151)--iw(180) are luparm parameters for LUSOL (some optional). */
/*     ------------------------------------------------------------------ */
/* solution file */
    nout = iw[151];
/* unit # for printed messages */
    maxcol = iw[153];
/*     ------------------------------------------------------------------ */
/*     Character  workspace. */
/*     cw(51)--cw(150): optional parameters */
/*     ------------------------------------------------------------------ */
/* lu1fac: max. # columns */
    s_copy(mprob, cw + 408, (ftnlen)8, (ftnlen)8);
/* Problem name */
    s_copy(mobj, cw + 416, (ftnlen)8, (ftnlen)8);
/* Objective name */
    s_copy(mrhs, cw + 424, (ftnlen)8, (ftnlen)8);
/* rhs name */
    s_copy(mrng, cw + 432, (ftnlen)8, (ftnlen)8);
/* range name */
    s_copy(mbnd, cw + 440, (ftnlen)8, (ftnlen)8);
/*     ------------------------------------------------------------------ */
/* bounds name */
    c4 = max(1e-4,eps3);
    c6 = max(1e-6,eps2);
/* Computing MAX */
    d__1 = 1. / (eps * 100.);
    chzbnd = max(d__1,1e6);
    never = 99999999;
/*     =============================================================== */
/*     Check the optional parameters. */
/*     =============================================================== */
    if (*nncon == 0) {
	*nnjac = 0;
    }
    if (*nnjac == 0) {
	*nncon = 0;
    }
    nnh = max(*nnjac,*nnobj);
    lincon = *nncon == 0;
    nonlinearcon = *nncon > 0;
    linear = nnh == 0;
    nonlinear = nnh > 0;
    if (iback == -11111) {
	iback = 0;
    }
    if (idump == -11111) {
	idump = 0;
    }
    if (iloadb == -11111) {
	iloadb = 0;
    }
    if (inewb == -11111) {
	inewb = 0;
    }
    if (iinsrt == -11111) {
	iinsrt = 0;
    }
    if (ioldb == -11111) {
	ioldb = 0;
    }
    if (ipnch == -11111) {
	ipnch = 0;
    }
    if (ireprt == -11111) {
	ireprt = 0;
    }
    if (isoln == -11111) {
	isoln = 0;
    }
/*     Set unspecified frequencies or silly values to defaults. */
    if (kchk == -11111) {
	kchk = 60;
    }
    if (kfac <= 0) {
	kfac = 100;
	if (nonlinear) {
	    kfac = 50;
	}
    }
    if (klog == -11111) {
	klog = 100;
    }
    if (ksumm == -11111) {
	ksumm = 100;
    }
    if (ksav == -11111) {
	ksav = 100;
    }
    if (kdegen == -11111) {
	kdegen = 10000;
    }
    if (mflush == -11111) {
	mflush = 0;
    }
/*     Sometimes, frequency 0 means "almost never". */
    if (kchk <= 0) {
	kchk = never;
    }
    if (mflush <= 0) {
	mflush = never;
    }
    if (klog <= 0) {
	klog = never;
    }
    if (ksav <= 0) {
	ksav = never;
    }
    if (ksumm <= 0) {
	ksumm = never;
    }
    if (kdegen <= 0) {
	kdegen = never;
    }
    if (kreset <= 0) {
	kreset = never;
    }
    if (icrash < 0) {
	icrash = 3;
    }
    if (minmax == -11111) {
	minmax = 1;
    }
    if (*nnobj == 0 && *iobj == 0 || objrow == 0) {
	minmax = 0;
    }
    if (minmax == -1) {
	minimize = -1;
    } else {
	minimize = 1;
    }
    if (objrow > 0 && minmax == 0) {
	objrow = 0;
    }
    if (mjrprint == -11111) {
	mjrprint = 1;
    }
    if (mnrprint == -11111) {
	mnrprint = 1;
    }
/*     if (mMinor   .lt. 0    ) mMinor   = max( 1000,5*max( n,m ) ) */
    if (mminor < 0) {
	mminor = 500;
    }
    if (mmajor < 0) {
/* Computing MAX */
	i__1 = 1000, i__2 = max(*n,*m) * 3;
	mmajor = max(i__1,i__2);
    }
    if (mskip < 0 && lincon) {
	mskip = never;
    }
    if (mskip < 0 && nonlinearcon) {
	mskip = 2;
    }
    if (mnewsb <= 0) {
	mnewsb = 99;
    }
    if (lprdbg < 0) {
	lprdbg = 0;
    }
    if (lprprm < 0) {
	lprprm = 1;
    }
    if (lprsch < 0) {
	lprsch = never;
    }
    if (lprscl < 0) {
	lprscl = 0;
    }
    if (lprsol < 0) {
	lprsol = 2;
    }
/*     lvlStart is checked in s3argA or s3argB */
/*     if (lvlStart.lt. 0     ) lvlStart =  0 */
    if (s_cmp(solver, "SNOPTA", (ftnlen)6, (ftnlen)6) == 0) {
	if (deropt != -11111) {
	    if (deropt < 0 || deropt > 1) {
		deropt = 1;
	    }
	    if (deropt == 0) {
		lvlder = 0;
	    }
	    if (deropt == 1) {
		lvlder = 3;
	    }
	} else {
	    deropt = 1;
	    lvlder = 3;
	}
	if (lvlsrch < 0) {
	    if (deropt == 0) {
		lvlsrch = 0;
	    }
	    if (deropt == 1) {
		lvlsrch = 1;
	    }
	}
    } else {
	if (lvlder != -11111) {
	    if (lvlder < 0 || lvlder > 3) {
		lvlder = 3;
	    }
	} else {
	    lvlder = 3;
	}
	if (lvlsrch < 0) {
	    if (lvlder != 3) {
		lvlsrch = 0;
	    }
	    if (lvlsrch < 0) {
		lvlsrch = 1;
	    }
	}
    }
    if (lvlver == -11111) {
	lvlver = 0;
    }
    if (lvlver < 0) {
	lvlver = -1;
    }
    if (lvlver > 3) {
	lvlver = 0;
    }
    lvlobje = 2;
    if (minmax == 0) {
	if (lvlder == 2) {
	    lvlder = 3;
	}
	if (lvlder == 0) {
	    lvlder = 1;
	}
	if (lvlver == 1) {
	    lvlver = 0;
	}
	if (lvlver == 3) {
	    lvlver = 2;
	}
    }
/*     Check  START and STOP  column numbers for derivative checking. */
    if (s_cmp(solver, "SNOPTA", (ftnlen)6, (ftnlen)6) == 0) {
	jverf1 = 1;
	jverf2 = *n;
	if (lvlver == 2 || lvlver == 0) {
	    jverf2 = 0;
	}
	jverf3 = 1;
	jverf4 = *n;
	if (lvlver == 1 || lvlver == 0) {
	    jverf4 = 0;
	}
	jverf5 = 1;
	jverf6 = *n;
	if (lvlver == 1 || lvlver == 0) {
	    jverf6 = 0;
	}
    } else {
	if (jverf1 <= 0) {
	    jverf1 = 1;
	}
	if (jverf2 < 0) {
	    jverf2 = *n;
	}
	if (lvlver == 2 || lvlver == 0) {
	    jverf2 = 0;
	}
	if (jverf3 <= 0) {
	    jverf3 = 1;
	}
	if (jverf4 < 0) {
	    jverf4 = *n;
	}
	if (lvlver == 1 || lvlver == 0) {
	    jverf4 = 0;
	}
	if (jverf5 <= 0) {
	    jverf5 = 1;
	}
	if (jverf6 < 0) {
	    jverf6 = *n;
	}
	if (lvlver == 1 || lvlver == 0) {
	    jverf6 = 0;
	}
    }
    if (lvlsys < 0) {
	lvlsys = 0;
    }
    if (lvlppm < 0) {
	lvlppm = 1;
    }
    emode = 1;
    if (stickyop < 0) {
	stickyop = 0;
    }
/*     Check superbasics limit maxS and Hessian dimension maxR. */
    if (nonlinear) {
	if (maxr < 0) {
/* Computing MIN */
	    i__1 = 2000, i__2 = nnh + 1;
	    maxr = min(i__1,i__2);
	}
	if (maxs < 0) {
	    maxs = nnh + 1;
	}
/* Computing MAX */
	i__1 = min(maxr,*n);
	maxr = max(i__1,0);
/* Computing MAX */
	i__1 = min(maxs,*n);
	maxs = max(i__1,1);
    } else {
/* Linear */
	if (maxs <= 0) {
	    maxs = 1;
	}
	if (maxr <= 0) {
	    maxr = 1;
	}
    }
    if (maxs < maxr) {
	maxr = maxs;
    }
    if (qpsolver < 0) {
	qpsolver = 0;
    }
    if (maxr == 0) {
	qpsolver = 1;
    }
    if (qpsolver == 2) {
	lvlpre = 0;
    }
    if (qpsolver == 0) {
	lvlpre = 0;
    }
    if (lvlpre > 1) {
	lvlpre = 1;
    }
    if (lvlpre < 0 && qpsolver == 1) {
	lvlpre = 0;
    }
    if (cgitmx < 0) {
	cgitmx = 100;
    }
    if (qpsolver == 1 || maxr < maxs) {
	if (lvlhess < 0) {
	    lvlhess = 0;
	}
	if (mqnmod < 0) {
	    mqnmod = 10;
	}
    } else {
	if (lvlhess < 0 && nnh > 75) {
	    lvlhess = 0;
	}
	if (lvlhess < 0 && nnh <= 75) {
	    lvlhess = 1;
	}
	if (lvlhess == 1) {
	    mqnmod = kreset;
	}
	if (mqnmod < 0) {
	    mqnmod = 10;
	}
    }
/*     --------------------------------- */
/*     CG QP optional parameters */
/*     --------------------------------- */
    if (etarg < 0. || etarg > 1.) {
	etarg = .1;
    }
/*     Check other options. */
/*     if (lvlScale .lt. 0  .and. */
/*    &    nnCon    .eq. 0      ) lvlScale   =  2 */
/*     if (lvlScale .lt. 0      ) lvlScale   =  1 */
/*                                lvlScale   =  min( lvlScale, 2 ) */
/*     if (lvlScale .eq. 1  .and.  nnJac .ge. n) */
/*    &                           lvlScale   = 0 */
    if (lvlscale < 0) {
	lvlscale = 0;
    }
    lvlscale = min(lvlscale,2);
    if (lvlscale == 1 && *nnjac >= *n) {
	lvlscale = 0;
    }
/*     Partial pricing parameters */
    minprc = 10;
    maxmn = max(*m,*n);
/*     Temporarily go back to what we had before. */
/*     Set nParPrU to previous value. */
/*     Then nParPrLP/nParPrQP will take that value. */
    if (nparpru <= 0) {
	nparpru = 10;
	if (nonlinear) {
	    nparpru = 1;
	}
    }
    minprc = 10;
    npr1 = *n / nparpru;
    npr2 = *m / nparpru;
    if (max(npr1,npr2) < minprc) {
	maxmn = max(*m,*n);
	nparpru = maxmn / min(maxmn,minprc);
    }
    nparprlp = nparpru;
    nparprqp = nparpru;
/* $$$!     The new part follows */
/* $$$      if (nParPrU .le. 0     ) then */
/* $$$         if (nParPrLP .le. 0 ) nParPrLP = 10 */
/* $$$                               nPr1     = n / nParPrLP */
/* $$$                               nPr2     = m / nParPrLP */
/* $$$         if (max( nPr1, nPr2 ) .lt. minPrc) then */
/* $$$                               nParPrLP = maxmn / min( maxmn, minPrc ) */
/* $$$         end if */
/* $$$         if (nParPrQP .le. 0 ) nParPrQP = 1 */
/* $$$                               nPr1     = n / nParPrQP */
/* $$$                               nPr2     = m / nParPrQP */
/* $$$         if (max( nPr1, nPr2 ) .lt. minPrc) then */
/* $$$                               nParPrQP = maxmn / min( maxmn, minPrc ) */
/* $$$         end if */
/* $$$      else */
/* $$$         nParPrLP = nParPrU */
/* $$$         nParPrQP = nParPrU */
/* $$$      end if */
    if (proxweight < 0.) {
	proxweight = 1e-4;
    }
    if (maxtime < 0.) {
	maxtime = 0.;
    }
    if (infbnd < 0.) {
	infbnd = 1e20;
    }
    if (bigfx <= 0.) {
	bigfx = 1e15;
    }
    if (bigdx <= 0.) {
	bigdx = infbnd;
    }
    if (hcondbnd <= 0.) {
	hcondbnd = chzbnd;
    }
    if (condumax <= 0.) {
	condumax = sqrt(chzbnd);
    }
    if (tcrash < 0. || tcrash >= 1.) {
	tcrash = .1;
    }
    if (vilim <= 0.) {
	vilim = 10.;
    }
    if (wolfeg < 0. || wolfeg > 1.) {
	wolfeg = .9;
    }
    if (wtmax < 0.) {
	wtmax = 1e10;
    }
    if (xdlim <= 0.) {
	xdlim = 2.;
    }
    if (xpen0 < 0.) {
	xpen0 = 0.;
    }
    if (zcndbd <= 0.) {
	if (qpsolver == 0) {
	    zcndbd = 1e4;
	} else {
	    zcndbd = 1e6;
	}
    }
/*     ---------------------------------------- */
/*     Set up the parameters for lu1fac. */
/*     LUprnt > 0 gives LU output on unit nout. */
/*     LUprnt is set in s2BLU. */
/*     ---------------------------------------- */
    if (maxcol < 0) {
	maxcol = 5;
    }
    if (nout == -11111) {
	nout = iprint;
    }
    if (lvlpiv <= 0) {
	lvlpiv = 0;
    }
    if (lvlpiv > 3) {
	lvlpiv = 0;
    }
    tpivot = lvlpiv;
    if (linear) {
	toldpp = 100.;
	toldrp = 10.;
	toldcp = 10.;
	tolddp = 10.;
	toldup = 10.;
    } else {
/* nonlinear */
	toldpp = 3.99;
	toldrp = 3.99;
	toldcp = 3.99;
	tolddp = 3.99;
	toldup = 3.99;
    }
    if (tolfac < 1.) {
	if (lvlpiv == 0) {
	    tolfac = toldpp;
	}
	if (lvlpiv == 1) {
	    tolfac = toldrp;
	}
	if (lvlpiv == 2) {
	    tolfac = toldcp;
	}
	if (lvlpiv == 3) {
	    tolfac = tolddp;
	}
    }
    if (tolupd < 1.) {
	tolupd = toldup;
    }
    lmax1 = tolfac;
    lmax2 = tolupd;
    if (utol1 <= 0.) {
	utol1 = eps1;
    }
    if (utol2 <= 0.) {
	utol2 = eps1;
    }
    if (utol1m <= 0.) {
	utol1m = eps1;
    }
    if (utol2m <= 0.) {
	utol2m = eps1;
    }
    if (dens2 < 0.) {
	dens2 = .6;
    }
    if (small <= 0.) {
	small = eps0;
    }
    if (uspace <= 0.) {
	uspace = 3.;
    }
    if (dens1 <= 0.) {
	dens1 = .3;
    }
/*     Set some SQP tolerances. */
/*     Set the minor and major optimality tolerances. */
/*     Solve the QP subproblems fairly accurately even if the */
/*     NLP Optimality Tolerance is big. */
    if (toloptnp <= 0.) {
	toloptnp = c6 * 2.;
	if (epsrf > 0.) {
/* Computing MAX */
	    d__1 = toloptnp, d__2 = sqrt(epsrf * 10.);
	    toloptnp = max(d__1,d__2);
	}
    }
    if (toloptqp <= 0.) {
/* Computing MIN */
	d__1 = c6, d__2 = toloptnp / 2.;
	toloptqp = min(d__1,d__2);
    }
    if (toloptfp < 0.) {
	toloptfp = c6;
    }
    if (tolcg <= 0.) {
	tolcg = .01;
    }
    if (tolrow <= 0.) {
	tolrow = c4;
    }
    if (tolswp <= 0.) {
	tolswp = eps4;
    }
    if (tolx <= 0.) {
	tolx = c6;
    }
    if (tolcon <= eps) {
	tolcon = c6;
    }
    toldj3 = toloptqp;
    if (scltol <= 0.) {
	scltol = .9;
    }
    if (scltol >= 1.) {
	scltol = .99;
    }
    if (tolpiv <= 0.) {
	tolpiv = eps1;
    }
    if (lincon) {
	if (wtinf0 < 0.) {
	    wtinf0 = 1.;
	}
    } else {
	if (wtinf0 < 0.) {
	    wtinf0 = 1e5;
	}
    }
    if (epsrf <= 0.) {
	epsrf = eps0;
    }
    if (fdint1 <= 0.) {
	fdint1 = sqrt(epsrf);
    }
    if (fdint2 <= 0.) {
	fdint2 = pow_dd(&epsrf, &c_b198);
    }
    if (iback == inewb) {
	iback = 0;
    }
    if (itnlim < 0) {
/* Computing MAX */
	i__1 = 10000, i__2 = max(*n,*m) * 10;
	itnlim = max(i__1,i__2);
    }
/*     Set default names (they may be printed by the basis routines). */
    if (s_cmp(mprob, "-1111111", (ftnlen)8, (ftnlen)8) == 0) {
	s_copy(mprob, "        ", (ftnlen)8, (ftnlen)8);
    }
    if (s_cmp(mobj, "-1111111", (ftnlen)8, (ftnlen)8) == 0) {
	s_copy(mobj, "        ", (ftnlen)8, (ftnlen)8);
    }
    if (s_cmp(mrhs, "-1111111", (ftnlen)8, (ftnlen)8) == 0) {
	s_copy(mrhs, "        ", (ftnlen)8, (ftnlen)8);
    }
    if (s_cmp(mrng, "-1111111", (ftnlen)8, (ftnlen)8) == 0) {
	s_copy(mrng, "        ", (ftnlen)8, (ftnlen)8);
    }
    if (s_cmp(mbnd, "-1111111", (ftnlen)8, (ftnlen)8) == 0) {
	s_copy(mbnd, "        ", (ftnlen)8, (ftnlen)8);
    }
/*     ------------------------------------------------------------------ */
/*     Done. */
/*     Re-assign the options to their respective work arrays. */
/*     ------------------------------------------------------------------ */
    rw[51] = toloptfp;
    rw[52] = toloptqp;
    rw[53] = toloptnp;
    rw[54] = tolcg;
    rw[56] = tolx;
    rw[57] = tolcon;
    rw[60] = tolpiv;
    rw[61] = tolrow;
    rw[62] = tcrash;
    rw[63] = utol1m;
    rw[64] = utol2m;
    rw[65] = tolswp;
    rw[66] = tolfac;
    rw[67] = tolupd;
    rw[70] = infbnd;
    rw[71] = bigfx;
    rw[72] = bigdx;
    rw[73] = epsrf;
    rw[76] = fdint1;
    rw[77] = fdint2;
    rw[79] = maxtime;
    rw[80] = xdlim;
    rw[81] = vilim;
    rw[83] = etarg;
    rw[84] = wolfeg;
    rw[85] = hcondbnd;
    rw[86] = zcndbd;
    rw[87] = condumax;
    rw[88] = wtinf0;
    rw[89] = xpen0;
    rw[90] = wtmax;
    rw[91] = proxweight;
    rw[92] = scltol;
    rw[151] = lmax1;
    rw[152] = lmax2;
    rw[153] = small;
    rw[154] = utol1;
    rw[155] = utol2;
    rw[156] = uspace;
    rw[157] = dens1;
    rw[158] = dens2;
/*     Dependent parameters set in s8defaults. */
    rw[181] = toldpp;
    rw[182] = toldcp;
    rw[183] = toldup;
    rw[186] = toldj3;
    rw[187] = toldrp;
/*     Addresses for integer quantities. */
    iw[52] = maxr;
    iw[53] = maxs;
    iw[54] = mqnmod;
    iw[55] = qpsolver;
    iw[56] = emode;
    iw[58] = kchk;
    iw[59] = kfac;
    iw[60] = ksav;
    iw[61] = klog;
    iw[62] = ksumm;
    iw[63] = kdegen;
    iw[64] = kreset;
    iw[66] = mflush;
    iw[67] = mskip;
/*     iw( 69)  =  lvlStart */
    iw[70] = lvlder;
    iw[71] = lvlsys;
    iw[72] = lvlhess;
    iw[73] = lvlobje;
    iw[75] = lvlscale;
    iw[76] = lvlsrch;
    iw[77] = lvlpre;
    iw[78] = lvlver;
    iw[79] = lvlppm;
    iw[80] = lvlpiv;
    iw[81] = lprprm;
    iw[82] = lprsch;
    iw[83] = lprscl;
    iw[84] = lprsol;
    iw[85] = lprdbg;
    iw[87] = minmax;
    iw[88] = icrash;
    iw[89] = itnlim;
    iw[90] = mmajor;
    iw[91] = mminor;
    iw[92] = mjrprint;
    iw[93] = mnrprint;
    iw[95] = mnewsb;
    iw[97] = cgitmx;
    iw[103] = objrow;
    iw[104] = deropt;
    iw[110] = jverf1;
    iw[111] = jverf2;
    iw[112] = jverf3;
    iw[113] = jverf4;
    iw[114] = jverf5;
    iw[115] = jverf6;
    iw[116] = stickyop;
    iw[120] = iback;
    iw[121] = idump;
    iw[122] = iloadb;
    iw[124] = inewb;
    iw[125] = iinsrt;
    iw[126] = ioldb;
    iw[127] = ipnch;
    iw[130] = ireprt;
    iw[131] = isoln;
    iw[151] = nout;
    iw[153] = maxcol;
    iw[156] = tpivot;
/*     Dependent parameters set in s8defaults. */
    iw[99] = nparprlp;
    iw[100] = nparprqp;
    iw[199] = minimize;
/*     Addresses for character quantities. */
    s_copy(cw + 408, mprob, (ftnlen)8, (ftnlen)8);
    s_copy(cw + 416, mobj, (ftnlen)8, (ftnlen)8);
    s_copy(cw + 424, mrhs, (ftnlen)8, (ftnlen)8);
    s_copy(cw + 432, mrng, (ftnlen)8, (ftnlen)8);
    s_copy(cw + 440, mbnd, (ftnlen)8, (ftnlen)8);
    return 0;
} /* s8defaults_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8defaults */
/* Subroutine */ int s8map_(integer *m, integer *n, integer *negcon, integer *
	nkx, integer *nncon, integer *nnjac, integer *nnobju, integer *nnobj, 
	integer *nnh, integer *lenr, integer *maxr, integer *maxs, integer *
	mqnmod, integer *lvlhess, integer *nextcw, integer *nextiw, integer *
	nextrw, integer *iw, integer *leniw)
{
    static integer lxscaled, lfeastype, nb, lr, ls, lv, ly, lr1, lr2, ls1, 
	    ls2, ls3, lu0, lx0, lx1, ly1, ly2, ly3, ldg, lhd, mbs, lfr, lrg, 
	    ldx, lfv, lfx, liy, lkx, lux, lrg2, liy1, liy2, lgbs, lkbs, lhdx, 
	    lpbs, lgqp, ngqp, lxbs, lxqp, lxqp0, lgobj, lblbs, llocg, lfcon, 
	    lgcon, nlocg, lbubs, lenfr, lblqp, lbuqp, lycon, lxpen, lgobj1, 
	    lgobj2, lfcon1, lfcon2, lgcon1, lgcon2, lycon1, lycon2, lgsave, 
	    lgobju, lgconu, ldycon, letype, lqprhs, lscales, lblsave, lbusave,
	     lestate;

/*     ================================================================== */
/*     s8Map   allocates all array storage for snopt, */
/*     using the values: */
/*        m    , n    , neJ */
/*        maxS                          Set in s8defaults. */
/*        nnObj, nnObjU, nnCon, nnJac   Set in specs file or arguments. */
/*        lenR , negCon                 Set in calling program */

/*     On exit, */
/*        nextcw, nextiw, nextrw are pointers to the next elements of */
/*                               free space in cw, iw, and rw. */

/*     29 Dec 2000: First version of s8Map. */
/*     24 Jan 2003: Added workspace for SYMMLQ */
/*     18 Jun 2008: Added space for iy2, pBS and rg2. */
/*     11 Sep 2014: Added nnObjU and nnObj for FP mode. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
/*     All dimensions are computed from */
/*        m     , n    , ne */
/*        lenR  , maxS , mQMmod */
/*        nnObjU, nnCon, nnJac */
/*        negCon */
    /* Parameter adjustments */
    --iw;

    /* Function Body */
    ngqp = *nnh;
    mbs = *m + *maxs;
    nb = *n + *m;
/*     Nonlinear constraints. */
    nlocg = *nnjac + 1;
/*     Addresses for the integer arrays. */
    lkx = *nextiw;
    lfeastype = lkx + *nkx;
    lkbs = lfeastype + mbs;
    lestate = lkbs + mbs;
    letype = lestate + nb;
    liy = letype + nb;
    liy1 = liy + nb;
    liy2 = liy1 + nb;
    *nextiw = liy2 + nb;
/*     Addresses for the double precision arrays. */
    lscales = *nextrw;
    ly = lscales + nb;
    ly1 = ly + nb;
    ly2 = ly1 + nb;
    ly3 = ly2 + nb;
/* SYMMLQ workspace */
    ls1 = ly3 + nb;
/* SYMMLQ workspace */
    ls2 = ls1 + *maxs;
/* SYMMLQ workspace */
    ls3 = ls2 + *maxs;
/* SYMMLQ workspace */
    lr1 = ls3 + *maxs;
/* SYMMLQ workspace */
    lr2 = lr1 + *maxs;
/* SYMMLQ workspace */
    lblqp = lr2 + *maxs;
    lbuqp = lblqp + nb;
    lblbs = lbuqp + nb;
    lbubs = lblbs + mbs;
    lxbs = lbubs + mbs;
    lxscaled = lxbs + mbs;
    lgbs = lxscaled + *nnh;
    lgqp = lgbs + mbs;
    lux = lgqp + ngqp;
    lhdx = lux + *nnh;
    lpbs = lhdx + *nnh;
    lu0 = lpbs + nb;
    ldg = lu0 + *nnh;
    lr = ldg + *nnh;
    lrg = lr + *lenr;
    lrg2 = lrg + *maxs;
    lblsave = lrg2 + *maxs;
    lbusave = lblsave + nb;
    *nextrw = lbusave + nb;
/*     Nonlinear Objective. */
    lgobj = *nextrw;
    lgobj1 = lgobj + *nnobj;
    lgobj2 = lgobj1 + *nnobj;
    lgobju = lgobj2 + *nnobj;
    lgsave = lgobju + *nnobj;
    *nextrw = lgsave + *nnobju;
/*     Nonlinear constraints. */
    llocg = *nextiw;
    *nextiw = llocg + nlocg;
    lfcon = *nextrw;
    lfcon1 = lfcon + *nncon;
    lfcon2 = lfcon1 + *nncon;
    lfx = lfcon2 + *nncon;
    lfv = lfx + *nncon;
    lycon = lfv + *nncon;
    lycon1 = lycon + *nncon;
    lycon2 = lycon1 + *nncon;
    ldycon = lycon2 + *nncon;
    lxpen = ldycon + *nncon;
    lgcon = lxpen + *nncon;
    lgcon1 = lgcon + *negcon;
    lgcon2 = lgcon1 + *negcon;
    lgconu = lgcon2 + *negcon;
    lqprhs = lgconu + *negcon;
    ldx = lqprhs + *m;
    lxqp = ldx + nb;
    lxqp0 = lxqp + nb;
    lx0 = lxqp0 + nb;
    lx1 = lx0 + nb;
    *nextrw = lx1 + nb;
/*     Store the addresses in iw. */
    iw[251] = lkx;
    iw[260] = llocg;
    iw[271] = lblqp;
    iw[272] = lbuqp;
    iw[273] = lblbs;
    iw[274] = lbubs;
    iw[275] = lblsave;
    iw[276] = lbusave;
    iw[277] = lpbs;
    iw[278] = lqprhs;
    iw[283] = letype;
    iw[284] = lfeastype;
    iw[285] = lestate;
    iw[287] = ldx;
    iw[288] = lhdx;
    iw[289] = ldg;
    iw[290] = lgqp;
    iw[291] = lgbs;
    iw[292] = lkbs;
    iw[293] = lrg;
    iw[294] = lrg2;
    iw[295] = lr;
    iw[296] = lscales;
    iw[297] = lgobj;
    iw[298] = lx0;
    iw[300] = lx1;
    iw[301] = lxbs;
    iw[302] = lxscaled;
    iw[304] = lxpen;
    iw[305] = lxqp;
    iw[306] = lxqp0;
    iw[308] = liy;
    iw[309] = liy1;
    iw[310] = liy2;
    iw[311] = ly;
    iw[312] = ly1;
    iw[313] = ly2;
    iw[314] = ly3;
    iw[316] = lfcon;
    iw[317] = lfcon1;
    iw[318] = lfcon2;
    iw[319] = lgconu;
    iw[320] = lgcon;
    iw[321] = lgcon1;
    iw[322] = lgcon2;
    iw[323] = lgobju;
    iw[324] = lgobj1;
    iw[325] = lgobj2;
    iw[336] = lfx;
    iw[337] = lfv;
    iw[339] = lgsave;
    iw[345] = lux;
    iw[346] = lu0;
    iw[348] = lycon;
    iw[349] = lycon1;
    iw[350] = lycon2;
    iw[351] = ldycon;
    iw[353] = lr1;
    iw[354] = lr2;
    iw[355] = ls1;
    iw[356] = ls2;
    iw[357] = ls3;
/*     Allocate space for an approximate Hessian. */
/*     The amount will depend on the method selected. */
    if (*lvlhess == 0) {
/*        --------------------------------------------------------------- */
/*        Compute the addresses of the limited-memory arrays. */
/*        These are saved and used for subsequent entries. */
/*        --------------------------------------------------------------- */
	lhd = *nextrw;
	ls = lhd + *nnh;
	lv = ls + *nnh * *mqnmod;
	*nextrw = lv + *nnh * *mqnmod;
	iw[347] = lhd;
	iw[401] = ls;
	iw[402] = lv;
    } else if (*lvlhess == 1) {
	lenfr = *nnh * (*nnh + 1) / 2;
	lhd = *nextrw;
	lfr = lhd + *nnh;
	*nextrw = lfr + lenfr;
	iw[347] = lhd;
	iw[391] = lfr;
	iw[392] = lenfr;
    }
    return 0;
} /* s8map_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8Map */
/* Subroutine */ int s8firstcall_(integer *iexit, logical *feasiblelc, 
	logical *gotfuns, logical *nonlinearcon, logical *nonlinearobj, 
	integer *n, integer *nb, integer *nncon0, integer *nncon, integer *
	nnjac, integer *nnh, integer *nnobj0, integer *nnobj, S_fp funwrapper,
	 U_fp funcon, U_fp funobj, U_fp userhv, doublereal *bl, doublereal *
	bu, doublereal *x, doublereal *x1, doublereal *ycon, integer *nej, 
	integer *nlocj, integer *locj, integer *indj, integer *neh, integer *
	nloch, integer *loch, integer *indh, doublereal *hcol, integer *
	negcon, integer *nlocg, integer *locg, doublereal *fobj, doublereal *
	fcon, doublereal *gcon, doublereal *gobj, doublereal *fcon2, 
	doublereal *gcon2, doublereal *gobj2, doublereal *y, doublereal *y1, 
	char *cu, integer *lencu, integer *iu, integer *leniu, doublereal *ru,
	 integer *lenru, char *cw, integer *lencw, integer *iw, integer *
	leniw, doublereal *rw, integer *lenrw, ftnlen cu_len, ftnlen cw_len)
{
    static integer savedlvlscl;
    extern /* Subroutine */ int s6getmissing_(integer *, integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     S_fp, U_fp, U_fp, U_fp, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, char *, integer *, integer *, integer *, doublereal 
	    *, integer *, char *, integer *, integer *, integer *, doublereal 
	    *, integer *, ftnlen, ftnlen);
    static integer modefg, inform__;
    extern /* Subroutine */ int s7checkg_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, S_fp, U_fp,
	     U_fp, U_fp, doublereal *, doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     doublereal *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, char *, integer *, integer *, integer 
	    *, doublereal *, integer *, char *, integer *, integer *, integer 
	    *, doublereal *, integer *, ftnlen, ftnlen);

/*     ================================================================== */
/*     s8firstCall  computes the first set of problem functions */

/*     Optionally, the derivatives are verified using finite differences. */

/*     17 Oct 2014: First version based on dnopt routine dnFirstCall. */
/*     19 Oct 2014: Added user-defined Hessian. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
/*     Scaling is turned off for now. */
/* scale option */
    /* Parameter adjustments */
    --x1;
    --x;
    --y1;
    --y;
    --bu;
    --bl;
    --fcon2;
    --fcon;
    --ycon;
    --gobj2;
    --gobj;
    --indj;
    --locj;
    --hcol;
    --indh;
    --loch;
    --gcon2;
    --gcon;
    --locg;
    cu -= 8;
    --iu;
    --ru;
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    savedlvlscl = iw[75];
    iw[75] = 0;
    modefg = 2;
    (*funwrapper)(&inform__, &modefg, nonlinearcon, nonlinearobj, n, negcon, 
	    nncon0, nncon, nnjac, nnh, nnobj0, nnobj, (U_fp)funcon, (U_fp)
	    funobj, (U_fp)userhv, &x[1], &ycon[1], nej, nlocj, &locj[1], &
	    indj[1], neh, nloch, &loch[1], &indh[1], &hcol[1], &fcon[1], fobj,
	     &gcon[1], &gobj[1], cu + 8, lencu, &iu[1], leniu, &ru[1], lenru, 
	    cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)8, (ftnlen)8)
	    ;
    *gotfuns = inform__ == 0;
    if (! (*gotfuns)) {
	if (inform__ < 0) {
	    if (*feasiblelc) {
		*iexit = 61;
/* Undefined fun at first feasible point */
	    } else {
		*iexit = 62;
/* Undefined fun at initial point */
	    }
	} else {
	    *iexit = inform__;
/* User wants to stop */
	}
	goto L999;
    }
/*     ------------------------------------------------------------------ */
/*     Check derivatives. */
/*     (One day, we will do this on the SCALED problem.) */

/*     Individual objective derivatives are not checked in FP mode. */
/*     ------------------------------------------------------------------ */
    s7checkg_(&inform__, n, nncon0, nncon, nnjac, nnh, nnobj0, nnobj, (S_fp)
	    funwrapper, (U_fp)funcon, (U_fp)funobj, (U_fp)userhv, &x[1], &x1[
	    1], &bl[1], &bu[1], fobj, &gobj[1], &ycon[1], nej, nlocj, &locj[1]
	    , &indj[1], neh, nloch, &loch[1], &indh[1], &hcol[1], negcon, 
	    nlocg, &locg[1], &fcon[1], &gcon[1], &gobj2[1], &fcon2[1], &gcon2[
	    1], &y[1], &y1[1], cu + 8, lencu, &iu[1], leniu, &ru[1], lenru, 
	    cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)8, (ftnlen)8)
	    ;
    if (inform__ != 0) {
	*iexit = inform__;
	goto L999;
    }
/*     ------------------------------------------------------------------ */
/*     Compute any missing derivatives. */
/*     ------------------------------------------------------------------ */
    s6getmissing_(&inform__, n, negcon, nncon0, nncon, nnjac, nnh, nnobj0, 
	    nnobj, (S_fp)funwrapper, (U_fp)funcon, (U_fp)funobj, (U_fp)userhv,
	     &bl[1], &bu[1], &x[1], &ycon[1], nej, nlocj, &locj[1], &indj[1], 
	    neh, nloch, &loch[1], &indh[1], &hcol[1], &fcon[1], fobj, &gcon[1]
	    , &gobj[1], &y1[1], cu + 8, lencu, &iu[1], leniu, &ru[1], lenru, 
	    cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)8, (ftnlen)8)
	    ;
    if (inform__ != 0) {
	*iexit = inform__;
	goto L999;
    }
L999:
    iw[75] = savedlvlscl;
    return 0;
} /* s8firstcall_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8firstCall */
/* Subroutine */ int s8sqp_(integer *iexit, S_fp funwrapper, U_fp funcon, 
	U_fp funobj, U_fp userhv, S_fp mjrlog, U_fp mnrlog, S_fp snstop, 
	logical *elastic, logical *gotr, integer *starttype, integer *itn, 
	integer *lenr, integer *m, integer *maxs, integer *mbs, integer *n, 
	integer *nb, integer *ns, integer *nncon0, integer *nncon, integer *
	nnobj0, integer *nnobj, integer *nnh0, integer *nnh, integer *majors, 
	integer *minors, integer *ndegen, doublereal *dualinf, integer *
	minimize, integer *iobj, doublereal *scaleobj, doublereal *objadd, 
	doublereal *fobj, doublereal *fmerit, doublereal *maxvi, doublereal *
	maxvirel, doublereal *supvi, integer *ninf, doublereal *sinf, integer 
	*ninfe, doublereal *sinfe, doublereal *wtinf0, doublereal *wtinf, 
	doublereal *penparm, doublereal *pinorm, doublereal *xnorm, integer *
	nej, integer *nlocj, integer *locj, integer *indj, doublereal *jcol, 
	integer *neh, integer *nloch, integer *loch, integer *indh, 
	doublereal *hcol, integer *negcon, integer *nlocg, integer *locg, 
	integer *etype, integer *estate, integer *feastype, integer *hs, 
	integer *kbs, doublereal *bl, doublereal *bu, doublereal *blqp, 
	doublereal *buqp, doublereal *blbs, doublereal *bubs, doublereal *fv, 
	doublereal *fx, doublereal *fcon, doublereal *gcon, doublereal *gobj, 
	doublereal *fcon1, doublereal *gcon1, doublereal *gobj1, doublereal *
	fcon2, doublereal *gcon2, doublereal *gobj2, doublereal *gbs, 
	doublereal *gqp, doublereal *dycon, doublereal *dx, doublereal *dg, 
	doublereal *udx, doublereal *hd, doublereal *hdx, doublereal *pbs, 
	doublereal *ycon, doublereal *ycon1, doublereal *ycon2, doublereal *
	pi, doublereal *qprhs, doublereal *r__, doublereal *rc, doublereal *
	rg, doublereal *rg2, doublereal *scales, doublereal *x, doublereal *
	x1, doublereal *xbs, doublereal *xqp0, doublereal *xqp, doublereal *
	xpen, integer *iy, integer *iy1, doublereal *y, doublereal *y1, 
	doublereal *y2, char *cu, integer *lencu, integer *iu, integer *leniu,
	 doublereal *ru, integer *lenru, char *cw, integer *lencw, integer *
	iw, integer *leniw, doublereal *rw, integer *lenrw, ftnlen cu_len, 
	ftnlen cw_len)
{
    /* Initialized data */

    static char line[4] = "----";

    /* Format strings */
    static char fmt_1000[] = "(1x,29a4)";
    static char fmt_1010[] = "(\002 Start of major itn\002,i6)";
    static char fmt_3020[] = "(\002 Itn\002,i7,\002 -- Central differences i"
	    "nvoked.\002,\002 Small step length.\002)";

    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    extern /* Subroutine */ int s8reseth_(integer *, doublereal *, doublereal 
	    *, integer *, integer *, doublereal *, integer *);
    static logical feasible, dualfeas;
    static integer jdualinf, elastics;
    static doublereal wtfactor;
    static logical printlog;
    static integer mjrprint, restarts;
    static doublereal toloptfp;
    static integer mnrprint, outeritn, qpsolver;
    extern /* Subroutine */ int s8optimizeslacks_(logical *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    static doublereal toloptnp;
    static logical printsum;
    static doublereal toloptqp;
    static integer j;
    extern /* Subroutine */ int s8initpen_(integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *), s8solveqp_(
	    integer *, U_fp, U_fp, U_fp, integer *, logical *, logical *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *, doublereal *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, ftnlen, ftnlen);
    static doublereal primalinf;
    static integer mrestarts, lurequest;
    static doublereal toloptqpk;
    extern /* Subroutine */ int s8fdswitch_(integer *, integer *, integer *, 
	    integer *, integer *, logical *, logical *, logical *, doublereal 
	    *, doublereal *, doublereal *, integer *, integer *, doublereal *,
	     integer *);
    extern /* Subroutine */ int s8hwrapper_();
    static logical primalfeas, needderivs;
    static doublereal ud0;
    static integer jprimalinf;
    static char str[132];
    static doublereal eps0, eps1, eps5;
    extern /* Subroutine */ int s8rc_(doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    extern /* Subroutine */ int s8hx_();
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static logical done, newb;
    static integer jobj, klog, maxr, itqp;
    static doublereal step;
    extern /* Subroutine */ int s6getmissing_(integer *, integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     S_fp, U_fp, U_fp, U_fp, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, char *, integer *, integer *, integer *, doublereal 
	    *, integer *, char *, integer *, integer *, integer *, doublereal 
	    *, integer *, ftnlen, ftnlen), s8getweights_(integer *, logical *,
	     integer *, doublereal *, doublereal *, doublereal *, doublereal *
	    , doublereal *, doublereal *, doublereal *, doublereal *, integer 
	    *, integer *), s6linesearch_(integer *, S_fp, U_fp, U_fp, U_fp, 
	    logical *, logical *, logical *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, integer *, integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, doublereal *, integer *, integer 
	    *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, char *, integer *, integer *, integer 
	    *, doublereal *, integer *, char *, integer *, integer *, integer 
	    *, doublereal *, integer *, ftnlen, ftnlen);
    static doublereal tolx, fobj1, fobj2;
    static logical nonlinearobj, nonlinearcon;
    static doublereal xpen0;
    extern /* Subroutine */ int s8hqn_(integer *, S_fp, U_fp, U_fp, U_fp, 
	    logical *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, doublereal *, integer *, doublereal *, integer *, logical *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, integer *, integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, doublereal *, integer *, integer 
	    *, integer *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, char *, integer *, integer *, integer *, doublereal 
	    *, integer *, char *, integer *, integer *, integer *, doublereal 
	    *, integer *, ftnlen, ftnlen), s8xhx_(integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, integer *), dload_(integer *, doublereal *, 
	    doublereal *, integer *);
    static logical fdobj;
    extern /* Subroutine */ int iload_(integer *, integer *, integer *, 
	    integer *);
    static logical fdcon;
    static integer nnjac;
    extern /* Subroutine */ int ddscl_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical goodg, usefd;
    static doublereal dxhdx;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static doublereal gnorm;
    static logical maxns, newlu;
    static integer isumm, ksumm, nswap;
    static doublereal wtmax;
    extern /* Subroutine */ int s2bfac_(integer *, integer *, logical *, 
	    logical *, logical *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, integer *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, integer *), daxpy_(integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *), 
	    s1page_(integer *, integer *, integer *);
    static doublereal sinfe1;
    extern /* Subroutine */ int s6rcnd_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    static doublereal gnorm0;
    extern doublereal dnrm1s_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int s8infs_(logical *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *), s2tols_(integer *, 
	    logical *, integer *, integer *, integer *, doublereal *, integer 
	    *);
    static doublereal utol1s, utol2s;
    static integer modefg;
    static logical needlu;
    static integer iabort;
    static logical ktcond[2];
    static integer cditns, mmajor;
    extern doublereal dnormi_(integer *, doublereal *, integer *);
    static integer inform__, lvlder, minmax, ninfqp;
    static logical maxits;
    static integer iprint;
    static logical fponly, feasibleslacks, newtol;
    static integer lvlpre, lvlpiv, typelu;
    static doublereal fobjqp, gmerit, drzmax, drzmin, phpmrt, sinfqp, tolcon, 
	    utolmn, weight, dxnorm;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen), s2aprod_(integer *, doublereal *, integer *, integer *
	    , integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *), s8inith_(
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *), 
	    s8merit_(integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *), s8gcopy_(integer *, integer *, integer *
	    , integer *, integer *, integer *, integer *, integer *, integer *
	    , doublereal *, integer *, integer *, integer *, doublereal *);
    static doublereal signobj;
    static logical boosted;
    static integer ninfeqp;
    static doublereal wtscale;
    static logical optimal;
    static doublereal sinfeqp;
    static logical nearopt;
    static integer lvlsrch;
    static doublereal condzhz;
    static integer lvlhess;
    static logical restart, fonlyls, firstqp;
    static integer qnskips, rtrmods;
    static logical gotnewx;

    /* Fortran I/O blocks */
    static icilist io___420 = { 0, str, 0, fmt_1000, 132, 1 };
    static icilist io___422 = { 0, str, 0, fmt_1010, 132, 1 };
    static icilist io___423 = { 0, str, 0, fmt_1000, 132, 1 };
    static icilist io___424 = { 0, str, 0, fmt_1010, 132, 1 };
    static icilist io___470 = { 0, str, 0, fmt_3020, 132, 1 };
    static icilist io___472 = { 0, str, 0, fmt_1000, 132, 1 };
    static icilist io___473 = { 0, str, 0, fmt_1010, 132, 1 };
    static icilist io___474 = { 0, str, 0, fmt_1000, 132, 1 };
    static icilist io___475 = { 0, str, 0, fmt_1010, 132, 1 };


/*     ================================================================== */
/*     s8SQP  solves a nonlinear programming problem. */
/*     A basis is assumed to be specified by nS, hs, x and the */
/*     superbasic parts of kBS. */
/*     In particular, there must be nS values hs(j) = 2, and the */
/*     corresponding j's must be listed in kBS(m+1) thru kBS(m+ns). */
/*     The ordering in kBS(m+1:m+nS) matches the reduced Hessian R. */

/*     On entry, if there are nonlinear constraints, Fx contains */
/*     the true nonlinear slacks (i.e., constraint values) */
/*     Fx  =  fCon + (linear A)*x,   excluding slacks. */

/*     On exit, if  iExit .lt. 30  it is safe to save the final */
/*     basis files and print the solution.  Otherwise, a fatal error */
/*     condition exists and numerous items will be undefined. */
/*     The last basis map saved (if any) retains the only useful */
/*     information. */

/*     30 Dec 1991: First version based on npsol routine npcore. */
/*     23 Oct 1993: Proximal point FP added. */
/*     29 Oct 1993: Crash on LG rows moved outside s5QP. */
/*     24 Apr 1994: Nx columns no longer in Q. */
/*     26 May 1995: Column order of R defined by kBS. */
/*     04 Aug 1995: Limited memory update */
/*     11 Aug 1995: tolg changed from 0.1 to 1.0d-4. */
/*     09 Nov 1995: Updated multipliers used to define Lagrangian. */
/*     19 Dec 1995: Finite-differences added. */
/*     09 Oct 1996: First Min Sum version. */
/*     16 Jul 1997: First thread-safe version. */
/*     09 Jul 1998: Quasi-Newton updates implemented correctly. */
/*     24 Aug 1998: Fixed bug in s8x1 found by Alan Brown at Nag. */
/*     06 Sep 1998: Pre- and post-QP diagonal Hessian scaling added. */
/*     16 Jan 1999: Name changed from s8core. */
/*     06 Apr 2001: For Hot Starts, don't mess with Utol. */
/*     27 Apr 2001: wtMax introduced as parameter to s8getWeights. */
/*     15 Jan 2003: CG and QN  QP solvers added. */
/*     03 Aug 2003: snEXIT and snPRNT adopted. */
/*     04 Jul 2005: Switched to vanilla CG for QN with nS > maxR. */
/*     16 Jun 2008: Call-status implemented correctly. */
/*     18 Jun 2008: pBS added as argument. */
/*     07 Oct 2014: infoTags added. */
/*     19 Oct 2014: Added user-defined Hessian. */
/*     29 Dec 2014: s8iQP and s8iQN merged to form s8solveQP. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/* abs tol for small diag of U. */
/* rel tol for small diag of U. */
/* =1(2) forwd (cntrl) diffs */
/* number of Hx products */
/* Current QP solver */
/* Current precon mode */
/* >0 => Mnr heading for iPrint */
/* >0 => Mjr heading for iPrint */
/* >0 => Mjr heading for iSumm */
/* infoTag(1:7) status info array */
/* (1): QN update type */
/* (2): QN mod type (A or A+B) */
/* (3): Line search result */
/* (4): QP Feasibility status */
/* (5)  QP Optimality  status */
/* (6) */
/*     ------------------------------------------------------------------ */
/* (7): Approx Hess type */
    /* Parameter adjustments */
    --r__;
    --qprhs;
    --pi;
    --rg2;
    --rg;
    --xbs;
    --gbs;
    --bubs;
    --blbs;
    --kbs;
    --feastype;
    --y2;
    --y1;
    --y;
    --iy1;
    --iy;
    --xqp;
    --xqp0;
    --x1;
    --x;
    --scales;
    --rc;
    --pbs;
    --dx;
    --buqp;
    --blqp;
    --bu;
    --bl;
    --hs;
    --estate;
    --etype;
    --xpen;
    --ycon2;
    --ycon1;
    --ycon;
    --dycon;
    --fcon2;
    --fcon1;
    --fcon;
    --fx;
    --fv;
    --gobj2;
    --gobj1;
    --gobj;
    --hdx;
    --hd;
    --udx;
    --dg;
    --gqp;
    --penparm;
    --jcol;
    --indj;
    --locj;
    --hcol;
    --indh;
    --loch;
    --gcon2;
    --gcon1;
    --gcon;
    --locg;
    cu -= 8;
    --iu;
    --ru;
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
/*     ------------------------------------------------------------------ */
    iprint = iw[12];
/* Print file */
    isumm = iw[13];
/* Summary file */
    nnjac = iw[21];
/* # nonlinear Jacobian variables */
    maxr = iw[52];
/* max columns of R. */
    qpsolver = iw[55];
/* = 0:1:2   => QPChol:CG:QN QP solver */
    klog = iw[61];
/* log/print frequency */
    ksumm = iw[62];
/* Summary print frequency */
    lvlhess = iw[72];
/* 0,1,2  => LM, FM, Exact Hessian */
    lvlsrch = iw[76];
/* >0     => use derivatives in the line search */
    lvlpre = iw[77];
/* >0     => QN preconditioned CG */
    lvlpiv = iw[80];
/* 0/1 Threshold partial/complete pivoting: use */
    minmax = iw[87];
/* 1, 0, -1  => MIN, FP, MAX */
    mmajor = iw[90];
/* limit on major iterations */
    mjrprint = iw[92];
/* Major print level */
    mnrprint = iw[93];
/* Minor print level */
    lvlder = iw[70];
/*     Constants */
/* = 0, 1, 2 or 3, the derivative level */
    eps0 = rw[2];
/* eps**(4/5)          IEEE DP  3.00e-13 */
    eps1 = rw[3];
/* eps**(2/3)          IEEE DP  3.67e-11 */
    eps5 = rw[7];
/* eps**(1/5)          IEEE DP  7.40e-04 */
    toloptfp = rw[51];
/* Minor Phase 1 Opt tol */
    toloptqp = rw[52];
/* Minor Phase 2 Opt tol */
    toloptnp = rw[53];
/* Major Optimality tolerance */
    tolx = rw[56];
/* Minor feasibility tolerance. */
    tolcon = rw[57];
/* Major feasibility tolerance. */
    xpen0 = rw[89];
/* initial penalty parameter. */
    wtmax = rw[90];
/* max     infeasibility weight */
    nonlinearcon = *nncon > 0;
    nonlinearobj = *nnobj > 0;
    fponly = minmax == 0;
    iw[188] = 0;
    iw[223] = 0;
    iw[224] = 0;
    iw[226] = 0;
    iw[208] = qpsolver;
/* Current QP solver */
    iw[209] = lvlpre;
/*     ------------------------------------------------------------------ */
/*     s8SQP  operates in either ``Normal'' or ``Elastic'' mode. */
/*     In elastic mode, the nonlinear slacks are allowed to be infeasible */
/*     while a weighted sum of the slack infeasibilities is minimized. */
/*     ------------------------------------------------------------------ */
/* Current precon mode */
    feasibleslacks = TRUE_;
    *elastic = FALSE_;
    iload_(nb, &c__0, &estate[1], &c__1);
    elastics = 0;
    *ninfe = 0;
    *sinfe = 0.;
    sinfe1 = 0.;
    dcopy_(nb, &bl[1], &c__1, &blqp[1], &c__1);
    dcopy_(nb, &bu[1], &c__1, &buqp[1], &c__1);
    *ninf = 0;
    *sinf = 0.;
    *iexit = 0;
    lurequest = 0;
    qnskips = 0;
    restarts = 0;
    if (*nnh > 0) {
	mrestarts = 2;
    } else {
	mrestarts = 0;
    }
    restart = FALSE_;
    rtrmods = 0;
    iload_(&c__7, &c__0, &iw[237], &c__1);
    signobj = (doublereal) (*minimize);
    gnorm = 0.;
    if (*iobj == 0) {
	gnorm0 = 0.;
    } else {
	gnorm0 = *scaleobj;
    }
/* $$$      gNorm0  = one */
/* $$$      if (NonlinearObj) then */
/* $$$!        gNorm0 = dnormi( nnObj, gObj, 1 ) */
/* $$$         gNorm0 = dnrm2 ( nnObj, gObj, 1 ) */
/* $$$      end if */
    primalinf = 0.;
    *dualinf = 0.;
    *wtinf = *wtinf0;
    toloptqpk = toloptqp * 10.;
    gmerit = 0.;
    step = 0.;
    ktcond[0] = FALSE_;
    ktcond[1] = FALSE_;
    firstqp = TRUE_;
    condzhz = 1.;
    fdobj = (lvlder == 0 || lvlder == 2) && *nnobj > 0;
    fdcon = (lvlder == 0 || lvlder == 1) && nnjac > 0;
    usefd = fdobj || fdcon;
    fonlyls = usefd || lvlsrch == 0;
    if (mjrprint >= 10 || mnrprint >= 10) {
	printlog = iprint > 0 && klog == 1;
	printsum = isumm > 0 && ksumm == 1;
	if (printlog) {
	    s_wsfi(&io___420);
	    for (j = 1; j <= 29; ++j) {
		do_fio(&c__1, line, (ftnlen)4);
	    }
	    e_wsfi();
	    snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)132);
	    s_wsfi(&io___422);
	    do_fio(&c__1, (char *)&(*majors), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)132);
	}
	if (printsum && mnrprint >= 10) {
	    s_wsfi(&io___423);
	    for (j = 1; j <= 19; ++j) {
		do_fio(&c__1, line, (ftnlen)4);
	    }
	    e_wsfi();
	    snprnt_(&c__2, str, &iw[1], leniw, (ftnlen)132);
	    s_wsfi(&io___424);
	    do_fio(&c__1, (char *)&(*majors), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__2, str, &iw[1], leniw, (ftnlen)132);
	}
    }
    jobj = *n + *iobj;
    if (nonlinearcon) {
	s8initpen_(nncon, &penparm[1], &xpen0, &xpen[1], &rw[1], lenrw);
    }
    if (*ns > maxr) {
	iw[208] = 1;
/* Use CG */
	iw[209] = 0;
/* with no preconditioning */
	*gotr = FALSE_;
    }
    dcopy_(nb, &x[1], &c__1, &xqp[1], &c__1);
    cditns = -1;
    needderivs = FALSE_;
/*     ======================Start of main loop========================== */
/*     Start of a Major Iteration. */
/*     ================================================================== */
    i__1 = mmajor;
    for (outeritn = 0; outeritn <= i__1; ++outeritn) {
	gotnewx = FALSE_;
/*        =============================================================== */
/*        Repeat                                          (until GotNewx) */
/*        =============================================================== */
L100:
	*minors = 0;
/*           ============================================================ */
/*           Repeat                 (until an accurate gradient is found) */
L200:
	if (needderivs) {
	    if (usefd) {
/*                    --------------------------------------------------- */
/*                    Compute any missing derivatives. */
/*                    --------------------------------------------------- */
		s6getmissing_(iexit, n, negcon, nncon0, nncon, &nnjac, nnh, 
			nnobj0, nnobj, (S_fp)funwrapper, (U_fp)funcon, (U_fp)
			funobj, (U_fp)userhv, &bl[1], &bu[1], &x[1], &ycon[1],
			 nej, nlocj, &locj[1], &indj[1], neh, nloch, &loch[1],
			 &indh[1], &hcol[1], &fcon[1], fobj, &gcon[1], &gobj[
			1], &y[1], cu + 8, lencu, &iu[1], leniu, &ru[1], 
			lenru, cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (
			ftnlen)8, (ftnlen)8);
		if (*iexit != 0) {
		    goto L999;
/* break major iteration loop */
		}
	    }
/* UseFD */
	    needderivs = FALSE_;
	}
	if (nonlinearobj) {
/*                 gNorm = dnormi( nnObj, gObj, 1 ) */
	    gnorm = dnrm1s_(nnobj, &gobj[1], &c__1);
/* Approximate 2-norm */
	}
	if (nonlinearcon) {
/*                 ------------------------------------------------------ */
/*                 Load the scaled Jacobian in J. */
/*                 Compute the QP right-hand side   QPrhs  =  Jx - fCon. */
/*                 Find Fx the nonlinear constraint values. */
/*                 ------------------------------------------------------ */
	    s8gcopy_(nncon, &nnjac, nej, nlocj, &locj[1], &indj[1], negcon, 
		    nlocg, &locg[1], &gcon[1], nej, nlocj, &locj[1], &jcol[1])
		    ;
	    dcopy_(nncon, &fcon[1], &c__1, &qprhs[1], &c__1);
	    s2aprod_(&c__0, &eps0, nej, nlocj, &locj[1], &indj[1], &jcol[1], &
		    c_b245, &x[1], &nnjac, &c_b246, &qprhs[1], nncon);
/*                 ------------------------------------------------------ */
/*                 s8sOpt  finds the nonlinear slacks  sN  that minimize */
/*                 the merit function with x(1:n) and yCon held fixed. */
/*                 The optimal slacks are loaded into  x(n+1:nb)  and the */
/*                 violations are calculated: */
/*                  Fv = fCon  + A(linear)x - nonlinear slacks */
/*                     = Fx                 - sN */
/*                 ------------------------------------------------------ */
	    if (! (*elastic)) {
		s8getweights_(&c__0, &boosted, itn, &gnorm0, &gnorm, wtinf0, 
			wtinf, &wtmax, &weight, &wtfactor, &wtscale, &iw[1], 
			leniw);
	    }
	    s8optimizeslacks_(elastic, n, nb, nncon, ninfe, sinfe, &tolx, 
		    wtinf, &bl[1], &bu[1], &fv[1], &x[1], &ycon[1], &xpen[1], 
		    &fx[1]);
	}
/*              --------------------------------------------------------- */
/*              Prepare to (re-)solve the QP subproblem (possibly after */
/*              the elastic weight has been increased). */
/*              --------------------------------------------------------- */
/*              Factorize the basis at x. */
/*              Compute xQP such that (J -I)*xQP = rhs. */
L300:
	if (firstqp) {
/*                 ------------------------------------------------------ */
/*                 First QP subproblem. */
/*                 ------------------------------------------------------ */
	    s8inith_(nnh0, nnh, &gnorm0, &gnorm, &ud0, &hd[1], &iw[1], leniw, 
		    &rw[1], lenrw);
	    needlu = TRUE_;
	    *gotr = FALSE_;
	    nswap = 0;
	    if (*ns == 0) {
		typelu = 0;
	    } else {
		typelu = 2;
	    }
	    utol1s = rw[154];
	    utol2s = rw[155];
/*                 To avoid an unnecessarily ill-conditioned starting */
/*                 basis for the first QP, use big singularity tols */
/*                 (except if it's a Hot Start!). */
	    if (*starttype == 3) {
		utolmn = eps1;
	    } else {
		utolmn = eps5;
	    }
	    rw[154] = max(utol1s,utolmn);
	    rw[155] = max(utol2s,utolmn);
	} else {
/*                 ------------------------------------------------------ */
/*                 Subsequent factorizations. */
/*                 ------------------------------------------------------ */
/*                 For linearly constrained problems, the factors L, U */
/*                 and R can be saved as long as a poor x does not force */
/*                 a new factorization. (Even in this case, R can be */
/*                 saved if there are no swaps.) */
	    needlu = nonlinearcon;
	    if (restarts == 0) {
		typelu = 3;
/*                    Reset the factors and update tols if changed */
/*                    during the previous major iteration. */
		s2tols_(&c__3, &newtol, itn, &iw[1], leniw, &rw[1], lenrw);
	    }
	}
	s2bfac_(iexit, &typelu, &needlu, &newlu, &newb, iobj, itn, &mjrprint, 
		&lurequest, m, mbs, n, nb, nnh, ns, &nswap, nej, nlocj, &locj[
		1], &indj[1], &jcol[1], &kbs[1], &hs[1], &blqp[1], &buqp[1], &
		blbs[1], &bubs[1], nncon0, nncon, &qprhs[1], &xqp[1], &xbs[1],
		 &iy[1], &iy1[1], &y[1], &y1[1], &iw[1], leniw, &rw[1], lenrw)
		;
	if (*iexit != 0) {
	    goto L999;
/* break major iteration loop */
	}
	if (iw[208] == 0) {
	    *gotr = *gotr && ! newlu;
	}
	needlu = FALSE_;
	if (mjrprint >= 10) {
	    iw[224] = 1;
	}
/*              xQP satisfies the general constraints. */
	if (firstqp) {
	    iw[80] = lvlpiv;
/* Reset original TPP or TCP */
	    rw[154] = utol1s;
	    rw[155] = utol2s;
	}
/*              --------------------------------------------------------- */
/*              Solve the QP subproblem to obtain kBS, xQP and pi. */
/*              The search direction will be dx = xQP - x. */
/*              Use x1 to store the first feasible point. */
/*              --------------------------------------------------------- */
	inform__ = 0;
	d__1 = *objadd + *fobj;
	s8solveqp_(&inform__, (U_fp)mnrlog, (U_fp)s8hwrapper_, (U_fp)s8hx_, &
		iw[188], elastic, gotr, itn, &itqp, lenr, m, maxs, mbs, n, nb,
		 nncon0, nncon, nnobj0, nnobj, nnh0, nnh, ns, ndegen, &
		mjrprint, &mnrprint, minimize, iobj, scaleobj, &d__1, &fobjqp,
		 &condzhz, &toloptfp, &toloptqpk, &tolx, &ninfqp, &sinfqp, &
		elastics, &ninfeqp, &sinfeqp, wtinf, &ud0, pinorm, nej, nlocj,
		 &locj[1], &indj[1], &jcol[1], neh, nloch, &loch[1], &indh[1],
		 &hcol[1], &etype[1], &estate[1], &feastype[1], &hs[1], &kbs[
		1], &bl[1], &bu[1], &blqp[1], &buqp[1], &blbs[1], &bubs[1], &
		gbs[1], &gqp[1], &gobj[1], &hd[1], &hdx[1], &pbs[1], &pi[1], &
		r__[1], &rc[1], &rg[1], &rg2[1], &qprhs[1], &scales[1], &x[1],
		 &xbs[1], &xqp0[1], &xqp[1], &iy[1], &iy1[1], &y[1], &y1[1], &
		y2[1], cu + 8, lencu, &iu[1], leniu, &ru[1], lenru, cw + 8, 
		lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)8, (ftnlen)8);
/*              inform    Status */
/*              ------    ------ */
/*               >0       Fatal error */
/*                0       QP solution found */
/*               -1       Too many iterations */
/*               -2       Too many superbasics */
	*minors += itqp;
	if (inform__ > 0) {
	    *iexit = inform__;
	    goto L999;
/* break major iteration loop */
	}
/*              QP inform values are saved until after printing */
	maxns = inform__ == -2;
	maxits = inform__ == -1;
	firstqp = FALSE_;
/*              --------------------------------------------------------- */
/*              In normal mode, the routine s8optimizeSlacks defines */
/*              feasible nonlinear slacks x(n+1:n+nnCon). */
/*              If the QP has just switched to elastic mode, the slacks */
/*              are no longer inside their bounds and x(n+1:n+nnCon) */
/*              can be defined as the true infeasibilities. */
/*              --------------------------------------------------------- */
	if (nonlinearcon) {
/*                 Compute the dual search direction. */
/*                 Set        dyCon = yCon - pi(qp) */
	    dcopy_(nncon, &pi[1], &c__1, &dycon[1], &c__1);
	    daxpy_(nncon, &c_b246, &ycon[1], &c__1, &dycon[1], &c__1);
	    if (*elastic && feasibleslacks) {
/*                    Initialize the elastic slacks */
/*                    This happens only once. */
		feasibleslacks = FALSE_;
		dcopy_(nncon, &fx[1], &c__1, &x[*n + 1], &c__1);
		dcopy_(nncon, &pi[1], &c__1, &ycon[1], &c__1);
		dload_(nncon, &c_b4, &dycon[1], &c__1);
/*                    yCon and  x  have changed, recompute Fv. */
		s8optimizeslacks_(elastic, n, nb, nncon, ninfe, sinfe, &tolx, 
			wtinf, &bl[1], &bu[1], &fv[1], &x[1], &ycon[1], &xpen[
			1], &fx[1]);
	    }
	}
/*              Compute the search direction dx. */
	dcopy_(nb, &xqp[1], &c__1, &dx[1], &c__1);
	daxpy_(nb, &c_b246, &x[1], &c__1, &dx[1], &c__1);
	*xnorm = dnormi_(n, &x[1], &c__1);
	dxnorm = dnormi_(n, &dx[1], &c__1);
/*              --------------------------------------------------------- */
/*              Test for convergence. */
/*              --------------------------------------------------------- */
/*              Compute all the QP reduced costs. */
/*              (We could use yCon for the nonlinear pi's). */
/*              Compute the maximum dual infeasibility. */
	s8rc_(scaleobj, minimize, iobj, m, n, nb, nnobj0, nnobj, nncon, &
		nnjac, negcon, nej, nlocj, &locj[1], &indj[1], &jcol[1], &
		gobj[1], &gcon[1], &pi[1], &rc[1]);
	s8infs_(elastic, n, nb, nncon0, nncon, &tolx, wtinf, &primalinf, 
		dualinf, &jprimalinf, &jdualinf, &bl[1], &bu[1], &fx[1], &rc[
		1], &x[1]);
	if (*gotr) {
	    s6rcnd_(&maxr, ns, lenr, &r__[1], &drzmax, &drzmin, &condzhz);
	}
	primalinf /= *xnorm + 1.;
	*dualinf /= *pinorm;
	primalfeas = primalinf <= tolcon;
	dualfeas = *dualinf <= toloptnp;
	ktcond[0] = primalfeas;
	ktcond[1] = dualfeas;
	optimal = dualfeas;
	feasible = primalfeas;
	done = feasible && optimal || optimal && *ninfe > 0 || feasible && 
		fponly;
	nearopt = ! done && (feasible && *dualinf < toloptnp * 10. || optimal 
		&& primalinf < tolcon * 10.);
	if (*elastic && done) {
	    s8getweights_(&c__1, &boosted, itn, &gnorm0, &gnorm, wtinf0, 
		    wtinf, &wtmax, &weight, &wtfactor, &wtscale, &iw[1], 
		    leniw);
	    if (boosted) {
/* Solve the QP again */
		needderivs = FALSE_;
		goto L300;
	    }
	}
	if (*majors == 0) {
	    toloptqpk = min(*dualinf,.001);
	}
	if (*minors == 0 && toloptqpk > toloptqp) {
	    toloptqpk *= .2;
	}
/* Computing MIN */
	d__1 = toloptqpk * .5, d__2 = *dualinf * .1;
	toloptqpk = min(d__1,d__2);
	toloptqpk = max(toloptqpk,toloptqp);
/*              --------------------------------------------------------- */
/*              Compute the current augmented Lagrangian merit function. */
/*              objAdd is added in the log routine. */
/*              --------------------------------------------------------- */
	if (*iobj == 0) {
	    *fmerit = 0.;
	} else {
	    *fmerit = signobj * x[jobj] * *scaleobj;
	}
	if (nonlinearobj) {
	    *fmerit += signobj * *fobj;
	}
	if (nonlinearcon) {
	    dcopy_(nncon, &fv[1], &c__1, &y[1], &c__1);
	    ddscl_(nncon, &xpen[1], &c__1, &y[1], &c__1);
	    *fmerit = *fmerit - ddot_(nncon, &ycon[1], &c__1, &fv[1], &c__1) 
		    + ddot_(nncon, &y[1], &c__1, &fv[1], &c__1) * .5;
	    if (*elastic) {
		*fmerit += *wtinf * *sinfe;
	    }
	}
/*              --------------------------------------------------------- */
/*              If the forward-difference estimate of the reduced gradient */
/*              of the Lagrangian is small,  prepare to: (i) switch to */
/*              central differences; (ii)  recompute the derivatives,  and */
/*              (iii) solve the QP again. */

/*              On the other hand, if central differences give a large */
/*              reduced-gradient norm, switch back to forward differences. */
/*              --------------------------------------------------------- */
	s8fdswitch_(nncon0, nncon, nnobj, itn, &cditns, &goodg, &needderivs, &
		usefd, dualinf, &fcon[1], fobj, &iw[1], leniw, &rw[1], lenrw);
/*              --------------------------------------------------------- */
/*              Print the details of this iteration. */
/*              --------------------------------------------------------- */
	(*mjrlog)(&iabort, ktcond, &mjrprint, minimize, n, nb, nncon0, nnobj, 
		ns, itn, majors, minors, &nswap, &condzhz, iobj, scaleobj, 
		objadd, fmerit, &penparm[1], &step, &primalinf, dualinf, 
		maxvi, maxvirel, &hs[1], nej, nlocj, &locj[1], &indj[1], &
		jcol[1], &scales[1], &bl[1], &bu[1], &fcon[1], &ycon[1], &x[1]
		, cu + 8, lencu, &iu[1], leniu, &ru[1], lenru, cw + 8, lencw, 
		&iw[1], leniw, &rw[1], lenrw, (ftnlen)8, (ftnlen)8);
	(*snstop)(&iabort, ktcond, &mjrprint, minimize, m, maxs, n, nb, 
		nncon0, nncon, nnobj0, nnobj, ns, itn, majors, minors, &nswap,
		 &condzhz, iobj, scaleobj, objadd, fmerit, &penparm[1], &step,
		 &primalinf, dualinf, maxvi, maxvirel, &hs[1], nej, nlocj, &
		locj[1], &indj[1], &jcol[1], negcon, &scales[1], &bl[1], &bu[
		1], &fcon[1], &gcon[1], &gobj[1], &ycon[1], &pi[1], &rc[1], &
		rg[1], &x[1], cu + 8, lencu, &iu[1], leniu, &ru[1], lenru, cw 
		+ 8, lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)8, (ftnlen)
		8);
	if (iabort != 0) {
	    *iexit = 74;
/* User has aborted the run via snSTOP */
	    goto L999;
/* break major iteration loop */
	}
	iw[239] = 0;
/* +          until (.not. (UseFD  .and.  .not.GoodG)) */
	if (usefd && ! goodg) {
	    goto L200;
	}
/*           ============================================================ */
	if (maxits) {
	    *iexit = 31;
	}
	if (*majors >= mmajor) {
	    *iexit = 32;
	}
	if (maxns && *majors >= mmajor) {
	    *iexit = 33;
	}
	if (feasible && iw[241] == 3) {
	    *iexit = 21;
	}
	if (nearopt && *iexit != 0 && (! usefd || iw[181] == 1)) {
	    *iexit = 3;
	}
	if (feasible && optimal) {
	    *iexit = 1;
	}
	if (feasible && fponly) {
	    *iexit = 2;
	}
	if (optimal && *ninfe > 0) {
	    *iexit = 13;
	}
	if (*iexit != 0) {
	    goto L999;
/* break major iteration loop */
	}
	step = 0.;
	nswap = 0;
/*           ============================================================ */
/*           Take a step in the right direction. */
/*           ============================================================ */
/*           Compute  dxHdx = s'Hs  and other directional derivatives. */
/*           Be prepared to fix up pHpMrt if there are linear variables. */
/*           ------------------------------------------------------------ */
	if (*nnh > 0) {
	    s8xhx_(nnh, &dx[1], &udx[1], &hdx[1], &dxhdx, &iw[1], leniw, &rw[
		    1], lenrw);
	} else {
	    dxhdx = 0.;
	}
	if (*nnh == *n) {
	    if (dxhdx == 0. && iw[243] != 2) {
		s8reseth_(nnh, &ud0, &hd[1], &iw[1], leniw, &rw[1], lenrw);
		needderivs = FALSE_;
		goto L100;
	    }
	    phpmrt = dxhdx;
	} else {
/* Computing MAX */
	    d__1 = eps1 * dxnorm * dxnorm;
	    phpmrt = max(d__1,dxhdx);
	}
/*           ------------------------------------------------------------ */
/*           Compute the contributions to the merit function and its */
/*           directional derivative from the nonlinear constraints. */
/*           The penalty parameters  xPen(j)  are increased if the */
/*           directional derivative is not sufficiently negative. */
/*           ------------------------------------------------------------ */
/*           First, compute the value and directional derivative of the */
/*           Lagrangian with respect to x and the multipliers. */
	if (*iobj == 0) {
	    *fmerit = 0.;
	    gmerit = 0.;
	} else {
	    *fmerit = signobj * x[jobj] * *scaleobj;
	    gmerit = signobj * dx[jobj] * *scaleobj;
	}
	if (nonlinearobj) {
	    *fmerit += signobj * *fobj;
	    gmerit += signobj * ddot_(nnobj, &gobj[1], &c__1, &dx[1], &c__1);
	}
	if (*elastic) {
	    *fmerit += *sinfe * *wtinf;
	    gmerit += (sinfeqp - *sinfe) * *wtinf;
	}
/*           ------------------------------------------------------------ */
/*           Compute the search direction for the multipliers and nonlinear */
/*           slacks, and the contributions to the merit function and its */
/*           directional derivative from the nonlinear constraints. */
/*           The penalty parameters  xPen(j)  are increased if the */
/*           directional derivative is not sufficiently negative. */
/*           ------------------------------------------------------------ */
	if (nonlinearcon) {
	    *fmerit -= ddot_(nncon, &ycon[1], &c__1, &fv[1], &c__1);
	    gmerit += ddot_(nncon, &ycon[1], &c__1, &fv[1], &c__1);
	    gmerit -= ddot_(nncon, &dycon[1], &c__1, &fv[1], &c__1);
	    s8merit_(nncon, fmerit, &gmerit, &phpmrt, &penparm[1], &fv[1], &
		    xpen[1], &y[1], &rw[1], lenrw);
	}
/*           ============================================================ */
/*           Compute the scalar step, the step length from x along dx. */
/*           Implicitly,  dx defines the change sInfEQP - sInfE in the */
/*           elastic variables. */

/*           x   is the base point, with associated  sInfE */
/*           x1  is x + step*dx,    with associated  sInfE1 */

/*           pBS == x2 is workspace. */
/*           ============================================================ */
	s6linesearch_(&inform__, (S_fp)funwrapper, (U_fp)funcon, (U_fp)funobj,
		 (U_fp)userhv, elastic, &fonlyls, &primalfeas, iobj, &signobj,
		 scaleobj, n, nb, nncon0, nncon, &nnjac, nnh, nnobj0, nnobj, 
		itn, majors, &qnskips, maxvi, maxvirel, supvi, &step, &dxnorm,
		 xnorm, fmerit, &gmerit, sinfe, &sinfeqp, &sinfe1, wtinf, &bl[
		1], &bu[1], &dx[1], &dycon[1], nej, nlocj, &locj[1], &indj[1],
		 &jcol[1], neh, nloch, &loch[1], &indh[1], &hcol[1], negcon, 
		nlocg, &locg[1], &fobj1, &fcon1[1], &gcon1[1], &gobj1[1], &
		fobj2, &fcon2[1], &gcon2[1], &gobj2[1], &fx[1], &x[1], &x1[1],
		 &pbs[1], &xqp[1], &pi[1], &ycon[1], &ycon1[1], &ycon2[1], &
		xpen[1], &y[1], &y1[1], &y2[1], cu + 8, lencu, &iu[1], leniu, 
		&ru[1], lenru, cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (
		ftnlen)8, (ftnlen)8);
	restart = restarts < mrestarts;
/* permission to restart */
	if (inform__ == 0) {
	    gotnewx = TRUE_;
	} else if (inform__ > 0) {
	    *iexit = inform__;
/* The user wants to stop or violn limit */
	    goto L999;
/* break major iteration loop */
	} else {
/*              ========================================================= */
/*              No acceptable step */
/*              ========================================================= */
	    if (nonlinearcon && *sinfe > 0.) {
		s8getweights_(&c__1, &boosted, itn, &gnorm0, &gnorm, wtinf0, 
			wtinf, &wtmax, &weight, &wtfactor, &wtscale, &iw[1], 
			leniw);
		if (boosted) {
		    *elastic = TRUE_;
		}
	    }
/*              ========================================================= */
/*              As yet, unknown reason for no acceptable step, */
/*              Deal with some obvious cases. */
/*              ========================================================= */
	    if (usefd && iw[181] != 2) {
/*                 ------------------------------------------------------ */
/*                 Switch to central differences and solve the QP again. */
/*                 ------------------------------------------------------ */
		cditns = 0;
		s_wsfi(&io___470);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)132);
		iw[181] = 2;
		iw[242] = 2;
	    } else {
		if (maxns) {
		    *iexit = 33;
/* Superbasics limit */
		} else if (nearopt) {
		    *iexit = 3;
/* Requested accuracy could not be ... */
		}
		if (*iexit > 0) {
		    goto L999;
/* break major iteration loop */
		}
		if (restart) {
		    s8inith_(nnh0, nnh, &gnorm0, &gnorm, &ud0, &hd[1], &iw[1],
			     leniw, &rw[1], lenrw);
		    if (nonlinearcon) {
			s8initpen_(nncon, &penparm[1], &xpen0, &xpen[1], &rw[
				1], lenrw);
		    }
		    ++restarts;
		} else {
/*                    --------------------------------------------------- */
/*                    We have run out of things to try. Bummer. */
/*                    --------------------------------------------------- */
		    *iexit = 41;
/* Current point cannot be improved... */
		    goto L999;
/* break major iteration loop */
		}
	    }
	}
	if (! gotnewx) {
	    goto L100;
	}
/* +       until       (GotNewx) */
/*        =============================================================== */
/*        The new point  x1  has been computed. */
/*        =============================================================== */
	restarts = 0;
	inform__ = 0;
/*        --------------------------------------------------------------- */
/*        Some unknown derivatives may need to be calculated at x1. */
/*        --------------------------------------------------------------- */
	if (fonlyls && *nnh > 0) {
	    modefg = 1;
	    (*funwrapper)(iexit, &modefg, &nonlinearcon, &nonlinearobj, n, 
		    negcon, nncon0, nncon, &nnjac, nnh, nnobj0, nnobj, (U_fp)
		    funcon, (U_fp)funobj, (U_fp)userhv, &x1[1], &ycon1[1], 
		    nej, nlocj, &locj[1], &indj[1], neh, nloch, &loch[1], &
		    indh[1], &hcol[1], &fcon2[1], &fobj2, &gcon1[1], &gobj1[1]
		    , cu + 8, lencu, &iu[1], leniu, &ru[1], lenru, cw + 8, 
		    lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)8, (ftnlen)8)
		    ;
	    if (*iexit != 0) {
		goto L999;
/* break major iteration loop */
	    }
	    if (usefd) {
		s6getmissing_(iexit, n, negcon, nncon0, nncon, &nnjac, nnh, 
			nnobj0, nnobj, (S_fp)funwrapper, (U_fp)funcon, (U_fp)
			funobj, (U_fp)userhv, &bl[1], &bu[1], &x1[1], &ycon1[
			1], nej, nlocj, &locj[1], &indj[1], neh, nloch, &loch[
			1], &indh[1], &hcol[1], &fcon1[1], &fobj1, &gcon1[1], 
			&gobj1[1], &y[1], cu + 8, lencu, &iu[1], leniu, &ru[1]
			, lenru, cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, 
			(ftnlen)8, (ftnlen)8);
		if (*iexit != 0) {
		    goto L999;
/* break major iteration loop */
		}
	    }
	}
	inform__ = 0;
	++(*majors);
	if (iw[181] == 2) {
	    ++cditns;
	}
	if (mjrprint >= 10 || mnrprint >= 10) {
	    printlog = iprint > 0 && klog == 1;
	    printsum = isumm > 0 && ksumm == 1;
	    if (printlog) {
		s1page_(&c__0, &iw[1], leniw);
		s_wsfi(&io___472);
		for (j = 1; j <= 29; ++j) {
		    do_fio(&c__1, line, (ftnlen)4);
		}
		e_wsfi();
		snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)132);
		s_wsfi(&io___473);
		do_fio(&c__1, (char *)&(*majors), (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)132);
	    }
	    if (printsum && mnrprint >= 10) {
		s_wsfi(&io___474);
		for (j = 1; j <= 19; ++j) {
		    do_fio(&c__1, line, (ftnlen)4);
		}
		e_wsfi();
		snprnt_(&c__2, str, &iw[1], leniw, (ftnlen)132);
		s_wsfi(&io___475);
		do_fio(&c__1, (char *)&(*majors), (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__2, str, &iw[1], leniw, (ftnlen)132);
	    }
	}
/*        =============================================================== */
/*        The problem functions have been defined at the new x. */
/*        =============================================================== */
	if (*nnh > 0 && (lvlhess == 0 || lvlhess == 1)) {
/* Update a QN approximate Hessian. */
	    s8hqn_(&inform__, (S_fp)funwrapper, (U_fp)funcon, (U_fp)funobj, (
		    U_fp)userhv, &usefd, &iw[208], starttype, lenr, m, mbs, n,
		     nb, nncon0, nncon, &nnjac, nnh, nnobj0, nnobj, ns, 
		    majors, &qnskips, &ud0, &step, minimize, &dxhdx, &rtrmods,
		     gotr, &penparm[1], fobj, &fcon[1], &gcon[1], &gobj[1], &
		    fcon1[1], &gcon1[1], &gobj1[1], nej, nlocj, &locj[1], &
		    indj[1], &jcol[1], neh, nloch, &loch[1], &indh[1], &hcol[
		    1], negcon, nlocg, &locg[1], &kbs[1], &bl[1], &bu[1], &dx[
		    1], &dg[1], &udx[1], &hdx[1], &hd[1], &ycon1[1], &r__[1], 
		    &x[1], &x1[1], &xqp0[1], &xpen[1], &y[1], &y1[1], &y2[1], 
		    cu + 8, lencu, &iu[1], leniu, &ru[1], lenru, cw + 8, 
		    lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)8, (ftnlen)8)
		    ;
	    if (inform__ != 0) {
		*iexit = inform__;
		goto L999;
/* break,  major iteration loop */
	    }
	}
/*        --------------------------------------------------------------- */
/*        Update the variables. */
/*        The QP solution, saved in xQP, is used to start the next QP. */
/*        (If a unit step was not taken last iteration, some more */
/*        nonbasics may be between their bounds. */
/*        Nov 10, 1994. Tried leaving the nonbasics between their */
/*        bounds after short step. In some cases, the number of minor */
/*        iterations increased dramatically with a very short step.) */
/*        --------------------------------------------------------------- */
	dcopy_(nb, &x1[1], &c__1, &x[1], &c__1);
	if (nonlinearcon) {
	    dcopy_(negcon, &gcon1[1], &c__1, &gcon[1], &c__1);
	    dcopy_(nncon, &ycon1[1], &c__1, &ycon[1], &c__1);
	    dcopy_(nncon, &fcon1[1], &c__1, &fcon[1], &c__1);
	}
	if (nonlinearobj) {
	    *fobj = fobj1;
	    dcopy_(nnobj, &gobj1[1], &c__1, &gobj[1], &c__1);
	}
	*sinfe = sinfe1;
    }
/*     ======================end of main loop============================ */
L999:
    return 0;
} /* s8sqp_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8SQP */
/* Subroutine */ int s8callstatus_(integer *status, integer *iw, integer *
	leniw)
{
    /* Format strings */
    static char fmt_9999[] = "(\002 XXX  user-function call-status not recog"
	    "nized.\002,\002 Requested status =\002,i6)";

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static char str[80];
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___477 = { 0, str, 0, fmt_9999, 80, 1 };


/*     ================================================================== */
/*     s8callStatus fetches the call-status for the snOpt user-defined */
/*     functions. */

/*     16 Jun 2008: First version of s8callStatus. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
/* NP user-routine call-status */
    /* Parameter adjustments */
    --iw;

    /* Function Body */
    if (iw[236] == 0) {
/* Standard call */
	*status = 0;
    } else if (iw[236] < 0) {
/* First call */
	*status = 1;
	iw[236] = 0;
    } else if (iw[236] >= 2) {
/* Last orders please */
	*status = iw[236];
	iw[236] = -1;
    } else {
	*status = iw[236];
	s_wsfi(&io___477);
	do_fio(&c__1, (char *)&(*status), (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)80);
    }
    return 0;
} /* s8callstatus_ */

