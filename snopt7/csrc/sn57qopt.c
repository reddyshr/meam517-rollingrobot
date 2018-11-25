/* ../snopt7/src/sn57qopt.f -- translated by f2c (version 20100827).
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
static integer c__0 = 0;
static doublereal c_b5 = 0.;
static integer c__1 = 1;
static integer c_n2 = -2;
static integer c__13 = 13;
static integer c__3 = 3;
static integer c_n3 = -3;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*     File  sn57qopt.f */

/*     s5solve */
/*     s5defaults  s5Map  s5CallStatus */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Subroutine */ int s5solve_(integer *iexit, char *solver, integer *
	starttype, U_fp qplog, U_fp hprod, U_fp hprod1, logical *gotr, 
	integer *m, integer *n, integer *nb, integer *nnh0, integer *nnh, 
	integer *nnames, integer *ngqp, integer *ngobj0, integer *ngobj, 
	integer *iobj, doublereal *objadd, doublereal *objqp, doublereal *
	objtrue, integer *ninf, doublereal *sinf, integer *ninfe, doublereal *
	sinfe, integer *nea, integer *nloca, integer *loca, integer *inda, 
	doublereal *acol, integer *neh, integer *nloch, integer *loch, 
	integer *indh, doublereal *hcol, doublereal *bl, doublereal *bu, 
	doublereal *gobj, char *names, integer *nrhs0, integer *nrhs, 
	doublereal *rhs, integer *lenx0, integer *nx0, doublereal *x0, 
	integer *etype, integer *hs, doublereal *x, doublereal *pi, 
	doublereal *rc, integer *ns, char *cu, integer *lencu, integer *iu, 
	integer *leniu, doublereal *ru, integer *lenru, char *cw, integer *
	lencw, integer *iw, integer *leniw, doublereal *rw, integer *lenrw, 
	ftnlen solver_len, ftnlen names_len, ftnlen cu_len, ftnlen cw_len)
{
    /* Format strings */
    static char fmt_1900[] = "(\002 No. of iterations\002,i20,2x,\002 Object"
	    "ive\002,6x,1p,e22.10)";
    static char fmt_1910[] = "(\002 No. of infeasibilities\002,i15,2x,\002 S"
	    "um of infeas\002,1p,e24.10)";
    static char fmt_1920[] = "(\002 No. of Hessian products\002,i14,2x,\002 "
	    "Linear    objective\002,1p,e18.10)";
    static char fmt_1930[] = "(40x,\002 Quadratic objective\002,1p,e18.10)";
    static char fmt_1940[] = "(\002 No. of superbasics\002,i19,2x,\002 No. o"
	    "f basic nonlinears\002,i14)";
    static char fmt_1950[] = "(\002 Elastic weight            \002,1p,e11.1,"
	    "2x,\002 Elastic objective\002,1p,e20.10)";
    static char fmt_1960[] = "(\002 No. of infeas elastics\002,i15,2x,\002 E"
	    "lastic infeas   \002,1p,e20.10)";
    static char fmt_1970[] = "(\002 No. of CG iterations\002,i17)";
    static char fmt_1980[] = "(\002 No. of degenerate steps\002,i14,2x,\002 "
	    "Percentage\002,f27.2)";

    /* System generated locals */
    address a__1[2];
    integer i__1, i__2[2];
    char ch__1[38];

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen), s_cat(char *,
	     char **, integer *, integer *, ftnlen);
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static doublereal hcondbnd, scaleobj, zcondbnd;
    static integer elastics, minimize, probtype, mjrprint;
    static doublereal toloptfp;
    static integer mnrprint, qpsolver;
    static doublereal toloptqp;
    static integer j, k;
    static logical unbounded;
    static integer lfeastype, lr;
    static doublereal fx[1];
    static integer ly;
    static logical infeasible, terminated;
    static integer ly1, ly2, itqptarget, printlevel;
    static logical switchtoqn;
    static integer nnb, mbs, lrg, itn, liy, lkx, nkx;
    static char str[132];
    static integer lrg2, liy1, liy2;
    extern /* Subroutine */ int s5lp_(integer *, integer *, char *, logical *,
	     integer *, U_fp, logical *, logical *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, doublereal *, doublereal *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, ftnlen, ftnlen), 
	    s5qn_(integer *, integer *, char *, integer *, U_fp, U_fp, U_fp, 
	    integer *, logical *, logical *, logical *, integer *, logical *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *, integer *, integer *, integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, ftnlen, ftnlen, 
	    ftnlen), s5qp_(integer *, integer *, char *, integer *, U_fp, 
	    U_fp, U_fp, integer *, integer *, logical *, logical *, logical *,
	     integer *, logical *, integer *, integer *, integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     integer *, doublereal *, doublereal *, doublereal *, doublereal *
	    , doublereal *, doublereal *, doublereal *, doublereal *, integer 
	    *, doublereal *, integer *, integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *, doublereal *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, ftnlen, ftnlen, 
	    ftnlen);
    static char str2[132];
    static integer suboptimize, eigh;
    static logical badz;
    static integer lgbs, lkbs, lenr, lhdx, lpbs, lgqp, maxr, lxbs, maxs, itqp;
    static doublereal tolx;
    static integer ngqp0;
    static logical needb;
    static doublereal degen;
    extern /* Subroutine */ int dload_(integer *, doublereal *, doublereal *, 
	    integer *);
    static logical bigsb, badzg;
    extern /* Subroutine */ int iload_(integer *, integer *, integer *, 
	    integer *);
    static integer emode, nnjac;
    extern /* Subroutine */ int ddscl_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dddiv_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static integer lblbs;
    static logical needx;
    static integer inewb;
    static char mprob[8];
    static integer lblqp, lbubs, lbuqp, nncon, nnobj, numlc;
    static doublereal objlp, vimax;
    static logical useqp;
    static doublereal wtinf;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static doublereal xnorm;
    extern /* Subroutine */ int s2amat_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, doublereal *, integer *), s5getb_(integer *, integer *,
	     U_fp, logical *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, integer *, doublereal *, doublereal *, doublereal *
	    , integer *, integer *, integer *, integer *, doublereal *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    char *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, char *, integer *, integer *, integer 
	    *, doublereal *, integer *, ftnlen, ftnlen);
    static integer nncon0;
    extern /* Subroutine */ int s1time_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *), s4newb_(integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     integer *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, char *, char *, integer *, integer *, integer *, 
	    ftnlen, ftnlen), s5qpfg_(U_fp, U_fp, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, char *, integer *, integer *, integer *, doublereal 
	    *, integer *, char *, integer *, integer *, integer *, doublereal 
	    *, integer *, ftnlen, ftnlen);
    static doublereal wtinf0, pnorm1, pnorm2;
    extern /* Subroutine */ int s5fixs_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *), s4stat_(integer *, char *, ftnlen);
    static integer ndegen;
    static doublereal objinf;
    static logical bigitn, needlu;
    static doublereal objmin;
    static char istate[4*3];
    static integer inform__, itnlim, lssave, minmax;
    extern doublereal dnormi_(integer *, doublereal *, integer *);
    static doublereal pinorm;
    static integer numliq;
    static doublereal rgnorm;
    static logical fponly;
    extern /* Subroutine */ int snwrap_(integer *, char *, char *, char *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen);
    static integer sqstat, typelu;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen), s4saveb_(integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, doublereal *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, char *, 
	    doublereal *, doublereal *, doublereal *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, ftnlen, ftnlen);
    static logical elastic, indefqp;
    static integer lscales, lblsave;
    static logical weakmin;
    static integer lvlobje;
    static char probtag[20];
    static doublereal signobj, targeth;
    static integer lbusave, lestate;
    static logical optimal;
    static doublereal condzhz, targetz;
    static logical subitns;

    /* Fortran I/O blocks */
    static icilist io___103 = { 0, str, 0, fmt_1900, 132, 1 };
    static icilist io___104 = { 0, str, 0, fmt_1910, 132, 1 };
    static icilist io___105 = { 0, str, 0, fmt_1920, 132, 1 };
    static icilist io___106 = { 0, str, 0, fmt_1930, 132, 1 };
    static icilist io___107 = { 0, str, 0, fmt_1940, 132, 1 };
    static icilist io___108 = { 0, str, 0, fmt_1950, 132, 1 };
    static icilist io___109 = { 0, str, 0, fmt_1960, 132, 1 };
    static icilist io___110 = { 0, str, 0, fmt_1970, 132, 1 };
    static icilist io___111 = { 0, str, 0, fmt_1980, 132, 1 };


/*     ================================================================== */
/*     s5solve solves the current problem. */

/*     On entry */
/*     --------- */
/*     the SPECS file has been read, */
/*     all data items have been loaded (including Acol, indA, locA, ...), */
/*     and workspace has been allocated within cw, iw and rw. */
/*     startType = lvlStart from s3argQ. */

/*     On exit, */
/*     -------- */
/*     iExit  =  0 if an optimal solution was found, */
/*            =  1 if the problem was infeasible, */
/*            =  2 if the problem was unbounded, */
/*            =  3 if the Iteration limit was exceeded, */
/*           ge  4 if iterations were terminated by some other */
/*                 error condition (see the SQOPT user's guide). */

/*     01 Oct 1994: First version of s5solve. */
/*     06 Aug 1996: Min Sum option added. */
/*     14 Jul 1997: Thread-safe version. */
/*     02 Aug 2003: snEXIT and snPRNT adopted. */
/*     08 Mar 2004: Hot starts implemented. */
/*     16 May 2006: Explicit target itQP added. */
/*     18 Jun 2008: Added space for iy2, pBS and rg2. */
/*     07 Mar 2013: mnrPrint changed to mjrPrint in call to s2Amat. */
/*     20 Sep 2014: Upper and lower bounds copied for elastic mode. */
/*     03 Nov 2014: neH, indH, locH, Hcol added as arguments. */
/*     27 Dec 2014: Implemented switch to QN when H is indefinite. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/* scale option */
/* Hessian- vector products */
/* Current QP solver */
/* Current precon mode */
/* # of LU factorizations */
/* # lines in log     file */
/* # lines in summary file */
/* >0 => Mnr heading for iPrint */
/* >0 => Minor heading for iSumm */
/* Save the LU factors */
/* Save the reduced Hessian */
/* QP user-routine call-status */
/* Number of symmlq iterations */
/* symmlq itns for last minor */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --pi;
    --rc;
    --x;
    --hs;
    --etype;
    --bu;
    --bl;
    names -= 8;
    --gobj;
    --acol;
    --inda;
    --loca;
    --hcol;
    --indh;
    --loch;
    --rhs;
    --x0;
    cu -= 8;
    --iu;
    --ru;
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    toloptfp = rw[51];
/* Minor Phase 1 Opt tol */
    toloptqp = rw[52];
/* Minor Phase 2 Opt tol */
    tolx = rw[56];
/* Minor feasibility tolerance. */
    hcondbnd = rw[85];
/* bound on the condition of Hz */
    zcondbnd = rw[86];
/* bound on the condition of Z */
    wtinf0 = rw[88];
/* infeasibility weight */
    nnobj = iw[22];
/* # of objective variables */
    lenr = iw[28];
/* R(lenR) is the reduced Hessian factor */
    maxr = iw[52];
/* max columns of R. */
    maxs = iw[53];
/* max # of superbasics */
    qpsolver = iw[55];
/* = 0:1:2   => QPChol:CG:QN QP solver */
    emode = iw[56];
/* >0    => use elastic mode */
    lvlobje = iw[73];
/* Elastic objective type */
    minmax = iw[87];
/* 1, 0, -1  => MIN, FP, MAX */
    itnlim = iw[89];
/* limit on total iterations */
    mnrprint = iw[93];
/* Minor print level */
    mjrprint = iw[92];
/* Major print level */
    inewb = iw[124];
/* new basis file */
    minimize = iw[199];
/* 1 (-1)    => minimize (maximize) */
    nkx = iw[247];
/* dimension of kx and its inverse, kxN */
    s_copy(mprob, cw + 408, (ftnlen)8, (ftnlen)8);
/* Addresses */
/* Problem name */
    lkx = iw[251];
/* j  = kx (jN) => col j of Jcol is variable jN */
    lfeastype = iw[284];
/* feasType(mBS) = feasibility types */
    lestate = iw[285];
/* eState(nb)    = status of elastics */
    lkbs = iw[292];
/* kBS(mBS)      = ( B  S ) list */
    lblqp = iw[271];
/* blQP(nb)      = working QP lower bounds */
    lbuqp = iw[272];
/* buQP(nb)      = working QP upper bounds */
    lblbs = iw[273];
/* blBS(mBS)     = lower bounds for xBS */
    lbubs = iw[274];
/* buBS(mBS)     = upper bounds for xBS */
    lpbs = iw[277];
/* pBS(nb)       = search direction */
    lxbs = iw[301];
/* xBS(mBS)      = basics, superbasics */
    lgqp = iw[290];
/* gQP(ngQP)     = QP gradient */
    lgbs = iw[291];
/* gBS(mBS)      = BS components of g */
    lrg = iw[293];
/* rg (maxS)     = reduced gradient */
    lrg2 = iw[294];
/* rg2(maxS)     = reduced gradient copy */
    lr = iw[295];
/* R(lenR)       = factor of Z'HZ */
    lscales = iw[296];
/* scales(nb)    = row and column scales */
    liy = iw[308];
/* iy (nb)       = integer work vector */
    liy1 = iw[309];
/* iy1(nb)       = integer work vector */
    liy2 = iw[310];
/* iy2(nb)       = integer work vector */
    ly = iw[311];
/*  y (nb)       = real work vector */
    ly1 = iw[312];
/*  y1(nb)       = real work vector */
    ly2 = iw[313];
/*  y2(nb)       = real work vector */
    lhdx = iw[288];
/* Hdx(nnH)      = product of H with  x - x0 */
    lblsave = iw[275];
/* blSave        = bl for the basis-finding LP */
    lbusave = iw[276];
/* buSave        = bu for the basis-finding LP */
    *iexit = 0;
    mbs = *m + maxs;
/* Figure out what type of objective we have. */
    if (minmax == 0 || emode == 2 && lvlobje == 2) {
	probtype = 0;
    } else if (*ngqp == 0) {
/* No explicit objective. Must be an LP. */
	if (*iobj == 0) {
	    probtype = 0;
	} else {
	    probtype = 1;
	}
    } else {
/*  Explicit objective. Check for quadratic term. */
	if (*nnh > 0) {
	    probtype = 2;
	} else {
	    probtype = 1;
	}
    }
    fponly = probtype == 0;
    s_copy(probtag, "linear constraints", (ftnlen)20, (ftnlen)18);
    iw[223] = 0;
/* Print the header for the Print   file */
    iw[225] = 0;
/* Print the header for the summary file */
    iw[220] = 0;
/* Line count for the print   file */
    iw[221] = 0;
/* Line count for the summary file */
    printlevel = mnrprint;
/*     Initialize counters based on gotHes and gotFac (set in s3prtQ) */
    if (iw[230] <= 0) {
	iw[210] = 0;
    }
    if (iw[231] <= 0) {
	iw[188] = 0;
    }
    iw[386] = 0;
    iw[387] = 0;
    itn = 0;
    itqp = 0;
    ndegen = 0;
    *ninf = 0;
    *ninfe = 0;
    nncon = 0;
    nncon0 = 1;
    nnjac = 0;
    numlc = *m;
    ngqp0 = max(*ngqp,1);
    iw[208] = qpsolver;
/* Local value of QPslvr */
    *objqp = 0.;
    *sinf = 0.;
    *sinfe = 0.;
    scaleobj = 1.;
    signobj = (doublereal) minimize;
    targeth = hcondbnd;
    targetz = zcondbnd;
    wtinf = wtinf0;
    suboptimize = 0;
/*     No suboptimization */
    itqptarget = itnlim;
/*     No suboptimization */
    switchtoqn = FALSE_;
/*     Start recording the solve time. */
    s1time_(&c__2, &c__0, &iw[1], leniw, &rw[1], lenrw);
/* Initialize quantities to avoid them being used before being set. */
    dload_(m, &c_b5, &pi[1], &c__1);
    dload_(&ngqp0, &c_b5, &rw[lgqp], &c__1);
    iload_(nb, &c__0, &iw[lestate], &c__1);
/*     Make copies of the upper and lower bounds. */
    dcopy_(nb, &bl[1], &c__1, &rw[lblqp], &c__1);
    dcopy_(nb, &bu[1], &c__1, &rw[lbuqp], &c__1);
/* ----------------------------------------------------------------- */
/* Print the matrix statistics. */
/* Find the rowtypes for use in s5getB (they are held in iy2). */
/* ----------------------------------------------------------------- */
    s2amat_(&c__1, &mjrprint, m, n, nb, &nncon, &nnjac, &nnobj, iobj, &numlc, 
	    &numliq, nea, nloca, &loca[1], &inda[1], &acol[1], &rw[lblqp], &
	    rw[lbuqp], &iw[liy2], &iw[1], leniw, &rw[1], lenrw);
/* ================================================================= */
/* Find a basis kBS(1:m) for the linear constraints and bounds. */
/* ================================================================= */
/* s5getB does the following. */
/*  1. The linear constraints are (optionally) scaled. */
/*  2. Elements x(n+1:n+m) of the initial x are assigned. */
/*  3. An LP is used to find a feasible x for the bounds and */
/*     linear equality constraints. */
/*  The base point x0 is not touched. */
    s5getb_(&inform__, starttype, (U_fp)qplog, &needb, m, &maxs, &mbs, n, nb, 
	    &nncon, &nnjac, &nnobj, nnames, ns, &itqp, &itnlim, &itn, &ndegen,
	     &numlc, &numliq, &toloptfp, &toloptqp, &tolx, ninf, sinf, &wtinf,
	     iobj, &scaleobj, &pinorm, &rgnorm, nea, nloca, &loca[1], &inda[1]
	    , &acol[1], &etype[1], &iw[lestate], &iw[liy2], &iw[lfeastype], &
	    hs[1], &iw[lkbs], names + 8, &bl[1], &bu[1], &rw[lblqp], &rw[
	    lbuqp], &rw[lblbs], &rw[lbubs], &rw[lblsave], &rw[lbusave], &rw[
	    lgbs], &pi[1], &rc[1], nrhs0, nrhs, &rhs[1], &rw[lscales], &c__1, 
	    &c__0, &x0[1], &x[1], &rw[lxbs], &iw[liy], &iw[liy1], &rw[ly], &
	    rw[ly1], &rw[ly2], cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (
	    ftnlen)8, (ftnlen)8);
/*     Potential inform values are: */
/*        0   basis found. */
/*            nInf = 0 => linear equalities are     satisfied. */
/*            nInf > 0 +> linear equalities are not satisfied */
/*                        unbounded FP problem */
/*       >0   fatal error. No feasible point for the equalities */
    if (inform__ > 0) {
	*iexit = inform__;
/* fatal error */
	goto L900;
    }
    if (*ngobj > 0 && iw[75] > 0) {
	ddscl_(ngobj, &rw[lscales], &c__1, &gobj[1], &c__1);
    }
    if (numliq > 0) {
	needlu = TRUE_;
    } else {
	needlu = FALSE_;
    }
    needx = needlu;
    typelu = 3;
/*     ------------------------------------------------------------------ */
/*     s5getB has already found a basis for the linear constraints that */
/*     may or may not be feasible. */
/*     ------------------------------------------------------------------ */
    elastic = emode == 2;
    elastics = 0;
    iload_(nb, &c__0, &iw[lestate], &c__1);
/*     ================================================================== */
/*     Solve the problem. */
/*     ================================================================== */
    useqp = probtype == 1 && *ngobj > 0 || probtype == 2;
    if (useqp) {
	iw[223] = 1;
/* Refresh print   heading. */
	iw[225] = 1;
/*        Unless CG was requested explicitly,  use maxS = maxR */
/*        with Cholesky and QN. Then switch to CG if necessary. */
/* Refresh summary heading */
	if (iw[208] == 0) {
/* ----------------------------------------------------------- */
/* Solve the QP. */
/* ----------------------------------------------------------- */
	    eigh = 0;
	    *gotr = FALSE_;
	    s5qp_(&inform__, &probtype, probtag, &suboptimize, (U_fp)qplog, (
		    U_fp)hprod, (U_fp)hprod1, &iw[188], &eigh, &elastic, gotr,
		     &needlu, &typelu, &needx, &lenr, m, &maxr, &mbs, n, nb, &
		    ndegen, &ngqp0, ngqp, ngobj0, ngobj, nnh0, nnh, ns, &itqp,
		     &itnlim, &itqptarget, &itn, &emode, &lvlobje, &
		    printlevel, &minimize, iobj, &scaleobj, objadd, objqp, &
		    targeth, &targetz, &toloptfp, &toloptqp, &tolx, ninf, 
		    sinf, &elastics, ninfe, sinfe, &wtinf, &pinorm, &rgnorm, 
		    nea, nloca, &loca[1], &inda[1], &acol[1], neh, nloch, &
		    loch[1], &indh[1], &hcol[1], &etype[1], &iw[lestate], &iw[
		    lfeastype], &hs[1], &iw[lkbs], &bl[1], &bu[1], &rw[lblqp],
		     &rw[lbuqp], &rw[lblbs], &rw[lbubs], &rw[lgbs], &gobj[1], 
		    &rw[lgqp], &rw[lhdx], &rw[lpbs], &pi[1], &rw[lr], &rc[1], 
		    &rw[lrg], nrhs0, nrhs, &rhs[1], &rw[lscales], lenx0, nx0, 
		    &x0[1], &x[1], &rw[lxbs], &x[1], &iw[liy], &iw[liy1], &rw[
		    ly], &rw[ly1], &rw[ly2], cu + 8, lencu, &iu[1], leniu, &
		    ru[1], lenru, cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw,
		     (ftnlen)20, (ftnlen)8, (ftnlen)8);
/* xFrozen = x, unused */
	} else if (iw[208] == 2) {
	    *gotr = iw[231] > 0;
	    s5qn_(&inform__, &probtype, probtag, &suboptimize, (U_fp)qplog, (
		    U_fp)hprod, (U_fp)hprod1, &iw[188], &elastic, gotr, &
		    needlu, &typelu, &needx, &lenr, m, &maxr, &mbs, n, nb, &
		    ndegen, &ngqp0, ngqp, ngobj0, ngobj, nnh0, nnh, ns, &itqp,
		     &itnlim, &itqptarget, &itn, &emode, &lvlobje, &
		    printlevel, &minimize, iobj, &scaleobj, objadd, objqp, &
		    condzhz, &targetz, &toloptfp, &toloptqp, &tolx, ninf, 
		    sinf, &elastics, ninfe, sinfe, &wtinf, &pinorm, nea, 
		    nloca, &loca[1], &inda[1], &acol[1], neh, nloch, &loch[1],
		     &indh[1], &hcol[1], &etype[1], &iw[lestate], &iw[
		    lfeastype], &hs[1], &iw[lkbs], &bl[1], &bu[1], &rw[lblqp],
		     &rw[lbuqp], &rw[lblbs], &rw[lbubs], &rw[lgbs], &gobj[1], 
		    &rw[lgqp], &rw[lhdx], &rw[lpbs], &pi[1], &rw[lr], &rc[1], 
		    &rw[lrg], &rw[lrg2], nrhs0, nrhs, &rhs[1], &rw[lscales], 
		    lenx0, nx0, &x0[1], &x[1], &rw[lxbs], &x[1], &iw[liy], &
		    iw[liy1], &rw[ly], &rw[ly1], &rw[ly2], cu + 8, lencu, &iu[
		    1], leniu, &ru[1], lenru, cw + 8, lencw, &iw[1], leniw, &
		    rw[1], lenrw, (ftnlen)20, (ftnlen)8, (ftnlen)8);
/* xFrozen = x, unused */
	}
	if (inform__ == -5 && maxr < maxs) {
/*           Too many superbasics. Switch to a CG solver */
	    if (iw[208] == 0) {
		iw[209] = 0;
/* with no preconditioning */
		*gotr = FALSE_;
	    } else if (iw[208] == 2) {
		iw[209] = 1;
/* with QN preconditioning */
	    }
	    iw[208] = 1;
/* Switch to CG */
	} else if (inform__ == -6) {
/*           Hessian Indefinite or ill-conditioned Switch to QN solver. */
	    iw[208] = 2;
/* Switch to QN */
	    switchtoqn = TRUE_;
	    *gotr = TRUE_;
	}
	if (iw[208] == 1 || switchtoqn) {
	    if (qpsolver == 1) {
		*gotr = FALSE_;
	    }
	    s5qn_(&inform__, &probtype, probtag, &suboptimize, (U_fp)qplog, (
		    U_fp)hprod, (U_fp)hprod1, &iw[188], &elastic, gotr, &
		    needlu, &typelu, &needx, &lenr, m, &maxs, &mbs, n, nb, &
		    ndegen, &ngqp0, ngqp, ngobj0, ngobj, nnh0, nnh, ns, &itqp,
		     &itnlim, &itqptarget, &itn, &emode, &lvlobje, &
		    printlevel, &minimize, iobj, &scaleobj, objadd, objqp, &
		    condzhz, &targetz, &toloptfp, &toloptqp, &tolx, ninf, 
		    sinf, &elastics, ninfe, sinfe, &wtinf, &pinorm, nea, 
		    nloca, &loca[1], &inda[1], &acol[1], neh, nloch, &loch[1],
		     &indh[1], &hcol[1], &etype[1], &iw[lestate], &iw[
		    lfeastype], &hs[1], &iw[lkbs], &bl[1], &bu[1], &rw[lblqp],
		     &rw[lbuqp], &rw[lblbs], &rw[lbubs], &rw[lgbs], &gobj[1], 
		    &rw[lgqp], &rw[lhdx], &rw[lpbs], &pi[1], &rw[lr], &rc[1], 
		    &rw[lrg], &rw[lrg2], nrhs0, nrhs, &rhs[1], &rw[lscales], 
		    lenx0, nx0, &x0[1], &x[1], &rw[lxbs], &x[1], &iw[liy], &
		    iw[liy1], &rw[ly], &rw[ly1], &rw[ly2], cu + 8, lencu, &iu[
		    1], leniu, &ru[1], lenru, cw + 8, lencw, &iw[1], leniw, &
		    rw[1], lenrw, (ftnlen)20, (ftnlen)8, (ftnlen)8);
/*           If the QP is unbounded and the reduced Hessian was flagged */
/*           as indefinite in s5QP, the problem is likely unbounded. */
/* xFrozen = x, unused */
	    if (inform__ == -2 && switchtoqn) {
		inform__ = -6;
	    }
	}
    } else {
/* -------------------------------------------------------------- */
/* LP with objective row in A. */
/* -------------------------------------------------------------- */
	*ns = 0;
/* Local value */
	*nnh = 0;
/* Local value */
	s5fixs_(&c__0, m, &maxs, &mbs, n, nb, ns, &hs[1], &iw[lkbs], &rw[
		lblqp], &rw[lbuqp], &rw[lblbs], &rw[lbubs], &x[1], &rw[lxbs]);
	iw[223] = 0;
	iw[225] = 0;
	s5lp_(iexit, &probtype, probtag, &elastic, &suboptimize, (U_fp)qplog, 
		&needlu, &needx, m, n, nb, &ndegen, &itqp, &itnlim, &itn, &
		emode, &lvlobje, &printlevel, &minimize, iobj, &scaleobj, 
		objadd, &toloptfp, &toloptqp, &tolx, ninf, sinf, &elastics, 
		ninfe, sinfe, &wtinf, &pinorm, &rgnorm, nea, nloca, &loca[1], 
		&inda[1], &acol[1], &etype[1], &iw[lestate], &iw[lfeastype], &
		hs[1], &iw[lkbs], &bl[1], &bu[1], &rw[lblqp], &rw[lbuqp], &rw[
		lblbs], &rw[lbubs], &rw[lgbs], &pi[1], &rc[1], nrhs0, nrhs, &
		rhs[1], &rw[lscales], &x[1], &rw[lxbs], &x[1], &iw[liy], &iw[
		liy1], &rw[ly], &rw[ly1], cw + 8, lencw, &iw[1], leniw, &rw[1]
		, lenrw, (ftnlen)20, (ftnlen)8);
/* xFrozen = x, unused */
    }
    terminated = inform__ > 0;
/* Fatal error (e.g., LU, time limit) */
    optimal = inform__ == 0;
/* Optimal multipliers */
    infeasible = inform__ == -1;
/* infeas nonelastics in elastic mode */
    unbounded = inform__ == -2;
/* LP is unbounded */
    bigitn = inform__ == -3;
/* Too many iterations */
    weakmin = inform__ == -4;
/* Weak QP minimizer */
    bigsb = inform__ == -5;
/* Too many superbasics */
    indefqp = inform__ == -6 || inform__ == -9;
/* QP Hessian not positive semidef */
    badzg = inform__ == -7;
/* Z'g could not be made small enough */
    badz = inform__ == -8;
/* Ill-conditioned Z */
    subitns = inform__ == -10;
/* Many possible outcomes! */
/* Too many subiterations. */
    if (terminated) {
/* LU error or time limit */
	*iexit = inform__;
    } else if (optimal && *ninf > 0) {
/* Optimal multipliers with infeasible constraints. */
	if (elastic) {
	    *iexit = 16;
/* infeasible nonelastics */
	} else {
	    *iexit = 11;
/* infeasible linear constraints */
	}
    } else if (optimal && *ninf == 0) {
/* Optimal multipliers with feasible nonelastic constraints. */
	if (! elastic) {
/* Optimal multipliers in normal mode. */
	    *iexit = 1;
/* Optimal */
	} else if (*ninfe > 0) {
/* Elastic mode with some infeasible elastic constraints. */
	    if (lvlobje == 1) {
		*iexit = 5;
/* elastic objective minimized */
	    } else if (lvlobje == 2) {
		if (fponly) {
		    *iexit = 6;
/* elastic infeasibilities minimized */
		} else {
		    *iexit = 14;
/* infeasibilities minimized */
		}
	    }
	} else {
/* Elastic mode with feasible elastics and nonelastics. */
	    if (fponly) {
		*iexit = 2;
/* Feasible point found. */
	    } else {
		*iexit = 1;
/* Optimal */
	    }
	}
    } else if (infeasible) {
	*iexit = 16;
/* infeasible nonelastics */
    } else if (unbounded) {
	*iexit = 21;
/* unbounded */
    } else if (bigitn) {
	*iexit = 31;
/* too many iterations */
    } else if (weakmin) {
	*iexit = 4;
/* weak minimizer */
    } else if (bigsb) {
	*iexit = 33;
/* too many superbasics */
    } else if (indefqp) {
	*iexit = 53;
/* QP Hessian is indefinite */
    } else {
	*iexit = 41;
/* Current point cannot be improved */
    }
/*     ================================================================== */
/*     Exit. */
/*     Set output variables and print a summary of the final solution. */
/*     objTrue is printed in s4newB */
/*     ================================================================== */
L900:
    snwrap_(iexit, solver, str, str2, &iw[1], leniw, (ftnlen)6, (ftnlen)132, (
	    ftnlen)132);
    s1time_(&c_n2, &c__0, &iw[1], leniw, &rw[1], lenrw);
    degen = ndegen * 100. / max(itn,1);
/*     Print statistics. */
    *objtrue = 0.;
    if (*iobj == 0) {
	objlp = *objadd;
    } else {
	objlp = *objadd + x[*n + *iobj] * scaleobj;
    }
    infeasible = *ninf > 0;
    xnorm = dnormi_(n, &x[1], &c__1);
    objmin = 0.;
    objinf = 0.;
    if (elastic) {
	*objtrue = objlp;
	if (*ngqp > 0) {
	    *objtrue += *objqp;
	}
	objmin = signobj * *objtrue + wtinf * *sinfe;
	objinf = *sinfe;
    } else {
/* Normal mode */
	*objtrue = objlp;
	if (*ngqp > 0) {
	    *objtrue += *objqp;
	}
	objmin = signobj * *objtrue;
    }
/* Count basic nonlinear variables (used only for printing). */
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
		sinf, objtrue, &iw[lkbs], &hs[1], &rw[lscales], &bl[1], &bu[1]
		, &x[1], &rw[lxbs], istate, cw + 8, lencw, &iw[1], leniw, (
		ftnlen)4, (ftnlen)8);
    }
/*     ------------------------------------------------------------------ */
/*     Print statistics. */
/*     ------------------------------------------------------------------ */
/* Writing concatenation */
    i__2[0] = 30, a__1[0] = " Problem name                 ";
    i__2[1] = 8, a__1[1] = mprob;
    s_cat(ch__1, a__1, i__2, &c__2, (ftnlen)38);
    snprnt_(&c__13, ch__1, &iw[1], leniw, (ftnlen)38);
    s_wsfi(&io___103);
    do_fio(&c__1, (char *)&itn, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*objtrue), (ftnlen)sizeof(doublereal));
    e_wsfi();
    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)132);
    if (infeasible) {
	s_wsfi(&io___104);
	do_fio(&c__1, (char *)&(*ninf), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&(*sinf), (ftnlen)sizeof(doublereal));
	e_wsfi();
	snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)132);
    }
    if (probtype == 2) {
	s_wsfi(&io___105);
	do_fio(&c__1, (char *)&iw[188], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&objlp, (ftnlen)sizeof(doublereal));
	e_wsfi();
	snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)132);
	s_wsfi(&io___106);
	do_fio(&c__1, (char *)&(*objqp), (ftnlen)sizeof(doublereal));
	e_wsfi();
	snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)132);
    }
    if (*ns > 0) {
	s_wsfi(&io___107);
	do_fio(&c__1, (char *)&(*ns), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&nnb, (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)132);
    }
    if (elastic) {
	if (*ninfe > 0) {
	    s_wsfi(&io___108);
	    do_fio(&c__1, (char *)&wtinf, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&objmin, (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)132);
	    s_wsfi(&io___109);
	    do_fio(&c__1, (char *)&(*ninfe), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&objinf, (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)132);
	}
    } else {
/* Normal mode */
	if (iw[386] > 0) {
	    s_wsfi(&io___110);
	    do_fio(&c__1, (char *)&iw[386], (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)132);
	}
    }
    s_wsfi(&io___111);
    do_fio(&c__1, (char *)&ndegen, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&degen, (ftnlen)sizeof(doublereal));
    e_wsfi();
    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)132);
/*     ------------------------------------------------------------------ */
/*     Unscale, save basis files and prepare to print the solution. */
/*     Clock 3 is "Output time". */
/*     ------------------------------------------------------------------ */
    s1time_(&c__3, &c__0, &iw[1], leniw, &rw[1], lenrw);
    s4saveb_(&inform__, &c__0, &minimize, m, n, nb, &nkx, &nncon0, &nncon, &
	    ngqp0, ngqp, nnames, ns, &itn, ninf, sinf, &wtinf, &vimax, iobj, &
	    scaleobj, objtrue, &pnorm1, &pnorm2, &pinorm, &xnorm, nea, nloca, 
	    &loca[1], &inda[1], &acol[1], &iw[lkx], &iw[lestate], &hs[1], &rw[
	    lscales], &bl[1], &bu[1], fx, &rw[lgqp], names + 8, &pi[1], &rc[1]
	    , &x[1], cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)8, (
	    ftnlen)8);
    if (*ngobj > 0 && iw[75] > 0) {
	dddiv_(ngobj, &rw[lscales], &c__1, &gobj[1], &c__1);
    }
/*     If task = PrintS, s4saveB prints the solution under the control */
/*     of lprSol (set by the  Solution  keyword in the SPECS file). */
/*     The printed solution may or may not be wanted, as follows: */

/*     lprSol = 0   means      No */
/*            = 1   means      If optimal, infeasible or unbounded */
/*            = 2   means      Yes */
/*            = 3   means      If error condition */
    s4saveb_(&inform__, &c__1, &minimize, m, n, nb, &nkx, &nncon0, &nncon, &
	    ngqp0, ngqp, nnames, ns, &itn, ninf, sinf, &wtinf, &vimax, iobj, &
	    scaleobj, objtrue, &pnorm1, &pnorm2, &pinorm, &xnorm, nea, nloca, 
	    &loca[1], &inda[1], &acol[1], &iw[lkx], &iw[lestate], &hs[1], &rw[
	    lscales], &bl[1], &bu[1], fx, &rw[lgqp], names + 8, &pi[1], &rc[1]
	    , &x[1], cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)8, (
	    ftnlen)8);
    s1time_(&c_n3, &c__0, &iw[1], leniw, &rw[1], lenrw);
/*     ------------------------------------------------------------------ */
/*     Set Obj for output. */
/*     Call  Hx  one last time with  nState .ge. 2. */
/*     Everything has been  unscaled, so we have to disable scaling. */
/*     ------------------------------------------------------------------ */
    lssave = iw[75];
    iw[75] = 0;
/* Computing MIN */
    i__1 = *iexit / 10;
    sqstat = min(i__1,4) + 2;
    iw[235] = sqstat;
    objlp = 0.;
    if (probtype == 0) {
	*objtrue = 0.;
    } else if (probtype == 1 || probtype == 2) {
	*objtrue = *objadd;
	if (*iobj > 0) {
	    objlp = x[*n + *iobj] * scaleobj;
	    *objtrue += objlp;
	}
	if (*ngqp > 0) {
	    s5qpfg_((U_fp)hprod, (U_fp)hprod1, &iw[188], ngqp, ngobj0, ngobj, 
		    nnh, neh, nloch, &loch[1], &indh[1], &hcol[1], &sqstat, 
		    objqp, &gobj[1], &rw[lgqp], lenx0, nx0, &x0[1], &x[1], &
		    rw[ly], cu + 8, lencu, &iu[1], leniu, &ru[1], lenru, cw + 
		    8, lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)8, (
		    ftnlen)8);
	    *objtrue += *objqp;
	}
    }
    iw[75] = lssave;
/*     Save some things needed by solvers calling s5Solve */
    rw[421] = *objtrue;
/* The true objective */
    rw[422] = pinorm;
/* Lagrange multiplier norm */
    rw[423] = xnorm;
/* Norm of the variables (for GAMS) */
    rw[424] = wtinf;
/* Infeasibility weight */
    rw[433] = *sinf + *sinfe;
/* Sum of infeasibilities */
    rw[434] = objlp;
/* Linear    objective term */
    iw[421] = itn;
/* Total iteration count */
    iw[423] = maxs;
/*     No nonlinear constraints. */
/* max # of superbasics */
    rw[432] = 0.;
/* Inf norm of the constraint violation */
    rw[435] = 0.;
/* Nonlinear objective term */
    rw[436] = 0.;
/* Norm of penalty parameters */
    iw[422] = 0;
/* Major iterations */
    return 0;
} /* s5solve_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5solve */
/* Subroutine */ int s5defaults_(integer *m, integer *n, integer *lencobj, 
	integer *ncolh, char *cw, integer *lencw, integer *iw, integer *leniw,
	 doublereal *rw, integer *lenrw, ftnlen cw_len)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static doublereal hcondbnd;
    static integer lvlscale, minimize, nparprlp, nparprqp, mjrprint, stickyop;
    static doublereal toloptfp;
    static integer mnrprint;
    static doublereal toloptqp, c4, c6;
    static logical qp;
    static doublereal eps, eps0, eps1, eps2, eps3, eps4;
    static integer npr1, npr2, kfac, kchk, klog, ksav, maxr, imps, maxs, nout;
    static doublereal tolx, dens1, dens2, lmax1, lmax2, utol1, utol2;
    static integer iback, ioldb, emode;
    static doublereal bigdx, bigfx;
    static integer ipnch;
    static doublereal etarg;
    static integer inewb;
    static doublereal small, tolcg, xdlim;
    static integer idump, maxmn, never, mskip, isoln;
    static doublereal rmaxs;
    static integer ksumm;
    static doublereal toldj3, wtinf0, utol1m, utol2m;
    static integer iloadb, kdegen;
    static doublereal infbnd, zcndbd, chzbnd;
    static integer icrash;
    static logical linear;
    static integer lprdbg;
    static doublereal tolfac, uspace;
    static integer maxcol;
    static doublereal tcrash, toldcp;
    static integer precon;
    static doublereal tolddp;
    static integer minprc, minmax, cgitmx, itnlim, kreset, mflush, lprscl, 
	    mnewsb, mminor;
    static doublereal scltol, tolcon;
    static integer lvlpre, iprint, ireprt;
    static doublereal toldpp;
    static integer iinsrt;
    static doublereal toldrp, toldup;
    static integer lprsol, lprprm, lvlpiv;
    static doublereal tolnlp, tolpiv, tolupd;
    static integer tpivot;
    static doublereal tolrow;
    static integer qpslvr;
    static doublereal tolswp;
    static integer lvlsys, lvlobje;
    static doublereal maxtime;
    static integer nparpru;

/*     ================================================================== */
/*     s5defaults checks and possibly prints the optional parameter values */
/*     for sqopt. */

/*     Optional parameters are checked and, if necessary,  changed to */
/*     reasonable values. */

/*     Note that parameters are checked before the amount of working */
/*     storage has been defined. */

/*     See  snworkspace.info  for full documentation of cw, iw and rw. */

/*     15 Nov 1991: first version. */
/*     02 Aug 2003: snPRNT adopted. */
/*     22 Jun 2004: Added default LU mod singularity tol */
/*     21 Dec 2004: Default LU tols fixed up. */
/*     02 May 2006: lvlTim removed. */
/*     01 Sep 2007: stickyOp added. */
/*     18 Jan 2010: PreCon initialized. */
/*     13 Jul 2013: Default lvldif set correctly. */
/*     25 Oct 2014: Bogus assignment of LUprnt removed (its set in s2BLU). */
/*     25 Oct 2014: nout set independently of lvlSys. */
/*     27 Dec 2014: tolOptFP set independently of tolOptQP. */
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
/*     ------------------------------------------------------------------ */
/*     rw(51)--rw(150): optional parameters set via the specs file. */
/*     ------------------------------------------------------------------ */
/* eps**(1/4)          IEEE DP  1.22e-04 */
    toloptfp = rw[51];
/* Minor Phase 1 Opt tol */
    toloptqp = rw[52];
/* Minor Phase 2 Opt tol */
    tolcg = rw[54];
/* cg tolerance */
    tolx = rw[56];
/* Minor feasibility tolerance */
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
/* LU factor tolerance */
    tolupd = rw[67];
/* LU update tolerance */
    infbnd = rw[70];
/* definition of plus infinity */
    bigfx = rw[71];
/* unbounded objective */
    bigdx = rw[72];
/* unbounded step */
    lvlpre = iw[77];
/* >0    => QN preconditioned CG */
    maxtime = rw[79];
/* time limit */
    xdlim = rw[80];
/* Step limit */
    etarg = rw[83];
/* Quasi-Newton QP rg tolerance */
    hcondbnd = rw[85];
/* bound on the condition of Hz */
    zcndbd = rw[86];
/* bound on the condition of Z */
    wtinf0 = rw[88];
/* infeasibility weight */
    scltol = rw[92];
/*     ------------------------------------------------------------------ */
/*     rw(151)--rw(180) contain  parmLU  parameters for LUSOL. */
/*     ------------------------------------------------------------------ */
/* scale tolerance. */
    lmax1 = rw[151];
/* max L-multiplier in factor */
    lmax2 = rw[152];
/* max L-multiplier in update */
    small = rw[153];
/* defn of small real */
    utol1 = rw[154];
/* abs tol for small diag of U */
    utol2 = rw[155];
/* rel tol for small diag of U */
    uspace = rw[156];
/* limit on waste space in U */
    dens1 = rw[157];
/* switch to search maxcol columns and no rows */
    dens2 = rw[158];
/*     ------------------------------------------------------------------ */
/*     rw(181)--rw(199) pass parameters into various routines. */
/*     ------------------------------------------------------------------ */
/*     toldj3    = rw(186) ! current optimality tol */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
/*     iw(1)--iw(50): I/O file numbers and dimensions. */
/*     ------------------------------------------------------------------ */
/* switch to dense LU */
    iprint = iw[12];
/*     ------------------------------------------------------------------ */
/*     iw(51)--iw(150): optional parameters set via the specs file. */
/*     ------------------------------------------------------------------ */
/* Print file */
    maxr = iw[52];
/* max columns of R. */
    maxs = iw[53];
/* max # of superbasics */
    qpslvr = iw[55];
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
/*     lvlStart  = iw( 69) ! = 0:1:2:3 => cold:warm:basis:hot start */
/* # largest value of nSkip */
    lvlsys = iw[71];
/* > 0   => print system info */
    lvlobje = iw[73];
/* Elastic option */
    lvlscale = iw[75];
/* scale option */
    lvlpiv = iw[80];
/* 0(1) LU threshold partial(complete) pivoting */
    lprprm = iw[81];
/* > 0    => parms are printed */
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
    mminor = iw[91];
/* limit on minor iterations */
    mnrprint = iw[93];
/* Minor print level */
    nparpru = iw[101];
/* # of partial pricing sections */
    mnewsb = iw[95];
/* # of working set changes */
    cgitmx = iw[97];
/* CG iteration limit */
    stickyop = iw[116];
/* > 0 => optional parameters are sticky */
    iback = iw[120];
/* backup file */
    idump = iw[121];
/* dump file */
    iloadb = iw[122];
/* load file */
    imps = iw[123];
/* MPS file */
    inewb = iw[124];
/* new basis file */
    iinsrt = iw[125];
/* insert file */
    ioldb = iw[126];
/* old basis file */
    ipnch = iw[127];
/* punch file */
    ireprt = iw[130];
/* Report file */
    isoln = iw[131];
/*     ------------------------------------------------------------------ */
/*     iw(151)--iw(180) contain luparm parameters for LUSOL. */
/*     ------------------------------------------------------------------ */
/* Solution file */
    nout = iw[151];
/* unit # for printed messages */
    maxcol = iw[153];
/*     ------------------------------------------------------------------ */
/* lu1fac: max. # columns */
    precon = iw[209];
/* Current precon mode (based on QPslvr) */
    c4 = max(1e-4,eps3);
    c6 = max(1e-6,eps2);
    never = 99999999;
    qp = *ncolh > 0;
    linear = ! qp;
/*     ================================================================== */
/*     Check the optional parameters. */
/*     ================================================================== */
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
	if (qp) {
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
/*     Sometimes, frequency 0 means "almost never". */
    if (kchk <= 0) {
	kchk = never;
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
    if (icrash < 0) {
	icrash = 3;
    }
    if (minmax == -11111) {
	minmax = 1;
    }
    if (minmax == -1) {
	minimize = -1;
    } else {
	minimize = 1;
    }
    if (mminor < 0) {
/* Computing MAX */
	i__1 = 1000, i__2 = max(*n,*m) * 5;
	mminor = max(i__1,i__2);
    }
    if (mnewsb <= 0) {
	mnewsb = never;
    }
    if (lprdbg < 0) {
	lprdbg = 0;
    }
    if (lprprm < 0) {
	lprprm = 1;
    }
    if (lprscl < 0) {
	lprscl = 0;
    }
    if (lprsol < 0) {
	lprsol = 2;
    }
/*     lvlStart is checked in s3argQ */
/*     if (lvlStart  .lt. 0     ) lvlStart = 0 */
    if (mnrprint < 0) {
	mnrprint = 1;
    }
    mjrprint = mnrprint;
    if (lvlobje < 0 || lvlobje > 2) {
	lvlobje = -11111;
    }
    if (lvlobje == -11111) {
	lvlobje = 2;
    }
    if (lvlsys < 0) {
	lvlsys = 0;
    }
    if (emode < 0 || emode > 2) {
	emode = -11111;
    }
    if (emode == -11111) {
	emode = 1;
    }
    if (stickyop < 0) {
	stickyop = 0;
    }
/*     Check superbasics limit and reduced Hessian size. */
    if (qp) {
	if (maxr < 0) {
/* Computing MIN */
	    i__1 = 2000, i__2 = *ncolh + 1;
	    maxr = min(i__1,i__2);
	}
	if (maxs < 0) {
	    maxs = *ncolh + 1;
	}
/* Computing MAX */
	i__1 = min(maxr,*n);
	maxr = max(i__1,0);
/* Computing MAX */
	i__1 = min(maxs,*n);
	maxs = max(i__1,1);
    } else {
/* linear */
	if (maxs <= 0) {
	    maxs = 1;
	}
	maxr = 1;
    }
    if (maxs < maxr) {
	maxs = maxr;
    }
    if (qpslvr < 0) {
	qpslvr = 0;
    }
    if (maxr == 0) {
	qpslvr = 1;
    }
    if (lvlpre < 0 || lvlpre > 1) {
	lvlpre = 0;
	precon = 0;
    } else {
	precon = 1;
    }
    if (cgitmx < 0) {
	cgitmx = 100;
    }
    if (etarg < 0. || etarg > 1.) {
	etarg = .5;
    }
/*     Check other options. */
    if (lvlscale < 0) {
	lvlscale = 0;
    }
    lvlscale = min(lvlscale,2);
/*     Partial pricing parameters */
    minprc = 10;
    maxmn = max(*m,*n);
/*     Temporarily go back to what we had before. */
/*     Set nParPrU to previous value. */
    if (nparpru <= 0) {
	if (qp) {
	    nparpru = 1;
	} else {
	    nparpru = 10;
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
/* $$$ */
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
    if (maxtime < 0.) {
	maxtime = 0.;
    }
    rmaxs = (doublereal) maxs;
/* Computing MAX */
    d__1 = 1. / (eps * 100. * rmaxs);
    chzbnd = max(d__1,1e6);
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
    if (xdlim <= 0.) {
	xdlim = 2.;
    }
    if (zcndbd <= 0.) {
	if (qpslvr == 0) {
	    zcndbd = 1e4;
	} else {
	    zcndbd = 1e6;
	}
    }
    if (tcrash < 0. || tcrash >= 1.) {
	tcrash = .1;
    }
/*     --------------------------------------- */
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
/* QP */
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
/*     Set some tolerances. */
/*     Set the optimality tolerance. */
/*     Solve the QP subproblems fairly accurately. */
    if (tolcg <= 0.) {
	tolcg = .01;
    }
    if (toloptqp <= 0.) {
	toloptqp = c6;
    }
    if (toloptfp < 0.) {
	toloptfp = c6;
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
    if (wtinf0 < 0.) {
	wtinf0 = 1.;
    }
    if (iback == inewb) {
	iback = 0;
    }
    if (itnlim < 0) {
/* Computing MAX */
	i__1 = 10000, i__2 = max(*n,*m) * 10;
	itnlim = max(i__1,i__2);
    }
/*     Load tolerances used to mark variables during printing in s4SavB. */
    tolnlp = toloptqp;
    tolcon = tolx;
/*     ------------------------------------------------------------------ */
/*     Re-assign the options to their respective work arrays. */
/*     ------------------------------------------------------------------ */
    rw[51] = toloptfp;
    rw[52] = toloptqp;
    rw[53] = tolnlp;
    rw[54] = tolcg;
    rw[56] = tolx;
    rw[57] = tolcon;
    rw[60] = tolpiv;
    rw[61] = tolrow;
    rw[62] = tcrash;
    rw[65] = tolswp;
    rw[66] = tolfac;
    rw[67] = tolupd;
    rw[70] = infbnd;
    rw[71] = bigfx;
    rw[72] = bigdx;
    rw[79] = maxtime;
    rw[80] = xdlim;
    rw[83] = etarg;
    rw[85] = hcondbnd;
    rw[86] = zcndbd;
    rw[88] = wtinf0;
    rw[92] = scltol;
    rw[151] = lmax1;
/* max L-multiplier in factor */
    rw[152] = lmax2;
/* max L-multiplier in update */
    rw[153] = small;
/* defn of small real */
    rw[154] = utol1;
/* abs tol for small diag of U */
    rw[155] = utol2;
/* rel tol for small diag of U */
    rw[156] = uspace;
/* limit on waste space in U */
    rw[157] = dens1;
/* switch to search maxcol columns and no rows */
    rw[158] = dens2;
/* switch to dense LU */
    rw[181] = toldpp;
    rw[182] = toldcp;
    rw[183] = toldup;
    rw[186] = toldj3;
    rw[187] = toldrp;
    iw[52] = maxr;
    iw[53] = maxs;
    iw[55] = qpslvr;
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
/*     iw( 69) = lvlStart */
    iw[71] = lvlsys;
    iw[73] = lvlobje;
    iw[75] = lvlscale;
    iw[77] = lvlpre;
    iw[80] = lvlpiv;
    iw[81] = lprprm;
    iw[83] = lprscl;
    iw[84] = lprsol;
    iw[85] = lprdbg;
    iw[87] = minmax;
    iw[88] = icrash;
    iw[89] = itnlim;
    iw[91] = mminor;
    iw[92] = mjrprint;
    iw[93] = mnrprint;
    iw[95] = mnewsb;
    iw[97] = cgitmx;
    iw[116] = stickyop;
    iw[120] = iback;
    iw[121] = idump;
    iw[122] = iloadb;
    iw[123] = imps;
    iw[124] = inewb;
    iw[125] = iinsrt;
    iw[126] = ioldb;
    iw[127] = ipnch;
    iw[130] = ireprt;
    iw[131] = isoln;
    iw[151] = nout;
    iw[153] = maxcol;
    iw[156] = tpivot;
    iw[199] = minimize;
    iw[209] = precon;
/* not optional parameters, but set here. */
    iw[99] = nparprlp;
    iw[100] = nparprqp;
    rw[186] = toldj3;
    return 0;
} /* s5defaults_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5defaults */
/* Subroutine */ int s5map_(integer *m, integer *n, integer *nkx, integer *
	ngobj, integer *nnh, integer *lenr, integer *maxr, integer *maxs, 
	integer *nextcw, integer *nextiw, integer *nextrw, integer *iw, 
	integer *leniw)
{
    static integer lxscaled, lfeastype, nb, lr, ly, lr1, lr2, ls1, ls2, ls3, 
	    ly1, ly2, ly3, mbs, lrg, ldx, liy, lkx, lrg2, liy1, liy2, lgbs, 
	    lkbs, lhdx, lpbs, lgqp, ngqp, lxbs, lblbs, lbubs, lblqp, lbuqp, 
	    lqprhs, lscales, lblsave, lbusave, lestate;

/*     ================================================================== */
/*     s5Map   allocates all array storage for sqopt, */
/*     using the values: */
/*        m    , n    , neA */
/*        maxS                    Set in s5defaults. */
/*        ngObj, nnH              Set from the argument list. */
/*        lenR                    Set in the calling program. */

/*     15 Nov 1991: First version based on Minos 5.4 routine m2core. */
/*     12 Nov 1994: Converted to integer and real storage. */
/*     06 Aug 1996: First min sum version. */
/*     14 Jul 1997: Thread-safe version. */
/*     01 May 1998: First version called by sqMem. This simplified */
/*                  version may slightly overestimate needed memory. */
/*     02 Aug 2003: snPRNT adopted. */
/*     13 May 2005: Bug fix: ly3 assigned to iw correctly */
/*     18 Jun 2008: Added space for iy2, pBS and rg2. */
/*     20 Sep 2014: Added space for bl and bu. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --iw;

    /* Function Body */
    ngqp = max(*ngobj,*nnh);
    mbs = *m + *maxs;
    nb = *n + *m;
/*     sqopt can use all of cw, iw and rw */
/*     except the first user workspace partitions. */
    lkx = *nextiw;
    lfeastype = lkx + *nkx;
    lkbs = lfeastype + mbs;
    lestate = lkbs + mbs;
    liy = lestate + nb;
    liy1 = liy + nb;
    liy2 = liy1 + nb;
    *nextiw = liy2 + nb;
/*     Addresses for the double precision arrays. */
    lscales = *nextrw;
    ly = lscales + nb;
    ly1 = ly + nb;
    ly2 = ly1 + nb;
    if (*maxr < *maxs) {
/* Define SYMMLQ workspace */
	ly3 = ly2 + nb;
	ls1 = ly3 + nb;
	ls2 = ls1 + *maxs;
	ls3 = ls2 + *maxs;
	lr1 = ls3 + *maxs;
	lr2 = lr1 + *maxs;
	lblqp = lr2 + *maxs;
    } else {
	ly3 = ly2 + nb;
	ls1 = ly3;
	ls2 = ls1;
	ls3 = ls2;
	lr1 = ls3;
	lr2 = lr1;
	lblqp = lr2;
    }
    lbuqp = lblqp + nb;
    lblbs = lbuqp + nb;
    lbubs = lblbs + mbs;
    lxbs = lbubs + mbs;
    lxscaled = lxbs + mbs;
    lhdx = lxscaled + *nnh;
    lpbs = lhdx + *nnh;
    lgqp = lpbs + nb;
    lgbs = lgqp + ngqp;
    lr = lgbs + mbs;
    lrg = lr + *lenr;
    lrg2 = lrg + *maxs;
    lblsave = lrg2 + *maxs;
    lbusave = lblsave + nb;
    lqprhs = lbusave + nb;
    ldx = lqprhs + *m;
    *nextrw = ldx + ngqp;
/*     --------------------------- */
/*     Store the addresses in iw. */
/*     --------------------------- */
    iw[251] = lkx;
    iw[271] = lblqp;
    iw[272] = lbuqp;
    iw[273] = lblbs;
    iw[274] = lbubs;
    iw[275] = lblsave;
    iw[276] = lbusave;
    iw[277] = lpbs;
    iw[278] = lqprhs;
    iw[284] = lfeastype;
    iw[285] = lestate;
    iw[287] = ldx;
    iw[288] = lhdx;
    iw[290] = lgqp;
    iw[291] = lgbs;
    iw[292] = lkbs;
    iw[293] = lrg;
    iw[294] = lrg2;
    iw[295] = lr;
    iw[296] = lscales;
    iw[301] = lxbs;
    iw[302] = lxscaled;
    iw[308] = liy;
    iw[309] = liy1;
    iw[310] = liy2;
    iw[311] = ly;
    iw[312] = ly1;
    iw[313] = ly2;
    iw[314] = ly3;
    iw[353] = lr1;
    iw[354] = lr2;
    iw[355] = ls1;
    iw[356] = ls2;
    iw[357] = ls3;
    return 0;
} /* s5map_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5Map */
/* Subroutine */ int s5callstatus_(integer *status, integer *iw, integer *
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
    static icilist io___261 = { 0, str, 0, fmt_9999, 80, 1 };


/* ================================================================= */
/* s5CallStatus fetches the call-status for the sqOpt user-defined */
/* matrix-vector product. */

/* 16 Jun 2008: First version of s5CallStatus. */
/* ================================================================= */
/* ----------------------------------------------------------------- */
/* ----------------------------------------------------------------- */
/* QP user-routine call-status */
    /* Parameter adjustments */
    --iw;

    /* Function Body */
    if (iw[235] == 0) {
/* Standard call */
	*status = 0;
    } else if (iw[235] < 0) {
/* First call */
	*status = 1;
	iw[235] = 0;
    } else if (iw[235] >= 2) {
/* Last orders please */
	*status = iw[235];
	iw[235] = -1;
    } else {
	*status = iw[235];
	s_wsfi(&io___261);
	do_fio(&c__1, (char *)&(*status), (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)80);
    }
    return 0;
} /* s5callstatus_ */

