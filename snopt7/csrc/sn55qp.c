/* ../snopt7/src/sn55qp.f -- translated by f2c (version 20100827).
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

static integer c__0 = 0;
static doublereal c_b5 = 0.;
static integer c__1 = 1;
static integer c__21 = 21;
static integer c__23 = 23;
static integer c__11 = 11;
static integer c__27 = 27;
static integer c__2 = 2;
static integer c_n2 = -2;
static doublereal c_b81 = -1.;
static integer c__31 = 31;
static doublereal c_b149 = 1.;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*     file  sn55qp.f */

/*     s5QP */
/*     s5checkp  s5chzq   s5getp   s5QPfg   s5QPitn   s5Rcheck   s5Rcol */
/*     s5rg      s5Rsng   s5Sdel   s5ZHZ    s5ZHZfac  s5Zp */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Subroutine */ int s5qp_(integer *iexit, integer *probtype, char *probtag, 
	integer *suboptimize, S_fp qplog, U_fp hprod, U_fp hprod1, integer *
	hvcalls, integer *eigh, logical *elastic, logical *gotr, logical *
	needlu, integer *typelu, logical *needx, integer *lenr, integer *m, 
	integer *maxs, integer *mbs, integer *n, integer *nb, integer *ndegen,
	 integer *ngqp0, integer *ngqp, integer *ngobj0, integer *ngobj, 
	integer *nnh0, integer *nnh, integer *ns, integer *itqp, integer *
	itqpmax, integer *itqptarget, integer *itn, integer *emode, integer *
	lvlobje, integer *printlevel, integer *minimize, integer *iobj, 
	doublereal *scaleobj, doublereal *objadd, doublereal *objqp, 
	doublereal *targeth, doublereal *targetz, doublereal *toloptfp, 
	doublereal *toloptqp, doublereal *tolx, integer *ninf, doublereal *
	sinf, integer *elastics, integer *ninfe, doublereal *sinfe, 
	doublereal *wtinf, doublereal *pinorm, doublereal *rgnorm, integer *
	nea, integer *nloca, integer *loca, integer *inda, doublereal *acol, 
	integer *neh, integer *nloch, integer *loch, integer *indh, 
	doublereal *hcol, integer *etype, integer *estate, integer *feastype, 
	integer *hs, integer *kbs, doublereal *bl, doublereal *bu, doublereal 
	*blqp, doublereal *buqp, doublereal *blbs, doublereal *bubs, 
	doublereal *gbs, doublereal *gobj, doublereal *gqp, doublereal *hdx, 
	doublereal *pbs, doublereal *pi, doublereal *r__, doublereal *rc, 
	doublereal *rg, integer *nrhs0, integer *nrhs, doublereal *rhs, 
	doublereal *scales, integer *lenx0, integer *nx0, doublereal *x0, 
	doublereal *x, doublereal *xbs, doublereal *xfrozen, integer *iy, 
	integer *iy1, doublereal *y, doublereal *y1, doublereal *y2, char *cu,
	 integer *lencu, integer *iu, integer *leniu, doublereal *ru, integer 
	*lenru, char *cw, integer *lencw, integer *iw, integer *leniw, 
	doublereal *rw, integer *lenrw, ftnlen probtag_len, ftnlen cu_len, 
	ftnlen cw_len)
{
    /* Format strings */
    static char fmt_1030[] = "(\002 Itn\002,i7,\002: Infeasible nonelastics."
	    "  Num =\002,i5,1p,\002  Sum of Infeasibilities =\002,e8.1)";
    static char fmt_1610[] = "(\002 Itn\002,i7,\002: Suboptimize: \002,i7"
	    ",\002 new superbasics\002)";
    static char fmt_1620[] = "(\002 Itn\002,i7,\002: Suboptimize: \002,i7"
	    ",\002 minor iterations\002)";
    static char fmt_1040[] = "(\002 Itn\002,i7,\002: Infinite pi-vector\002)";
    static char fmt_1600[] = "(\002 Itn\002,i7,\002: Singularity after a "
	    "\002,\002bound swap.  Basis refactorized\002)";
    static char fmt_8050[] = "(\002 Itn\002,i7,\002: Infeasible \002,a)";
    static char fmt_8060[] = "(\002 Itn\002,i7,\002: Elastic Phase 1 -- maki"
	    "ng \002,\002non-elastic variables feasible\002)";
    static char fmt_1010[] = "(\002 Biggest dj =\002,1p,e11.3,\002 (variabl"
	    "e\002,i7,\002)\002,\002    norm rg =\002,e11.3,\002   norm pi "
	    "=\002,e11.3)";
    static char fmt_1020[] = "(\002 Norm rg =\002,1p,e11.3,\002   norm pi "
	    "=\002,e11.3)";
    static char fmt_1000[] = "(\002 ==> LU file has increased by a factor o"
	    "f\002,f6.1)";

    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    extern /* Subroutine */ int s2unpack_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, doublereal *), s5zhzfac_(integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *), s5ereset_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    static logical feasible;
    static doublereal objslack;
    static logical increase;
    static doublereal objprint;
    static logical printlog, printsum;
    static doublereal rowerror;
    static logical checkfeas, deadpoint;
    static doublereal c6;
    static logical firstfeas;
    static integer lurequest, jq, kp;
    static logical justphase1, stationary;
    static integer jbq, jbr;
    static doublereal djq;
    static integer nbs, jsq, jsr;
    static char str[132];
    static doublereal djq0, eps0, eps2;
    extern /* Subroutine */ int s5rg_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *), s5hs_(integer *, integer *, doublereal *
	    , doublereal *, integer *, doublereal *);
    static integer kfac, nfac, kchk;
    static doublereal bold;
    static logical newb;
    static integer klog;
    static logical goth;
    static integer kprc, ksav, maxr;
    static logical luok;
    static integer nfix[2];
    static doublereal step;
    static logical newx;
    extern /* Subroutine */ int s5inf_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal tolx0;
    static logical needf;
    extern /* Subroutine */ int dload_(integer *, doublereal *, doublereal *, 
	    integer *), s5zhz_(integer *, U_fp, U_fp, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, char *, integer *, integer *, integer *, 
	    doublereal *, integer *, char *, integer *, integer *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen);
    static logical needv;
    static doublereal drmin;
    static logical newsb;
    static doublereal condz, drmax;
    static integer lumax;
    static logical newlu;
    static integer ksumm, nsmax, nswap;
    static doublereal norma, normg, pivot, rgtol[2], normz;
    extern /* Subroutine */ int s2bfac_(integer *, integer *, logical *, 
	    logical *, logical *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, integer *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, integer *), dcopy_(integer *, 
	    doublereal *, integer *, doublereal *, integer *), s5einf_(
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *, doublereal *, doublereal *, doublereal *, doublereal *), 
	    s6rcnd_(integer *, integer *, integer *, doublereal *, doublereal 
	    *, doublereal *, doublereal *), s1time_(integer *, integer *, 
	    integer *, integer *, doublereal *, integer *), s2bsol_(integer *,
	     integer *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, doublereal *, integer *), s5qpfg_(U_fp, U_fp, integer *
	    , integer *, integer *, integer *, integer *, integer *, integer *
	    , integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, char *, integer *, integer *, integer 
	    *, doublereal *, integer *, char *, integer *, integer *, integer 
	    *, doublereal *, integer *, ftnlen, ftnlen);
    static integer lusiz0;
    extern /* Subroutine */ int s5setx_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *, integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *);
    static integer kdegen;
    static doublereal infbnd;
    static logical needlm, checkx, needpi;
    static doublereal featol;
    static logical maxref, qpdone;
    static doublereal weight;
    extern doublereal dnormi_(integer *, doublereal *, integer *);
    static integer jqsave, inform__, itnlim, eigzhz;
    static logical gotgqp;
    static integer itnfix, mnewsb;
    static logical usegqp;
    static integer nfmove, frozen, nuncon;
    static doublereal djqprt;
    static integer lvltol, nonopt, sqstat;
    static doublereal zhzmin;
    static integer kprprt;
    static doublereal zhzmax, rgtest, normzq, tolinc;
    extern /* Subroutine */ int s5degen_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *), snprnt_(integer *, char *, integer *, 
	    integer *, ftnlen), s5egrad_(integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *), s5price_(logical *, logical *
	    , logical *, logical *, integer *, integer *, integer *, integer *
	    , integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, doublereal *, integer *), s4ksave_(integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, doublereal *, doublereal *, integer *, integer *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, char *, 
	    integer *, integer *, integer *, ftnlen), s5setpi_(integer *, 
	    integer *, logical *, doublereal *, doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, integer *), s5qpitn_(integer *
	    , U_fp, U_fp, integer *, integer *, integer *, logical *, logical 
	    *, logical *, logical *, logical *, logical *, logical *, logical 
	    *, logical *, logical *, logical *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, integer *, integer 
	    *, integer *, doublereal *, integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, char *, integer *, integer *, integer 
	    *, doublereal *, integer *, char *, integer *, integer *, integer 
	    *, doublereal *, integer *, ftnlen, ftnlen);
    static logical checkpi;
    extern /* Subroutine */ int s2trylu_(integer *, integer *, integer *, 
	    integer *, logical *, integer *, integer *, integer *, doublereal 
	    *, integer *);
    static doublereal signobj;
    static logical bndswap;
    static doublereal maxtime;
    static logical optimal;
    static doublereal bgrowth, condzhz;
    static logical prtlvl10;
    static integer rankzhz;
    extern /* Subroutine */ int s2gather_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *);

    /* Fortran I/O blocks */
    static icilist io___77 = { 0, str, 0, fmt_1030, 132, 1 };
    static icilist io___83 = { 0, str, 0, fmt_1610, 132, 1 };
    static icilist io___84 = { 0, str, 0, fmt_1620, 132, 1 };
    static icilist io___85 = { 0, str, 0, fmt_1040, 132, 1 };
    static icilist io___91 = { 0, str, 0, fmt_1600, 132, 1 };
    static icilist io___98 = { 0, str, 0, fmt_8050, 132, 1 };
    static icilist io___99 = { 0, str, 0, fmt_8060, 132, 1 };
    static icilist io___100 = { 0, str, 0, fmt_1010, 132, 1 };
    static icilist io___101 = { 0, str, 0, fmt_1020, 132, 1 };
    static icilist io___106 = { 0, str, 0, fmt_1000, 132, 1 };


/* ================================================================= */
/* s5QP   solves a linear or quadratic program. */

/* The optimization can pass through the following phases: */

/*   Phase 1               find a feasible point for all variables */

/*   Elastic Phase 1       make the non-elastic variables feasible */
/*                         while allowing infeasible elastics */

/*   Phase 2               minimize the objective */

/*   Elastic Phase 2       minimize a composite objective while */
/*                         keeping the non-elastics feasible */

/*                         In this phase, lvlObjE means: */

/*            lvlObjE = 0  zero     weight on the infeasibilities */
/*                                  (infeasibillities are ignored) */
/*                      1  finite   weight on the infeasibilities */
/*                      2  infinite weight on the infeasibilities */
/*                                  (the objective is ignored) */

/* On entry: */
/* --------- */
/* elastics        is the number of elastic variables allowed to */
/*                 move outside their bounds. */
/* nInfE           is the number of free elastics. */
/* subOptimize =-1 means no suboptimization. */
/* subOptimize = 0 allows optimization if necessary. */
/*                 During suboptimization, a subset of the */
/*                 variables are frozen at their initial values. */
/*                 Suboptimization is triggered if: */
/*                 (1) the number of new superbasics exceeds a */
/*                     preassigned limit. */
/*                 (2) the number of iterations exceeds itQPTarget. */

/* problemType     the type of QP being solved. */
/*                 Val */
/*                 --- */
/*                  0  FP   feasible point only */
/*                  1  LP   LP problem */
/*                  2  QP   QP problem */
/*                  3  FPE  feasible point for equalities only */
/*                  4  FPS  feasible point for QP subProblem */
/*                  5  QPS  QP subproblem */
/*                  6  QPP  FP subproblem with proximal point obj */

/* ngQP = max( nnH, ngObj ) */

/* The array kBS is a permutation on the column indices. */
/* kBS(1  :m )    holds the col. indices of the basic variables. */
/* kBS(m+1:m+nS)  holds the col. indices of the superbasic variable */
/*                These nS columns indices must have hs(j) = 2. */

/* On exit: */
/* --------- */
/* subOptimize > 0 means the QP was suboptimized. */
/*               1 new superbasic limit exceeded. */
/*               2 itQP > itQPtarget */

/*  iExit       Result */
/*  -----       ------ */
/*   >0         Fatal LU error */
/*    0         QP solution found */
/*   -1         QP is infeasible */
/*   -2         QP is unbounded */
/*   -3         Too many iterations */
/*   -4         Weak QP minimizer */
/*   -5         Too many superbasics */
/*   -6         QP Hessian not positive semidefinite after pricing */
/*   -7         Z'g could not be made sufficiently small */
/*   -8         Ill-conditioned Z */
/*   -9         Z'HZ indefinite after s5ZHZfac (Cholesky on Z'HZ) */

/* 30 Sep 1991: First version of s5QP  based on Qpsol routine qpcor */
/* 29 Oct 1993: QP objective computed separately. */
/* 19 May 1995: Bordered Hessian updated. */
/* 30 Jul 1995: Border updates removed. */
/* 04 Jan 1996: Positive semi-definite H treated correctly. */
/* 20 Jul 1996: Slacks changed to be the row value. */
/* 09 Aug 1996: First Min Sum version. */
/* 15 Jul 1997: Thread-safe version. */
/* 02 Feb 1998: Piecewise linear line search added. */
/* 07 Nov 1998: Explicit Hessian option added. */
/* 24 Dec 1999: Sub-optimization option added. */
/* 25 Jul 2001: Exit on nS > maxR activated. */
/* 30 Jul 2003: Superbasic slacks allowed. */
/* 02 Aug 2003: snEXIT and snPRNT adopted. */
/* 24 Dec 2003: pi checked for NaN and Inf entries. */
/* 07 May 2006: s4ksave handles negative values of hs. */
/* 16 May 2006: Explicit target itQP added */
/* 26 May 2013: infBnd used to identify infinite bounds. */
/* 02 Nov 2014: Refactor B and ZHZ before iExit = -6. */
/* 03 Nov 2014: neH, indH, locH, Hcol added as arguments. */
/* 22 Nov 2014: iExit set correctly when eMode = 0. */
/* ================================================================= */
/*     ------------------------------------------------------------------ */
/* partial pricing in use */
/* partial pricing for LPs */
/* partial pricing for QPs */
/* size of L0 */
/* size of initial  U */
/* size of current  L */
/* size of current  U */
/* phase 1 dj tol for p.p. */
/* phase 2 dj tol for p.p. */
/* current optimality tol */
/* xBS(kObj) is the obj. slack */
/* itns since last factorize */
/* number of LU mods */
/* (on/off) log     status */
/* (on/off) summary status */
/* # lines in log     file */
/* # lines in summary file */
/* >0 => Minor head in log file */
/* >0 => Minor head for iSumm */
/* Solve time */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --r__;
    --pi;
    --rg;
    --xbs;
    --pbs;
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
    --xfrozen;
    --x;
    --scales;
    --rc;
    --buqp;
    --blqp;
    --bu;
    --bl;
    --hs;
    --estate;
    --etype;
    --gqp;
    --hdx;
    --acol;
    --inda;
    --loca;
    --hcol;
    --indh;
    --loch;
    --gobj;
    --rhs;
    --x0;
    cu -= 8;
    --iu;
    --ru;
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    eps0 = rw[2];
/* eps**(4/5)       IEEE DP  3.00e-13 */
    eps2 = rw[4];
/* eps**(1/2)       IEEE DP  1.49e-08 */
    infbnd = rw[70];
/* definition of an infinite bound */
    maxtime = rw[79];
/* max time allowed */
    maxr = iw[52];
/* max columns of R. */
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
    itnlim = iw[89];
/* limit on total iterations */
    mnewsb = iw[95];
/* max # of new superbasics */
    nfac = iw[210];
/* # of LU factorizations */
    *iexit = 0;
    sqstat = 0;
    c6 = max(1e-6,eps2);
    if (nfac > 0) {
	lusiz0 = iw[171] + iw[172];
	lumax = lusiz0 << 1;
    }
    prtlvl10 = *printlevel >= 10;
    printlog = *printlevel >= 1 && (*itqp % klog == 0 && *itqp != 0 || klog ==
	     1);
    printsum = *printlevel >= 1 && (*itqp % ksumm == 0 && *itqp != 0 || ksumm 
	    == 1);
    iw[218] = 0;
    iw[219] = 0;
    if (printlog) {
	iw[218] = 1;
    }
    if (printsum) {
	iw[219] = 1;
    }
/*     ------------------------------------------------------------------ */
/*     s5QP operates in either ``Normal'' or ``Elastic'' mode. */
/*     Everything is normal unless a weighted sum is being minimized or */
/*     the constraints are infeasible. */
/*     The logical Feasible refers to the feasibility of the nonelastics. */
/*     wtInf  is the optional parameter Infeasibility Weight. */
/*     ------------------------------------------------------------------ */
    feasible = FALSE_;
    goth = *nnh > 0;
    gotgqp = *ngqp > 0;
/*     JustPhase1 = stop at the end of phase1 (either regular or elastic) */
    justphase1 = *probtype == 0 || *probtype == 3 || *probtype == 4;
/*     The phase 2 objective is F1 + wtInf*F2. */
    if (*elastic) {
	needf = *lvlobje != 2;
/* F1 required in phase 2 */
	needv = *lvlobje != 0;
/* F2 required in phase 2 */
    } else {
	needf = TRUE_;
	needv = FALSE_;
    }
    *objqp = 0.;
    objslack = 0.;
    pivot = 0.;
    step = 0.;
    *sinfe = 0.;
    *ninfe = 0;
    nonopt = -1;
    norma = 1.;
    normz = 1.;
    condz = 1.;
    jq = 0;
    djq = 0.;
    djq0 = 0.;
    djqprt = 0.;
    jbq = 0;
/* x(jBq) is the incoming  BS */
    jbr = 0;
/* x(jBr) is the outgoing  BS */
    jsq = 0;
/* x(jSq) is the incoming SBS */
    jsr = 0;
/* x(jSr) is the outgoing SBS */
    jqsave = 0;
    kprprt = 0;
    signobj = (doublereal) (*minimize);
    rgtol[0] = min(*toloptqp,c6);
/* relaxed   rgTol */
    rgtol[1] = eps0;
/* stringent rgTol */
    lvltol = 2;
/* working   rgTol */
    rw[184] = 100.;
    rw[185] = 100.;
    if (justphase1) {
	iw[94] = iw[99];
    } else {
	iw[94] = iw[100];
    }
    kprc = 0;
/* last sec scanned in part. prc */
    lurequest = 0;
    eigzhz = 0;
    bndswap = FALSE_;
    checkfeas = TRUE_;
/* Check that x is feasible. */
    checkpi = TRUE_;
/* Get norm pi in next call of s5setpi */
    deadpoint = FALSE_;
    needpi = TRUE_;
    newlu = TRUE_;
    newx = FALSE_;
    qpdone = FALSE_;
/*     nUncon  counts the number of unconstrained (i.e., Newton) steps. */
/*             If the test for a minimizer were scale-independent, */
/*             Uncon would never be larger than 1. */
/*     nfmove  counts the number of times that the QP obj is decreased, */
    nfmove = 0;
    nuncon = 0;
/*     subOptimize ne 0 forces optimization with a subset of the */
/*     variables frozen at their initial values. */
    frozen = 0;
    nsmax = *ns + mnewsb;
    s5hs_(&c__0, nb, &blqp[1], &buqp[1], &hs[1], &x[1]);
    s5degen_(&inform__, &c__0, printlevel, nb, ninf, itn, &featol, tolx, &
	    tolinc, &hs[1], &blqp[1], &buqp[1], &x[1], &itnfix, nfix, &tolx0, 
	    &iw[1], leniw, &rw[1], lenrw);
/*     ======================Start of main loop========================== */
/* +    do while (.not. QPdone  .and.  iExit .eq. 0) */
L100:
    if (! qpdone && *iexit == 0) {
/* ============================================================== */
/* Check the initial  x  and move it onto  ( A  -I )*x = b. */
/* If NeedLU is true, this will require a basis factorization. */
/* ============================================================== */
/* If necessary,  factorize the basis  ( B = LU ) and set x. */
/* If NeedLU is false on entry to s5QP, the first call to s2Bfac */
/* will try to use existing factors. */
/* If NeedLU is true on entry to s5QP, an LU factorization of */
/* type typeLU is computed. */

/* The reason for the requested LU is as follows. */

/* LUrequest =  0  First LU for a given subproblem */
/* LUrequest =  1  Frequency */
/* LUrequest =  2  LU nonzeros increased */
/* LUrequest =  3 */
/* LUrequest =  4 */
/* LUrequest =  5  Singular after LU mod */
/* LUrequest =  6  Unstable LU mod (growth in new column of U) */
/* LUrequest =  7  Not enough memory */
/* LUrequest =  8 */
/* LUrequest =  9 */
/* LUrequest = 10  Row error in setx */
/* LUrequest = 11  Big  dx   in setx */

/* LUrequest = 20 */
/* LUrequest = 21  Iterative refinement failed in QP */
/* LUrequest = 22  Unbounded QP */
/* LUrequest = 23  Infeasibility after refactorization */
/* LUrequest = 24  Small directional derivative in QP */
/* LUrequest = 25  Ill-conditioned Z in QP */
/* LUrequest = 26  Indefinite Z'HZ in QP */
/* LUrequest = 27  R singular after bound swap in QP */
/* -------------------------------------------------------------- */
	firstfeas = FALSE_;
	if (lurequest > 0) {
	    *needlu = TRUE_;
	}
	if (*needx || *needlu) {
	    s2bfac_(iexit, typelu, needlu, &newlu, &newb, iobj, itn, 
		    printlevel, &lurequest, m, mbs, n, nb, nnh, ns, &nswap, 
		    nea, nloca, &loca[1], &inda[1], &acol[1], &kbs[1], &hs[1],
		     &blqp[1], &buqp[1], &blbs[1], &bubs[1], nrhs0, nrhs, &
		    rhs[1], &x[1], &xbs[1], &iy[1], &iy1[1], &y[1], &y1[1], &
		    iw[1], leniw, &rw[1], lenrw);
	    if (newlu) {
		lusiz0 = iw[171] + iw[172];
		lumax = lusiz0 << 1;
		*gotr = FALSE_;
/* Reset R. */
		if (prtlvl10) {
		    iw[223] = 1;
		}
	    }
	    if (*iexit != 0) {
		goto L100;
	    }
	    needpi = TRUE_;
/* Recalculate the pi's. */
	    *needlu = FALSE_;
	    *needx = FALSE_;
	    newx = TRUE_;
	    checkfeas = TRUE_;
	    checkpi = TRUE_;
/* Check for NaNs when getting piNorm */
	    pivot = 0.;
	    jqsave = 0;
	    nuncon = 0;
	}
	newsb = FALSE_;
	optimal = FALSE_;
	nbs = *m + *ns;
	*ninf = 0;
	*sinf = 0.;
	dload_(&nbs, &c_b5, &gbs[1], &c__1);
	normg = 1.;
	if (checkfeas) {
/*           In Phase 1 or just after a factorize, check the feasibility */
/*           of the basic and superbasic non-elastics. */
	    if (*elastic) {
/*              Check that the elastic variables satisfy blQP and buQP. */
		s5ereset_(&nbs, nb, elastics, &featol, &infbnd, &etype[1], &
			estate[1], &kbs[1], &bl[1], &bu[1], &blqp[1], &buqp[1]
			, &blbs[1], &bubs[1], &x[1]);
	    }
/*           FirstFeas  indicates that we have just become feasible. */
/*           FirstFeas  is turned off once a step is taken. */
	    s5inf_(&nbs, &featol, &infbnd, ninf, sinf, &feastype[1], &blbs[1],
		     &bubs[1], &gbs[1], &xbs[1]);
	    if (*ninf > 0) {
/*              Non-elastics are infeasible. */
/*              If necessary, switch back to the feasibility phase, after */
/*              refactorization (possibly with tighter tols). */
/*              Print something if the basis has just been refactorized. */
		if (prtlvl10 && iw[215] == 0) {
		    s_wsfi(&io___77);
		    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&(*ninf), (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&(*sinf), (ftnlen)sizeof(doublereal)
			    );
		    e_wsfi();
		    snprnt_(&c__21, str, &iw[1], leniw, (ftnlen)132);
		}
		if (feasible) {
		    s2trylu_(itn, &c__23, ns, &lurequest, &luok, typelu, &iw[
			    1], leniw, &rw[1], lenrw);
		    if (! luok) {
			*iexit = 11;
		    }
		    feasible = FALSE_;
		    goto L100;
		}
		*gotr = FALSE_;
/* Is this needed? */
	    }
/*           Feasible => the nonelastics are feasible. */
/*           Feasible => normal or elastic Phase 2 */
	    if (! feasible) {
		firstfeas = *ninf == 0;
/* Feasible for the first time */
	    }
	    feasible = *ninf == 0;
	    checkfeas = *ninf > 0;
	}
/* if CheckFeas */
	if (*elastic) {
/* ----------------------------------------------------------- */
/* Compute the sum of infeasibilities of the elastics. */
/* If  nInfE .ne. elastics, then there is at least one */
/* nonbasic elastic fixed at its current value. */
/* ----------------------------------------------------------- */
	    s5einf_(nb, &nbs, &estate[1], &kbs[1], &featol, ninfe, sinfe, &bl[
		    1], &bu[1], &x[1]);
	}
	if (feasible && justphase1) {
/* The non-elastics are feasible, prepare to exit. */
	    condzhz = 0.;
	    djqprt = 0.;
	    *rgnorm = 0.;
	    dload_(m, &c_b5, &pi[1], &c__1);
	    *pinorm = 1.;
/* pinorm = max(norm(pi), 1.0) */
	} else {
	    if (feasible) {
/*              --------------------------------------------------------- */
/*              Feasible for the nonelastics. */
/*              (Elastic = false means no elastics.) */
/*              --------------------------------------------------------- */
/*              If just feasible, compute the QP objective (and gradient) */
/*              and R. */
		if (firstfeas || newx) {
		    if (needf) {
/*                    =================================================== */
/*                    Initialize the QP objective and gradient. */
/*                    objQP is the explicit linear plus quadratic obj. */
/*                    objQP is not scaled by signObj. */
/*                    objQP is updated after each QP step. */
/*                    =================================================== */
			if (gotgqp) {
			    if (*hvcalls == 0) {
				sqstat = 1;
			    }
			    s5qpfg_((U_fp)hprod, (U_fp)hprod1, hvcalls, ngqp, 
				    ngobj0, ngobj, nnh, neh, nloch, &loch[1], 
				    &indh[1], &hcol[1], &sqstat, objqp, &gobj[
				    1], &gqp[1], lenx0, nx0, &x0[1], &x[1], &
				    y[1], cu + 8, lencu, &iu[1], leniu, &ru[1]
				    , lenru, cw + 8, lencw, &iw[1], leniw, &
				    rw[1], lenrw, (ftnlen)8, (ftnlen)8);
			    sqstat = 0;
			}
			if (goth && ! (*gotr)) {
/*                       ------------------------------------------------ */
/*                       Load and factor the reduced Hessian. */
/*                       This happens after every LU factorize. */
/*                       If the reduced Hessian is not positive definite, */
/*                       reduce the LU factor tolerances to get a better */
/*                       conditioned Z. */
/*                       ------------------------------------------------ */
			    if (*ns > 0) {
				s5zhz_(&inform__, (U_fp)hprod, (U_fp)hprod1, 
					hvcalls, &maxr, lenr, minimize, m, 
					mbs, n, nb, nnh, ns, nea, nloca, &
					loca[1], &inda[1], &acol[1], neh, 
					nloch, &loch[1], &indh[1], &hcol[1], &
					zhzmin, &zhzmax, &condz, &normz, &kbs[
					1], &r__[1], &y[1], &y1[1], &y2[1], 
					cu + 8, lencu, &iu[1], leniu, &ru[1], 
					lenru, cw + 8, lencw, &iw[1], leniw, &
					rw[1], lenrw, (ftnlen)8, (ftnlen)8);
/*                          inform = -1, 0, >0 */
				if (inform__ != 0) {
				    if (inform__ == -1) {
					*iexit = -8;
/* Ill-conditioned Z */
				    } else {
					*iexit = inform__;
/* Fatal error in LU */
				    }
				    goto L100;
				}
				s5zhzfac_(&inform__, eigh, &c__0, itn, lenr, 
					m, &maxr, mbs, nb, ns, targeth, &
					zhzmin, &zhzmax, &normz, &rankzhz, &
					hs[1], &kbs[1], &iy[1], &blqp[1], &
					buqp[1], &blbs[1], &bubs[1], &x[1], &
					xbs[1], &r__[1], &y1[1], &iw[1], 
					leniw);
/*                          inform = -2, -1, 0 */
				if (inform__ != 0) {
				    *iexit = -9;
/* Z'HZ not positive definite */
				    goto L100;
				}
			    }
/* nS > 0 */
			    *gotr = TRUE_;
			}
/* GotH and not GotR */
		    }
/* Needf */
		    newx = FALSE_;
/* objQP and gQP are now updated */
		    nbs = *m + *ns;
		    if (*gotr || *ns == 0) {
			eigzhz = 1;
		    } else {
			eigzhz = 0;
		    }
		}
/*              --------------------------------------------------------- */
/*              Gather the QP gradient in BS order. */
/*              Assign the nonzero components of gBS. */
/*              --------------------------------------------------------- */
/* FirstFeas .or. Newx */
		if (needf) {
		    if (gotgqp) {
			s2gather_(ngqp, &nbs, &kbs[1], &signobj, &gqp[1], &
				gbs[1]);
		    }
		    if (*iobj > 0) {
			gbs[iw[205]] = signobj * *scaleobj;
		    }
		}
		if (*elastic && *ninfe > 0 && needv) {
		    s5egrad_(nb, &nbs, wtinf, &estate[1], &kbs[1], &gbs[1]);
		}
		normg = dnormi_(&nbs, &gbs[1], &c__1);
/*              --------------------------------------------------------- */
/*              See if it's time to suboptimize. */
/*              NOTE: We must not suboptimize if all steps have been */
/*              degenerate. */
/*              --------------------------------------------------------- */
		if (*suboptimize != 0 || nfmove == 0) {
/*                 Relax */
		} else {
		    if (*ns >= nsmax) {
			*suboptimize = 1;
			if (prtlvl10) {
			    s_wsfi(&io___83);
			    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(
				    integer));
			    do_fio(&c__1, (char *)&mnewsb, (ftnlen)sizeof(
				    integer));
			    e_wsfi();
			    snprnt_(&c__21, str, &iw[1], leniw, (ftnlen)132);
			}
		    } else if (*itqp >= *itqptarget) {
			*suboptimize = 2;
			if (prtlvl10) {
			    s_wsfi(&io___84);
			    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(
				    integer));
			    do_fio(&c__1, (char *)&(*itqptarget), (ftnlen)
				    sizeof(integer));
			    e_wsfi();
			    snprnt_(&c__21, str, &iw[1], leniw, (ftnlen)132);
			}
		    }
		}
	    }
/* Feasible */
	    if (needpi) {
		dcopy_(m, &gbs[1], &c__1, &y[1], &c__1);
		s5setpi_(&inform__, m, &checkpi, pinorm, &y[1], &pi[1], &iw[1]
			, leniw, &rw[1], lenrw);
		checkpi = FALSE_;
		if (inform__ != 0) {
		    if (inform__ > 0) {
			*iexit = inform__;
		    } else {
/* pi is infinite or contains a NaN/Inf. */
			s_wsfi(&io___85);
			do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer)
				);
			e_wsfi();
			snprnt_(&c__21, str, &iw[1], leniw, (ftnlen)132);
			s2trylu_(itn, &c__11, ns, &lurequest, &luok, typelu, &
				iw[1], leniw, &rw[1], lenrw);
			if (! luok) {
			    *iexit = 43;
			}
		    }
		    goto L100;
		}
		needpi = FALSE_;
	    }
	    *rgnorm = 0.;
	    if (*ns > 0) {
		s5rg_(m, &nbs, n, ns, &eps0, nea, nloca, &loca[1], &inda[1], &
			acol[1], &gbs[1], &pi[1], &rg[1], rgnorm, &kbs[1]);
	    }
/*           ============================================================ */
/*           Determine if the reduced Hessian is positive definite. */
/*           ============================================================ */
	    condzhz = 0.;
	    if (*gotr) {
		s6rcnd_(&maxr, ns, lenr, &r__[1], &drmax, &drmin, &condzhz);
	    }
/*           ============================================================ */
/*           Check for a stationary point.  Use a stringent rgTol after */
/*           a constrained step to help avoid false stationary points. */
/*           In theory, the reduced gradient is zero and the reduced */
/*           Hessian is positive definite after a bound swap. */

/*           If x is a minimizer,  reduced costs are calculated. */
/*           ============================================================ */
	    if (feasible) {
		rw[186] = *toloptqp;
	    } else {
		rw[186] = *toloptfp;
	    }
	    rgtest = max(*pinorm,normg);
	    if (! feasible) {
		stationary = *rgnorm <= rgtol[0] * rgtest;
	    } else if (nuncon >= 1) {
		stationary = *rgnorm <= rgtol[lvltol - 1] * rgtest;
	    } else {
		stationary = *rgnorm <= rgtol[1] * rgtest;
	    }
	    if (feasible) {
		maxref = nuncon > 1;
		if ((maxref || bndswap) && ! stationary) {
/*                 If this point should be stationary but isn't. */
/*                 If possible, relax the reduced-gradient tolerance. */
		    if (lvltol == 2) {
			lvltol = 1;
			stationary = *rgnorm <= rgtol[lvltol - 1] * rgtest;
		    }
		    if (! stationary) {
			s2trylu_(itn, &c__21, ns, &lurequest, &luok, typelu, &
				iw[1], leniw, &rw[1], lenrw);
			if (! luok) {
			    *iexit = -7;
/* Large Z'g */
			}
			goto L100;
		    }
		}
		deadpoint = stationary && needf && eigzhz == 0;
	    }
	    if (stationary || bndswap) {
		jqsave = 0;
		nuncon = 0;
		bndswap = FALSE_;
		if (*gotr && eigzhz == 0) {
		    s_wsfi(&io___91);
		    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		    e_wsfi();
		    snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)132);
		    s2trylu_(itn, &c__27, ns, &lurequest, &luok, typelu, &iw[
			    1], leniw, &rw[1], lenrw);
		    if (! luok) {
			*iexit = 44;
		    }
/* Ill-conditioned Z */
		    goto L100;
		}
	    }
	    needlm = stationary;
	    kprprt = kprc;
	    jq = 0;
	    if (needlm) {
/*              --------------------------------------------------------- */
/*              Compute Lagrange multipliers. */
/*              --------------------------------------------------------- */
		djq0 = djq;
/* save djq in case of bad statpt */
		djq = 0.;
		nuncon = 0;
		usegqp = feasible && needf && gotgqp;
		weight = 0.;
		if (*elastic && feasible) {
		    weight = *wtinf;
		}
		s5price_(elastic, &feasible, &increase, &usegqp, suboptimize, 
			itn, m, n, nb, ngqp0, ngqp, &frozen, &nonopt, &weight,
			 &signobj, pinorm, &jq, &djq, &kprc, &rw[184], nea, 
			nloca, &loca[1], &inda[1], &acol[1], &etype[1], &hs[1]
			, &gqp[1], &pi[1], &rc[1], &x[1], &xfrozen[1], &iw[1],
			 leniw, &rw[1], lenrw);
		optimal = nonopt == 0;
		newsb = nonopt > 0;
	    }
/* NeedLM */
	}
/* Feasible and JustPhase1 */
	qpdone = optimal || deadpoint || feasible && justphase1;
	if (qpdone) {
/*           ------------------------------------------------------------ */
/*           Apparently we are finished. */
/*           See if any nonbasics have to be set back on their bounds. */
/*           ------------------------------------------------------------ */
	    s5degen_(&inform__, &c__1, printlevel, nb, ninf, itn, &featol, 
		    tolx, &tolinc, &hs[1], &blqp[1], &buqp[1], &x[1], &itnfix,
		     nfix, &tolx0, &iw[1], leniw, &rw[1], lenrw);
	    qpdone = inform__ == 0;
	    if (qpdone) {
/* -------------------------------------------------------- */
/* So far so good.  Now check the row residuals. */
/* -------------------------------------------------------- */
		if (iw[215] > 0) {
		    s5setx_(&inform__, &c__1, itn, m, n, nb, &nbs, &rowerror, 
			    nea, nloca, &loca[1], &inda[1], &acol[1], &kbs[1],
			     &xbs[1], nrhs0, nrhs, &rhs[1], &x[1], &y[1], &y1[
			    1], &iw[1], leniw, &rw[1], lenrw);
		    qpdone = inform__ == 0;
		    lurequest = inform__;
		    if (lurequest > 0) {
			*typelu = 2;
		    }
		}
	    }
	    if (qpdone) {
		if (deadpoint) {
		    *iexit = -4;
		}
	    } else {
		*needx = TRUE_;
		needpi = TRUE_;
		goto L100;
	    }
	}
/* done */
	if (firstfeas && justphase1) {
/*           Relax, we are about to exit without printing anything. */
	} else {
/*           ============================================================ */
/*           Print the details of this iteration. */
/*           ============================================================ */
	    objprint = 0.;
	    if (feasible) {
		if (needf) {
		    if (*iobj != 0) {
			objslack = xbs[iw[205]] * *scaleobj;
		    }
		    objprint = *objadd + objslack + *objqp;
		}
		if (needv) {
		    objprint += signobj * *wtinf * *sinfe;
		}
	    }
	    (*qplog)(probtype, probtag, elastic, gotr, &firstfeas, &feasible, 
		    m, mbs, nnh, ns, &jsq, &jbr, &jsr, &iw[220], &iw[221], 
		    itn, itqp, &kprprt, lvlobje, &pivot, &step, ninf, sinf, 
		    ninfe, sinfe, wtinf, &nonopt, &objprint, &condzhz, &
		    djqprt, rgnorm, &kbs[1], &xbs[1], &iw[1], leniw, (ftnlen)
		    20);
	}
	jbq = 0;
	jbr = 0;
	jsq = 0;
	jsr = 0;
	kprprt = 0;
	djqprt = 0.;
	if (qpdone) {
/*           ------------------------------------------------------------ */
/*           Optimal .or. Deadpoint  .or. (Feasible .and. JustPhase1) */
/*           ------------------------------------------------------------ */
	    if (*ninf > 0) {
/* No feasible point for the nonelastic variables */
/* Stop or continue in elastic mode, depending on the */
/* specified level of infeasibility. */
		if (*emode == 0 || *elastic) {
/*                 ------------------------------------------------------ */
/*                 The nonelastic constraints cannot be satisfied. */
/*                 ------------------------------------------------------ */
		    *iexit = -1;
/* Infeasible nonelastics */
		} else {
/*                 ------------------------------------------------------ */
/*                 Infeasible nonelastics in Normal mode. */
/*                 Print a message and start elastic Phase 1. */
/*                 ------------------------------------------------------ */
		    if (prtlvl10) {
			s_wsfi(&io___98);
			do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer)
				);
			do_fio(&c__1, probtag, (ftnlen)20);
			e_wsfi();
			snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)132);
			s_wsfi(&io___99);
			do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer)
				);
			e_wsfi();
			snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)132);
			iw[223] = 1;
			iw[225] = 1;
		    }
		    *elastic = TRUE_;
		    checkfeas = TRUE_;
/* call s5eReset */
		    qpdone = FALSE_;
		    needf = *lvlobje != 2;
/* Need F1 in e-phase 2 */
		    needv = *lvlobje != 0;
/* Need F2 in e-phase 2 */
		    needpi = TRUE_;
		    djq = 0.;
		    step = 0.;
		}
		goto L100;
	    }
	    if (prtlvl10 && ! justphase1) {
		if (jq != 0) {
		    djqprt = signobj * djq;
		    if (prtlvl10 && klog == 1) {
			s_wsfi(&io___100);
			do_fio(&c__1, (char *)&djq, (ftnlen)sizeof(doublereal)
				);
			do_fio(&c__1, (char *)&jq, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&(*rgnorm), (ftnlen)sizeof(
				doublereal));
			do_fio(&c__1, (char *)&(*pinorm), (ftnlen)sizeof(
				doublereal));
			e_wsfi();
			snprnt_(&c__11, str, &iw[1], leniw, (ftnlen)132);
		    }
		} else {
		    if (prtlvl10 && klog == 1) {
			s_wsfi(&io___101);
			do_fio(&c__1, (char *)&(*rgnorm), (ftnlen)sizeof(
				doublereal));
			do_fio(&c__1, (char *)&(*pinorm), (ftnlen)sizeof(
				doublereal));
			e_wsfi();
			snprnt_(&c__11, str, &iw[1], leniw, (ftnlen)132);
		    }
		}
	    }
	} else {
/* ---------------------------------------------------------- */
/* Do another QP iteration. */
/* A nonbasic has been selected to become superbasic. */
/* Compute the vector y such that B y = column jq. */
/* ---------------------------------------------------------- */
	    if (newsb) {
/* -------------------------------------------------------- */
/* The price has selected a nonbasic to become superbasic. */
/* -------------------------------------------------------- */
		if (*ns + 1 > maxr) {
		    *iexit = -5;
		    goto L100;
		}
		lvltol = 2;
		djqprt = djq;
/*              --------------------------------------------------------- */
/*              Compute the vector pBS such that B pB = column jq. */
/*              pBS is a multiple of part of the new column of  Z  and */
/*              is used to define the QP search direction and update R. */
/*              --------------------------------------------------------- */
/*              Unpack column jq into  y1  and solve  B*y = y1. */
/*              The solve computes  y1  such that  L*y1 = ajq. */
/*              It is used below to modify L and U in s5QPitn. */
		s2unpack_(&jq, m, n, nea, &norma, nloca, &loca[1], &inda[1], &
			acol[1], &y1[1]);
		s2bsol_(iexit, &c__1, m, &y1[1], &pbs[1], &iw[1], leniw, &rw[
			1], lenrw);
		if (*iexit != 0) {
		    return 0;
		}
		normzq = dnormi_(m, &y1[1], &c__1);
		normz = max(normzq,normz);
/* Computing MAX */
		d__1 = normzq / norma;
		condz = max(d__1,condz);
	    }
/*           ============================================================ */
/*           Take a step. */
/*           ============================================================ */
	    if (*itn >= itnlim || *itqp >= *itqpmax) {
		*iexit = -3;
		goto L100;
	    }
	    ++(*itqp);
	    ++(*itn);
/*           Decide if we want to print something this iteration. */
	    printlog = *printlevel >= 1 && *itqp % klog == 0;
	    printsum = *printlevel >= 1 && *itqp % ksumm == 0;
	    iw[218] = 0;
	    iw[219] = 0;
	    if (printlog) {
		iw[218] = 1;
	    }
	    if (printsum) {
		iw[219] = 1;
	    }
/*           ------------------------------------------------------------ */
/*           Take a reduced gradient step. */
/*           The new  x  will either minimize the objective on the */
/*           working set or lie on the boundary of a new constraint. */
/*           ------------------------------------------------------------ */
	    s5qpitn_(&inform__, (U_fp)hprod, (U_fp)hprod1, hvcalls, eigh, &
		    eigzhz, &bndswap, elastic, &feasible, &gotgqp, &goth, 
		    gotr, &increase, &needf, &needv, &needpi, &newsb, itn, 
		    lenr, m, mbs, &maxr, maxs, n, nb, nnh0, nnh, ns, ngqp0, 
		    ngqp, ndegen, elastics, &lurequest, &kp, &jbq, &jsq, &jbr,
		     &jsr, &jq, &jqsave, &nfmove, &nuncon, &djq0, &djq, 
		    minimize, objqp, &featol, &pivot, &step, &tolinc, wtinf, 
		    nea, nloca, &loca[1], &inda[1], &acol[1], neh, nloch, &
		    loch[1], &indh[1], &hcol[1], &etype[1], &estate[1], &
		    feastype[1], &hs[1], &kbs[1], &bl[1], &bu[1], &blqp[1], &
		    buqp[1], &blbs[1], &bubs[1], &gbs[1], &gqp[1], &hdx[1], &
		    pbs[1], &rg[1], &r__[1], &x[1], &xbs[1], &y[1], &y1[1], &
		    y2[1], cu + 8, lencu, &iu[1], leniu, &ru[1], lenru, cw + 
		    8, lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)8, (
		    ftnlen)8);
/*           Check for trouble in s5QPitn. */
/*           inform values are -2, -1, 0, >0 */
	    if (inform__ != 0) {
		if (inform__ > 0) {
		    *iexit = inform__;
/* Fatal LU error */
		} else if (inform__ == -1) {
		    *iexit = -2;
/* unbounded */
		} else if (inform__ == -2) {
/*                 Exit if B was just factored or Z is well-conditioned. */
/*                 Otherwise, refactor B, possibly with a BS factorize. */
/*                 The value of targetZ is increased to avoid multiple */
/*                 factorizations. */
		    if (iw[215] == 0) {
			*iexit = -6;
/* Z'HZ is indefinite */
		    }
		}
		if (*iexit != 0) {
		    goto L100;
		}
	    }
	    if (lurequest > 0) {
		s2trylu_(itn, &lurequest, ns, &lurequest, &luok, typelu, &iw[
			1], leniw, &rw[1], lenrw);
		if (! luok) {
		    *iexit = 43;
		    goto L100;
		}
	    }
	    ++iw[215];
/*           Increment featol every iteration. */
	    featol += tolinc;
/*           ============================================================ */
/*           Test for error condition and/or frequency interrupts. */
/*           ============================================================ */
/*           (1) Save a basis map (frequency controlled). */
/*           (2) Every kdegen iterations, reset featol and move nonbasic */
/*               variables onto their bounds if they are very close. */
/*           (3) Refactorize the basis if it has been modified too many */
/*               times. */
/*           (4) Update the LU factors of the basis if requested. */
/*           (5) Check row error (frequency controlled). */
	    if (*itn % ksav == 0) {
		s4ksave_(minimize, m, n, nb, ns, mbs, itn, ninf, sinf, objqp, 
			&kbs[1], &hs[1], &scales[1], &blqp[1], &buqp[1], &x[1]
			, &xbs[1], cw + 8, lencw, &iw[1], leniw, (ftnlen)8);
	    }
	    if (*itn % kdegen == 0) {
		s5degen_(&inform__, &c__2, printlevel, nb, ninf, itn, &featol,
			 tolx, &tolinc, &hs[1], &blqp[1], &buqp[1], &x[1], &
			itnfix, nfix, &tolx0, &iw[1], leniw, &rw[1], lenrw);
		*needx = inform__ > 0;
	    }
	    if (lurequest == 0) {
		if (iw[216] >= kfac - 1) {
		    lurequest = 1;
		} else if (iw[216] >= 20 && iw[173] + iw[174] > lumax) {
		    bgrowth = (doublereal) (iw[173] + iw[174]);
		    bold = (doublereal) lusiz0;
		    bgrowth /= bold;
		    if (prtlvl10) {
			s_wsfi(&io___106);
			do_fio(&c__1, (char *)&bgrowth, (ftnlen)sizeof(
				doublereal));
			e_wsfi();
			snprnt_(&c__21, str, &iw[1], leniw, (ftnlen)132);
		    }
		    lurequest = 2;
		} else {
		    checkx = iw[215] % kchk == 0;
		    if (checkx && ! (*needx)) {
			s5setx_(&inform__, &c__1, itn, m, n, nb, &nbs, &
				rowerror, nea, nloca, &loca[1], &inda[1], &
				acol[1], &kbs[1], &xbs[1], nrhs0, nrhs, &rhs[
				1], &x[1], &y[1], &y1[1], &iw[1], leniw, &rw[
				1], lenrw);
			lurequest = inform__;
		    }
		}
		if (lurequest > 0) {
		    *typelu = 3;
		}
	    }
/* Check time limit */
	    if (maxtime > 0. && *itn % 20 == 0) {
		s1time_(&c_n2, &c__0, &iw[1], leniw, &rw[1], lenrw);
		s1time_(&c__2, &c__0, &iw[1], leniw, &rw[1], lenrw);
		if (rw[462] > maxtime) {
		    *iexit = 34;
/* time limit reached */
		}
	    }
	}
/* not done */
	goto L100;
/* +    end while */
    }
/*     ======================end of main loop============================ */

    s5hs_(&c__1, nb, &blqp[1], &buqp[1], &hs[1], &x[1]);
    if (*suboptimize > 0) {
	if (frozen > 0) {
/*           Relax */
	} else {
	    *suboptimize = 0;
	}
    }
    return 0;
/* L1800: */
} /* s5qp_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5QP */
/* Subroutine */ int s5checkp_(integer *iexit, integer *itn, integer *nbs, 
	integer *jqsave, integer *kbs, doublereal *gtp, doublereal *pbs, 
	integer *iw, integer *leniw)
{
    /* Format strings */
    static char fmt_1000[] = "(\002 Itn\002,i7,\002: Bad directional derivat"
	    "ive \002,1p,e9.1)";
    static char fmt_9000[] = "(\002 XXX  s5checkp.  kSave not found. jqSave "
	    "= \002,i5)";

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static integer j, k, jq;
    static char str[80];
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer ksave;
    static doublereal psave;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___114 = { 0, str, 0, fmt_1000, 80, 1 };
    static icilist io___115 = { 0, str, 0, fmt_9000, 80, 1 };


/*     ================================================================== */
/*     s5checkp  makes  pBS  a feasible direction. */

/*     16 Jun 1995: First version of s5checkp. */
/*     02 Aug 2003: snPRNT adopted. */
/*     02 Aug 2003: Current version of s5checkp. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --pbs;
    --kbs;
    --iw;

    /* Function Body */
    *iexit = 0;
/*     ------------------------------------------------------------------ */
/*     Find the element of  pBS  corresponding to the most recently freed */
/*     variable. Usually, it will be pBS(nBS). */
/*     ------------------------------------------------------------------ */
    jq = abs(*jqsave);
    ksave = 0;
    for (k = *nbs; k >= 1; --k) {
	j = kbs[k];
	if (j == jq) {
	    ksave = k;
	    goto L100;
	}
    }
/*     ------------------------------------------------------------------ */
/*     Choose the sign of  pBS  so that the most recently freed */
/*     variable continues to increase or decrease. */
/*     ------------------------------------------------------------------ */
L100:
    if (ksave > 0) {
	psave = pbs[ksave];
	if (*jqsave < 0 && psave > 0. || *jqsave > 0 && psave < 0.) {
	    dscal_(nbs, &c_b81, &pbs[1], &c__1);
	    *gtp = -(*gtp);
	}
	if (*gtp > 0.) {
/*           ------------------------------------------------------------ */
/*           Looks as though the sign of gtp cannot be relied upon. */
/*           In later versions we'll fix this variable. */
/*           For now, we just print a warning and stop. */
/*           ------------------------------------------------------------ */
	    s_wsfi(&io___114);
	    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*gtp), (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
	    *iexit = 1;
/* Bad directional derivative */
	}
    } else {
/*        --------------------------------------------------------------- */
/*        Couldn't find the index of the most recently freed variable. */
/*        This should never happen! */
/*        --------------------------------------------------------------- */
	s_wsfi(&io___115);
	do_fio(&c__1, (char *)&(*jqsave), (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
    }
    return 0;
} /* s5checkp_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5checkp */
/* Subroutine */ int s5chzq_(integer *m, integer *mbs, integer *n, integer *
	nb, integer *ns, integer *kbsq, doublereal *pivot, doublereal *tolpiv,
	 integer *nea, integer *nloca, integer *loca, integer *inda, 
	doublereal *acol, integer *kbs, doublereal *bl, doublereal *bu, 
	doublereal *xbs, doublereal *y, integer *iw, integer *leniw, 
	doublereal *rw, integer *lenrw)
{
    /* Format strings */
    static char fmt_1000[] = "(\002 XXX  s5chzq.  Max pivot is too small:"
	    "\002,1p,e11.1)";

    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static integer j, k;
    static doublereal d1, d2;
    static integer m1;
    static doublereal xj, tol;
    static char str[80];
    static doublereal eps0, dpiv;
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen), s2bprod_(integer *, doublereal *, integer *, integer *
	    , integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, integer *);

    /* Fortran I/O blocks */
    static icilist io___119 = { 0, str, 0, fmt_1000, 80, 1 };


/*     ================================================================== */
/*     s5chzq  selects a superbasic to replace the kp-th basic variable. */
/*     On entry,  y  contains the kp-th row of B(inverse). */
/*     On exit, pivot and  y(m+1), ..., y(m+nS) define the S-part of */
/*     the modifying vector w. */

/*     01 Dec 1991: First version based on Minos routine m7chzq. */
/*     02 Aug 2003: snPRNT adopted. */
/*     30 Jun 2005: Current version of s5chzq. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --y;
    --xbs;
    --kbs;
    --bu;
    --bl;
    --acol;
    --inda;
    --loca;
    --iw;
    --rw;

    /* Function Body */
    eps0 = rw[2];
/*     Set yS = 0 -  S'*y. */
/* eps**(4/5) */
    m1 = *m + 1;
    s2bprod_(&c__1, &eps0, n, ns, &kbs[m1], nea, nloca, &loca[1], &inda[1], &
	    acol[1], &c_b81, &y[1], m, &c_b5, &y[m1], ns);
    *kbsq = *m + idamax_(ns, &y[m1], &c__1);
    *pivot = (d__1 = y[*kbsq], abs(d__1));
/*     Exit if the pivot is too small. */
    if (*pivot < *tolpiv) {
	s_wsfi(&io___119);
	do_fio(&c__1, (char *)&(*pivot), (ftnlen)sizeof(doublereal));
	e_wsfi();
	snprnt_(&c__31, str, &iw[1], leniw, (ftnlen)80);
	*kbsq = -(*m + *ns);
    } else {
/*        Choose one away from its bounds if possible. */
	tol = *pivot * .1;
	dpiv = -1.;
	i__1 = *m + *ns;
	for (k = m1; k <= i__1; ++k) {
	    if ((d__1 = y[k], abs(d__1)) >= tol) {
		j = kbs[k];
		xj = xbs[k];
		d1 = xj - bl[j];
		d2 = bu[j] - xj;
/* Computing MIN */
		d__1 = abs(d1), d__2 = abs(d2);
		d1 = min(d__1,d__2);
		if (dpiv <= d1) {
		    dpiv = d1;
		    *kbsq = k;
		}
	    }
	}
	*pivot = -y[*kbsq];
    }
/* pivot .ge. tolpiv */
    return 0;
} /* s5chzq_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5chzq */
/* Subroutine */ int s5getp_(logical *feasible, logical *gotr, integer *
	eigzhz, integer *maxr, integer *lenr, integer *ns, doublereal *r__, 
	doublereal *rg, doublereal *p, doublereal *gp, doublereal *php)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dcopy_(integer *, doublereal *, integer *, doublereal 
	    *, integer *);
    static doublereal rlast;
    extern /* Subroutine */ int s6rcol_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *), 
	    s6rsol_(integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *);
    static integer lrlast;

/*     ================================================================== */
/*     s5getp  computes a search direction  p  for the superbasic */
/*     variables, using the current reduced gradient  g. */

/*     29 Mar 2001: R stored by rows. */
/*     20 May 2001: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --r__;
    --p;
    --rg;

    /* Function Body */
    if (*feasible && *gotr) {
	if (*eigzhz == 1) {
/* ----------------------------------------------------------- */
/* The Newton direction. */
/* ----------------------------------------------------------- */
	    dcopy_(ns, &rg[1], &c__1, &p[1], &c__1);
	    s6rsol_(&c__1, maxr, ns, lenr, &r__[1], &p[1]);
	    *php = ddot_(ns, &p[1], &c__1, &p[1], &c__1);
	    s6rsol_(&c__0, maxr, ns, lenr, &r__[1], &p[1]);
	    dscal_(ns, &c_b81, &p[1], &c__1);
	    *gp = ddot_(ns, &rg[1], &c__1, &p[1], &c__1);
	} else {
/* ----------------------------------------------------------- */
/* A direction of zero or negative curvature. */
/* ----------------------------------------------------------- */
	    s6rcol_(&c__1, ns, maxr, ns, lenr, &r__[1], &p[1], &lrlast);
	    if (*ns > 1) {
		i__1 = *ns - 1;
		s6rsol_(&c__0, maxr, &i__1, lenr, &r__[1], &p[1]);
	    }
	    p[*ns] = -1.;
	    rlast = r__[lrlast];
/* the last diag of R. */
	    if (rlast >= 0.) {
/* Computing 2nd power */
		d__1 = rlast;
		*php = d__1 * d__1;
	    } else {
/* Computing 2nd power */
		d__1 = rlast;
		*php = -(d__1 * d__1);
	    }
	    *gp = ddot_(ns, &rg[1], &c__1, &p[1], &c__1);
	    if (*gp > 0.) {
		dscal_(ns, &c_b81, &p[1], &c__1);
		*gp = -(*gp);
	    }
	}
    } else {
/* -------------------------------------------------------------- */
/* A direction of steepest-descent. */
/* -------------------------------------------------------------- */
	dcopy_(ns, &rg[1], &c__1, &p[1], &c__1);
	dscal_(ns, &c_b81, &p[1], &c__1);
	*gp = ddot_(ns, &rg[1], &c__1, &p[1], &c__1);
	*php = 0.;
    }
/* Feasible and GotR */
    return 0;
} /* s5getp_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5getp */
/* Subroutine */ int s5zhzeig_(logical *testrd, integer *eigh, integer *
	eigzhz, integer *itn, integer *maxr, integer *lenr, integer *ns, 
	doublereal *rlast, doublereal *r__, integer *iw, integer *leniw, 
	doublereal *rw, integer *lenrw)
{
    /* Format strings */
    static char fmt_1000[] = "(\002 Itn\002,i7,\002: Reduced Hessian is inde"
	    "finite.\002,\002 Square of diag, min diag = \002,1p,2e9.1)";
    static char fmt_2000[] = "(\002 Itn\002,i7,\002: Reduced Hessian is semi"
	    "definite.\002,\002 Square of diag, min diag = \002,1p,2e9.1)";

    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static doublereal hcondbnd;
    static char str[110];
    static doublereal drsq, condh, drmin, drmax;
    extern /* Subroutine */ int s6rcnd_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), snprnt_(
	    integer *, char *, integer *, integer *, ftnlen);
    static doublereal drsqmin;

    /* Fortran I/O blocks */
    static icilist io___136 = { 0, str, 0, fmt_1000, 110, 1 };
    static icilist io___137 = { 0, str, 0, fmt_2000, 110, 1 };


/* ================================================================= */
/* s5ZHZeig  estimates the inertia of the current reduced Hessian. */

/* On entry, */
/*   eigH   encodes the inertia of the Hessian H */
/*              eigH   Hessian */
/*             -----   --------------------- */
/*              -1                indefinite */
/*               0     positive semidefinite */
/*               1     positive semidefinite */

/*   Rlast   is the last diagonal element of R */
/*   TestRd  indicates if the sign of Rlast is to be tested. */
/*           If TestRd = true, then */

/*                  (>)                       (Pos Def ) */
/*            Rlast (=) 0 implies that ZHZ is (singular) */
/*                  (<)                       (indef   ) */

/* On exit, */
/*   eigZHZ  encodes the inertia of the reduced Hessian. */


/* 15 Jul 1995: First version of s5ZHZeig. */
/* 02 Aug 2003: snPRNT adopted. */
/* 16 Apr 2005: Signed square root of SC stored in Rlast. */
/* 26 Dec 2014: Fixed bug in the test for a singular R. */
/* ================================================================= */
/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
    /* Parameter adjustments */
    --r__;
    --iw;
    --rw;

    /* Function Body */
    hcondbnd = rw[85];
/* bound on the condition of ZHZ */
    if (*ns == 0) {
/* ZHZ is positive definite at a vertex. */
	*eigzhz = 1;
    } else {
/* Compute dRsqmin, the square of the */
/* smallest allowable diagonal of a */
/* positive-definite ZHZ. */
	i__1 = *ns - 1;
	s6rcnd_(maxr, &i__1, lenr, &r__[1], &drmax, &drmin, &condh);
	drsqmin = drmax * (drmax / hcondbnd);
/* Computing 2nd power */
	d__1 = *rlast;
	drsq = d__1 * d__1;
	if (*testrd && *rlast < 0.) {
	    drsq = -drsq;
	}
	if (drsq >= drsqmin) {
	    *eigzhz = 1;
	} else if (drsq >= -drsqmin) {
	    *eigzhz = 0;
	} else {
	    *eigzhz = -1;
	}
	if (*eigh == 0 || *eigh == 1) {
/* should be a positive semidefinite ZHZ */
	    if (*eigzhz == -1) {
		s_wsfi(&io___136);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&drsq, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&drsqmin, (ftnlen)sizeof(doublereal));
		e_wsfi();
		snprnt_(&c__21, str, &iw[1], leniw, (ftnlen)110);
	    }
	    if (*eigh == 1 && *eigzhz == 0) {
/* should be a positive-definite ZHZ */
		s_wsfi(&io___137);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&drsq, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&drsqmin, (ftnlen)sizeof(doublereal));
		e_wsfi();
		snprnt_(&c__21, str, &iw[1], leniw, (ftnlen)110);
	    }
	}
    }
    return 0;
} /* s5zhzeig_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5ZHZeig */
/* Subroutine */ int s5qpfg_(S_fp hprod, U_fp hprod1, integer *hvcalls, 
	integer *ngqp, integer *ngobj0, integer *ngobj, integer *nnh, integer 
	*neh, integer *nloch, integer *loch, integer *indh, doublereal *hcol, 
	integer *sqstat, doublereal *fqp, doublereal *gobj, doublereal *gqp, 
	integer *lenx0, integer *nx0, doublereal *x0, doublereal *x, 
	doublereal *dx, char *cu, integer *lencu, integer *iu, integer *leniu,
	 doublereal *ru, integer *lenru, char *cw, integer *lencw, integer *
	iw, integer *leniw, doublereal *rw, integer *lenrw, ftnlen cu_len, 
	ftnlen cw_len)
{
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern /* Subroutine */ int dload_(integer *, doublereal *, doublereal *, 
	    integer *), dcopy_(integer *, doublereal *, integer *, doublereal 
	    *, integer *), daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static integer nzero;

/*     ================================================================== */
/*     s5QPfg  computes various quantities associated with the LP/QP. */

/*       1.  fQP =  gObj'*(x-x0)  + half*(x - x0)'*H*(x - x0) */
/*       2.  gQP =  gradient of fQP */

/*     On entry, */
/*     ngQP         is max( ngObj, nnH ) */
/*     x(ngQP)      are the nonlinear variables */
/*     x0(ngQP)     is the base point x0 (scaled if necessary) */
/*     gObj(ngObj)  defines the explicit QP linear term */

/*     On exit, */
/*     fQP          is the QP quadratic term (1) above */
/*     gQP(ngQP)    is the gradient of fQP */
/*     dx(ngQP)     is  x-x0 */

/*     02 May 1992: First version of s5QPfg. */
/*     23 Oct 1993: Hx added as an argument. */
/*     29 Oct 1993: Modified to compute only the QP objective. */
/*     07 Oct 1994: gQP added as an argument. */
/*     09 Dec 2004: f77+ version. */
/*     02 Mar 2013: HvCalls added. */
/*     03 Nov 2014: neH, indH, locH, Hcol added as arguments. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --dx;
    --x;
    --gqp;
    --gobj;
    --hcol;
    --indh;
    --loch;
    --x0;
    cu -= 8;
    --iu;
    --ru;
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    if (*ngqp <= 0) {
	return 0;
    }
    dcopy_(ngqp, &x[1], &c__1, &dx[1], &c__1);
    if (*nx0 > 0) {
	daxpy_(ngqp, &c_b81, &x0[1], &c__1, &dx[1], &c__1);
    }
    *fqp = 0.;
    if (*nnh > 0) {
	(*hprod)((U_fp)hprod1, nnh, neh, nloch, &loch[1], &indh[1], &hcol[1], 
		&dx[1], &gqp[1], sqstat, cu + 8, lencu, &iu[1], leniu, &ru[1],
		 lenru, cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)
		8, (ftnlen)8);
	++(*hvcalls);
	*fqp = ddot_(nnh, &dx[1], &c__1, &gqp[1], &c__1) * .5;
    }
    nzero = *ngqp - *nnh;
    if (nzero > 0) {
	dload_(&nzero, &c_b5, &gqp[*nnh + 1], &c__1);
    }
    if (*ngobj > 0) {
	*fqp += ddot_(ngobj, &gobj[1], &c__1, &dx[1], &c__1);
	daxpy_(ngobj, &c_b149, &gobj[1], &c__1, &gqp[1], &c__1);
    }
    return 0;
} /* s5qpfg_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5QPfg */
/* Subroutine */ int s5qpitn_(integer *iexit, S_fp hprod, U_fp hprod1, 
	integer *hvcalls, integer *eigh, integer *eigzhz, logical *bndswap, 
	logical *elastic, logical *feasible, logical *gotgqp, logical *goth, 
	logical *gotr, logical *increase, logical *needf, logical *needv, 
	logical *needpi, logical *newsb, integer *itn, integer *lenr, integer 
	*m, integer *mbs, integer *maxr, integer *maxs, integer *n, integer *
	nb, integer *nnh0, integer *nnh, integer *ns, integer *ngqp0, integer 
	*ngqp, integer *ndegen, integer *elastics, integer *lurequest, 
	integer *kp, integer *jbq, integer *jsq, integer *jbr, integer *jsr, 
	integer *jq, integer *jqsave, integer *nfmove, integer *nuncon, 
	doublereal *djq0, doublereal *djq, integer *minimize, doublereal *
	objqp, doublereal *featol, doublereal *pivot, doublereal *step, 
	doublereal *tolinc, doublereal *wtinf, integer *nea, integer *nloca, 
	integer *loca, integer *inda, doublereal *acol, integer *neh, integer 
	*nloch, integer *loch, integer *indh, doublereal *hcol, integer *
	etype, integer *estate, integer *feastype, integer *hs, integer *kbs, 
	doublereal *bl, doublereal *bu, doublereal *blqp, doublereal *buqp, 
	doublereal *blbs, doublereal *bubs, doublereal *gbs, doublereal *gqp, 
	doublereal *hdx, doublereal *pbs, doublereal *rg, doublereal *r__, 
	doublereal *x, doublereal *xbs, doublereal *y, doublereal *y1, 
	doublereal *y2, char *cu, integer *lencu, integer *iu, integer *leniu,
	 doublereal *ru, integer *lenru, char *cw, integer *lencw, integer *
	iw, integer *leniw, doublereal *rw, integer *lenrw, ftnlen cu_len, 
	ftnlen cw_len)
{
    /* Format strings */
    static char fmt_1000[] = "(\002 Itn\002,i7,\002: Bad direction after add"
	    "ing a superbasic.\002)";
    static char fmt_9999[] = "(\002 Itn\002,i7,\002: Choose q failed in s5QP"
	    "itn!!\002)";

    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    extern /* Subroutine */ int s2unpack_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, doublereal *), s5zhzeig_(logical *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, doublereal *, integer *, integer *, doublereal *, integer *);
    static integer jqestate, jrestate;
    extern /* Subroutine */ int s2scatter_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *);
    static logical unbounded;
    static integer jr;
    static doublereal gp;
    static integer ns1, nbs;
    static doublereal eps, php;
    static integer ksq, lrs;
    static char str[80];
    static integer nbs1;
    static doublereal pbs1, eps0;
    extern /* Subroutine */ int s5zp_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *, doublereal *, doublereal *, integer *, integer *, doublereal *,
	     integer *);
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer kbsq;
    static logical move;
    static doublereal gpqp, tolp;
    static integer ntry;
    static doublereal tolp0;
    extern /* Subroutine */ int s5bsx_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *), dload_(integer *, 
	    doublereal *, doublereal *, integer *), dscal_(integer *, 
	    doublereal *, doublereal *, integer *);
    static doublereal bigdx, exact, bound, norma;
    static logical uncon;
    static doublereal stepb, phpqp, rlast, pnorm, stepp;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *), s5sdel_(integer *, integer *,
	     integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *), s6rdel_(integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *), 
	    s2bsol_(integer *, integer *, integer *, doublereal *, doublereal 
	    *, integer *, integer *, doublereal *, integer *), s5rcol_(
	    integer *, S_fp, U_fp, integer *, integer *, integer *, integer *,
	     doublereal *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, integer *, integer 
	    *, integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, ftnlen, ftnlen), 
	    s5getp_(logical *, logical *, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *), s5chzq_(integer *, integer *, integer *, integer *
	    , integer *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, doublereal *, integer *), s5step_(integer *, integer *
	    , doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, logical *, logical *, logical *, logical *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    static logical rcheck;
    static doublereal infbnd;
    extern /* Subroutine */ int s2bmod2_(integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *);
    static logical hitcon;
    static integer inform__, infpiv, lrlast;
    static doublereal sclpiv;
    static logical testrd, hitlow;
    static doublereal tolpiv;
    static integer sqstat;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen), s6rswap_(integer *, integer *, integer *, doublereal *
	    , doublereal *, doublereal *, integer *, doublereal *);
    static doublereal signobj;
    static logical onbound;
    static integer jqstate, jrstate;
    static doublereal stepmax;
    extern /* Subroutine */ int s5checkp_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    integer *), s5rcheck_(integer *, S_fp, U_fp, integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     integer *, integer *, doublereal *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___162 = { 0, str, 0, fmt_1000, 80, 1 };
    static icilist io___184 = { 0, str, 0, fmt_9999, 80, 1 };


/*     ================================================================== */
/*     s5QPitn performs a QP step. */

/*     On entry, */
/*        NewSB = true implies that variable jq just went superbasic. */
/*                In this case: */
/*                pBS  satisfies B pBS = a(jq). */
/*                y1   satisfies L  y1 = a(jq). */

/*     On exit, */
/*        pBS contains the most recent QP search direction. */

/*      iExit       Result */
/*      -----       ------ */
/*       -2         reduced Hessian is not positive semidefinite */
/*       -1         unbounded */
/*        0         normal exit */
/*       >0         Fatal LU error */

/*     25 Nov 1991: First version of s5QPitn. */
/*     05 Jan 1996: Positive semidefinite R treated correctly. */
/*     29 Aug 1996: First min sum version added. */
/*     27 Jul 1997: Thread-safe version. */
/*     02 Feb 1998: Piecewise linear line search added. */
/*     23 Mar 2000: gQP  and  H  scaled. */
/*     16 Oct 2000: Reverted to non-bordered version of s5QPitn. */
/*     04 Dec 2000: R converted to row-wise storage. */
/*     02 Aug 2003: snPRNT adopted. */
/*     07 May 2006: s5Zp added to compute Z*p. */
/*     08 Apr 2008: eState accessed in elastic mode only. */
/*     04 Jul 2008: Modify both bl and bu in elastic phase 1. */
/*     02 Mar 2013: HvCalls added. */
/*     02 Nov 2014: Set LUrequest = 25 when Z'HZ is indefinite. */
/*     03 Nov 2014: neH, indH, locH, Hcol added as arguments. */
/*     10 Dec 2014: Removed unused argument obj. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/* number of LU mods */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --r__;
    --xbs;
    --pbs;
    --gbs;
    --bubs;
    --blbs;
    --kbs;
    --feastype;
    --rg;
    --y2;
    --y1;
    --y;
    --x;
    --buqp;
    --blqp;
    --bu;
    --bl;
    --hs;
    --estate;
    --etype;
    --hdx;
    --gqp;
    --acol;
    --inda;
    --loca;
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
    eps = rw[1];
/* machine precision.  IEEE DP  2.22e-16 */
    eps0 = rw[2];
/* eps**(4/5)          IEEE DP  3.00e-13 */
    tolpiv = rw[60];
/* excludes small elements of pBS. */
    infbnd = rw[70];
/* definition of an infinite bound */
    bigdx = rw[72];
/* unbounded step. */
    *iexit = 0;
    sqstat = 0;
    unbounded = FALSE_;
    signobj = (doublereal) (*minimize);
    nbs = *m + *ns;
    if (*newsb) {
/*        --------------------------------------------------------------- */
/*        New superbasic. */
/*        Definite must be true if there is a new superbasic. */
/*        --------------------------------------------------------------- */
	ns1 = *ns + 1;
	nbs1 = nbs + 1;
	kbs[nbs1] = *jq;
	xbs[nbs1] = x[*jq];
	blbs[nbs1] = blqp[*jq];
	bubs[nbs1] = buqp[*jq];
	jqstate = hs[*jq];
	*eigzhz = 0;
	if (*gotr) {
/*           ------------------------------------------------------------ */
/*           Add the new column to R at position nS+1. */
/*           Check for a singular or indefinite reduced Hessian. */
/*           ------------------------------------------------------------ */
	    kbs[nbs1] = *jq;
	    s5rcol_(iexit, (S_fp)hprod, (U_fp)hprod1, hvcalls, minimize, jq, &
		    ns1, &rlast, &lrlast, maxr, lenr, m, mbs, n, nb, nnh, &
		    ns1, nea, nloca, &loca[1], &inda[1], &acol[1], neh, nloch,
		     &loch[1], &indh[1], &hcol[1], &kbs[1], &r__[1], &y[1], &
		    y2[1], &pbs[1], cu + 8, lencu, &iu[1], leniu, &ru[1], 
		    lenru, cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (
		    ftnlen)8, (ftnlen)8);
	    if (*iexit != 0) {
		goto L900;
	    }
/* Test sign(Rlast) in s5ZHZeig */
	    testrd = TRUE_;
	    s5zhzeig_(&testrd, eigh, eigzhz, itn, maxr, lenr, &ns1, &rlast, &
		    r__[1], &iw[1], leniw, &rw[1], lenrw);
	    if (*feasible) {
		if (*eigzhz == -1) {
		    *iexit = -2;
		    *lurequest = 25;
		    goto L900;
		}
	    } else {
/*              ZHZ needs to be pos def in phase 1 */
/*              If not, stop updating it. */
		*gotr = *eigzhz == 1;
		*eigzhz = 0;
/* ZHZ isn't used to define pBS */
	    }
	}
/* -->  R may be checked here by setting  Rcheck = .true. */
/* GotR */
	rcheck = FALSE_;
	if (rcheck && *gotr) {
	    s5rcheck_(iexit, (S_fp)hprod, (U_fp)hprod1, hvcalls, eigh, itn, 
		    minimize, maxr, lenr, m, mbs, n, nb, nnh, &ns1, nea, 
		    nloca, &loca[1], &inda[1], &acol[1], neh, nloch, &loch[1],
		     &indh[1], &hcol[1], &kbs[1], &r__[1], &y[1], &y2[1], &
		    pbs[1], cu + 8, lencu, &iu[1], leniu, &ru[1], lenru, cw + 
		    8, lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)8, (
		    ftnlen)8);
	    if (*iexit != 0) {
		goto L900;
	    }
	}
	if (*increase) {
	    *jqsave = *jq;
	} else {
	    *jqsave = -(*jq);
	}
	*ns = ns1;
	nbs = nbs1;
	feastype[nbs] = 0;
	if (*feasible && *needf && *gotgqp && *jq <= *ngqp) {
	    gbs[nbs] = signobj * gqp[*jq];
	} else {
	    gbs[nbs] = 0.;
	}
/*        =============================================================== */
/*        Set eState(jq) and the elastic parts of blBS and buBS. */
/*        =============================================================== */
	if (*elastic) {
	    jqestate = estate[*jq];
/*           Check if the new superbasic is an elastic variable that */
/*           wants to move infeasible. If so, set its elastic state and */
/*           modify the working bounds. */
/* Saved value */
	    if (etype[*jq] > 0 && estate[*jq] == 0) {
		if (*increase) {
		    if (jqstate == 1 || jqstate == 4) {
			++(*elastics);
			estate[*jq] = 2;
			blqp[*jq] = buqp[*jq];
			buqp[*jq] = infbnd;
			blbs[nbs] = blqp[*jq];
			bubs[nbs] = buqp[*jq];
			if (*feasible && *needv) {
			    gbs[nbs] += *wtinf;
			}
		    }
		} else {
		    if (jqstate == 0 || jqstate == 4) {
			++(*elastics);
			estate[*jq] = 1;
			buqp[*jq] = blqp[*jq];
			blqp[*jq] = -infbnd;
			blbs[nbs] = blqp[*jq];
			bubs[nbs] = buqp[*jq];
			if (*feasible && *needv) {
			    gbs[nbs] -= *wtinf;
			}
		    }
		}
	    }
	}
/*        --------------------------------------------------------------- */
/*        In phase 1, or phase 2 for an LP, price can select nonbasics */
/*        floating free between their bounds with zero reduced cost. */
/*        We have to check that dqj is not zero. */
/*        --------------------------------------------------------------- */
/* Elastic */
	rg[*ns] = *djq;
	if (! (*feasible) || *needf && ! (*goth)) {
	    if (hs[*jq] == -1) {
		if (*increase) {
		    rg[*ns] = -1.;
		} else {
		    rg[*ns] = 1.;
		}
	    }
	}
	*jsq = *jq;
	hs[*jq] = 2;
    }
/*     ------------------------------------------------------------------ */
/*     Compute the search direction for the superbasics. */
/*     Store the free components of the search direction in pBS(1:nBS). */
/*     First, find the search direction pS for the superbasics, store */
/*     it in  pBS(m+1:nBS), and find its norm.  Store the search */
/*     direction for the basic variables in pBS(1)  ,...,pBS(m). */
/*     ------------------------------------------------------------------ */
/* newSB */
L100:
    s5getp_(feasible, gotr, eigzhz, maxr, lenr, ns, &r__[1], &rg[1], &pbs[*m 
	    + 1], &gp, &php);
    pbs1 = pbs[*m + *ns];
    s5zp_(iexit, m, mbs, n, nb, ns, &eps0, &pnorm, nea, nloca, &loca[1], &
	    inda[1], &acol[1], &kbs[1], &pbs[1], &y2[1], &iw[1], leniw, &rw[1]
	    , lenrw);
    if (*iexit != 0) {
	goto L900;
    }
    if (*feasible) {
/* -------------------------------------------------------------- */
/* If R is singular, ensure that pBS is a feasible direction. */
/* A nonzero exit value of inform implies that the directional */
/* derivative is too small to be relied upon. */
/* -------------------------------------------------------------- */
	if (*gotr && *eigzhz == 0) {
	    s5checkp_(&inform__, itn, &nbs, jqsave, &kbs[1], &gp, &pbs[1], &
		    iw[1], leniw);
	    if (inform__ > 0) {
		*lurequest = 24;
		goto L900;
	    }
	}
	if (*newsb && *gotr) {
/* Check for a feasible direction. */
/* A large  rgTol  may give a pBS(nBS) with the wrong sign. */
/* If so, continue minimizing with the old superbasic set. */
	    if (*djq * pbs1 > 0.) {
		s_wsfi(&io___162);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
		--(*ns);
		--nbs;
		hs[*jq] = jqstate;
		if (*elastic) {
/*                 If a variable went elastic, reset it as not elastic. */
		    if (estate[*jq] != jqestate) {
			--(*elastics);
			blqp[*jq] = bl[*jq];
			buqp[*jq] = bu[*jq];
			estate[*jq] = jqestate;
		    }
		}
		*jq = 0;
		*djq = *djq0;
		*jsq = -(*jsq);
		*jqsave = 0;
		*eigzhz = 1;
		*newsb = FALSE_;
		goto L100;
	    }
	    *bndswap = FALSE_;
	}
/*        --------------------------------------------------------------- */
/*        Compute y = pBS(scattered) and Hdx(scattered). */
/*        The vector Hdx is used to update the objective and gradient of */
/*        the QP.  Form  gpQP  and  pHpQP  for the quadratic. */
/*        gp = gpQP - pBS(kObj) + terms from the elastic gradient. */
/*        --------------------------------------------------------------- */
	if (*needf && (*gotgqp || *goth)) {
	    s2scatter_(ngqp, &nbs, &kbs[1], &c_b149, &pbs[1], &y[1]);
	    if (*gotgqp) {
		gpqp = ddot_(ngqp, &gqp[1], &c__1, &y[1], &c__1);
	    }
	    if (*goth) {
		phpqp = 0.;
		(*hprod)((U_fp)hprod1, nnh, neh, nloch, &loch[1], &indh[1], &
			hcol[1], &y[1], &hdx[1], &sqstat, cu + 8, lencu, &iu[
			1], leniu, &ru[1], lenru, cw + 8, lencw, &iw[1], 
			leniw, &rw[1], lenrw, (ftnlen)8, (ftnlen)8);
		++(*hvcalls);
		phpqp += ddot_(nnh, &y[1], &c__1, &hdx[1], &c__1);
	    }
	}
    }
/*     ------------------------------------------------------------------ */
/*     Find the nearest constraint in direction  x + step*pBS (step > 0). */
/*     Exact  is the step that takes xBS(kp) exactly onto bound. It may */
/*     be positive or slightly negative. (Not defined if Unbounded.) */

/*     If OnBound  is true, step is a step that reaches a bound exactly. */
/*     xBS(kp) reaches the value bound.  If we take a constrained step, */
/*     bound is used to put the new nonbasic variable x(jr) exactly on */
/*     its bound. */

/*     If Unbounded is true, step = stepMax. */
/*     ------------------------------------------------------------------ */
/* Feasible */
    stepmax = bigdx / pnorm;
    sclpiv = 1.;
    tolp0 = tolpiv;
    tolp = tolpiv * pnorm;
    ntry = 0;
/* +    Repeat */
L200:
    tolp /= sclpiv;
    tolp0 /= sclpiv;
    s5step_(&nbs, ndegen, featol, &infbnd, &stepmax, tolinc, &tolp, &feastype[
	    1], &blbs[1], &bubs[1], &xbs[1], &pbs[1], &hitlow, &move, &
	    onbound, &unbounded, &infpiv, kp, &bound, &exact, &stepb, &stepp);
/*        Find if the step is constrained or unconstrained. */
/*        If R has been flagged as singular, we double check by trying */
/*        to compute the QP minimizer along pBS.  If the minimizer */
/*        exists,  the singularity tolerance must be too large. */
    if (*feasible) {
	if (*eigzhz == 1) {
	    uncon = stepp > 1.;
	} else if (*eigzhz == 0) {
	    if (unbounded) {
		uncon = FALSE_;
	    } else {
		uncon = php > 0. && stepp * php > -gp;
	    }
	} else {
/*  eigZHZ .eq. INDEF */
	    if (*jq <= *nnh) {
/*                 Nonlinear x(jq), limit the step size */
		uncon = stepp > 1.;
	    } else {
/*                 Linear    x(jq). step to the boundary */
		uncon = unbounded;
	    }
	}
	unbounded = unbounded && ! uncon || stepmax <= 1.;
    } else {
/* infeasible */
	uncon = FALSE_;
    }
    sclpiv = 10.;
    ++ntry;
/* +    until    ( infpiv .eq. 0 .and. (.not.Unbounded .or. Feasible) .or. */
/* +                 ntry .ge. mtry) */
    if (! (infpiv == 0 && (! unbounded || *feasible) || ntry >= 6)) {
	goto L200;
    }
    if (unbounded) {
	*iexit = -1;
	goto L900;
    }
    hitcon = ! uncon;
    *needpi = TRUE_;
    if (hitcon) {
	*nuncon = 0;
	*step = stepb;
    } else {
	*pivot = 0.;
	if (*eigzhz == 1) {
	    ++(*nuncon);
	    *step = 1.;
	} else if (*eigzhz == 0) {
	    ++(*nuncon);
	    *step = -gp / php;
	    *eigzhz = 1;
	} else {
	    *step = 1.;
	}
    }
/* ----------------------------------------------------------------- */
/* Compute the new objective function. */
/* Note: pHp = signObj*pHpQP */
/* ----------------------------------------------------------------- */
    if (*feasible && *needf && *step > 0.) {
	if (*gotgqp) {
	    *objqp += *step * gpqp;
	}
	if (*goth) {
/* Computing 2nd power */
	    d__1 = *step;
	    *objqp += phpqp * .5 * (d__1 * d__1);
	    daxpy_(nnh, step, &hdx[1], &c__1, &gqp[1], &c__1);
	}
    }
    if (*feasible && move) {
	++(*nfmove);
    }
/*     ------------------------------------------------------------------ */
/*     Update the basic variables xBS. */
/*     ------------------------------------------------------------------ */
    daxpy_(&nbs, step, &pbs[1], &c__1, &xbs[1], &c__1);
    s5bsx_(&c__1, &nbs, nb, &kbs[1], &x[1], &xbs[1]);
    if (hitcon) {
/*        =============================================================== */
/*        There is a blocking variable. */
/*        It could be a fixed variable, whose new state must be 4. */
/*        =============================================================== */
	*pivot = -pbs[*kp];
	jr = kbs[*kp];
	*bndswap = jr == abs(*jqsave);
/* $$$!        10 Mar 2004: Care is needed to prevent the */
/* $$$!        new nonbasic variable jr from ending up slightly inside */
/* $$$!        its bound.  EXPAND normally ensures that x(jr) will be */
/* $$$!        ON or slightly OUTSIDE its bound, but now we realise that */
/* $$$!        rounding error might make it slightly INSIDE. */
/* $$$ */
/* $$$         if (OnBound) then */
/* $$$            x(jr) = bound */
/* $$$         else if (HitLow) then */
/* $$$            x(jr) = min( x(jr), bl(jr) ) */
/* $$$         else */
/* $$$            x(jr) = max( x(jr), bu(jr) ) */
/* $$$         end if */
	if (onbound) {
	    x[jr] = bound;
	}
	if (*elastic) {
	    jrestate = estate[jr];
	} else {
	    jrestate = 0;
	}
	if (jrestate == 0) {
	    if (blbs[*kp] == bubs[*kp]) {
		jrstate = 4;
	    } else if (hitlow) {
		jrstate = 0;
	    } else {
		jrstate = 1;
	    }
	} else {
/*           Elastic x(jr) hits its bound. */
/*           Reset the true upper and lower bounds. */
	    estate[jr] = 0;
	    --(*elastics);
	    blqp[jr] = bl[jr];
	    buqp[jr] = bu[jr];
	    if (jrestate == 1) {
		if (blqp[jr] == buqp[jr]) {
		    jrstate = 4;
		} else if (onbound) {
		    jrstate = 0;
		} else if (x[jr] < buqp[jr]) {
		    jrstate = -1;
		} else {
		    jrstate = 1;
		}
	    } else {
/* jreState =  2 */
		if (blqp[jr] == buqp[jr]) {
		    jrstate = 4;
		} else if (onbound) {
		    jrstate = 1;
		} else if (x[jr] > blqp[jr]) {
		    jrstate = -1;
		} else {
		    jrstate = 0;
		}
	    }
	}
	if (*kp <= *m) {
/*           ============================================================ */
/*           A variable in B hit a bound. */
/*           Find column kSq = kBSq-m  of S to replace column kp of B. */
/*           If nS = 1 there is no choice. */
/*           ============================================================ */
	    if (*ns == 1) {
		kbsq = nbs;
		*pivot /= pbs1;
	    } else {
		dload_(m, &c_b5, &y2[1], &c__1);
		y2[*kp] = 1.;
		s2bsol_(iexit, &c__2, m, &y2[1], &y[1], &iw[1], leniw, &rw[1],
			 lenrw);
		if (*iexit != 0) {
		    return 0;
		}
		s5chzq_(m, mbs, n, nb, ns, &kbsq, pivot, &tolp0, nea, nloca, &
			loca[1], &inda[1], &acol[1], &kbs[1], &blqp[1], &buqp[
			1], &xbs[1], &y[1], &iw[1], leniw, &rw[1], lenrw);
		if (kbsq <= 0) {
		    s_wsfi(&io___184);
		    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		    e_wsfi();
		    snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
		    kbsq = nbs;
		}
	    }
	    ksq = kbsq - *m;
	    hs[jr] = jrstate;
	    *jbr = jr;
/* Outgoing basic */
	    *jsr = kbs[kbsq];
/* Outgoing superbasic */
	    kbs[kbsq] = *jbr;
	    *jbq = *jsr;
/* Incoming basic */
	    kbs[*kp] = *jsr;
	    blbs[*kp] = blbs[kbsq];
	    bubs[*kp] = bubs[kbsq];
	    xbs[*kp] = xbs[kbsq];
	    hs[*jbq] = 3;
	    if (*ns > 1 && *gotr) {
/*              Finish computing y(m+1), ..., y(m+nS). */
		y[kbsq] = -(*pivot + 1.);
		d__1 = 1. / *pivot;
		dscal_(ns, &d__1, &y[*m + 1], &c__1);
		s6rswap_(maxr, ns, lenr, &r__[1], &y2[1], &y[*m + 1], &ksq, &
			eps0);
	    }
/*           ------------------------------------------------------------ */
/*           Get a new  y1, used to modify L and U.  If the outgoing */
/*           superbasic just came in, we already have it. */
/*           ------------------------------------------------------------ */
	    if (*jsr != *jq) {
		s2unpack_(jbq, m, n, nea, &norma, nloca, &loca[1], &inda[1], &
			acol[1], &y1[1]);
		s2bsol_(iexit, &c__0, m, &y1[1], &y[1], &iw[1], leniw, &rw[1],
			 lenrw);
		if (*iexit != 0) {
		    return 0;
		}
	    }
/*           Update the LU factors. */
	    ++iw[216];
	    s2bmod2_(&inform__, kp, m, &y1[1], &iw[1], leniw, &rw[1], lenrw);
/*           If necessary, complete the update but refactor on exit. */
	    if (inform__ == -1) {
		*lurequest = 5;
	    }
/* Singular after LU mod */
	    if (inform__ == 2) {
		*lurequest = 6;
	    }
/* Unstable LU mod */
	    if (inform__ == 7) {
		*lurequest = 7;
	    }
/* Insufficient free memory */
	} else {
/*           ============================================================ */
/*           A variable in S hit a bound. */
/*           ============================================================ */
	    hs[jr] = jrstate;
	    *jsr = jr;
	    kbsq = *kp;
	    ksq = kbsq - *m;
	}
/*        Delete the kSq-th superbasic and adjust all arrays in BS order. */
	s5sdel_(&ksq, m, ns, &nbs, &kbs[1], &blbs[1], &bubs[1], &gbs[1], &rg[
		1], &xbs[1]);
	if (*gotr) {
/*           ------------------------------------------------------------ */
/*           Cyclically demote column kSq of R to position nS. */
/*           ------------------------------------------------------------ */
	    if (ksq < *ns) {
		s6rdel_(&ksq, maxr, ns, lenr, &r__[1], &eps);
	    }
	}
/* Feasible and GotH */
	--(*ns);
	--nbs;
	if (*eigzhz == 0) {
/*           Recheck the last diagonal of R */
/*           It can only increase in magnitude */
	    if (*feasible && *gotr) {
		if (*ns > 0) {
		    lrs = (*ns - 1) * *maxr + (3 - *ns) * *ns / 2;
/* Magic formula! */
		    rlast = r__[lrs];
		}
		testrd = FALSE_;
/* ignore the sign of Rlast */
		s5zhzeig_(&testrd, eigh, eigzhz, itn, maxr, lenr, ns, &rlast, 
			&r__[1], &iw[1], leniw, &rw[1], lenrw);
	    } else if (*ns == 0) {
		*eigzhz = 1;
	    } else {
		*eigzhz = 0;
	    }
	}
/* -->     R can be checked here. */
    }
/* HitCon */
L900:
    return 0;
} /* s5qpitn_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5QPitn */
/* Subroutine */ int s5rcheck_(integer *iexit, S_fp hprod, U_fp hprod1, 
	integer *hvcalls, integer *eigh, integer *itn, integer *minimize, 
	integer *maxr, integer *lenr, integer *m, integer *mbs, integer *n, 
	integer *nb, integer *nnh, integer *ns, integer *nea, integer *nloca, 
	integer *loca, integer *inda, doublereal *acol, integer *neh, integer 
	*nloch, integer *loch, integer *indh, doublereal *hcol, integer *kbs, 
	doublereal *r__, doublereal *v, doublereal *w, doublereal *y, char *
	cu, integer *lencu, integer *iu, integer *leniu, doublereal *ru, 
	integer *lenru, char *cw, integer *lencw, integer *iw, integer *leniw,
	 doublereal *rw, integer *lenrw, ftnlen cu_len, ftnlen cw_len)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    extern /* Subroutine */ int s2unpack_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, doublereal *), s5zhzeig_(logical *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, doublereal *, integer *, integer *, doublereal *, integer *), 
	    s2scatter_(integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *);
    static integer jq, js, nbs;
    static doublereal whw, eps0;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal drsq;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal norma, rlast;
    extern /* Subroutine */ int s2bsol_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *), s6rcol_(integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *), s6rsol_(
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *);
    static integer lencol, eigzhz, lrlast;
    static logical testrd;
    static integer status;
    extern /* Subroutine */ int s2bprod_(integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, integer *);
    static doublereal signobj, rnormsq;
    extern /* Subroutine */ int s2gather_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *);

/*     ================================================================== */
/*     s5Rcheck  computes the Cholesky factor R such that */
/*     R'R = Z'HZ.  The update corresponds to the addition of a new */
/*     column to Z. */

/*     On entry, */
/*        R     holds the columns of the factor associated with the */
/*              first jRadd-1 columns of Q. */

/*        nS    is the number of columns in R. */

/*     14 Mar 2001: First version (s6Rchk) based on SNOPT routine s5Rcol. */
/*     20 Jun 2008: Renamed s5Rcheck. */
/*     29 Jun 2008: Revamped. */
/*     02 Mar 2013: HvCalls added. */
/*     03 Nov 2014: neH, indH, locH, Hcol added as arguments. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --r__;
    --y;
    --kbs;
    --w;
    --v;
    --acol;
    --inda;
    --loca;
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
    eps0 = rw[2];
    nbs = *m + *ns;
    signobj = (doublereal) (*minimize);
    *iexit = 0;
/*     ------------------------------------------------------------------ */
/*     Main loop to find a column of Z'HZ. */
/*     ------------------------------------------------------------------ */
    i__1 = *ns;
    for (js = 1; js <= i__1; ++js) {
/* Computing MIN */
	i__2 = js - 1;
	lencol = min(i__2,*nnh);
/*        --------------------------------------------------------------- */
/*        Get the nonlinear elements of the column of Z. */
/*        Find y such that B y = column jq. */
/*        Scatter the nonlinear part of y into w. */
/*        --------------------------------------------------------------- */
	jq = kbs[*m + js];
	s2unpack_(&jq, m, n, nea, &norma, nloca, &loca[1], &inda[1], &acol[1],
		 &w[1]);
	s2bsol_(iexit, &c__1, m, &w[1], &y[1], &iw[1], leniw, &rw[1], lenrw);
	if (*iexit != 0) {
	    return 0;
	}
	s2scatter_(nnh, m, &kbs[1], &c_b81, &y[1], &w[1]);
	if (jq <= *nnh) {
	    w[jq] = 1.;
	}
/*        --------------------------------------------------------------- */
/*        Compute  H*w  and  w'*H*w. */
/*        --------------------------------------------------------------- */
	whw = 0.;
	if (*nnh > 0) {
	    status = 0;
	    (*hprod)((U_fp)hprod1, nnh, neh, nloch, &loch[1], &indh[1], &hcol[
		    1], &w[1], &v[1], &status, cu + 8, lencu, &iu[1], leniu, &
		    ru[1], lenru, cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw,
		     (ftnlen)8, (ftnlen)8);
	    ++(*hvcalls);
	    whw += ddot_(nnh, &w[1], &c__1, &v[1], &c__1);
	    if (*minimize < 0) {
		dscal_(nnh, &signobj, &v[1], &c__1);
		whw = signobj * whw;
	    }
	}
	rnormsq = 0.;
	if (js > 1) {
/*           ------------------------------------------------------------ */
/*           Gather the nonlinear elements of v in w (= vBS). */
/*           Compute Z'w  (solve  B'vB = wB and form  wS = wS - S'vB). */
/*           ------------------------------------------------------------ */
	    s2gather_(nnh, &nbs, &kbs[1], &c_b149, &v[1], &w[1]);
	    s2bsol_(iexit, &c__2, m, &w[1], &v[1], &iw[1], leniw, &rw[1], 
		    lenrw);
	    if (*iexit != 0) {
		return 0;
	    }
	    if (*ns > 0) {
		s2bprod_(&c__1, &eps0, n, ns, &kbs[*m + 1], nea, nloca, &loca[
			1], &inda[1], &acol[1], &c_b81, &v[1], m, &c_b149, &w[
			*m + 1], ns);
	    }
/*           ------------------------------------------------------------ */
/*           Solve  R'v = Z(j)'Hw.  Store v in w(m+1:m+jS). */
/*           ------------------------------------------------------------ */
	    s6rsol_(&c__1, maxr, &lencol, lenr, &r__[1], &w[*m + 1]);
	    rnormsq = ddot_(&lencol, &w[*m + 1], &c__1, &w[*m + 1], &c__1);
	}
	if (js <= *nnh) {
	    drsq = whw - rnormsq;
/* square of new diagonal of R. */
	    if (drsq >= 0.) {
		rlast = sqrt(drsq);
	    } else {
		rlast = -sqrt((abs(drsq)));
	    }
	} else {
	    drsq = 0.;
	    rlast = 0.;
	}
	w[*m + js] = rlast;
/*        Insert w(m+1:m+jS) as column jS of R. */
	s6rcol_(&c__0, &js, maxr, &js, lenr, &r__[1], &w[*m + 1], &lrlast);
	testrd = TRUE_;
/* the sign of Rlast is relevant */
	s5zhzeig_(&testrd, eigh, &eigzhz, itn, maxr, lenr, ns, &rlast, &r__[1]
		, &iw[1], leniw, &rw[1], lenrw);
	if (eigzhz == -1) {
	    *iexit = 6;
	    return 0;
	}
    }
    return 0;
} /* s5rcheck_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5Rcheck */
/* Subroutine */ int s5rcol_(integer *iexit, S_fp hprod, U_fp hprod1, integer 
	*hvcalls, integer *minimize, integer *jq, integer *jradd, doublereal *
	rlast, integer *lrlast, integer *maxr, integer *lenr, integer *m, 
	integer *mbs, integer *n, integer *nb, integer *nnh, integer *ns, 
	integer *nea, integer *nloca, integer *loca, integer *inda, 
	doublereal *acol, integer *neh, integer *nloch, integer *loch, 
	integer *indh, doublereal *hcol, integer *kbs, doublereal *r__, 
	doublereal *v, doublereal *w, doublereal *y, char *cu, integer *lencu,
	 integer *iu, integer *leniu, doublereal *ru, integer *lenru, char *
	cw, integer *lencw, integer *iw, integer *leniw, doublereal *rw, 
	integer *lenrw, ftnlen cu_len, ftnlen cw_len)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    extern /* Subroutine */ int s2scatter_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *);
    static integer nbs;
    static doublereal whw, eps0;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal drsq;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), s2bsol_(integer *, integer *, integer *, doublereal *,
	     doublereal *, integer *, integer *, doublereal *, integer *), 
	    s6rcol_(integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *), s6rsol_(integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *);
    static integer lencol, sqstat;
    extern /* Subroutine */ int s2bprod_(integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, integer *);
    static doublereal signobj, rnormsq;
    extern /* Subroutine */ int s2gather_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *);

/*     ================================================================== */
/*     s5Rcol  computes column jRadd of the Cholesky factor R such that */
/*     R'R = Z'HZ.  The update corresponds to the addition of a new */
/*     column to Z. */

/*     On entry, */
/*        R     holds the columns of the factor associated with the */
/*              first jRadd-1 columns of Q. */

/*        y     is the vector such that B y = a(jq). */

/*        nS    is the number of columns in R. */

/*     11 Dec 1991: First version based on Qpsol routine Qpcolr. */
/*     24 Apr 1994: Columns of Nx no longer in Q. */
/*     27 Oct 2000: Previous version of s5Rcol. */
/*     04 Dec 2000: R converted to row-wise storage. */
/*     09 Dec 2004: Current version of s5Rcol. */
/*     02 Mar 2013: HvCalls added. */
/*     03 Nov 2014: neH, indH, locH, Hcol added as arguments. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --r__;
    --y;
    --kbs;
    --w;
    --v;
    --acol;
    --inda;
    --loca;
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
    eps0 = rw[2];
    *iexit = 0;
    nbs = *m + *ns;
    signobj = (doublereal) (*minimize);
/* Computing MIN */
    i__1 = *jradd - 1;
    lencol = min(i__1,*nnh);
/*     ------------------------------------------------------------------ */
/*     Get w, the vector of nonlinear components of the new column of Z. */
/*     ------------------------------------------------------------------ */
/*     The input vector y satisfies B y = column jq. */
/*     Scatter the nonlinear components of y into w. */
    s2scatter_(nnh, m, &kbs[1], &c_b81, &y[1], &w[1]);
    if (*jq <= *nnh) {
	w[*jq] = 1.;
    }
/*     ------------------------------------------------------------------ */
/*     Compute  H*w  and  w'*H*w. */
/*     ------------------------------------------------------------------ */
    whw = 0.;
    if (*nnh > 0) {
	sqstat = 0;
	(*hprod)((U_fp)hprod1, nnh, neh, nloch, &loch[1], &indh[1], &hcol[1], 
		&w[1], &v[1], &sqstat, cu + 8, lencu, &iu[1], leniu, &ru[1], 
		lenru, cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)8,
		 (ftnlen)8);
	++(*hvcalls);
	whw += ddot_(nnh, &w[1], &c__1, &v[1], &c__1);
	if (*minimize < 0) {
	    dscal_(nnh, &signobj, &v[1], &c__1);
	    whw = signobj * whw;
	}
    }
    rnormsq = 0.;
    if (*jradd > 1) {
/*        -------------------------------------------------------------- */
/*        Gather the nonlinear elements of v in w (= vBS). */
/*        Compute Z'w  (solve  B'vB = wB and form  wS = wS - S'vB). */
/*        -------------------------------------------------------------- */
	s2gather_(nnh, &nbs, &kbs[1], &c_b149, &v[1], &w[1]);
	s2bsol_(iexit, &c__2, m, &w[1], &v[1], &iw[1], leniw, &rw[1], lenrw);
	if (*iexit != 0) {
	    return 0;
	}
	if (*ns > 0) {
	    s2bprod_(&c__1, &eps0, n, ns, &kbs[*m + 1], nea, nloca, &loca[1], 
		    &inda[1], &acol[1], &c_b81, &v[1], m, &c_b149, &w[*m + 1],
		     ns);
	}
/*        -------------------------------------------------------------- */
/*        Solve  R'v = Z(j)'Hw.  Store v in w(m+1:m+jRadd). */
/*        -------------------------------------------------------------- */
	s6rsol_(&c__1, maxr, &lencol, lenr, &r__[1], &w[*m + 1]);
	rnormsq = ddot_(&lencol, &w[*m + 1], &c__1, &w[*m + 1], &c__1);
    }
    if (*jradd <= *nnh) {
	drsq = whw - rnormsq;
/* square of the new diagonal of R. */
	if (drsq >= 0.) {
	    *rlast = sqrt(drsq);
	} else {
	    *rlast = -sqrt((abs(drsq)));
	}
    } else {
	drsq = 0.;
	*rlast = 0.;
    }
    w[*m + *jradd] = *rlast;
/*     Insert w(m+1:m+jRadd) as column jRadd of R. */
    s6rcol_(&c__0, jradd, maxr, jradd, lenr, &r__[1], &w[*m + 1], lrlast);
    return 0;
} /* s5rcol_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5Rcol */
/* Subroutine */ int s5rg_(integer *m, integer *nbs, integer *n, integer *ns, 
	doublereal *tolz, integer *nea, integer *nloca, integer *loca, 
	integer *inda, doublereal *acol, doublereal *gbs, doublereal *pi, 
	doublereal *rg, doublereal *rgnorm, integer *kbs)
{
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    extern doublereal dnormi_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int s2bprod_(integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, integer *);

/*     ================================================================== */
/*     s5rg    calculates the reduced gradient  rg = gS - S'*pi. */

/*     23 Nov 1991: First version based on Minos routine m7rg. */
/*     16 Nov 2001: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --pi;
    --kbs;
    --gbs;
    --rg;
    --acol;
    --inda;
    --loca;

    /* Function Body */
    dcopy_(ns, &gbs[*m + 1], &c__1, &rg[1], &c__1);
    s2bprod_(&c__1, tolz, n, ns, &kbs[*m + 1], nea, nloca, &loca[1], &inda[1],
	     &acol[1], &c_b81, &pi[1], m, &c_b149, &rg[1], ns);
    *rgnorm = dnormi_(ns, &rg[1], &c__1);
    return 0;
} /* s5rg_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5rg */
/* Subroutine */ int s5sdel_(integer *ksq, integer *m, integer *ns, integer *
	nbs, integer *kbs, doublereal *blbs, doublereal *bubs, doublereal *
	gbs, doublereal *rg, doublereal *xbs)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j, k;

/*     ================================================================== */
/*     s5Sdel  deletes the kSqth superbasic variable from the arrays */
/*     kBS, blBS, blBS, gBS, rg and xBS. */

/*     16 Jun 2001: First version of s5Bswp. */
/*     16 Jun 2001: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     Shift all the arrays one place to the left. */
    /* Parameter adjustments */
    --rg;
    --xbs;
    --gbs;
    --bubs;
    --blbs;
    --kbs;

    /* Function Body */
    i__1 = *ns - 1;
    for (j = *ksq; j <= i__1; ++j) {
	k = *m + j;
	kbs[k] = kbs[k + 1];
	blbs[k] = blbs[k + 1];
	bubs[k] = bubs[k + 1];
	gbs[k] = gbs[k + 1];
	xbs[k] = xbs[k + 1];
	rg[j] = rg[j + 1];
    }
    return 0;
} /* s5sdel_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5Sdel */
/* Subroutine */ int s5zhz_(integer *iexit, S_fp hprod, U_fp hprod1, integer *
	hvcalls, integer *maxr, integer *lenr, integer *minimize, integer *m, 
	integer *mbs, integer *n, integer *nb, integer *nnh, integer *ns, 
	integer *nea, integer *nloca, integer *loca, integer *inda, 
	doublereal *acol, integer *neh, integer *nloch, integer *loch, 
	integer *indh, doublereal *hcol, doublereal *zhzmin, doublereal *
	zhzmax, doublereal *condz, doublereal *normz, integer *kbs, 
	doublereal *r__, doublereal *v, doublereal *w, doublereal *y, char *
	cu, integer *lencu, integer *iu, integer *leniu, doublereal *ru, 
	integer *lenru, char *cw, integer *lencw, integer *iw, integer *leniw,
	 doublereal *rw, integer *lenrw, ftnlen cu_len, ftnlen cw_len)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    extern /* Subroutine */ int s2unpack_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, doublereal *), s2scatter_(integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *);
    static integer jq, js, nbs;
    static doublereal eps0, diag;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal flmax;
    extern /* Subroutine */ int s2bsol_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *), s6rrow_(integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *);
    static doublereal normaj;
    extern doublereal dnormi_(integer *, doublereal *, integer *);
    static integer lrlast, sqstat;
    static doublereal normzj;
    extern /* Subroutine */ int s2bprod_(integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, integer *);
    static doublereal signobj;
    extern /* Subroutine */ int s2gather_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *);

/*     ================================================================== */
/*     s5ZHZ computes the reduced Hessian and loads it by columns into */
/*     the upper triangle R. */

/*      iExit       Status */
/*      -----       ------ */
/*        0         reduced Hessian computed successfully */
/*       >0         Fatal error in LU solve */

/*     13 Oct 1992: First version based on QPSOL routine Qpcrsh. */
/*     15 Oct 1994: Dependent columns fixed at their current value. */
/*     04 Dec 2000: R converted to row-wise storage. */
/*     02 Mar 2013: HvCalls added. */
/*     01 Nov 2014: normZ and ZHZmin added as arguments. */
/*     01 Nov 2014: Removed exit based on large condZ. */
/*     03 Nov 2014: neH, indH, locH, Hcol added as arguments. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --r__;
    --kbs;
    --y;
    --w;
    --v;
    --acol;
    --inda;
    --loca;
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
    eps0 = rw[2];
/* eps**(4/5)       IEEE DP  3.00e-13 */
    flmax = rw[8];
/* est. of the largest pos. real */
    *iexit = 0;
    *zhzmin = flmax;
    *zhzmax = 0.;
    *condz = 1.;
    *normz = 0.;
    if (*ns == 0) {
	return 0;
    }
    signobj = (doublereal) (*minimize);
    nbs = *m + *ns;
    sqstat = 0;
/*     ------------------------------------------------------------------ */
/*     Main loop to find a column of Z'HZ. */
/*     ------------------------------------------------------------------ */
    i__1 = *ns;
    for (js = 1; js <= i__1; ++js) {
/* -------------------------------------------------------------- */
/* Get the nonlinear elements of the column of Z. */
/* Find y such that B y = column jq. */
/* Scatter the nonlinear part of y into w. */
/* -------------------------------------------------------------- */
	jq = kbs[*m + js];
	s2unpack_(&jq, m, n, nea, &normaj, nloca, &loca[1], &inda[1], &acol[1]
		, &w[1]);
	s2bsol_(iexit, &c__1, m, &w[1], &y[1], &iw[1], leniw, &rw[1], lenrw);
	if (*iexit > 0) {
	    return 0;
	}
	normzj = dnormi_(m, &y[1], &c__1);
	*normz = max(normzj,*normz);
/* Computing MAX */
	d__1 = normzj / normaj;
	*condz = max(d__1,*condz);
	s2scatter_(nnh, m, &kbs[1], &c_b81, &y[1], &w[1]);
	if (jq <= *nnh) {
	    w[jq] = 1.;
	}
/* Set v = H w. */
	if (*nnh > 0) {
	    (*hprod)((U_fp)hprod1, nnh, neh, nloch, &loch[1], &indh[1], &hcol[
		    1], &w[1], &v[1], &sqstat, cu + 8, lencu, &iu[1], leniu, &
		    ru[1], lenru, cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw,
		     (ftnlen)8, (ftnlen)8);
	    ++(*hvcalls);
	    if (*minimize < 0) {
		dscal_(nnh, &signobj, &v[1], &c__1);
	    }
	}
/* -------------------------------------------------------------- */
/* Gather w = vBS and compute v = Z' w. */
/* Solve  B' vB = wB  and  form  wS = wS - S' vB. */
/* -------------------------------------------------------------- */
	s2gather_(nnh, &nbs, &kbs[1], &c_b149, &v[1], &w[1]);
	s2bsol_(iexit, &c__2, m, &w[1], &v[1], &iw[1], leniw, &rw[1], lenrw);
	if (*iexit > 0) {
	    return 0;
	}
	s2bprod_(&c__1, &eps0, n, ns, &kbs[*m + 1], nea, nloca, &loca[1], &
		inda[1], &acol[1], &c_b81, &v[1], m, &c_b149, &w[*m + 1], ns);
/* -------------------------------------------------------------- */
/* Store w(1:nS) in the jS-th row of R. */
/* R is NO LONGER SYMMETRIZED. */
/* -------------------------------------------------------------- */
	s6rrow_(&js, maxr, ns, lenr, &r__[1], &w[*m + 1], &lrlast);
	diag = r__[lrlast];
/* Computing MIN */
	d__1 = *zhzmin, d__2 = abs(diag);
	*zhzmin = min(d__1,d__2);
/* Computing MAX */
	d__1 = *zhzmax, d__2 = abs(diag);
	*zhzmax = max(d__1,d__2);
    }
    return 0;
} /* s5zhz_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5ZHZ */
/* Subroutine */ int s5zhzfac_(integer *iexit, integer *eigh, integer *
	typefac, integer *itn, integer *lenr, integer *m, integer *maxr, 
	integer *mbs, integer *nb, integer *ns, doublereal *hcondbnd, 
	doublereal *zhzmin, doublereal *zhzmax, doublereal *normz, integer *
	rankzhz, integer *hs, integer *kbs, integer *perm, doublereal *bl, 
	doublereal *bu, doublereal *blbs, doublereal *bubs, doublereal *x, 
	doublereal *xbs, doublereal *r__, doublereal *e, integer *iw, integer 
	*leniw)
{
    /* Format strings */
    static char fmt_9060[] = "(\002 Itn\002,i7,\002: Reduced Hessian appears"
	    " to be indefinite.\002,\002 dpiv, Hdmin = \002,1p,e9.2,\002,\002"
	    ",e9.2)";
    static char fmt_9000[] = "(\002 Itn\002,i7,\002: Reduced Hessian appears"
	    " to have \002,i6,\002 small eigenvalues.  PD tol = \002,1p,e9.2)";

    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static integer j, k;
    static doublereal s;
    static integer js;
    static doublereal eps;
    static char str[90];
    static integer jmax, kmax;
    static doublereal dpiv, hdmin;
    static integer nmodh, ksave, pivot;
    extern /* Subroutine */ int s6mchl_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, doublereal *), 
	    s6chol_(integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *);
    static integer inform__, nssave;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___237 = { 0, str, 0, fmt_9060, 90, 1 };
    static icilist io___238 = { 0, str, 0, fmt_9000, 90, 1 };


/*     ================================================================== */
/*     s5ZHZfac  factors the reduced Hessian Z'HZ. */

/*       On entry, */
/*         eigH     defines the class of QP Hessian:  positive definite, */
/*                  semidefinite or indefinite. */
/*         R(lenR)  holds the upper-triangular part of Z'HZ. */

/*       On exit, */
/*          R(lenR) holds the factor of the largest positive-definite */
/*                  subset of the rows and columns of ZHZ.  Excluded */
/*                  superbasics are made nonbasic at their current value. */

/*         iExit    Result */
/*         -----    ------ */
/*          -2      H  singular when it should be positive definite */
/*          -1      H  indefinite when it should be positive definite */
/*           0      positive-definite factor computed successfully */

/*     13 Oct 1992: First version based on Qpsol routine Qpcrsh. */
/*     15 Oct 1994: Dependent columns fixed at their current value. */
/*     02 Aug 2003: snPRNT adopted. */
/*     01 Nov 2014: ZHZmin and normZ added as arguments. */
/*     03 Nov 2014: Modified Cholesky added as an option. */
/*     03 Nov 2014: rankZHZ added to the argument list. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --r__;
    --e;
    --perm;
    --xbs;
    --bubs;
    --blbs;
    --kbs;
    --x;
    --bu;
    --bl;
    --hs;
    --iw;

    /* Function Body */
    *iexit = 0;
    *rankzhz = 0;
    eps = max(*zhzmax,1.) / *hcondbnd;
/* Computing MAX */
    d__1 = *zhzmax / *hcondbnd;
    hdmin = max(d__1,eps);
/*     eps    = max ( ZHZmax/normZ/normZ, one )/Hcondbnd */
/*     Hdmin  = max ( ZHZmax/normZ/normZ/Hcondbnd, eps ) */
    if (*eigh == 1) {
	pivot = 0;
    } else if (*eigh == -1 || *eigh == 0) {
	pivot = 1;
    } else {
/* Relax -- there are no other options. */
    }
    if (*typefac == 0) {
	s6chol_(&inform__, &pivot, maxr, ns, lenr, &r__[1], &hdmin, &dpiv, 
		rankzhz, &perm[1]);
    } else if (*typefac == 1) {
	s6mchl_(&inform__, &pivot, maxr, ns, lenr, &r__[1], &hdmin, &eps, &
		dpiv, rankzhz, &nmodh, &perm[1], &e[1]);
    } else {
/* Relax -- there are no other options. */
    }
/*     inform > 0 implies rankZHZ < nS. */
    if (pivot == 1) {
/*        ----------------------- */
/*        Apply any interchanges. */
/*        ----------------------- */
	i__1 = min(*rankzhz,*ns);
	for (j = 1; j <= i__1; ++j) {
	    jmax = perm[j];
	    if (jmax > j) {
		kmax = *m + jmax;
		k = *m + j;
		ksave = kbs[kmax];
		kbs[kmax] = kbs[k];
		kbs[k] = ksave;
		s = xbs[kmax];
		xbs[kmax] = xbs[k];
		xbs[k] = s;
		s = blbs[kmax];
		blbs[kmax] = blbs[k];
		blbs[k] = s;
		s = bubs[kmax];
		bubs[kmax] = bubs[k];
		bubs[k] = s;
	    }
	}
    }
    if (dpiv < hdmin) {
/*        --------------------------------------- */
/*        H is not positive definite. */
/*        rankZHZ < nS */
/*        --------------------------------------- */
	if (*eigh == 1) {
/*           H is meant to be positive definite. Exit */
	    *iexit = -2;
	    s_wsfi(&io___237);
	    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&dpiv, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&hdmin, (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    snprnt_(&c__21, str, &iw[1], leniw, (ftnlen)90);
	} else if (dpiv < -hdmin && *eigh != -1) {
/* H is indefinite but should be positive semidefinite. */
	    *iexit = -1;
	} else {
/*           H is singular or an indefinite H is possible */
/*           This is okay, but we have to set hs to match the part */
/*           of Rz that we have. */
	    if (*eigh == 0) {
		s_wsfi(&io___238);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		i__1 = *ns - *rankzhz;
		do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&hdmin, (ftnlen)sizeof(doublereal));
		e_wsfi();
		snprnt_(&c__21, str, &iw[1], leniw, (ftnlen)90);
	    }
	    nssave = *ns;
	    i__1 = nssave;
	    for (js = *rankzhz + 1; js <= i__1; ++js) {
		k = *m + js;
		j = kbs[k];
/* Make variable  j  nonbasic (it is already feasible). */
/* hs(j) = -1 means x(j) is strictly between its bounds. */
		if (x[j] <= bl[j]) {
		    x[j] = bl[j];
		    hs[j] = 0;
		} else if (x[j] >= bu[j]) {
		    x[j] = bu[j];
		    hs[j] = 1;
		} else {
		    hs[j] = -1;
		}
		if (bl[j] == bu[j]) {
		    hs[j] = 4;
		}
		--(*ns);
	    }
	    *ns = min(*ns,*rankzhz);
	}
/* rankZHZ < nS */
    }
    return 0;
} /* s5zhzfac_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5ZHZfac */
/* Subroutine */ int s5zp_(integer *iexit, integer *m, integer *mbs, integer *
	n, integer *nb, integer *ns, doublereal *eps0, doublereal *pnorm, 
	integer *nea, integer *nloca, integer *loca, integer *inda, 
	doublereal *acol, integer *kbs, doublereal *pbs, doublereal *y, 
	integer *iw, integer *leniw, doublereal *rw, integer *lenrw)
{
    /* System generated locals */
    doublereal d__1;

    /* Local variables */
    static integer nbs;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), s2bsol_(integer *, integer *, integer *, doublereal *,
	     doublereal *, integer *, integer *, doublereal *, integer *);
    extern doublereal dnormi_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int s2bprod_(integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, integer *);

/*     ================================================================== */
/*     s5Zp computes the free components of the search direction */
/*     p = Z pS, where pS is the search direction for the superbasics, */
/*     stored in  pBS(m+1:nBS) */

/*     On exit, the  free components of the search direction are stored */
/*     in pBS(1:nBS). The search direction for the basic variables is */
/*     stored in pBS(1),...,pBS(m). */

/*     20 Dec 2005: First version of s5Zp. */
/*     01 Dec 2012: Added pNorm <= 0.0 test. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --pbs;
    --kbs;
    --y;
    --acol;
    --inda;
    --loca;
    --iw;
    --rw;

    /* Function Body */
    *iexit = 0;
    nbs = *m + *ns;
    *pnorm = dnormi_(ns, &pbs[*m + 1], &c__1);
    if (*pnorm <= 0.) {
	*pnorm = 1.;
    }
/* First, compute  y = - S*pS and prepare to solve  B*pB = y */
/* for pB, the search direction for the basic variables. */
/* We first normalize y so the LU solver won't ignore */
/* too many "small" elements while computing pB. */
    d__1 = 1. / *pnorm;
    dscal_(ns, &d__1, &pbs[*m + 1], &c__1);
    s2bprod_(&c__0, eps0, n, ns, &kbs[*m + 1], nea, nloca, &loca[1], &inda[1],
	     &acol[1], &c_b81, &pbs[*m + 1], ns, &c_b5, &y[1], m);
/* Solve  B*pBS = y  and unnormalize all of pBS. */
    s2bsol_(iexit, &c__1, m, &y[1], &pbs[1], &iw[1], leniw, &rw[1], lenrw);
    if (*iexit != 0) {
	return 0;
    }
    dscal_(&nbs, pnorm, &pbs[1], &c__1);
    *pnorm = dnormi_(&nbs, &pbs[1], &c__1);
    return 0;
} /* s5zp_ */

