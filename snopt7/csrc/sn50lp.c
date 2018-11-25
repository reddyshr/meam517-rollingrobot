/* ../snopt7/src/sn50lp.f -- translated by f2c (version 20100827).
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
static integer c__23 = 23;
static integer c__21 = 21;
static integer c__11 = 11;
static integer c__31 = 31;
static integer c__2 = 2;
static integer c_n2 = -2;
static doublereal c_b69 = 1.;
static integer c__32 = 32;
static integer c__3 = 3;
static doublereal c_b159 = -1.;
static integer c__22 = 22;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*     File  sn50lp.f */

/*     s5LP     s5BSx    s5degen   s5eGrad  s5eInf  s5eReset */
/*     s5erc    s5FixS   s5FixX    s5getB   s5hs    s5Inf */
/*     s5LG     s5LPit   s5price   s5rc     s5setpi s5setx   s5step */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Subroutine */ int s5lp_(integer *iexit, integer *probtype, char *probtag, 
	logical *elastic, integer *suboptimize, S_fp lplog, logical *needlu, 
	logical *needx, integer *m, integer *n, integer *nb, integer *ndegen, 
	integer *itlp, integer *itlpmax, integer *itn, integer *emode, 
	integer *lvlobje, integer *printlevel, integer *minimize, integer *
	iobj, doublereal *scaleobj, doublereal *objadd, doublereal *toloptfp, 
	doublereal *toloptlp, doublereal *tolx, integer *ninf, doublereal *
	sinf, integer *elastics, integer *ninfe, doublereal *sinfe, 
	doublereal *wtinf, doublereal *pinorm, doublereal *rgnorm, integer *
	nea, integer *nloca, integer *loca, integer *inda, doublereal *acol, 
	integer *etype, integer *estate, integer *feastype, integer *hs, 
	integer *kbs, doublereal *bl, doublereal *bu, doublereal *blqp, 
	doublereal *buqp, doublereal *blbs, doublereal *bubs, doublereal *gbs,
	 doublereal *pi, doublereal *rc, integer *nrhs0, integer *nrhs, 
	doublereal *rhs, doublereal *scales, doublereal *x, doublereal *xbs, 
	doublereal *xfrozen, integer *iy, integer *iy1, doublereal *y, 
	doublereal *y1, char *cw, integer *lencw, integer *iw, integer *leniw,
	 doublereal *rw, integer *lenrw, ftnlen probtag_len, ftnlen cw_len)
{
    /* Format strings */
    static char fmt_1030[] = "(\002 Itn\002,i7,\002: Infeasible nonelastics."
	    "  Num =\002,i5,1p,\002   Sum of Infeasibilities =\002,e8.1)";
    static char fmt_8050[] = "(\002 Itn\002,i7,\002: Infeasible \002,a)";
    static char fmt_8060[] = "(\002 Itn\002,i7,\002: Elastic Phase 1 -- maki"
	    "ng\002,\002 nonelastic variables feasible\002)";
    static char fmt_1010[] = "(\002 Biggest dj =\002,1p,e11.3,\002 (variabl"
	    "e\002,i7,\002)\002,\002    norm rg =\002,e11.3,\002   norm pi "
	    "=\002,e11.3)";
    static char fmt_1020[] = "(\002 Norm rg =\002,1p,e11.3,\002   norm pi "
	    "=\002,e11.3)";
    static char fmt_1000[] = "(\002 ==> LU file has increased by a factor o"
	    "f\002,f6.1)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    extern /* Subroutine */ int s2unpack_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, doublereal *), s5ereset_(integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *);
    static logical feasible, increase, printlog, printsum;
    static doublereal rowerror;
    static logical checkfeas, firstfeas;
    static integer lurequest, jq, kp, ns;
    static logical justphase1;
    static integer neg, jbq, jbr;
    static doublereal djq;
    static integer mbs, nbs, nnh, jsq, jsr;
    static char str[115];
    extern /* Subroutine */ int s5hs_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *);
    static integer kfac, nfac, kchk;
    static doublereal bold;
    static integer leng;
    static logical newb;
    static integer klog;
    static logical gotg;
    static integer kprc, ksav, nfix[2];
    static logical gotr, luok;
    static doublereal step;
    extern /* Subroutine */ int s5inf_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal tolx0;
    static logical needf;
    extern /* Subroutine */ int dload_(integer *, doublereal *, doublereal *, 
	    integer *);
    static logical needv;
    static doublereal objlp, anorm;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer lumax;
    static logical newlu;
    static integer ksumm, nswap;
    static doublereal dummy[1];
    extern /* Subroutine */ int s2bfac_(integer *, integer *, logical *, 
	    logical *, logical *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, integer *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, integer *);
    static doublereal pivot;
    extern /* Subroutine */ int s5einf_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *), s1time_(integer *, integer *, 
	    integer *, integer *, doublereal *, integer *), s2bsol_(integer *,
	     integer *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, doublereal *, integer *);
    static integer lusiz0;
    extern /* Subroutine */ int s5setx_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *, integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *);
    static integer kdegen;
    static doublereal infbnd;
    static logical needpi, checkx;
    static doublereal featol;
    static logical lpdone;
    static doublereal weight, tolinc;
    static integer inform__, itnlim;
    static doublereal objprt;
    static integer itnfix, frozen;
    static doublereal djqprt;
    static integer nonopt, kprprt, typelu;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen), s5degen_(integer *, integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, integer *, doublereal *, 
	    integer *), s5egrad_(integer *, integer *, doublereal *, integer *
	    , integer *, doublereal *), s5price_(logical *, logical *, 
	    logical *, logical *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, doublereal *, doublereal *, integer *, doublereal *, integer *,
	     doublereal *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, integer *), s4ksave_(integer *, integer *, integer *
	    , integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, char *, 
	    integer *, integer *, integer *, ftnlen), s5setpi_(integer *, 
	    integer *, logical *, doublereal *, doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, integer *), s5lpitn_(integer *
	    , logical *, logical *, logical *, logical *, integer *, integer *
	    , integer *, integer *, integer *, integer *, integer *, integer *
	    , integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, doublereal *, integer *);
    static logical checkpi;
    extern /* Subroutine */ int s2trylu_(integer *, integer *, integer *, 
	    integer *, logical *, integer *, integer *, integer *, doublereal 
	    *, integer *);
    static doublereal signobj, maxtime;
    static logical optimal;
    static doublereal bgrowth, condzhz;
    static logical prtlvl10;

    /* Fortran I/O blocks */
    static icilist io___62 = { 0, str, 0, fmt_1030, 115, 1 };
    static icilist io___69 = { 0, str, 0, fmt_8050, 115, 1 };
    static icilist io___70 = { 0, str, 0, fmt_8060, 115, 1 };
    static icilist io___71 = { 0, str, 0, fmt_1010, 115, 1 };
    static icilist io___72 = { 0, str, 0, fmt_1020, 115, 1 };
    static icilist io___77 = { 0, str, 0, fmt_1000, 115, 1 };


/* ================================================================= */
/* s5LP   solves a linear program. */

/* The optimization can pass through the following phases: */

/* Phase 1               find a feasible point for all variables */

/* Elastic Phase 1       make the nonelastic variables feasible */
/*                       while allowing infeasible elastics */

/* Phase 2               minimize the objective */

/* Elastic Phase 2       minimize a composite objective while */
/*                       keeping the nonelastics feasible */

/*                       In this phase, lvlObjE means: */

/*          lvlObjE = 0  zero     weight on the infeasibilities */
/*                                (infeasibillities are ignored) */
/*                    1  finite   weight on the infeasibilities */
/*                    2  infinite weight on the infeasibilities */
/*                                (the objective is ignored) */

/* On entry: */
/* --------- */
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
/*                  3  FPE  feasible point for equalities only */
/*                  4  FPS  feasible point for QP subProblem */


/* On exit: */
/* -------- */
/*   iExit         Status */
/*   -----         ------ */
/*     -3          Too many iterations */
/*     -2          LP is unbounded */
/*     -1          Nonelastic variables are infeasible */
/*      0          LP solution found */
/*     >0          Fatal error */

/*  The array kBS is a permutation on the column indices. */
/*  kBS(1  :m )  holds the column indices of the basic variables. */
/*  Superbasics have been temporarily fixed at their current value. */

/* 30 Sep 1991: First version of s5LP based on Qpsol's lpcore. */
/* 20 Jul 1996: Slacks changed to be the row value. */
/* 06 Aug 1996: First Min Sum version. */
/* 14 Jul 1997: Thread-safe version. */
/* 24 Dec 1999: Suboptimization option added. */
/* 01 Aug 2003: snEXIT and snPRNT adopted. */
/* 24 Dec 2003: pi checked for NaN and Inf entries. */
/* 07 May 2006: s4ksave handles negative values of hs. */
/* 17 Jun 2008: Real workspace reorganized. */
/* 23 Oct 2010: pinorm initialized to one instead of zero. */
/* 26 May 2013: infBnd used to identify infinite bounds. */
/* 20 Sep 2014: bl and bu added for elastic mode. */
/* 17 Nov 2014: Time limit added. */
/* ================================================================= */
/*     ------------------------------------------------------------------ */
/* partial pricing in use */
/* partial pricing for LPs */
/* phase 1 dj tol for p.p. */
/* phase 2 dj tol for p.p. */
/* current optimality tol */
/* size of L0 */
/* size of initial  U */
/* size of current  L */
/* size of current  U */
/* xBS(kObj) is the obj. slack */
/* itns since last factorize */
/* number of LU mods */
/* (on/off) log     status */
/* (on/off) summary status */
/* # lines in log     file */
/* # lines in summary file */
/* >0 => Minor heading for iPrint */
/* >0 => Minor heading for iSumm */
/* Solve time */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --xbs;
    --pi;
    --gbs;
    --bubs;
    --blbs;
    --kbs;
    --feastype;
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
    --acol;
    --inda;
    --loca;
    --rhs;
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    infbnd = rw[70];
/* definition of an infinite bound */
    maxtime = rw[79];
/* max time allowed */
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
    nfac = iw[210];
/* # of LU factorizations */
    if (nfac > 0) {
	lusiz0 = iw[171] + iw[172];
	lumax = lusiz0 << 1;
    }
    nnh = 0;
/* local value of nnH */
    neg = 0;
    leng = 1;
    ns = 0;
/* local value of nS */
    nbs = *m;
    mbs = *m + 1;
/* simplex steps */
    prtlvl10 = *printlevel >= 10;
    printlog = *printlevel >= 1 && (*itlp % klog == 0 && *itlp != 0 || klog ==
	     1);
    printsum = *printlevel >= 1 && (*itlp % ksumm == 0 && *itlp != 0 || ksumm 
	    == 1);
    iw[218] = 0;
    iw[219] = 0;
    if (printlog) {
	iw[218] = 1;
    }
    if (printsum) {
	iw[219] = 1;
    }
    kprc = 0;
/* last section scanned in part. pricing */
    *iexit = 0;
    lurequest = 0;
    checkfeas = TRUE_;
/* Check that x is feasible. */
    checkpi = TRUE_;
/* Compute norm pi in s5setpi */
    feasible = FALSE_;
    gotr = FALSE_;
    firstfeas = FALSE_;
/* ----------------------------------------------------------------- */
/* s5LP operates in either ``Normal'' or ``Elastic'' mode. */
/* Everything is normal unless a weighted sum is being minimized or */
/* the constraints are infeasible. */
/* The logical Feasible refers to the nonelastic variables. */
/* wtInf  is the optional parameter Infeasibility Weight. */
/* ----------------------------------------------------------------- */
/* JustPhase1 = stop after phase1 (either normal or elastic) */
    justphase1 = *probtype == 0 || *probtype == 3 || *probtype == 4;
/* The phase 2 objective is F1 + wtInf*F2. */
    if (*elastic) {
	needf = *lvlobje != 2;
/* F1 required in phase 2 */
	needv = *lvlobje != 0;
/* F2 required in phase 2 */
    } else {
	needf = TRUE_;
	needv = FALSE_;
    }
    needpi = TRUE_;
    newlu = TRUE_;
    optimal = FALSE_;
    lpdone = FALSE_;
    condzhz = 0.;
    objlp = 0.;
    pivot = 0.;
    *rgnorm = 0.;
    step = 0.;
    *sinfe = 0.;
    *ninfe = 0;
    nonopt = -1;
    jq = 0;
    djq = 0.;
    jbq = 0;
/* x(jBq) is the incoming  BS */
    jbr = 0;
/* x(jBr) is the outgoing  BS */
    jsq = 0;
/* x(jSq) is the incoming SBS */
    jsr = 0;
/* x(jSr) is the outgoing SBS */
    kprprt = 0;
    signobj = (doublereal) (*minimize);
    typelu = 0;
    frozen = 0;
    dummy[0] = 0.;
    rw[184] = 100.;
/* Used only for LP partial pricing */
    rw[185] = 100.;

    iw[94] = iw[99];
    s5hs_(&c__0, nb, &blqp[1], &buqp[1], &hs[1], &x[1]);
    s5degen_(&inform__, &c__0, printlevel, nb, ninf, itn, &featol, tolx, &
	    tolinc, &hs[1], &blqp[1], &buqp[1], &x[1], &itnfix, nfix, &tolx0, 
	    &iw[1], leniw, &rw[1], lenrw);
/* !    ======================Start of main loop========================== */
/* +    do while (.not. LPdone  .and.  iExit .eq. 0) */
L100:
    if (! lpdone && *iexit == 0) {
/* ============================================================== */
/* Check the initial  x  and move it onto  ( A -I )*x = b. */
/* If NeedLU is true, this will require a basis factorization. */
/* ============================================================== */
/* If necessary,  factorize the basis  ( B = LU ) and set x. */
/* If NeedLU is false on entry to s5LP, the first call to s2Bfac */
/* will try to use existing factors. */
/* If NeedLU is true on entry to s5LP, an LU factorization of */
/* type typeLU is computed. */

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
/* LUrequest = 11  Big  dx   in setx or setpi */
/* LUrequest = 23  Infeasibility after refactorization */

/* In elastic mode, reset eState, blQP, buQP, blBS and buBS for */
/* the elastics. elastics gives the number of elastic variables. */
/* -------------------------------------------------------------- */
	firstfeas = FALSE_;
	if (lurequest > 0) {
	    *needlu = TRUE_;
	}
	if (*needx || *needlu) {
	    s2bfac_(iexit, &typelu, needlu, &newlu, &newb, iobj, itn, 
		    printlevel, &lurequest, m, &mbs, n, nb, &nnh, &ns, &nswap,
		     nea, nloca, &loca[1], &inda[1], &acol[1], &kbs[1], &hs[1]
		    , &blqp[1], &buqp[1], &blbs[1], &bubs[1], nrhs0, nrhs, &
		    rhs[1], &x[1], &xbs[1], &iy[1], &iy1[1], &y[1], &y1[1], &
		    iw[1], leniw, &rw[1], lenrw);
	    if (newlu) {
		lusiz0 = iw[171] + iw[172];
		lumax = lusiz0 << 1;
		if (prtlvl10) {
		    iw[223] = 1;
		}
	    }
	    if (*iexit > 0) {
		goto L100;
	    }
	    needpi = TRUE_;
/* Recalculate the pi's. */
	    *needx = FALSE_;
	    checkpi = TRUE_;
	    checkfeas = TRUE_;
	}
	optimal = FALSE_;
	nbs = *m + ns;
	*ninf = 0;
	*sinf = 0.;
	dload_(&nbs, &c_b5, &gbs[1], &c__1);
	if (checkfeas) {
/* In Phase 1 or just after a factorize, check that the basic */
/* and superbasic nonelastics are feasible. */
/* Any elastic variables are checked first */
	    if (*elastic) {
/*              Check that the elastic variables satisfy blQP and buQP. */
/*              Change eState if necessary. */
		s5ereset_(&nbs, nb, elastics, &featol, &infbnd, &etype[1], &
			estate[1], &kbs[1], &bl[1], &bu[1], &blqp[1], &buqp[1]
			, &blbs[1], &bubs[1], &x[1]);
	    }
/* FirstFeas = true means that we have just become feasible. */
/* FirstFeas is turned off once a step is taken. */
	    s5inf_(&nbs, &featol, &infbnd, ninf, sinf, &feastype[1], &blbs[1],
		     &bubs[1], &gbs[1], &xbs[1]);
	    if (*ninf > 0) {
/* The nonelastics are infeasible.  If necessary, switch */
/* back to the feasibility phase,  after refactorization. */
/* Print something if the basis has been refactorized. */
		if (feasible) {
		    s2trylu_(itn, &c__23, &ns, &lurequest, &luok, &typelu, &
			    iw[1], leniw, &rw[1], lenrw);
		    if (! luok) {
			*iexit = 11;
		    }
		    feasible = FALSE_;
		    goto L100;
		}
		if (prtlvl10 && iw[215] == 0) {
		    s_wsfi(&io___62);
		    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&(*ninf), (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&(*sinf), (ftnlen)sizeof(doublereal)
			    );
		    e_wsfi();
		    snprnt_(&c__21, str, &iw[1], leniw, (ftnlen)115);
		}
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
/* CheckFeas */
	if (*elastic) {
/* ----------------------------------------------------------- */
/* Find the sum of infeasibilities of the elastic variables. */
/* If  nInfE .ne. elastics, then there is at least one */
/* nonbasic elastic fixed at its current value. */
/* ----------------------------------------------------------- */
	    s5einf_(nb, &nbs, &estate[1], &kbs[1], &featol, ninfe, sinfe, &bl[
		    1], &bu[1], &x[1]);
	}
	objlp = 0.;
	if (*iobj > 0) {
	    objlp = xbs[iw[205]] * *scaleobj;
	}
	if (feasible && justphase1) {
/* The non-elastics are feasible.  Prepare to exit. */
	    djqprt = 0.;
	    dload_(m, &c_b5, &pi[1], &c__1);
	    *pinorm = 1.;
/* pinorm = max(norm(pi), 1.0) */
	} else {
	    if (feasible) {
/* -------------------------------------------------------- */
/* Feasible for the nonelastics. */
/* (Elastic = false means no elastics.) */
/* -------------------------------------------------------- */
		if (needf) {
		    if (*iobj != 0) {
			gbs[iw[205]] = signobj * *scaleobj;
		    }
		}
		if (*elastic && *elastics > 0 && needv) {
		    s5egrad_(nb, &nbs, wtinf, &estate[1], &kbs[1], &gbs[1]);
		}
	    }
	    if (needpi) {
/* -------------------------------------------------------- */
/* Compute pi, the multipliers for Ax - s = b. */
/* -------------------------------------------------------- */
		dcopy_(m, &gbs[1], &c__1, &y[1], &c__1);
		s5setpi_(&inform__, m, &checkpi, pinorm, &y[1], &pi[1], &iw[1]
			, leniw, &rw[1], lenrw);
		if (inform__ != 0) {
		    if (inform__ > 0) {
			*iexit = inform__;
		    } else {
/* pi is infinite or contains a NaN/Inf. */
			s2trylu_(itn, &c__11, &ns, &lurequest, &luok, &typelu,
				 &iw[1], leniw, &rw[1], lenrw);
			if (! luok) {
			    *iexit = 43;
			}
		    }
		    goto L100;
		}
		needpi = FALSE_;
	    }
/* =========================================================== */
/* Check for optimality. */
/* Find the reduced costs. */
/* =========================================================== */
	    if (feasible) {
		rw[186] = *toloptlp;
	    } else {
		rw[186] = *toloptfp;
	    }
	    kprprt = kprc;
	    jq = 0;
	    djqprt = djq;
	    djq = 0.;
	    gotg = FALSE_;
	    weight = 0.;
	    if (*elastic && feasible) {
		weight = *wtinf;
	    }
	    s5price_(elastic, &feasible, &increase, &gotg, suboptimize, itn, 
		    m, n, nb, &leng, &neg, &frozen, &nonopt, &weight, &
		    signobj, pinorm, &jq, &djq, &kprc, &rw[184], nea, nloca, &
		    loca[1], &inda[1], &acol[1], &etype[1], &hs[1], dummy, &
		    pi[1], &rc[1], &x[1], &xfrozen[1], &iw[1], leniw, &rw[1], 
		    lenrw);
	    optimal = nonopt == 0;
	}
/* Feasible and JustPhase1 */
	lpdone = optimal || feasible && justphase1;
	if (lpdone) {
/* ----------------------------------------------------------- */
/* Apparently we are optimal. */
/* See if any nonbasics have to be set back on their bounds. */
/* ----------------------------------------------------------- */
	    s5degen_(&inform__, &c__1, printlevel, nb, ninf, itn, &featol, 
		    tolx, &tolinc, &hs[1], &blqp[1], &buqp[1], &x[1], &itnfix,
		     nfix, &tolx0, &iw[1], leniw, &rw[1], lenrw);
	    lpdone = inform__ == 0;
	    if (lpdone) {
/* -------------------------------------------------------- */
/* So far so good.  Now check the row residuals. */
/* -------------------------------------------------------- */
		if (iw[215] > 0) {
		    s5setx_(&inform__, &c__1, itn, m, n, nb, &nbs, &rowerror, 
			    nea, nloca, &loca[1], &inda[1], &acol[1], &kbs[1],
			     &xbs[1], nrhs0, nrhs, &rhs[1], &x[1], &y[1], &y1[
			    1], &iw[1], leniw, &rw[1], lenrw);
		    lpdone = inform__ == 0;
		    lurequest = inform__;
		}
	    }
/* If x is not optimal, set  x  so that ( A  -I )*x = b */
/* and check feasibility. */
	    if (! lpdone) {
		*needx = TRUE_;
		goto L100;
	    }
	}
/*        ============================================================ */
/*        Print the details of this iteration. */
/*        ============================================================ */
/* LPdone */
	if (firstfeas && justphase1) {
/* Relax, we are about to exit without printing */
	} else {
	    objprt = 0.;
	    if (feasible) {
		if (needf) {
		    objprt = *objadd + objlp;
		}
		if (needv) {
		    objprt += signobj * *wtinf * *sinfe;
		}
	    }
	    (*lplog)(probtype, probtag, elastic, &gotr, &firstfeas, &feasible,
		     m, &mbs, &nnh, &ns, &jsq, &jbr, &jsr, &iw[220], &iw[221],
		     itn, itlp, &kprprt, lvlobje, &pivot, &step, ninf, sinf, 
		    ninfe, sinfe, wtinf, &nonopt, &objprt, &condzhz, &djqprt, 
		    rgnorm, &kbs[1], &xbs[1], &iw[1], leniw, (ftnlen)20);
	}
	jbq = 0;
	jbr = 0;
	jsq = 0;
	jsr = 0;
	kprprt = 0;
	if (lpdone) {
/* ----------------------------------------------------------- */
/* Optimal  .or. (Feasible .and. JustPhase1) */
/* ----------------------------------------------------------- */
	    if (*ninf > 0) {
/* No feasible point. */
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
			s_wsfi(&io___69);
			do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer)
				);
			do_fio(&c__1, probtag, (ftnlen)20);
			e_wsfi();
			snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)115);
			s_wsfi(&io___70);
			do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer)
				);
			e_wsfi();
			snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)115);
			iw[223] = 1;
			iw[225] = 1;
		    }
		    *elastic = TRUE_;
		    checkfeas = TRUE_;
/* call s5eReset */
		    lpdone = FALSE_;
		    needf = *lvlobje != 2;
/* Use F1 in e-phase 2 */
		    needv = *lvlobje != 0;
/* Use F2 in e-phase 2 */
		    needpi = TRUE_;
		    djq = 0.;
		    step = 0.;
		}
		goto L100;
	    }
	    if (prtlvl10 && ! justphase1) {
		if (jq != 0) {
		    djq = signobj * djq;
		    if (klog == 1) {
			s_wsfi(&io___71);
			do_fio(&c__1, (char *)&djq, (ftnlen)sizeof(doublereal)
				);
			do_fio(&c__1, (char *)&jq, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&(*rgnorm), (ftnlen)sizeof(
				doublereal));
			do_fio(&c__1, (char *)&(*pinorm), (ftnlen)sizeof(
				doublereal));
			e_wsfi();
			snprnt_(&c__31, str, &iw[1], leniw, (ftnlen)115);
		    }
		} else {
		    if (klog == 1) {
			s_wsfi(&io___72);
			do_fio(&c__1, (char *)&(*rgnorm), (ftnlen)sizeof(
				doublereal));
			do_fio(&c__1, (char *)&(*pinorm), (ftnlen)sizeof(
				doublereal));
			e_wsfi();
			snprnt_(&c__31, str, &iw[1], leniw, (ftnlen)115);
		    }
		}
	    }
	} else {
/* ----------------------------------------------------------- */
/* Do another LP iteration. */
/* A nonbasic has been selected to become superbasic. */
/* Compute the vector y such that B y = column jq. */
/* ----------------------------------------------------------- */
/* Unpack column jq into  y1  and solve  B*y = y1. */
/* The solve computes  y1  such that  L*y1 = ajq. */
/* It is used below to modify L and U in s5LPit. */
	    s2unpack_(&jq, m, n, nea, &anorm, nloca, &loca[1], &inda[1], &
		    acol[1], &y1[1]);
	    s2bsol_(iexit, &c__1, m, &y1[1], &y[1], &iw[1], leniw, &rw[1], 
		    lenrw);
	    if (*iexit > 0) {
		return 0;
	    }
/* =========================================================== */
/* Take a simplex step.  A variable will become nonbasic */
/* at the new x. */
/* =========================================================== */
	    if (*itn >= itnlim || *itlp >= *itlpmax) {
		*iexit = -3;
/* Excess iterations */
		goto L100;
	    }
	    ++(*itlp);
	    ++(*itn);
	    ++iw[215];
	    newlu = FALSE_;
	    checkpi = FALSE_;
/* Decide if we want to print something this iteration. */
	    printlog = *printlevel >= 1 && *itlp % klog == 0;
	    printsum = *printlevel >= 1 && *itlp % ksumm == 0;
	    iw[218] = 0;
	    iw[219] = 0;
	    if (printlog) {
		iw[218] = 1;
	    }
	    if (printsum) {
		iw[219] = 1;
	    }
/* ----------------------------------------------------------- */
/* Take a simplex step. */
/* The new x will still be at a vertex (possibly temporary) */
/* Check for unboundedness (inform = 1). */
/* ----------------------------------------------------------- */
	    i__1 = *m + 1;
	    s5lpitn_(&inform__, &feasible, &increase, &needpi, elastic, &i__1,
		     m, nb, ndegen, elastics, &lurequest, &kp, &jbq, &jsq, &
		    jbr, &jsr, &jq, &featol, &pivot, &step, &tolinc, &etype[1]
		    , &estate[1], &feastype[1], &hs[1], &kbs[1], &bl[1], &bu[
		    1], &blqp[1], &buqp[1], &blbs[1], &bubs[1], &x[1], &xbs[1]
		    , &y[1], &y1[1], &iw[1], leniw, &rw[1], lenrw);
	    if (inform__ == 1) {
		*iexit = -2;
/* Unbounded direction */
		goto L100;
	    }
/* Increment featol every iteration. */
	    featol += tolinc;
/* =========================================================== */
/* Test for error condition and/or frequency interrupts. */
/* =========================================================== */
/* (1) Save a basis map (frequency controlled). */
/* (2) Every kdegen iterations, reset featol and move nonbasic */
/*    variables onto their bounds if they are very close. */
/* (3) Refactorize the basis if it has been modified */
/*    too many times. */
/* (4) Update the LU factors of the basis if requested. */
/* (5) Check row error (frequency controlled). */
	    if (*itn % ksav == 0) {
		s4ksave_(minimize, m, n, nb, &ns, &mbs, itn, ninf, sinf, &
			objlp, &kbs[1], &hs[1], &scales[1], &blqp[1], &buqp[1]
			, &x[1], &xbs[1], cw + 8, lencw, &iw[1], leniw, (
			ftnlen)8);
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
			s_wsfi(&io___77);
			do_fio(&c__1, (char *)&bgrowth, (ftnlen)sizeof(
				doublereal));
			e_wsfi();
			snprnt_(&c__21, str, &iw[1], leniw, (ftnlen)115);
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
		    typelu = 0;
		}
	    }
/* Check time limit */
	    if (maxtime > 0. && *itn % 50 == 0) {
		s1time_(&c_n2, &c__0, &iw[1], leniw, &rw[1], lenrw);
		s1time_(&c__2, &c__0, &iw[1], leniw, &rw[1], lenrw);
		if (rw[462] > maxtime) {
		    *iexit = 34;
/* time limit reached */
		}
	    }
	}
/* not Optimal */
	goto L100;
/* +    end while */
    }
/*     ======================end of main loop============================ */

    s5hs_(&c__1, nb, &blqp[1], &buqp[1], &hs[1], &x[1]);
    return 0;
} /* s5lp_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5LP */
/* Subroutine */ int s5bsx_(integer *task, integer *nbs, integer *nb, integer 
	*kbs, doublereal *x, doublereal *xbs)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j, k;

/*     ================================================================= */
/*     s5BSx   copies free variables from  xBS  into  x  or vice versa, */
/*             depending on whether  Task is 'xBS to x' or 'x to xBS'. */

/*     07 Nov 1991: First version based on Minos routine m5bsx. */
/*     21 Aug 1999: Current version of s5BSx. */
/*     ================================================================= */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --xbs;
    --kbs;
    --x;

    /* Function Body */
    if (*task == 1) {
	i__1 = *nbs;
	for (k = 1; k <= i__1; ++k) {
	    j = kbs[k];
	    x[j] = xbs[k];
	}
    } else if (*task == 0) {
	i__1 = *nbs;
	for (k = 1; k <= i__1; ++k) {
	    j = kbs[k];
	    xbs[k] = x[j];
	}
    }
    return 0;
} /* s5bsx_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5BSx */
/* Subroutine */ int s5degen_(integer *inform__, integer *task, integer *
	printlevel, integer *nb, integer *ninf, integer *itn, doublereal *
	featol, doublereal *tolx, doublereal *tolinc, integer *hs, doublereal 
	*bl, doublereal *bu, doublereal *x, integer *itnfix, integer *nfix, 
	doublereal *tolx0, integer *iw, integer *leniw, doublereal *rw, 
	integer *lenrw)
{
    /* Format strings */
    static char fmt_1000[] = "(\002 Itn\002,i7,\002: Basics recomputed aft"
	    "er \002,i7,\002  nonbasics set on bound\002)";

    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static integer j;
    static doublereal b1, b2, d1, d2;
    static char str[80];
    static doublereal eps1, tolz, tolx1;
    static integer kdegen, maxfix;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___92 = { 0, str, 0, fmt_1000, 80, 1 };


/*     ================================================================== */
/*     s5degen performs most of the tasks associated with degeneracy. */
/*     The degeneracy-resolving strategy operates in the following way. */

/*     Over a cycle of iterations, the feasibility tolerance featol */
/*     increases slightly (from tolx0 to tolx1 in steps of tolinc). */
/*     This ensures that all steps taken will be positive. */

/*     After kdegen consecutive iterations, nonbasic variables within */
/*     featol of their bounds are set exactly on their bounds and the */
/*     basic variables are recomputed to satisfy ( A  -I )*x = b. */
/*     featol is then reduced to tolx0 for the next cycle of iterations. */


/*     If Task = Init, s5degen initializes the parameters: */

/*     featol  is the current feasibility tolerance. */
/*     tolx0   is the minimum feasibility tolerance. */
/*     tolx1   is the maximum feasibility tolerance. */
/*     tolinc  is the increment to featol. */
/*     kdegen  is the expand frequency (specified by the user). */
/*             it is the frequency of resetting featol to tolx0. */
/*     nDegen  counts the number of degenerate steps (not used here, but */
/*             incremented by s5step). */
/*     itnfix  is the last iteration at which an 'Optimal' or 'Cycle' */
/*             entry set nonbasics onto their bound. */
/*     nfix(j) counts the number of times an 'Optimal' entry has */
/*             set nonbasics onto their bound, */
/*             where j=1 if infeasible, j=2 if feasible. */

/*     tolx0 and tolx1 are both close to the feasibility tolerance tolx */
/*     specified by the user.  (They must both be less than tolx.) */


/*     If Task = Cycle,  s5degen has been called after a cycle of */
/*     kdegen iterations.  Nonbasic x(j)s are examined to see if any are */
/*     off their bounds by an amount approaching featol.  inform returns */
/*     how many.  Deviations as small as tolz (e.g. 1.0d-11) are not */
/*     counted. If inform is positive, the basic variables are */
/*     recomputed.  It is assumed that s5LP or s5QP will then continue */
/*     iterations. */

/*     itnfix, nfix, tolx0 could be treated as SAVED variables. */
/*     They are included as arguments to prevent threading problems in a */
/*     multiprocessing environment. */

/*     If Task = Optml,  s5degen is being called after a subproblem */
/*     has been judged optimal, infeasible or unbounded. */
/*     Nonbasic x(j)s are examined as above. */

/*     07 Nov 1991: First version based on Minos routine m5dgen. */
/*     12 Jul 1997: Thread-safe version. */
/*     11 May 2001: Removed summary file printing. */
/*     01 Aug 2003: snPRNT adopted. */
/*     01 Aug 2003: Current version of s5degen. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --x;
    --bu;
    --bl;
    --hs;
    --nfix;
    --iw;
    --rw;

    /* Function Body */
    eps1 = rw[3];
/* eps**(2/3) */
    kdegen = iw[63];
/* max. expansions of featol */
    *inform__ = 0;
    if (*task == 0) {
/*        Task = Initialize. */
/*        Initialize at the start of each major iteration. */
/*        kdegen is the expand frequency      and */
/*        tolx   is the feasibility tolerance */
/*        (specified by the user).  They are not changed. */
/*        nDegen counts the total number of degenerate steps, and is */
/*        initialized by s5solve or s8solve. */
	*itnfix = 0;
	nfix[1] = 0;
	nfix[2] = 0;
	*tolx0 = *tolx * .5;
	tolx1 = *tolx * .99;
	if (kdegen < 99999999) {
	    *tolinc = (tolx1 - *tolx0) / kdegen;
	} else {
	    *tolinc = 0.;
	}
	*featol = *tolx0;
    } else if (*task == 2 || *task == 1) {
/*        --------------------------------------------------------------- */
/*        Task = 'E'nd of cycle or 'O'ptimal. */
/*        initialize local variables maxfix and tolz. */
/*        --------------------------------------------------------------- */
	maxfix = 2;
	tolz = eps1;
	if (*task == 1) {
/*           Task = Optimal. */
/*           Return with inform = 0 if the last call was at the */
/*           same itn, or if there have already been maxfix calls */
/*           with the same state of feasibility. */
	    if (*itnfix == *itn) {
		return 0;
	    }
	    if (*ninf > 0) {
		j = 1;
	    } else {
		j = 2;
	    }
	    if (nfix[j] >= maxfix) {
		return 0;
	    }
	    ++nfix[j];
	}
/*        Set nonbasics on their nearest bound if they are within */
/*        the current featol of that bound. */
	*itnfix = *itn;
	i__1 = *nb;
	for (j = 1; j <= i__1; ++j) {
	    if (hs[j] <= 1 || hs[j] == 4) {
		b1 = bl[j];
		b2 = bu[j];
		d1 = (d__1 = x[j] - b1, abs(d__1));
		d2 = (d__1 = x[j] - b2, abs(d__1));
		if (d1 > d2) {
		    b1 = b2;
		    d1 = d2;
		}
		if (d1 <= *featol) {
		    if (d1 > tolz) {
			++(*inform__);
		    }
		    x[j] = b1;
		}
	    }
	}
/*        Reset featol to its minimum value. */
	*featol = *tolx0;
	if (*inform__ > 0) {
/*           The basic variables will be reset. */
	    if (*printlevel >= 10) {
		s_wsfi(&io___92);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&(*inform__), (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__21, str, &iw[1], leniw, (ftnlen)80);
	    }
	}
    }
    return 0;
} /* s5degen_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5degen */
/* Subroutine */ int s5egrad_(integer *nb, integer *nbs, doublereal *wtinf, 
	integer *estate, integer *kbs, doublereal *gbs)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j, k, estatej;

/*     ================================================================== */
/*     s5eGrad  is called when elastic variables are allowed to violate */
/*     their true bounds bl and bu.  It updates the gradient gBS */
/*     to include the gradient of the elastic penalty term. */

/*     On exit, */
/*       gBS(nBS)    is the rhs for the equations for pi. */

/*     08 Oct 1996: First version of s5Egrd. */
/*     21 Apr 1999: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --estate;
    --gbs;
    --kbs;

    /* Function Body */
    i__1 = *nbs;
    for (k = 1; k <= i__1; ++k) {
	j = kbs[k];
	estatej = estate[j];
	if (estatej == 0) {
/* Relax */
	} else if (estatej == 1) {
	    gbs[k] -= *wtinf;
	} else {
/*  (eStatej .eq. 2) */
	    gbs[k] += *wtinf;
	}
    }
    return 0;
} /* s5egrad_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5eGrad */
/* Subroutine */ int s5einf_(integer *nb, integer *nbs, integer *estate, 
	integer *kbs, doublereal *featol, integer *ninfe, doublereal *sinfe, 
	doublereal *bl, doublereal *bu, doublereal *x)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j, k;
    static doublereal res;
    static integer numinf;
    static doublereal suminf;
    static integer estatej;

/*     ================================================================== */
/*     s5eInf  computes the sum and number of elastic infeasibilities. */

/*     On exit, */
/*     nInfE is the number of elastics allowed to go infeasible. */

/*     20 Aug 1996: First version of s5eInf. */
/*     15 Oct 2014: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/* ----------------------------------------------------------------- */
/* Find the number and sum of the elastic infeasibilities. */
/* ----------------------------------------------------------------- */
    /* Parameter adjustments */
    --x;
    --bu;
    --bl;
    --estate;
    --kbs;

    /* Function Body */
    numinf = 0;
    suminf = 0.;
    i__1 = *nbs;
    for (k = 1; k <= i__1; ++k) {
	j = kbs[k];
	estatej = estate[j];
	if (estatej == 0) {
/* Relax, this variable is not elastic */
	} else if (estatej == 1) {
	    ++numinf;
	    res = bl[j] - x[j];
	    if (res > *featol) {
		suminf += res;
	    }
	} else if (estatej == 2) {
	    ++numinf;
	    res = x[j] - bu[j];
	    if (res > *featol) {
		suminf += res;
	    }
	}
    }
    *sinfe = suminf;
    *ninfe = numinf;
    return 0;
} /* s5einf_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5eInf */
/* Subroutine */ int s5ereset_(integer *nbs, integer *nb, integer *elastics, 
	doublereal *featol, doublereal *infbnd, integer *etype, integer *
	estate, integer *kbs, doublereal *bl, doublereal *bu, doublereal *
	blqp, doublereal *buqp, doublereal *blbs, doublereal *bubs, 
	doublereal *x)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j, k;
    static doublereal xj, blj, buj, violl, violu;
    static integer etypej, estatej;

/* In */
/* Out */
/* InOut */
/* In */
/*     ================================================================== */
/*     s5eReset is called after a factorize in elastic mode. */

/*     The upper and lower bounds on the elastic variables are reset for */
/*     Elastic Phase 1 and 2. */

/*     x = x(feasible) + w - v */

/*     s5eReset is called in elastic mode after a factorize. */
/*     The bounds blBS and buBS are redefined for the basic elastics. */

/*     In the case of a singular basis, variables associated with */
/*     dependent basic columns are fixed inside their bounds. If one of */
/*     these basics is elastic, it will be fixed inside its working */
/*     bounds but outside bl or bu.  These variables are not counted in */
/*     elastic infeasibilities because they will be made basic in */
/*     s5price. */

/*     23 Aug 1996: First version of e5Reset (s5Eset). */
/*     12 Jun 2000: Elastic mode cleaned up. */
/*     15 Oct 2014: blQP and buQP added to handle fixed elastics. */
/*     ================================================================== */
/* InOut */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --bubs;
    --blbs;
    --kbs;
    --x;
    --buqp;
    --blqp;
    --bu;
    --bl;
    --estate;
    --etype;

    /* Function Body */
    *elastics = 0;
    i__1 = *nbs;
    for (k = 1; k <= i__1; ++k) {
	j = kbs[k];
	etypej = etype[j];
	estatej = estate[j];
	blj = bl[j];
	buj = bu[j];
	xj = x[j];
	if (estatej > 0) {
/*           ------------------------------------------------------------ */
/*           x(j) is elastic. */

/*           EXPAND may give a slightly feasible elastic variable, i.e., */
/*               x(j) le bl(j) + featol  or x(j) ge bu(j) - featol. */
/*           x(j) is kept elastic regardless. */
/*           ------------------------------------------------------------ */
	    if (estatej == 1) {
/*              eStatej predicts that xj violates its true lower bound, */
/*              i.e., eTypej must be either 1 or 3. */

/*              In this case  x(j) = bl(j) - v(j), with  v(j) ge -featol */
		violl = blj - xj;
		if (violl > -(*featol)) {
/*                 x(j) violates its true lower bound as predicted. */
/*                 Keep eState and the QP bounds as they are. */
		    ++(*elastics);
		} else {
/*                 x(j) satisfies its true lower bound. */
/*                 Reset eState and the QP bounds to their true values. */
		    estatej = 0;
		    blqp[j] = blj;
		    buqp[j] = buj;
/*                 x(j) is made nonelastic unless it violates its */
/*                 (elastic) bu(j) by more than +featol. */
		    if (etypej == 3) {

/*                    x(j) = bu(j) + w(j),  with  w(j) ge +featol */
			violu = xj - buj;
			if (violu > *featol) {
			    ++(*elastics);
			    estatej = 2;
			    blqp[j] = buj;
			    buqp[j] = *infbnd;
			}
		    }
		    blbs[k] = blqp[j];
		    bubs[k] = buqp[j];
		}
	    } else if (estatej == 2) {
/*              x(j) is predicted to violate its true upper bound. */
/*              eTypej must be either 2 or 3. */
/*              x(j) = bu(j) + w(j),  with  w(j) ge -featol */
		violu = xj - buj;
		if (violu > -(*featol)) {
/*                 x(j) violates its true upper bound as predicted. */
/*                 Keep eState and the QP bounds as they are. */
		    ++(*elastics);
		} else {
/*                 x(j) satisfies its true upper bound. */
/*                 Reset the QP bounds to their true values. */
		    estatej = 0;
		    blqp[j] = blj;
		    buqp[j] = buj;
/*                 x(j) is made nonelastic unless it violates its */
/*                 (elastic) bl(j) by more than +featol. */
		    if (etypej == 3) {
			violl = blj - xj;
			if (violl > *featol) {
			    ++(*elastics);
			    estatej = 1;
			    blqp[j] = -(*infbnd);
			    buqp[j] = blj;
			}
		    }
		    blbs[k] = blqp[j];
		    bubs[k] = buqp[j];
		}
	    }
	} else {
/*           ----------------------------------------------------------- */
/*           Check if any basic or superbasic elastic variables that are */
/*           predicted to be feasible violate their elastic bounds. */
/*           ----------------------------------------------------------- */
/* eStatej = 0 */
	    if (etypej == 1) {
/*              If xj violates its lower bound by more than featol, */
/*              make it elastic. */
		violl = blj - xj;
		if (violl > *featol) {
		    ++(*elastics);
		    estatej = 1;
		    blqp[j] = -(*infbnd);
		    buqp[j] = blj;
		    blbs[k] = blqp[j];
		    bubs[k] = buqp[j];
		}
	    } else if (etypej == 2) {
		violu = xj - buj;
		if (violu > *featol) {
		    ++(*elastics);
		    estatej = 2;
		    blqp[j] = buj;
		    buqp[j] = *infbnd;
		    blbs[k] = blqp[j];
		    bubs[k] = buqp[j];
		}
	    } else if (etypej == 3) {
		violl = blj - xj;
		violu = xj - buj;
		if (violl > *featol) {
		    ++(*elastics);
		    estatej = 1;
		    blqp[j] = -(*infbnd);
		    buqp[j] = blj;
		    blbs[k] = blqp[j];
		    bubs[k] = buqp[j];
		} else if (violu > *featol) {
		    ++(*elastics);
		    estatej = 2;
		    blqp[j] = buj;
		    buqp[j] = *infbnd;
		    blbs[k] = blqp[j];
		    bubs[k] = buqp[j];
		}
	    }
	}
/*        Mark the basic and superbasic elastics. */
	estate[j] = estatej;
    }
    return 0;
} /* s5ereset_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5eReset */
/* Subroutine */ int s5erc_(integer *j1, integer *j2, logical *gotg, integer *
	m, integer *n, integer *leng, integer *neg, doublereal *signobj, 
	integer *nea, integer *nloca, integer *loca, integer *inda, 
	doublereal *acol, integer *etype, integer *hs, doublereal *g, 
	doublereal *pi, doublereal *rc)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, l;
    static doublereal dj;

/*     ================================================================== */
/*     s5erc  computes reduced costs rc(j) in the range j = j1 to j2 */
/*     for fixed nonbasic columns that have one or more elastic bounds. */
/*     It is called by s5price. */

/*     07 Feb 1998: First version based on s5rc. */
/*     20 Sep 2014: Compute reduced costs for  hs(j) = -1. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --pi;
    --rc;
    --hs;
    --etype;
    --g;
    --acol;
    --inda;
    --loca;

    /* Function Body */
    i__1 = *j2;
    for (j = *j1; j <= i__1; ++j) {
	if ((hs[j] == 4 || hs[j] == -1) && etype[j] > 0) {
	    dj = 0.;
	    i__2 = loca[j + 1] - 1;
	    for (l = loca[j]; l <= i__2; ++l) {
		i__ = inda[l];
		dj += pi[i__] * acol[l];
	    }
	    rc[j] = -dj;
	}
    }
/* Include the nonlinear gradient term if present. */
    if (*gotg) {
	if (*j1 <= *neg) {
	    i__1 = min(*j2,*neg);
	    for (j = *j1; j <= i__1; ++j) {
		if ((hs[j] == 4 || hs[j] == -1) && etype[j] > 0) {
		    rc[j] += *signobj * g[j];
		}
	    }
	}
    }
    return 0;
} /* s5erc_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5erc */
/* Subroutine */ int s5fixs_(integer *task, integer *m, integer *maxs, 
	integer *mbs, integer *n, integer *nb, integer *ns, integer *hs, 
	integer *kbs, doublereal *bl, doublereal *bu, doublereal *blbs, 
	doublereal *bubs, doublereal *x, doublereal *xbs)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j, k;

/*     ================================================================== */
/*     s5fixS   concerns temporary bounds on superbasic variables. */
/*     If Task = Fix,  s5fixS sets hs(j) = -1, 0, 1 or 4 for certain */
/*     superbasic variables. */

/*     If Task = Free, s5fixS changes -1 values to hs(j) = 2. */

/*     30 May 1995: First version of s5fixS. */
/*     12 Jul 2001: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --xbs;
    --bubs;
    --blbs;
    --kbs;
    --x;
    --bu;
    --bl;
    --hs;

    /* Function Body */
    if (*task == 0) {
	if (*ns > 0) {
/*           ------------------------------------------------------------ */
/*           Change superbasic hs(j) to be temporarily fixed. */
/*           ------------------------------------------------------------ */
	    *ns = 0;
	    i__1 = *nb;
	    for (j = 1; j <= i__1; ++j) {
		if (hs[j] == 2) {
		    if (bl[j] == bu[j]) {
			hs[j] = 4;
		    } else if (x[j] <= bl[j]) {
			hs[j] = 0;
		    } else if (x[j] >= bu[j]) {
			hs[j] = 1;
		    } else {
			hs[j] = -1;
		    }
		}
	    }
	}
    } else if (*task == 1) {
/*        --------------------------------------------------------------- */
/*        Free the temporarily fixed structurals. */
/*        Load the superbasic variables/bounds into xBS, blBS, buBS. */
/*        --------------------------------------------------------------- */
	j = 1;
/* +       while (j .le. n  .and.  nS .lt. maxS) do */
L100:
	if (j <= *n && *ns < *maxs) {
	    if (hs[j] == -1) {
		++(*ns);
		k = *m + *ns;
		hs[j] = 2;
		xbs[k] = x[j];
		blbs[k] = bl[j];
		bubs[k] = bu[j];
		kbs[k] = j;
	    }
	    ++j;
	    goto L100;
/* +       end while */
	}
    }
    return 0;
} /* s5fixs_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5fixS */
/* Subroutine */ int s5fixx_(integer *task, integer *b1, integer *b2, 
	doublereal *tolx, integer *hs, doublereal *bl, doublereal *bu, 
	doublereal *x)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static integer j;
    static doublereal xj;

/*     ================================================================== */
/*     s5FixX  ensures that variables satisfy their simple bounds. */

/*     If Task = xBound, variables x(b1) through x(b2) are made to */
/*                       satisfy their bounds. */
/*     If Task = xMove , variables x(b1) through x(b2) are made to */
/*                       satisfy their bounds. In addition, any nonbasic */
/*                       variable close to its bound is moved onto it. */

/*     29 Apr 1999: First version of s5FixX. */
/*     29 Apr 1999: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --x;
    --bu;
    --bl;
    --hs;

    /* Function Body */
    if (*task == 0) {
	i__1 = *b2;
	for (j = *b1; j <= i__1; ++j) {
/* Computing MAX */
	    d__1 = x[j], d__2 = bl[j];
	    x[j] = max(d__1,d__2);
/* Computing MIN */
	    d__1 = x[j], d__2 = bu[j];
	    x[j] = min(d__1,d__2);
	}
    } else if (*task == 1) {
	i__1 = *b2;
	for (j = *b1; j <= i__1; ++j) {
/* Computing MAX */
	    d__1 = x[j], d__2 = bl[j];
	    xj = max(d__1,d__2);
/* Computing MIN */
	    d__1 = xj, d__2 = bu[j];
	    xj = min(d__1,d__2);
	    if (hs[j] <= 1) {
		if (xj <= bl[j] + *tolx) {
		    xj = bl[j];
		}
		if (xj >= bu[j] - *tolx) {
		    xj = bu[j];
		}
	    }
	    x[j] = xj;
	}
    }
    return 0;
} /* s5fixx_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5FixX */
/* Subroutine */ int s5getb_(integer *iexit, integer *start, S_fp lplog, 
	logical *needb, integer *m, integer *maxs, integer *mbs, integer *n, 
	integer *nb, integer *nncon, integer *nnjac, integer *nnobj, integer *
	nnames, integer *ns, integer *itqp, integer *itqpmax, integer *itn, 
	integer *ndegen, integer *numlc, integer *numliq, doublereal *
	toloptfp, doublereal *toloptqp, doublereal *tolx, integer *ninf, 
	doublereal *sinf, doublereal *wtinf, integer *iobj, doublereal *
	scaleobj, doublereal *pinorm, doublereal *rgnorm, integer *nea, 
	integer *nloca, integer *loca, integer *inda, doublereal *acol, 
	integer *etype, integer *estate, integer *rowtypes, integer *feastype,
	 integer *hs, integer *kbs, char *names, doublereal *bl, doublereal *
	bu, doublereal *blqp, doublereal *buqp, doublereal *blbs, doublereal *
	bubs, doublereal *blsave, doublereal *busave, doublereal *gbs, 
	doublereal *pi, doublereal *rc, integer *nrhs0, integer *nrhs, 
	doublereal *rhs, doublereal *scales, integer *lenx0, integer *nx0, 
	doublereal *x0, doublereal *x, doublereal *xbs, integer *iy, integer *
	iy1, doublereal *y, doublereal *y1, doublereal *y2, char *cw, integer 
	*lencw, integer *iw, integer *leniw, doublereal *rw, integer *lenrw, 
	ftnlen names_len, ftnlen cw_len)
{
    /* Initialized data */

    static char line[4] = "----";

    /* Format strings */
    static char fmt_1332[] = "(1x,28a4)";
    static char fmt_1316[] = "(1x,16a4)";
    static char fmt_2100[] = "(\002 Itn\002,i7,\002: Making linear equality "
	    "rows feasible\002)";
    static char fmt_2200[] = "(\002 Itn\002,i7,\002: Feasible linear equalit"
	    "y rows\002)";
    static char fmt_2300[] = "(\002 Itn\002,i7,\002: Infeasible linear equal"
	    "ity rows\002)";
    static char fmt_2400[] = "(\002 Itn\002,i7,\002: Making all linear rows "
	    "feasible\002)";
    static char fmt_2410[] = "(\002 Itn\002,i7,\002: Making the linear rows "
	    "feasible\002)";

    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer emodeleq, elastics, minimize;
    extern /* Subroutine */ int s4checkhs_(integer *, integer *, integer *, 
	    integer *, integer *, logical *, integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, integer *);
    static integer j, mjrprtlvl, mnrprtlvl, lc1;
    extern /* Subroutine */ int s2getscales_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, doublereal *, integer *);
    static integer nnl;
    static char str[120];
    extern /* Subroutine */ int s5lg_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *), s5hs_(integer *, integer *, doublereal *
	    , doublereal *, integer *, doublereal *), s5lp_(integer *, 
	    integer *, char *, logical *, integer *, S_fp, logical *, logical 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, integer *, integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, doublereal *, doublereal *, doublereal *, doublereal *
	    , doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    char *, integer *, integer *, integer *, doublereal *, integer *, 
	    ftnlen, ftnlen);
    static integer suboptimize;
    extern /* Subroutine */ int dload_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer emode, ioldb, ninfe;
    static logical needx;
    static doublereal sinfe;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), s2applyscales_(integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), s1page_(integer *, integer *, integer *), s2amat_(
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, doublereal *, 
	    integer *), s4getb_(integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, char *, integer *, integer *, 
	    doublereal *, integer *, ftnlen), s5fixs_(integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *);
    static doublereal objadd;
    extern /* Subroutine */ int s5fixx_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *)
	    ;
    static integer iloadb;
    static doublereal infbnd;
    static integer iobjfp;
    static logical needlu;
    static integer lcrash;
    static doublereal tcrash;
    static integer inform__, lssave, numleq, istart, iinsrt, nrhslp;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen), s2crash_(integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, integer *, integer 
	    *, integer *, integer *, doublereal *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, doublereal *, integer *);
    static integer nrhslp0;
    static logical elastic;
    static integer lvlobje;
    static char probtag[20];

    /* Fortran I/O blocks */
    static icilist io___140 = { 0, str, 0, fmt_1332, 120, 1 };
    static icilist io___142 = { 0, str, 0, fmt_1316, 120, 1 };
    static icilist io___143 = { 0, str, 0, fmt_2100, 120, 1 };
    static icilist io___155 = { 0, str, 0, fmt_2200, 120, 1 };
    static icilist io___156 = { 0, str, 0, fmt_2300, 120, 1 };
    static icilist io___157 = { 0, str, 0, fmt_2400, 120, 1 };
    static icilist io___158 = { 0, str, 0, fmt_2410, 120, 1 };


/*     ================================================================== */
/*     s5getB   finds an initial basis kBS(1:m) for the linear */
/*     constraints and bounds.  First, we attempt to find a  feasible */
/*     point for the bounds and general linear equalities. This may fail, */
/*     so there is the chance that the initial basis is not feasible. */
/*     This difficulty is taken care of in subsequent calls to s5LP from */
/*     s5solveLP, s5solveQP, s5solveQN or s8feasLC. */

/*     1. The linear constraints are (optionally) scaled. */

/*     2. Elements x(n+1:n+m) of the initial x are assigned. */
/*        A scaled copy of the initial point is stored in x0. */

/*     3. An LP is used to find a feasible x for the bounds and */
/*        linear equality constraints. */


/*     On entry, */
/*     --------- */
/*     Start        is 0,1,2,3: an integer version of the solver's Start */
/*                  (character for sqopt, snoptb, snoptc, npopt */
/*                  integer    for snopta). */

/*     RowTypes(nb) is the vector of row types defined in s2Amat. */

/*     bl, bu     contain the user-defined lower and upper bounds */

/*     x(nb)        contains the initial x. Only the first n elements */
/*                  have been defined (by the user).  Here, elements */
/*                  x(n+1:n+m) are assigned, (optionally) scaled and */
/*                  saved in x0. */


/*     blSave, buSave are used to save the scaled lower and upper bounds. */

/*     iExit  Status */
/*     -----  ------ */
/*        0   basis found. */
/*            nInf = 0 => linear equalities are     satisfied. */
/*            nInf > 0 +> linear equalities are not satisfied */
/*                        (unbounded problem with emode > 0) */

/*       12   no feasible point and emode = 0. */
/*       31   Iterations limit. */
/*       >0   fatal error */

/*     31 Jul 1996: First version of s5getB. */
/*     12 Jul 1997: Thread-safe version. */
/*     01 Aug 2003: snEXIT and snPRNT adopted. */
/*     08 Mar 2004: Implemented gotHes, gotScl for Hot starts. */
/*     17 Jun 2008: Real workspace reorganized. */
/*     20 Nov 2014: Normal exit if infeasible and eMode > 0. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/* scale option */
/* Crash option */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --pi;
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
    --x;
    --scales;
    --rc;
    --busave;
    --blsave;
    --buqp;
    --blqp;
    --bu;
    --bl;
    --hs;
    --rowtypes;
    --estate;
    --etype;
    names -= 8;
    --acol;
    --inda;
    --loca;
    --rhs;
    --x0;
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
/*     ------------------------------------------------------------------ */
    iloadb = iw[122];
/* load file */
    iinsrt = iw[125];
/* insert file */
    ioldb = iw[126];
/* old basis file */
    emode = iw[56];
/* >0   => use elastic mode */
    mjrprtlvl = iw[92];
/* Major print level */
    mnrprtlvl = iw[93];
/* Minor print level */
    infbnd = rw[70];
/* definition of an infinite bound */
    tcrash = rw[62];
/* Initialize a few things. */
/* crash tolerance. */
    *iexit = 0;
    inform__ = 0;
    minimize = 1;
    nnl = max(*nnobj,*nnjac);
    *ninf = 0;
    ninfe = 0;
/* Local to s5getB */
    sinfe = 0.;
/* Local to s5getB */
    objadd = 0.;
/* Local to s5getB */
    iobjfp = 0;
/* Used for FP calculation */
    *scaleobj = 1.;
    numleq = *numlc - *numliq;
/* Initialize the row and column scales to "unit scaling" */
    if (iw[75] > 0) {
	dload_(nb, &c_b69, &scales[1], &c__1);
    }
/* ================================================================= */
/* Decode Start. */
/* ================================================================= */
    *needb = TRUE_;
    if (*start == 0 || *start == 1) {
/* ------------------------------- */
/* Cold start  or  Basis file. */
/* ------------------------------- */
	istart = 0;
/* Computing MAX */
	i__1 = max(ioldb,iinsrt);
	*needb = max(i__1,iloadb) <= 0;
	*ns = 0;
    } else if (*start == 2) {
/* ------------------------------- */
/* Warm start. */
/* ------------------------------- */
	istart = 1;
	*needb = FALSE_;
    } else if (*start == 3) {
/* ------------------------------- */
/* Hot start. */
/* ------------------------------- */
	istart = 1;
	*needb = FALSE_;
    }
    s1page_(&c__1, &iw[1], leniw);
    if (istart == 0) {
/* -------------------------------------------------------------- */
/* Cold start, or Basis file provided. */
/* Input a basis file if one exists, thereby defining hs and x. */
/* (Otherwise, s2crash will be called later to define hs.) */
/* -------------------------------------------------------------- */
/* Initialize x(n+1:nb) and pi(1:m) before scaling the problem. */
/* The basis files initialize all of x. */
/* One day they may load pi for nonlinear problems. */
	dload_(m, &c_b5, &x[*n + 1], &c__1);
	dload_(m, &c_b5, &pi[1], &c__1);
	snprnt_(&c__21, " Initial basis", &iw[1], leniw, (ftnlen)14);
	snprnt_(&c__21, " -------------", &iw[1], leniw, (ftnlen)14);
	if (*needb) {
	    snprnt_(&c__31, " No basis file supplied", &iw[1], leniw, (ftnlen)
		    23);
	    if (iw[88] == 0) {
		*needb = FALSE_;
		lcrash = 0;
		s2crash_(&lcrash, &mnrprtlvl, m, n, nb, nncon, &iw[88], &
			tcrash, nea, nloca, &loca[1], &inda[1], &acol[1], &
			kbs[1], &hs[1], &rowtypes[1], &bl[1], &bu[1], &x[1], &
			iw[1], leniw, &rw[1], lenrw);
	    }
	} else {
	    s4getb_(iexit, m, n, nb, nnames, ns, iobj, &hs[1], &bl[1], &bu[1],
		     &x[1], names + 8, &iw[1], leniw, &rw[1], lenrw, (ftnlen)
		    8);
	    if (*iexit > 0) {
/* Set blSave and buSave to avoid uninitialized copying. */
		dcopy_(nb, &bl[1], &c__1, &blsave[1], &c__1);
		dcopy_(nb, &bu[1], &c__1, &busave[1], &c__1);
		goto L900;
	    }
	}
    }
/* ----------------------------------------------------------------- */
/* Move x inside its bounds. */
/* ----------------------------------------------------------------- */
/* iStart = 0 */
    s5fixx_(&c__0, &c__1, n, tolx, &hs[1], &bl[1], &bu[1], &x[1]);
/* ----------------------------------------------------------------- */
/* Scale the linear part of the constraints. */
/* (Any nonlinear elements in A contain fake nonzeros.) */
/* ----------------------------------------------------------------- */
    if (iw[75] > 0 && *numlc > 0) {
	if (iw[232] == 0) {
	    iw[232] = 1;
	    lssave = iw[75];
	    iw[75] = 2;
	    s2getscales_(&mjrprtlvl, m, n, nb, &nnl, nncon, nnjac, &rowtypes[
		    1], nea, nloca, &loca[1], &inda[1], &acol[1], &scales[1], 
		    &bl[1], &bu[1], &y[1], &y1[1], &iw[1], leniw, &rw[1], 
		    lenrw);
	    iw[75] = lssave;
	}
	s2applyscales_(&c__0, m, n, nb, iobj, &infbnd, scaleobj, nea, nloca, &
		loca[1], &inda[1], &acol[1], &scales[1], &bl[1], &bu[1], &pi[
		1], &x[1]);
    }
/* Save the (scaled) bounds.  Save the scaled initial point in x0. */
    dcopy_(nb, &bl[1], &c__1, &blsave[1], &c__1);
    dcopy_(nb, &bu[1], &c__1, &busave[1], &c__1);
    if (*nx0 > 0) {
	dcopy_(nx0, &x[1], &c__1, &x0[1], &c__1);
    }
/* ----------------------------------------------------------------- */
/* Prepare to get feasible for the linear constraints. */
/* Relax any nonlinear rows. */
/* ----------------------------------------------------------------- */
    if (*nncon > 0 && *numlc > 0) {
	d__1 = -infbnd;
	dload_(nncon, &d__1, &bl[*n + 1], &c__1);
	dload_(nncon, &infbnd, &bu[*n + 1], &c__1);
    }
/* ----------------------------------------------------------------- */
/* Compute a starting basis for the linear constraints. */
/* This means attempting to get feasible for the linear equalities. */
/* If the E rows are infeasible, s5LP  will take care of it. */
/* ----------------------------------------------------------------- */
    if (mnrprtlvl > 10) {
	s_wsfi(&io___140);
	for (j = 1; j <= 28; ++j) {
	    do_fio(&c__1, line, (ftnlen)4);
	}
	e_wsfi();
	snprnt_(&c__31, str, &iw[1], leniw, (ftnlen)120);
	s_wsfi(&io___142);
	for (j = 1; j <= 16; ++j) {
	    do_fio(&c__1, line, (ftnlen)4);
	}
	e_wsfi();
	snprnt_(&c__32, str, &iw[1], leniw, (ftnlen)120);
    }
/* ----------------------------------------------------------------- */
/* Define the working bounds. */
/* ----------------------------------------------------------------- */
    dcopy_(nb, &bl[1], &c__1, &blqp[1], &c__1);
    dcopy_(nb, &bu[1], &c__1, &buqp[1], &c__1);
    if (*needb) {
/* -------------------------------------------------------------- */
/* Crash is needed to find a basis. */
/* -------------------------------------------------------------- */
/* Treat Crash 1 the same as Crash 2. */
	iw[88] = max(iw[88],2);
	iw[88] = min(iw[88],3);
	if (*numlc > 0) {
/* =========================================================== */
/* Find a feasible point for the linear constraints. */
/* =========================================================== */
/* Crash 2 finds a basis for ALL the LINEAR constraints */
/*         (equalities and inequalities). */
/* Crash 3 treats linear EQUALITIES separately. */
/* ----------------------------------------------------------- */
	    if (iw[88] == 2) {
/* Relax. */
	    } else if (iw[88] == 3) {
		if (numleq > 0) {
/* Find a basis for the linear EQUALITIES. */
/* Linear inequality rows are made to appear to be free */
		    if (mnrprtlvl >= 10) {
			s_wsfi(&io___143);
			do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer)
				);
			e_wsfi();
			snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)120);
		    }
		}
/* numLEQ > 0 */
	    }
/* iCrash = 2 or 3 */
	}
/* Call Crash even if there are no linear constraints. */
/* kBS(1:m) is used as workspace. */
/* We haven't done any solve yet, so we should NOT need */
/* call s5hs  ( Extern, nb, bl, bu, hs, x ) . */
/* numLC > 0 */
	lcrash = iw[88];
	s2crash_(&lcrash, &mnrprtlvl, m, n, nb, nncon, &iw[88], &tcrash, nea, 
		nloca, &loca[1], &inda[1], &acol[1], &kbs[1], &hs[1], &
		rowtypes[1], &bl[1], &bu[1], &x[1], &iw[1], leniw, &rw[1], 
		lenrw);
    }
/* ================================================================= */
/* 1. Set nS to match hs(*). */
/* 2. Set kBS(m+1:m+nS) to define the initial superbasics. */
/* 3. Set all nonbasic x to be within bounds, which may change some */
/*    hs values from 0 to 1. */
/* 4. Set nonbasic x to be exactly on nearly satisfied bounds. */
/*    (Some nonbasics may still be between bounds.) */
/* ================================================================= */
/* NeedB */
    s4checkhs_(m, maxs, mbs, n, nb, needb, &iw[231], ns, iobj, &hs[1], &kbs[1]
	    , &bl[1], &bu[1], &x[1], &iw[1], leniw, &rw[1], lenrw);
    if (*needb && *numlc > 0) {
/* Fix SBs. This forces simplex steps. */
	s5fixs_(&c__0, m, maxs, mbs, n, nb, ns, &hs[1], &kbs[1], &bl[1], &bu[
		1], &blbs[1], &bubs[1], &x[1], &xbs[1]);
	if (numleq > 0) {
	    if (*numliq == 0) {
		nrhslp = 0;
	    } else {
		lc1 = *n + *nncon + 1;
/* Relax the bounds bl and bu associated with the */
/* linear INEQUALITIES. */
/* The rhs stord in y2  makes the relaxed rows satisfied */
/* at the current x. */
		s5lg_(m, n, nb, nncon, &nrhslp, nea, nloca, &loca[1], &inda[1]
			, &acol[1], &bl[1], &bu[1], nrhs0, nrhs, &rhs[1], &y2[
			1], &x[1], &y[1], &rw[1], lenrw);
/* Initialize the working upper and lower bounds. */
/* bl and bu are unaltered in this call because eMode = 0. */
		dcopy_(numlc, &bl[lc1], &c__1, &blqp[lc1], &c__1);
		dcopy_(numlc, &bu[lc1], &c__1, &buqp[lc1], &c__1);
	    }
	    s_copy(probtag, "linear equalities", (ftnlen)20, (ftnlen)17);
	    nrhslp0 = max(nrhslp,1);
	    elastic = FALSE_;
	    lvlobje = 0;
	    emodeleq = 0;
/* No Elastic mode for E-row phase 1 */
	    needlu = TRUE_;
	    needx = needlu;
	    suboptimize = -1;
/* Set hs = 4, -1 for fixed variables and nonbasics */
/* between their bounds. */
	    s5hs_(&c__0, nb, &bl[1], &bu[1], &hs[1], &x[1]);
	    s5lp_(&inform__, &c__3, probtag, &elastic, &suboptimize, (S_fp)
		    lplog, &needlu, &needx, m, n, nb, ndegen, itqp, itqpmax, 
		    itn, &emodeleq, &lvlobje, &mnrprtlvl, &minimize, &iobjfp, 
		    scaleobj, &objadd, toloptfp, toloptqp, tolx, ninf, sinf, &
		    elastics, &ninfe, &sinfe, wtinf, pinorm, rgnorm, nea, 
		    nloca, &loca[1], &inda[1], &acol[1], &etype[1], &estate[1]
		    , &feastype[1], &hs[1], &kbs[1], &bl[1], &bu[1], &blqp[1],
		     &buqp[1], &blbs[1], &bubs[1], &gbs[1], &pi[1], &rc[1], &
		    nrhslp0, &nrhslp, &y2[1], &scales[1], &x[1], &xbs[1], &x[
		    1], &iy[1], &iy1[1], &y[1], &y1[1], cw + 8, lencw, &iw[1],
		     leniw, &rw[1], lenrw, (ftnlen)20, (ftnlen)8);
/* Reset the original bounds on the linear inequality rows. */
	    if (*numliq > 0) {
		lc1 = *n + *nncon + 1;
		dcopy_(numlc, &blsave[lc1], &c__1, &bl[lc1], &c__1);
		dcopy_(numlc, &busave[lc1], &c__1, &bu[lc1], &c__1);
		dcopy_(numlc, &blsave[lc1], &c__1, &blqp[lc1], &c__1);
		dcopy_(numlc, &busave[lc1], &c__1, &buqp[lc1], &c__1);
		s5hs_(&c__0, nb, &bl[1], &bu[1], &hs[1], &x[1]);
	    }
/* numLIQ > 0 */
	    if (mnrprtlvl >= 10) {
		if (*ninf == 0) {
/* The linear E rows are now satisfied. */
		    s_wsfi(&io___155);
		    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		    e_wsfi();
		    snprnt_(&c__21, " ", &iw[1], leniw, (ftnlen)1);
		    snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)120);
		} else {
/* The linear E rows are infeasible. */
		    s_wsfi(&io___156);
		    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		    e_wsfi();
		    snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)120);
		}
	    }
	}
/* Potential inform values from s5LP are: */
/*   -3    Too many iterations */
/*   -2    variable is unbounded (this should not happen) */
/*   -1    linear constraints are infeasible */
/*    0    feasible point found */
/*   >0    Fatal error */
/* numLEQ > 0 */
	if (inform__ > 0) {
	    *iexit = inform__;
/*  fatal error. Exit */
	} else if (inform__ == -3) {
	    *iexit = 31;
/* Too many iterations */
	} else if ((inform__ == -1 || inform__ == -2) && emode == 0) {
/* Infeasible linear equalities and no elastic mode. */
	    *iexit = 12;
/* infeasible linear equalities */
	} else {
/* A feasible point has been found. */

/* Define a basis that includes the linear INEQUALITIES. */
/* If the linear equalities are infeasible, the feasibility */
/* phase for the inequalities (with an appropriate value of */
/* eMode) will determine what happens next. */
	    *iexit = 0;
	    if (*numliq > 0) {
		if (mnrprtlvl >= 10) {
		    if (numleq > 0) {
			s_wsfi(&io___157);
			do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer)
				);
			e_wsfi();
			snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)120);
		    } else {
			s_wsfi(&io___158);
			do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer)
				);
			e_wsfi();
			snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)120);
		    }
		}
		s5hs_(&c__1, nb, &bl[1], &bu[1], &hs[1], &x[1]);
		s2amat_(&c__0, &mnrprtlvl, m, n, nb, nncon, nnjac, nnobj, 
			iobj, numlc, numliq, nea, nloca, &loca[1], &inda[1], &
			acol[1], &bl[1], &bu[1], &rowtypes[1], &iw[1], leniw, 
			&rw[1], lenrw);
		lcrash = 4;
		s2crash_(&lcrash, &mnrprtlvl, m, n, nb, nncon, &iw[88], &
			tcrash, nea, nloca, &loca[1], &inda[1], &acol[1], &
			kbs[1], &hs[1], &rowtypes[1], &bl[1], &bu[1], &x[1], &
			iw[1], leniw, &rw[1], lenrw);
	    }
/* numLIQ > 0 */
	}
    }
/* NeedB and numLC > 0 */
L900:
    return 0;
} /* s5getb_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5getB */
/* Subroutine */ int s5hs_(integer *mode, integer *nb, doublereal *bl, 
	doublereal *bu, integer *hs, doublereal *x)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j;

/*     ================================================================== */
/*     s5hs   sets the state vector hs. */
/*     if mode = 'Internal', s5hs sets hs(j) = -1 or 4 for certain */
/*        nonbasic variables.  This allows s5price to operate more */
/*        efficiently.  The internal values of hs are now as follows: */

/*        hs(j) = -1  Nonbasic between bounds (bl     <  x <  bu    ) */
/*        hs(j) =  0  Nonbasic on lower bound (bl-tol <  x <= bl    ) */
/*        hs(j) =  1  Nonbasic on upper bound (bu     <= x <  bu+tol) */
/*        hs(j) =  2  Superbasic */
/*        hs(j) =  3  Basic */
/*        hs(j) =  4  Nonbasic and fixed      (bl     = x  =  bu    ) */

/*        where 0 <= tol < the feasibility tolerance. */

/*     if mode = 'External', s5hs changes -1 or 4 values to hs(j) = 0, */
/*        ready for basis saving and the outside world. */

/*     08 Apr 1992: First version of s5hs. */
/*     21 Aug 1999: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --x;
    --hs;
    --bu;
    --bl;

    /* Function Body */
    if (*mode == 0) {
/*        --------------------------------------------------------------- */
/*        Change nonbasic hs(j) to internal values (including 4 and -1). */
/*        This may change existing internal values if bl and bu have been */
/*        changed -- e.g. at the start of each major iteration. */
/*        --------------------------------------------------------------- */
	i__1 = *nb;
	for (j = 1; j <= i__1; ++j) {
	    if (hs[j] <= 1) {
		if (bl[j] == bu[j]) {
		    hs[j] = 4;
		} else if (x[j] <= bl[j]) {
		    hs[j] = 0;
		} else if (x[j] >= bu[j]) {
		    hs[j] = 1;
		} else {
		    hs[j] = -1;
		}
	    }
	}
    } else if (*mode == 1) {
/*        --------------------------------------------------------------- */
/*        Change hs to external values. */
/*        Some nonbasic hs(j) may be 4 or -1.  Change them to 0. */
/*        --------------------------------------------------------------- */
	i__1 = *nb;
	for (j = 1; j <= i__1; ++j) {
	    if (hs[j] == 4 || hs[j] == -1) {
		hs[j] = 0;
	    }
	}
    }
    return 0;
} /* s5hs_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5hs */
/* Subroutine */ int s5inf_(integer *nbs, doublereal *featol, doublereal *
	infbnd, integer *ninf, doublereal *sinf, integer *feastype, 
	doublereal *blbs, doublereal *bubs, doublereal *gbs, doublereal *xbs)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer k;
    static doublereal xk, res, inflow, infupp;
    static integer inftype;

/*     ================================================================== */
/*     s5Inf  computes the sum and number of infeasibilities. */

/*     On entry: */
/*     -------- */
/*     featol        is the current feasibility tolerance. */
/*     infBnd        is the magnitude of an infinite bound. */
/*     gBB(nBS)      is the zero vector. */
/*     xBS(nBS)      is the vector of basics and superbasics */

/*     On exit: */
/*     -------- */
/*     nInf          is the number violated nonelastic bounds. */
/*     gBS(nBS)      is the rhs for the equations for pi. */
/*     feasType(nBS) defines the type of infeasibility. */

/*     feasType     x(j)                          Meaning */
/*     --------     -----                         ------- */
/*      -2      infeasible                           x(j) .le. bl(j)-tol */
/*       0        feasible            bl(j)-tol .le. x(j) .le. bu(j)+tol */
/*      +2      infeasible                           x(j) .ge. bu(j)+tol */

/*     feasType is used in s5step. */

/*     29 Oct 1993: First version of s5Inf. */
/*     26 May 2013: infBnd used to identify infinite bounds. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --xbs;
    --gbs;
    --bubs;
    --blbs;
    --feastype;

    /* Function Body */
    infupp = *infbnd * .99;
    inflow = *infbnd * -.99;
    *ninf = 0;
    *sinf = 0.;
    i__1 = *nbs;
    for (k = 1; k <= i__1; ++k) {
	inftype = 0;
	xk = xbs[k];
/* Check if the lower bound (if any) is violated. */
	if (blbs[k] > inflow) {
	    res = blbs[k] - xk;
	    if (res > *featol) {
		gbs[k] = -1.;
		inftype = -2;
	    }
	}
/* Check if the upper bound (if any) is violated. */
	if (inftype == 0 && bubs[k] < infupp) {
	    res = xk - bubs[k];
	    if (res > *featol) {
		gbs[k] = 1.;
		inftype = 2;
	    }
	}
	if (inftype != 0) {
	    ++(*ninf);
	    *sinf += res;
	}
	feastype[k] = inftype;
    }
    return 0;
} /* s5inf_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5Inf */
/* Subroutine */ int s5lg_(integer *m, integer *n, integer *nb, integer *
	nncon, integer *nrhslp, integer *nea, integer *nloca, integer *loca, 
	integer *inda, doublereal *acol, doublereal *bl, doublereal *bu, 
	integer *nrhs0, integer *nrhs, doublereal *rhs, doublereal *rhslp, 
	doublereal *x, doublereal *y, doublereal *rw, integer *lenrw)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j;
    static doublereal eps0;
    extern /* Subroutine */ int dload_(integer *, doublereal *, doublereal *, 
	    integer *), daxpy_(integer *, doublereal *, doublereal *, integer 
	    *, doublereal *, integer *);
    static doublereal infbnd;
    extern /* Subroutine */ int s2aprod_(integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *);

/*     ================================================================== */
/*     s5LG  relaxes the linear inequality constraints for an LP that */
/*     gets feasible for the linear equality constraints. A right-hand */
/*     side is computed that makes the relaxed rows satisfied at the */
/*     current x.   Then x will not be disturbed more than */
/*     necessary during a Warm start. */

/*     02 Apr 2005: First version of s5LG. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --y;
    --rhslp;
    --x;
    --bu;
    --bl;
    --acol;
    --inda;
    --loca;
    --rhs;
    --rw;

    /* Function Body */
    eps0 = rw[2];
/* eps**(4/5) */
    infbnd = rw[70];
/* definition of an infinite bound */
    *nrhslp = *m;
    s2aprod_(&c__0, &eps0, nea, nloca, &loca[1], &inda[1], &acol[1], &c_b69, &
	    x[1], n, &c_b5, &y[1], m);
    daxpy_(nrhslp, &c_b159, &x[*n + 1], &c__1, &y[1], &c__1);
    if (*nrhs > 0) {
	daxpy_(nrhs, &c_b159, &rhs[1], &c__1, &y[1], &c__1);
    }
    dload_(nrhslp, &c_b5, &rhslp[1], &c__1);
    i__1 = *m;
    for (i__ = *nncon + 1; i__ <= i__1; ++i__) {
	j = *n + i__;
	if (bl[j] < bu[j]) {
	    bl[j] = -infbnd;
	    bu[j] = infbnd;
	    rhslp[i__] = y[i__];
	}
    }
    return 0;
} /* s5lg_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5LG */
/* Subroutine */ int s5lpitn_(integer *iexit, logical *feasible, logical *
	increase, logical *needpi, logical *elastic, integer *m1, integer *m, 
	integer *nb, integer *ndegen, integer *elastics, integer *lurequest, 
	integer *kp, integer *jbq, integer *jsq, integer *jbr, integer *jsr, 
	integer *jq, doublereal *featol, doublereal *pivot, doublereal *step, 
	doublereal *tolinc, integer *etype, integer *estate, integer *
	feastype, integer *hs, integer *kbs, doublereal *bl, doublereal *bu, 
	doublereal *blqp, doublereal *buqp, doublereal *blbs, doublereal *
	bubs, doublereal *x, doublereal *xbs, doublereal *pbs, doublereal *y1,
	 integer *iw, integer *leniw, doublereal *rw, integer *lenrw)
{
    /* Format strings */
    static char fmt_1000[] = "(\002 Variable\002,i6,\002  can increase indef"
	    "initely\002)";
    static char fmt_1100[] = "(\002 Variable\002,i6,\002  can decrease indef"
	    "initely\002)";

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static logical unbounded;
    static integer jr;
    static char str[50];
    static logical move;
    static doublereal tolp;
    static integer ntry;
    static doublereal tolp0;
    extern /* Subroutine */ int s5bsx_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *), dscal_(integer *, 
	    doublereal *, doublereal *, integer *);
    static doublereal bigdx, exact, bound;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static doublereal pnorm, stepp;
    extern /* Subroutine */ int s2bmod_(integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *), 
	    s5step_(integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, doublereal *, doublereal *, logical *, logical *, 
	    logical *, logical *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    static doublereal infbnd;
    extern doublereal dnormi_(integer *, doublereal *, integer *);
    static integer inform__, infpiv;
    static doublereal sclpiv;
    static logical hitlow;
    static doublereal tolpiv;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);
    static integer jestate;
    static logical onbound;
    static integer jqstate, jrstate;
    static doublereal stepmax;

    /* Fortran I/O blocks */
    static icilist io___193 = { 0, str, 0, fmt_1000, 50, 1 };
    static icilist io___194 = { 0, str, 0, fmt_1100, 50, 1 };


/*     ================================================================== */
/*     s5LPitn computes one step of the primal simplex method. */
/*     jq is the variable entering the basis and djq is its reduced cost. */

/*      iExit       Status */
/*      -----       ------ */
/*        0         Normal exit */
/*        1         LP is unbounded */

/*     10 Sep 1991: First version based on Minos routine m5lpit. */
/*     02 Aug 1996: First min sum version added by PEG. */
/*     12 Jul 1997: Thread-safe version. */
/*     01 Aug 2003: snPRNT adopted. */
/*     08 Apr 2008: eState accessed in elastic mode only. */
/*     04 Jul 2008: Modify both bl and bu in elastic phase 1. */
/*     ================================================================= */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
/* number of LU mods */
    /* Parameter adjustments */
    --y1;
    --pbs;
    --xbs;
    --bubs;
    --blbs;
    --kbs;
    --feastype;
    --x;
    --buqp;
    --blqp;
    --bu;
    --bl;
    --hs;
    --estate;
    --etype;
    --iw;
    --rw;

    /* Function Body */
    tolpiv = rw[60];
/* excludes small elements of pBS. */
    infbnd = rw[70];
/* definition of an infinite bound */
    bigdx = rw[72];
/* s5step assumes that the first m1 components of xBS can Move. */
/* unbounded step. */
    *jsq = *jq;
/* Entering superbasic */
    *jsr = *jq;
/* leaving  superbasic */
    feastype[*m1] = 0;
    kbs[*m1] = *jq;
    xbs[*m1] = x[*jq];
    blbs[*m1] = blqp[*jq];
    bubs[*m1] = buqp[*jq];
/* ================================================================= */
/* Set eState(jq) and the elastic parts of blBS and buBS. */
/* ================================================================= */
    if (*elastic) {
/* If the new superbasic is an elastic variable */
/* and it wants to move infeasible, set its elastic state. */
	if (etype[*jq] > 0 && estate[*jq] == 0) {
	    jqstate = hs[*jq];
	    if (*increase) {
		if (jqstate == 1 || jqstate == 4) {
		    ++(*elastics);
		    estate[*jq] = 2;
		    blqp[*jq] = buqp[*jq];
		    buqp[*jq] = infbnd;
		    blbs[*m1] = blqp[*jq];
		    bubs[*m1] = buqp[*jq];
		}
	    } else {
		if (jqstate == 0 || jqstate == 4) {
		    ++(*elastics);
		    estate[*jq] = 1;
		    buqp[*jq] = blqp[*jq];
		    blqp[*jq] = -infbnd;
		    blbs[*m1] = blqp[*jq];
		    bubs[*m1] = buqp[*jq];
		}
	    }
	}
    }
/*     ================================================================== */
/*     Select a variable to be dropped from B. */
/*     s5step  uses the (m+1)th element of  blBS, buBS, xBS and pBS. */
/*     ================================================================== */
/* Elastic mode */
    if (*increase) {
	dscal_(m, &c_b159, &pbs[1], &c__1);
	pbs[*m1] = 1.;
    } else {
	pbs[*m1] = -1.;
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
    pnorm = dnormi_(m1, &pbs[1], &c__1);
    stepmax = bigdx / pnorm;
    sclpiv = 1.;
    tolp0 = tolpiv;
    tolp = tolpiv * pnorm;
    ntry = 0;
/* +    Repeat */
L200:
    tolp /= sclpiv;
    tolp0 /= sclpiv;
    s5step_(m1, ndegen, featol, &infbnd, &stepmax, tolinc, &tolp, &feastype[1]
	    , &blbs[1], &bubs[1], &xbs[1], &pbs[1], &hitlow, &move, &onbound, 
	    &unbounded, &infpiv, kp, &bound, &exact, step, &stepp);
    sclpiv = 10.;
    ++ntry;
/* +    until    ( infpiv .eq. 0 .and. (.not.Unbounded .or. Feasible) .or. */
/* +                 ntry .ge. mtry) */
    if (! (infpiv == 0 && (! unbounded || *feasible) || ntry >= 6)) {
	goto L200;
    }
    if (! unbounded) {
/* -------------------------------------------------------------- */
/* Update the basic variables xBS and copy them into x. */
/* -------------------------------------------------------------- */
	jr = kbs[*kp];
	daxpy_(m1, step, &pbs[1], &c__1, &xbs[1], &c__1);
	s5bsx_(&c__1, m1, nb, &kbs[1], &x[1], &xbs[1]);
/* $$$!        10 Mar 2004: Care is needed to prevent the */
/* $$$!        new nonbasic variable jr from ending up slightly inside */
/* $$$!        its bound.  EXPAND normally ensures that x(jr) will be */
/* $$$!        ON or slightly OUTSIDE its bound, but now we realise that */
/* $$$!        rounding error might make it slightly INSIDE. */
/* $$$ */
/* $$$         if (OnBound) then */
/* $$$            x(jr) = bound */
/* $$$         else if (HitLow) then */
/* $$$            x(jr) = min( x(jr), blQP(jr) ) */
/* $$$         else */
/* $$$            x(jr) = max( x(jr), buQP(jr) ) */
/* $$$         end if */
	if (onbound) {
	    x[jr] = bound;
	}
	if (*kp == *m1) {
/* ----------------------------------------------------------- */
/* Variable jq reaches its opposite bound. */
/* ----------------------------------------------------------- */
	    if (*increase) {
		hs[*jq] = 1;
	    } else {
		hs[*jq] = 0;
	    }
	    feastype[*kp] = 0;
	    *pivot = 0.;
	    if (! (*feasible) || *elastic) {
		*needpi = TRUE_;
	    }
	} else {
/* ----------------------------------------------------------- */
/* Variable jq replaces the kp-th variable of  B. */
/* It could be a fixed variable, whose new state must be 4. */
/* ----------------------------------------------------------- */
	    *needpi = TRUE_;
	    *jbq = *jq;
	    *jbr = jr;
	    hs[*jq] = 3;
	    *pivot = -pbs[*kp];
	    if (*elastic) {
		jestate = estate[jr];
	    } else {
		jestate = 0;
	    }
	    if (jestate == 0) {
/* Normal variable hits its bound. */
		if (blbs[*kp] == bubs[*kp]) {
		    jrstate = 4;
		} else if (hitlow) {
		    jrstate = 0;
		} else {
		    jrstate = 1;
		}
	    } else {
/* Elastic x(jr) hits its bound. */
/* Reset the true upper and lower bounds. */
		estate[jr] = 0;
		--(*elastics);
		blqp[jr] = bl[jr];
		buqp[jr] = bu[jr];
		if (jestate == 1) {
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
/*   jeState .eq. 2 */
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
	    hs[jr] = jrstate;
	    kbs[*kp] = *jq;
	    xbs[*kp] = xbs[*m1];
	    blbs[*kp] = blbs[*m1];
	    bubs[*kp] = bubs[*m1];
/* Update the LU factors. */
	    ++iw[216];
	    s2bmod_(&inform__, kp, m, &y1[1], &iw[1], leniw, &rw[1], lenrw);
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
	}
/* kp ne m1 */
    } else {
/* Apparently the solution is unbounded. */
	if (*increase) {
	    s_wsfi(&io___193);
	    do_fio(&c__1, (char *)&(*jq), (ftnlen)sizeof(integer));
	    e_wsfi();
	} else {
	    s_wsfi(&io___194);
	    do_fio(&c__1, (char *)&(*jq), (ftnlen)sizeof(integer));
	    e_wsfi();
	}
	snprnt_(&c__21, str, &iw[1], leniw, (ftnlen)50);
	*iexit = 1;
/* Unbounded direction */
    }
    return 0;
} /* s5lpitn_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5LPitn */
/* Subroutine */ int s5price_(logical *elastic, logical *feasible, logical *
	increase, logical *gotg, integer *suboptimize, integer *itn, integer *
	m, integer *n, integer *nb, integer *leng, integer *neg, integer *
	frozen, integer *nonopt, doublereal *weight, doublereal *signobj, 
	doublereal *pinorm, integer *jq, doublereal *djq, integer *kprc, 
	doublereal *toldj, integer *nea, integer *nloca, integer *loca, 
	integer *inda, doublereal *acol, integer *etype, integer *hs, 
	doublereal *g, doublereal *pi, doublereal *rc, doublereal *x, 
	doublereal *xfrozen, integer *iw, integer *leniw, doublereal *rw, 
	integer *lenrw)
{
    /* Format strings */
    static char fmt_1000[] = "(\002 Itn\002,i7,\002: toldj =\002,1p,e8.1,"
	    "\002    Norm pi =\002,e8.1,\002    weight = \002,e8.1)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static doublereal d__;
    static integer j;
    static doublereal d1, d2;
    static integer j1, j2, k1, k2;
    static doublereal dj;
    static integer jj, js, np;
    static doublereal dj1, dj2;
    static char str[100];
    extern /* Subroutine */ int s5rc_(integer *, integer *, logical *, 
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *, integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *);
    static integer npr1, npr2, nprc, jslk;
    static doublereal told;
    extern /* Subroutine */ int s5erc_(integer *, integer *, logical *, 
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *, integer *, integer *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *);
    static integer jfree;
    static doublereal djmax;
    static integer lvldj, nparp, kpsav;
    static doublereal infbnd;
    static integer lprdbg, etypej, nparpr;
    static doublereal tolmin;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___225 = { 0, str, 0, fmt_1000, 100, 1 };


/*     ================================================================== */
/*     s5price  selects a nonbasic variable to enter the basis, */
/*     using the reduced gradients  dj = g(j) - pi'*a(j). */

/*     This version does partial pricing on both structurals and slacks. */
/*     Dynamic tolerances are used if partial price is in effect. */

/*     Partial pricing here means sectional pricing, because the three */
/*     blocks of  (A  -I)  are sliced up into nParPr sections */
/*     of equal size.  (The last section of each may be a little bigger, */
/*     since nParPr is unlikely to divide evenly into  n  or  m.) */

/*     input    g      = gradient for nonlinear variables. */
/*              pi     = pricing vector. */
/*              kPrc   = the no. of the section where s5price last found */
/*                       a useful dj. */
/*                       (kPrc = 0 at the start of each major iteration.) */
/*              toldj(1:2) hold the current told if partial pricing, for */
/*                       phase 1 and 2 respectively.  told is used to */
/*                       determine if a dj is significant. */
/*              toldj(3) holds the specified optimality tolerance. */
/*              biggst   keeps track of the biggest dj found during the */
/*                       last scan of all sections of ( A  -I ). */

/*     output   kPrc   = the last section scanned. */
/*              nonOpt = the no. of useful djs found in that section. */
/*              jq     = best column found. */
/*              djq    = best dj. */
/*              toldj(1:2) save the current told if partial pricing. */
/*              Increase says if variable jq should increase or decrease. */

/*     In the code below, */
/*     the next section of  A  contains nPr1 structurals (j1+1 thru k1), */
/*     the next section of -I  contains nPr2 slacks      (j2+1 thru k2). */
/*     If  nParPr  is rather large, either nPr1 or nPr2 could be zero, */
/*     but not both. */

/*     If subOptimize > 0,  variables that haven't moved are not */
/*     priced, thereby limiting the number of superbasic variables. */

/*     ------------------------------------------------------------------ */
/*     09 Aug 1992: First version of s5pric based on Minos 5.4 m5pric. */
/*     29 Jul 1996: Multiple pricing removed. */
/*     05 Aug 1996: First version with Elastic mode. */
/*     12 Jul 1997: Thread-safe version. */
/*     22 Dec 1999: subOptimize implemented. */
/*     31 Mar 2001: Free variables with hs(j)=-1 eligible for the basis. */
/*     19 Sep 2014: Name changed from s5pric to s5price. */
/*     19 Sep 2014: Added case for a temporarily fixed elastic. */
/*     21 Feb 2015: Deleted code that removes nonbasics floating free */
/*                  between their bounds. This code can cause a bogus */
/*                  unbounded indication in s5step. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/* -->  parameter         (zero = 0.0d+0, reduce = 0.25d+0) */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --pi;
    --xfrozen;
    --x;
    --rc;
    --hs;
    --etype;
    --g;
    --toldj;
    --acol;
    --inda;
    --loca;
    --iw;
    --rw;

    /* Function Body */
    infbnd = rw[70];
/* definition of an infinite bound */
    lprdbg = iw[85];
/* > 0    => private debug print */
    nparpr = iw[94];
/* # of partial pricing sections */
    djmax = -infbnd;
    *djq = 0.;
    *jq = 0;
    jfree = 0;
    *frozen = 0;
    *nonopt = 0;
    nprc = 0;
    nparp = nparpr;
    npr1 = *n / nparp;
    npr2 = *m / nparp;
    if (max(npr1,npr2) <= 0) {
	nparp = 1;
    }
/*     Set the tolerance for a significant dj. */
    tolmin = toldj[3] * *pinorm;
    if (*feasible) {
	lvldj = 2;
    } else {
	lvldj = 1;
    }
    told = toldj[lvldj];
    if (nparp == 1) {
	told = tolmin;
    }
/*     Set pointers to the next section of  A  and  -I. */
/*     nPrc counts how many sections have been scanned in this call. */
/*     kPrc keeps track of which one to start with. */
L100:
    ++nprc;
    ++(*kprc);
    if (*kprc > nparp) {
	*kprc = 1;
    }
    npr1 = *n / nparp;
    j1 = (*kprc - 1) * npr1;
    k1 = j1 + npr1;
    if (*kprc == nparp) {
	k1 = *n;
    }
/* Computing MAX */
    i__1 = 0, i__2 = k1 - j1;
    npr1 = max(i__1,i__2);
    npr2 = *m / nparp;
    j2 = *n + (*kprc - 1) * npr2;
    k2 = j2 + npr2;
    if (*kprc == nparp) {
	k2 = *nb;
    }
/* Computing MAX */
    i__1 = 0, i__2 = k2 - j2;
    npr2 = max(i__1,i__2);
/*     ------------------------------------------------------------------ */
/*     Main loops for partial pricing (or full pricing). */
/*     Compute reduced costs rc(*) */
/*     for the kPrc-th section of structurals */
/*     and the kPrc-th section of slacks. */
/*     ------------------------------------------------------------------ */
    i__1 = j1 + 1;
    s5rc_(&i__1, &k1, gotg, m, n, leng, neg, signobj, nea, nloca, &loca[1], &
	    inda[1], &acol[1], &hs[1], &g[1], &pi[1], &rc[1]);
    i__1 = k2;
    for (j = j2 + 1; j <= i__1; ++j) {
	rc[j] = pi[j - *n];
    }
/*     ------------------------------------------------------------------ */
/*     Main loop for pricing structural and slack reduced costs. */
/*     dj is rc(j), the reduced cost. */
/*     d  is -dj or +dj, depending on which way x(j) can move. */
/*     We are looking for the largest d (which will be positive). */
/*     ------------------------------------------------------------------ */
    np = npr1 + npr2;
    j = j1;
    jslk = npr1 + 1;
    i__1 = np;
    for (jj = 1; jj <= i__1; ++jj) {
	if (jj == jslk) {
	    j = j2;
	}
	++j;
	js = hs[j];
	if (js <= 1) {
	    dj = rc[j];
	    if (js == 0) {
/*              xj  is allowed to increase. */
		d__ = -dj;
	    } else if (js == 1) {
/*              xj  is allowed to decrease. */
		d__ = dj;
	    } else {
/*              js is -1. */
/*              xj  is free to move either way. */
/*              Save the index as jfree in case it is the only one. */
		d__ = abs(dj);
		jfree = j;
	    }
	    if (*suboptimize > 0) {
		if (x[j] == xfrozen[j]) {
		    if (d__ > told) {
			++(*frozen);
			d__ = 0.;
		    }
		}
	    }
/*           See if this dj is significant. */
/*           Also see if it is the biggest dj so far. */
	    if (d__ > told) {
		++(*nonopt);
	    }
	    if (djmax < d__) {
		djmax = d__;
		*djq = dj;
		*jq = j;
		kpsav = *kprc;
	    }
	}
    }
    if (*elastic) {
/*        --------------------------------------------------------------- */
/*        Scan this section again, looking for nonbasic elastics. */
/*        --------------------------------------------------------------- */
/*        Compute reduced costs rc(j) for fixed nonbasic columns. */
/*        (These columns are skipped in s5rc) */
	i__1 = j1 + 1;
	s5erc_(&i__1, &k1, gotg, m, n, leng, neg, signobj, nea, nloca, &loca[
		1], &inda[1], &acol[1], &etype[1], &hs[1], &g[1], &pi[1], &rc[
		1]);
	j = j1;
	i__1 = np;
	for (jj = 1; jj <= i__1; ++jj) {
	    if (jj == jslk) {
		j = j2;
	    }
	    ++j;
	    etypej = etype[j];
	    if (etypej > 0) {
		js = hs[j];
		dj = rc[j];
		if (js == 0) {
/*                 ------------------------------------------------------ */
/*                 Nonbasic at its lower bound. */
/*                 An elastic xj can decrease through the bound. */
/*                 ------------------------------------------------------ */
		    if (etypej == 1 || etypej == 3) {
			dj -= *weight;
			d__ = dj;
		    }
		} else if (js == 1) {
/*                 ------------------------------------------------------ */
/*                 Nonbasic at its upper bound. */
/*                 The default is to allow xj to decrease. */
/*                 However, an elastic xj can increase through the bound. */
/*                 ------------------------------------------------------ */
		    if (etypej == 2 || etypej == 3) {
			dj += *weight;
			d__ = -dj;
		    }
		} else if (js == 4 || js == -1) {
/*                 ------------------------------------------------------ */
/*                 Fixed elastic variable. */
/*                 xj is free to move either way. */
/*                 ------------------------------------------------------ */
		    if (etypej == 2) {
			dj1 = 0.;
			d1 = 0.;
		    } else {
			dj1 = dj - *weight;
			d1 = dj1;
		    }
		    if (etypej == 1) {
			dj2 = 0.;
			d2 = 0.;
		    } else {
			dj2 = dj + *weight;
			d2 = -dj2;
		    }
		    if (d1 >= d2) {
/*                    xj  is allowed to decrease. */
			dj = dj1;
			d__ = d1;
		    } else {
/*                    xj  is allowed to increase. */
			dj = dj2;
			d__ = d2;
		    }
		} else {
		    d__ = 0.;
		    dj = 0.;
		}
		if (*suboptimize > 0) {
		    if (x[j] == xfrozen[j]) {
			if (d__ > told) {
			    ++(*frozen);
			    d__ = 0.;
			}
		    }
		}
/*              See if this dj is significant. */
/*              Also see if it is the biggest dj so far. */
		if (d__ > told) {
		    ++(*nonopt);
		}
		if (djmax < d__) {
		    djmax = d__;
		    *djq = dj;
		    *jq = j;
		    kpsav = *kprc;
		}
	    }
	}
    }
/*     ------------------------------------------------------------------ */
/*     End of loop looking for biggest dj in the kPrc-th section. */
/*     ------------------------------------------------------------------ */
    if (*nonopt == 0) {
	if (nparp > 1) {
/*           ============================================================ */
/*           No significant dj has been found.  (All are less than told.) */
/*           Price the next section, if any remain. */
/*           ============================================================ */
	    if (nprc < nparp) {
		goto L100;
	    }
/*           ============================================================ */
/*           All sections have been scanned.  Reduce told */
/*           and grab the best dj if it is bigger than tolmin. */
/*           ============================================================ */
	    if (djmax > tolmin) {
		*nonopt = 1;
		*kprc = kpsav;
/* Computing MAX */
		d__1 = djmax * .2;
		told = max(d__1,tolmin);
		toldj[lvldj] = told;
		if (lprdbg >= 1) {
		    s_wsfi(&io___225);
		    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&told, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&(*pinorm), (ftnlen)sizeof(
			    doublereal));
		    do_fio(&c__1, (char *)&(*weight), (ftnlen)sizeof(
			    doublereal));
		    e_wsfi();
		    snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)100);
		}
	    }
	}
    }
    *increase = *djq < 0.;
    return 0;
} /* s5price_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5price */
/* Subroutine */ int s5rc_(integer *j1, integer *j2, logical *gotg, integer *
	m, integer *n, integer *leng, integer *neg, doublereal *signobj, 
	integer *nea, integer *nloca, integer *loca, integer *inda, 
	doublereal *acol, integer *hs, doublereal *g, doublereal *pi, 
	doublereal *rc)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, l;
    static doublereal dj;

/*     ================================================================== */
/*     s5rc   computes reduced costs rc(j) for nonbasic columns of A */
/*     in the range j = j1 to j2.  It is called by s5price. */

/*     The loop for computing dj for each column could conceivably be */
/*     optimized on some machines.  However, there are seldom more than */
/*     5 or 10 entries in a column. */

/*     Note that we could skip fixed variables by passing in the bounds */
/*     and testing if bl(j) .eq. bu(j), but these are relatively rare. */
/*     But see comment for 08 Apr 1992 in m5pric. */

/*     31 Jan 1992: First version of s5rc. */
/*     08 Apr 1992: Internal values of hs(j) are now used, so fixed */
/*                  variables (hs(j) = 4) are skipped as we would like. */
/*     03 Apr 1999: Linear objective stored as row 0 of A. */
/*     11 Apr 1999: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --pi;
    --rc;
    --hs;
    --g;
    --acol;
    --inda;
    --loca;

    /* Function Body */
    i__1 = *j2;
    for (j = *j1; j <= i__1; ++j) {
	if (hs[j] <= 1) {
	    dj = 0.;
	    i__2 = loca[j + 1] - 1;
	    for (l = loca[j]; l <= i__2; ++l) {
		i__ = inda[l];
		dj += pi[i__] * acol[l];
	    }
	    rc[j] = -dj;
	}
    }
/*     Include the nonlinear objective gradient if relevant. */
    if (*gotg) {
	if (*j1 <= *neg) {
	    i__1 = min(*j2,*neg);
	    for (j = *j1; j <= i__1; ++j) {
		if (hs[j] <= 1) {
		    rc[j] += *signobj * g[j];
		}
	    }
	}
    }
    return 0;
} /* s5rc_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5rc */
/* Subroutine */ int s5setpi_(integer *iexit, integer *m, logical *checkpi, 
	doublereal *pinorm, doublereal *rhs, doublereal *pi, integer *iw, 
	integer *leniw, doublereal *rw, integer *lenrw)
{
    /* System generated locals */
    doublereal d__1;

    /* Local variables */
    static doublereal flmax;
    extern /* Subroutine */ int s2bsol_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *);
    extern doublereal dnormi_(integer *, doublereal *, integer *), dnormj_(
	    integer *, doublereal *, integer *);

/*     ================================================================== */
/*     s5setpi solves  B' pi = rhs.  Beware -- rhs is altered by s2Bsol. */
/*     If a new x has just been computed, the norm is computed by dnormj. */

/*     On Exit: */
/*      iExit  = -1 if pi contains a NaN or Inf. */
/*      iExit  =  0 if pi was computed successfully. */
/*      iExit  >  0 if there was an unexpected error in the solver. */

/*     08 Aug 1996: First version of s5setp. */
/*     16 Nov 2001: Infinity norm used for piNorm (no longer dnrm1s). */
/*     19 Sep 2014: Name changed to s5setpi. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --pi;
    --rhs;
    --iw;
    --rw;

    /* Function Body */
    flmax = rw[8];
/* est. of the largest pos. real */
    *iexit = 0;
    s2bsol_(iexit, &c__2, m, &rhs[1], &pi[1], &iw[1], leniw, &rw[1], lenrw);
    if (*iexit != 0) {
	return 0;
    }
    if (*checkpi) {
	*pinorm = dnormj_(m, &pi[1], &c__1);
	if (*pinorm < flmax) {
/* false if pi = inf, nan */
	    *pinorm = max(*pinorm,1.);
	} else {
	    *iexit = -1;
	}
    } else {
/* Computing MAX */
	d__1 = dnormi_(m, &pi[1], &c__1);
	*pinorm = max(d__1,1.);
    }
    return 0;
} /* s5setpi_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5setpi */
/* Subroutine */ int s5setx_(integer *iexit, integer *task, integer *itn, 
	integer *m, integer *n, integer *nb, integer *nbs, doublereal *
	rowerror, integer *nea, integer *nloca, integer *loca, integer *inda, 
	doublereal *acol, integer *kbs, doublereal *xbs, integer *nrhs0, 
	integer *nrhs, doublereal *rhs, doublereal *x, doublereal *y, 
	doublereal *y1, integer *iw, integer *leniw, doublereal *rw, integer *
	lenrw)
{
    /* Format strings */
    static char fmt_1000[] = "(\002 Itn\002,i7,\002: Row check\002,\002.  Ma"
	    "x residual =\002,1p,e8.1,\002 on row\002,i5,\002.  Norm x =\002,"
	    "e8.1,\002.  Norm dx =\002,e8.1)";
    static char fmt_1001[] = "(\002 Itn\002,i7,\002: Row check\002,\002.  Ma"
	    "x residual =\002,1p,e8.1,\002 on row\002,i5)";

    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static char str[110];
    static doublereal eps0;
    static integer imax;
    static doublereal rmax;
    extern /* Subroutine */ int s5bsx_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *), dload_(integer *, 
	    doublereal *, doublereal *, integer *), dcopy_(integer *, 
	    doublereal *, integer *, doublereal *, integer *);
    static logical goodx;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static doublereal xnorm;
    extern /* Subroutine */ int s2bsol_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *);
    extern integer jdamax_(integer *, doublereal *, integer *);
    static integer lprdbg;
    static logical bigres;
    extern doublereal dnormi_(integer *, doublereal *, integer *), dnormj_(
	    integer *, doublereal *, integer *);
    static doublereal dxnorm;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);
    static doublereal tolrow;
    extern /* Subroutine */ int s2aprod_(integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *);

    /* Fortran I/O blocks */
    static icilist io___241 = { 0, str, 0, fmt_1000, 110, 1 };
    static icilist io___242 = { 0, str, 0, fmt_1001, 110, 1 };


/*     ================================================================== */
/*     s5setx performs the following functions: */
/*      Task            Function */
/*      ====           ======== */
/*       0 (Resetx)    the basic components of x are computed to satisfy */
/*                     Ax - s = b; that is  (A -I)*x = b. Then a row */
/*                     check is performed to see how well  (A -I)*x = b */
/*                     is satisfied.  y is set to be the row residuals, */
/*                     y = b - Ax + s,  and the row error is norm(y). */

/*       1 (GetRes)    just get the row error. */

/*     The row error is a measure of how well x satisfies (A -I)*x = b. */

/*     18 Nov 1991: First version of s5setx based on Minos routine m5setx. */
/*     12 Jul 1997: Thread-safe version. */
/*     25 Jul 2003: Realized dx can sometimes be a NAN (or INF) but */
/*                  norm(NAN) > 1/eps is always false. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --y1;
    --y;
    --x;
    --xbs;
    --kbs;
    --acol;
    --inda;
    --loca;
    --rhs;
    --iw;
    --rw;

    /* Function Body */
    eps0 = rw[2];
/* eps**(4/5)       IEEE DP  3.00e-13 */
    tolrow = rw[61];
/* tolerance for the row error. */
    lprdbg = iw[85];
/* > 0    => private debug print */
    *iexit = 0;
    s5bsx_(&c__0, nbs, nb, &kbs[1], &x[1], &xbs[1]);
    xnorm = dnormi_(nbs, &xbs[1], &c__1);
    dxnorm = 0.;
    goodx = TRUE_;
/*     ------------------------------------------------------------------ */
/*     Compute row residuals  y  =  rhs - (A -I)*x. */
/*     The slack columns are done separately. */
/*     ------------------------------------------------------------------ */
    if (*nrhs > 0) {
	dcopy_(nrhs, &rhs[1], &c__1, &y[1], &c__1);
    }
    if (*nrhs < *m) {
	i__1 = *m - *nrhs;
	dload_(&i__1, &c_b5, &y[*nrhs + 1], &c__1);
    }
    s2aprod_(&c__0, &eps0, nea, nloca, &loca[1], &inda[1], &acol[1], &c_b159, 
	    &x[1], n, &c_b69, &y[1], m);
    daxpy_(m, &c_b69, &x[*n + 1], &c__1, &y[1], &c__1);
/*     ------------------------------------------------------------------ */
/*     Do a row check, perhaps after recomputing the basic x. */
/*     ----------------------------------------------------------------- */
    if (*task == 0) {
/*        ================================================ */
/*        Extract xBS, the basics and superbasics, from x. */
/*        See if iterative refinement is worth doing. */
/*        ================================================ */
	*rowerror = dnormj_(m, &y[1], &c__1);
	bigres = *rowerror > eps0;
	if (bigres) {
/*           ------------------------------------------------------------ */
/*           Compute a correction to basic x from  B*y1 = y. */
/*           Extract the basic and superbasic variables from x. */
/*           Set basic x = x + y1. */
/*           Store the new basic variables in x. */
/*           ------------------------------------------------------------ */
	    s2bsol_(iexit, &c__1, m, &y[1], &y1[1], &iw[1], leniw, &rw[1], 
		    lenrw);
	    if (*iexit != 0) {
		return 0;
	    }
	    dxnorm = dnormj_(m, &y1[1], &c__1);
	    goodx = dxnorm * eps0 <= 1.;
/* false if dxNorm = inf, nan */
	    if (goodx) {
		daxpy_(m, &c_b69, &y1[1], &c__1, &xbs[1], &c__1);
		s5bsx_(&c__1, m, nb, &kbs[1], &x[1], &xbs[1]);
/*              Compute  y  =  rhs  -  (A -I)*x  again for the new x. */
		if (*nrhs > 0) {
		    dcopy_(nrhs, &rhs[1], &c__1, &y[1], &c__1);
		}
		if (*nrhs < *m) {
		    i__1 = *m - *nrhs;
		    dload_(&i__1, &c_b5, &y[*nrhs + 1], &c__1);
		}
		s2aprod_(&c__0, &eps0, nea, nloca, &loca[1], &inda[1], &acol[
			1], &c_b159, &x[1], n, &c_b69, &y[1], m);
		daxpy_(m, &c_b69, &x[*n + 1], &c__1, &y[1], &c__1);
	    } else {
		*iexit = 11;
/* big dx (may be nan or inf) */
	    }
	}
/* BigRes */
    }
/*     Find the norm of xBS, the basic and superbasic x. */
/*     Find the maximum row residual. */
/* Task .eq. Reset */
    imax = jdamax_(m, &y[1], &c__1);
    if (imax > 0) {
	rmax = (d__1 = y[imax], abs(d__1));
	*rowerror = rmax / (xnorm + 1.);
    } else {
	imax = -imax;
	rmax = dnormj_(m, &y[1], &c__1);
/* = flmax! */
	*rowerror = rmax;
    }
    bigres = *rowerror > tolrow;
    if (bigres) {
	*iexit = 10;
    }
    if (*iexit > 0 || lprdbg >= 2) {
	s_wsfi(&io___241);
	do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&imax, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&xnorm, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&dxnorm, (ftnlen)sizeof(doublereal));
	e_wsfi();
	snprnt_(&c__21, str, &iw[1], leniw, (ftnlen)110);
	s_wsfi(&io___242);
	do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&imax, (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__22, str, &iw[1], leniw, (ftnlen)110);
    }
    return 0;
} /* s5setx_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5setx */
/* Subroutine */ int s5step_(integer *nbs, integer *ndegen, doublereal *
	featol, doublereal *infbnd, doublereal *stepmax, doublereal *tolinc, 
	doublereal *tolpiv, integer *feastype, doublereal *blbs, doublereal *
	bubs, doublereal *xbs, doublereal *pbs, logical *hitlow, logical *
	move, logical *onbound, logical *unbounded, integer *infpiv, integer *
	kp, doublereal *bound, doublereal *exact, doublereal *stepb, 
	doublereal *stepp)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j;
    static doublereal res, delta;
    static integer jhitf, jhiti;
    static doublereal stepi;
    static integer jtype;
    static doublereal pivot;
    static logical blockf, blocki;
    static doublereal pivabs, inflow, infupp, pivmaxf, pivmaxi, stepmin;

/*     ================================================================== */
/*     s5step  finds a steplength  stepB  such that  xBS + stepB*pBS */
/*     reaches one of the bounds on  xBS. */

/*     In this  version of s5step,  when x  is infeasible, the  number of */
/*     infeasibilities  will never  increase.   If the  number stays  the */
/*     same, the  sum of  infeasibilities will  decrease.  If  the number */
/*     decreases by one or more,  the sum of infeasibilities will usually */
/*     decrease also, but  occasionally it will increase  after the step- */
/*     length stepB is taken.  (Convergence  is still assured because the */
/*     number has decreased.) */

/*     Two possible steps are computed as follows: */

/*     stepF = the maximum step that can be taken without violating */
/*             one of the bounds that are currently satisfied. */

/*     stepI = the maximum (nonzero) step that has the property of */
/*             reaching a  bound that  is currently violated,  subject to */
/*             the  pivot being  reasonably  close to  the maximum  pivot */
/*             among infeasible variables.  (stepI is not defined if x is */
/*             feasible.) */

/*     stepI  is needed  occasionally when  infeasible, to  prevent going */
/*     unnecessarily far when stepF is  quite large.  It will always come */
/*     into  effect when  x  is about  to become  feasible.   The sum  of */
/*     infeasibilities will  decrease initially  as stepB  increases from */
/*     zero, but may start increasing for larger steps.  choosing a large */
/*     stepI allows several elements of x  to become feasible at the same */
/*     time. */

/*     In the end, we take  stepB = stepF  if x is feasible, or if */
/*     stepI > stepP (where stepP is the perturbed step from pass 1). */
/*     Otherwise,  we take  stepB = stepI. */

/*     Input parameters */
/*     ---------------- */
/*     nBS      is  m + 1  for s5lpit,  m + nBS  for s5QPit. */
/*     stepMax  defines what should be treated as an unbounded step. */
/*     infBnd   provides insurance for detecting unboundedness. */
/*              if stepB reaches a bound as large as infBnd, it is */
/*              classed as an unbounded step. */
/*     tolpiv   is a tolerance to exclude negligible elements of pBS. */
/*     featol   is the current feasibility tolerance used by s5QP. */
/*              Typically in the range 0.5*tolx to 0.99*tolx, */
/*              where tolx is the featol specified by the user. */
/*     tolinc   is used to determine stepMin (see below), the minimum */
/*              positive step. */
/*     feasType is set by  s5Inf  as follows: */
/*              feasType(j) = -2  if x(j) .lt. bl(j) - featol */
/*                          =  0  if x(j)  is feasible */
/*                          = +2  if x(j) .gt. bu(j) + featol */
/*     blBS     the lower bounds on the basic and superbasic variables. */
/*     buBS     the upper bounds on ditto. */
/*     xBS      the values of       ditto. */
/*     pBS      the search direction for the basics and superbasics. */


/*     Output parameters */
/*     ----------------- */
/*     HitLow = true  if a lower bound restricted stepB. */
/*            = false otherwise. */
/*     Move   = true  if exact ge stepMin (defined at end of code). */
/*     OnBound = true  if stepB =  exact. */
/*                    this means that the step stepB moves xBS(kp) */
/*                    exactly onto one of its bounds, namely bound. */
/*            = false if the exact step would be too small */
/*                    ( exact lt stepMin ). */
/*              (with these definitions,  Move = OnBound). */
/*     Unbounded = true  if stepB = stepMax.  kp may possibly be zero. */
/*              the parameters HitLow, Move, OnBound, bound and exact */
/*              should not be used. */
/*     infPiv = the number of  indices such  that |pBS(i)| <  tolpiv and */
/*              xBS(i) + stepP*pBS(i) is infeasible by more than featol. */
/*     kp     = the index (if any) such that xBS(kp) reaches a bound. */
/*     bound  = the bound value blBS(kp) or buBS(kp) corresponding */
/*              to HitLow. */
/*     exact  = the step that would take xBS(kp) exactly onto bound. */
/*     stepB  = an allowable, positive steplength. */
/*              if Unbounded is true,  stepB = stepMax. */
/*              otherwise,          stepB = max(stepMin, exact). */
/*     stepP  = the perturbed steplength from pass 1. */

/*     07 Nov 1991: First version based on Minos routine m5chzr. */
/*     27 Dec 2003: infPiv added to monitor unwanted infeasibilities. */
/*     12 Dec 2012: Recoded to remove goto statements. */
/*     13 Dec 2012: Fixed definition of infPiv. */
/*     25 May 2013: Fixed missing else block giving undefined values of */
/*                  kp and HitLow. */
/*     26 May 2013: infBnd used to identify infinite bounds. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --pbs;
    --xbs;
    --bubs;
    --blbs;
    --feastype;

    /* Function Body */
    infupp = *infbnd * .99;
    inflow = *infbnd * -.99;
/* ----------------------------------------------------------------- */
/* First pass. */
/* For feasible variables, find the steplength stepP that reaches */
/* the nearest perturbed (expanded) bound.  stepP will be slightly */
/* larger than the step to the nearest true bound. */
/* For infeasible variables, find the maximum pivot pivmaxI. */
/* ----------------------------------------------------------------- */
    delta = *featol;
    *stepp = *stepmax;
    pivmaxi = 0.;
    jhitf = 0;
    i__1 = *nbs;
    for (j = 1; j <= i__1; ++j) {
	pivot = pbs[j];
	pivabs = abs(pivot);
	if (pivabs > *tolpiv) {
	    jtype = feastype[j];
	    if (pivot < 0.) {
/* x is decreasing. */
		if (jtype >= 0) {
/* The lower bound (if any) is satisfied. */
/* Test for a smaller stepP. */
		    if (blbs[j] > inflow) {
			res = xbs[j] - blbs[j] + delta;
			if (*stepp * pivabs > res) {
			    *stepp = res / pivabs;
			    jhitf = j;
			}
		    }
/* If the upper bound is violated, */
/* test if this variable has a bigger pivot. */
		    if (jtype > 0) {
			pivmaxi = max(pivmaxi,pivabs);
		    }
		}
	    } else {
/* x is increasing. */
		if (jtype <= 0) {
/* The upper bound (if any) is satisfied. */
/* Test for a smaller stepP. */
		    if (bubs[j] < infupp) {
			res = bubs[j] - xbs[j] + delta;
			if (*stepp * pivabs > res) {
			    *stepp = res / pivabs;
			    jhitf = j;
			}
		    }
/* If the lower bound is violated, */
/* test if this variable has a bigger pivot. */
		    if (jtype < 0) {
			pivmaxi = max(pivmaxi,pivabs);
		    }
		}
	    }
	}
    }
/* ----------------------------------------------------------------- */
/* Second pass. */
/* For feasible variables, compute the steps without perturbation. */
/* Choose the largest pivot element subject to the step being no */
/* greater than stepP. */
/* For infeasible variables, find the largest step subject to */
/* the pivot element being no smaller than gamma*pivmaxI. */
/* ----------------------------------------------------------------- */
    stepi = 0.;
    pivmaxf = 0.;
    pivmaxi *= .001;
    jhiti = 0;
    *infpiv = 0;
    i__1 = *nbs;
    for (j = 1; j <= i__1; ++j) {
	pivot = pbs[j];
	pivabs = abs(pivot);
	jtype = feastype[j];
	if (pivabs > *tolpiv) {
	    if (pivot < 0.) {
/* x is decreasing. */
		if (jtype >= 0) {
/* The lower bound (if any) is satisfied. */
/* Test for a bigger pivot. */
		    if (pivabs > pivmaxf && blbs[j] > inflow) {
			res = xbs[j] - blbs[j];
			if (pivabs * *stepp >= res) {
			    pivmaxf = pivabs;
			    jhitf = j;
			}
		    }
		    if (jtype > 0) {
/* An upper bound is present and violated. */
/* Test for a bigger stepI. */
			if (pivabs >= pivmaxi) {
			    res = xbs[j] - bubs[j];
			    if (pivabs * stepi < res) {
				stepi = res / pivabs;
				jhiti = j;
			    }
			}
		    }
		}
	    } else {
/* x is increasing. */
		if (jtype <= 0) {
/* The upper bound (if any) is satisfied. */
/* Test for a bigger pivot */
		    if (pivabs > pivmaxf && bubs[j] < infupp) {
			res = bubs[j] - xbs[j];
			if (pivabs * *stepp >= res) {
			    pivmaxf = pivabs;
			    jhitf = j;
			}
		    }
		    if (jtype < 0) {
/* A lower bound is present and violated. */
/* Test for a bigger stepI. */
			if (pivabs >= pivmaxi) {
			    res = blbs[j] - xbs[j];
			    if (pivabs * stepi < res) {
				stepi = res / pivabs;
				jhiti = j;
			    }
			}
		    }
		}
	    }
	} else if (jtype == 0 && pivabs > 0.) {
/* Feasible variable with a negligible pivot. */
/* Check the step makes us infeasible by no more than stepP. */
	    if (pivot < 0.) {
/* x is decreasing. */
		if (blbs[j] > inflow) {
		    res = xbs[j] - blbs[j] + delta;
		    if (*stepp * pivabs > res) {
			++(*infpiv);
		    }
		}
	    } else {
/* x is increasing. */
		if (bubs[j] < infupp) {
		    res = bubs[j] - xbs[j] + delta;
		    if (*stepp * pivabs > res) {
			++(*infpiv);
		    }
		}
	    }
	}
    }
/* ----------------------------------------------------------------- */
/* See if a feasible and/or infeasible variable blocks. */
/* ----------------------------------------------------------------- */
    blockf = jhitf > 0;
    blocki = jhiti > 0;
    *unbounded = ! (blockf || blocki);
    if (*unbounded) {
	*stepb = *stepmax;
	*move = TRUE_;
	*onbound = FALSE_;
    } else {
	if (blockf) {
/* ----------------------------------------------------------- */
/* Variable hits a bound for which it is currently feasible. */
/* The step length stepF is not used, so no need to get it, */
/* but we know that stepF .le. stepP, the step from pass 1. */
/* ----------------------------------------------------------- */
	    *kp = jhitf;
	    pivot = pbs[*kp];
	    *hitlow = pivot < 0.;
	} else {
/* BlockI */
	    *kp = jhiti;
	    pivot = pbs[*kp];
	    *hitlow = pivot > 0.;
	}
/* If there is a choice between stepF and stepI, it is probably */
/* best to take stepI (so that the infeasible variable jhitI can */
/* be kicked out of the basis). */
/* However, we can't if stepI is bigger than stepP. */
	if (blocki && stepi <= *stepp) {
	    *kp = jhiti;
	    pivot = pbs[*kp];
	    *hitlow = pivot > 0.;
	}
/* -------------------------------------------------------------- */
/* Try to step exactly onto a bound, but make sure the exact */
/* stepB is sufficiently positive (exact will be either stepF or */
/* stepI). As featol increases by tolinc each iteration, we know */
/* that a step as large as stepMin (below) will not cause  any */
/* feasible variables to become infeasible (where feasibility */
/* is measured by the current featol). */
/* -------------------------------------------------------------- */
	if (*hitlow) {
	    *bound = blbs[*kp];
	} else {
	    *bound = bubs[*kp];
	}
	stepmin = *tolinc / abs(pivot);
	*exact = (*bound - xbs[*kp]) / pivot;
	*stepb = max(stepmin,*exact);
	*onbound = *stepb == *exact;
	*move = *exact >= stepmin;
	if (! (*move)) {
	    ++(*ndegen);
	}
    }
    return 0;
} /* s5step_ */

