/* ../snopt7/src/sn80ncon.f -- translated by f2c (version 20100827).
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
static integer c__23 = 23;
static integer c__0 = 0;
static doublereal c_b12 = 1.;
static doublereal c_b21 = -1.;
static integer c__3 = 3;
static integer c__31 = 31;
static integer c__22 = 22;
static integer c__33 = 33;
static integer c__6 = 6;
static integer c__4 = 4;
static integer c__21 = 21;
static integer c__5 = 5;
static integer c__20 = 20;
static integer c__25 = 25;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*     File  sn80ncon.f */

/*     s8fdSwitch        s8Fv         s8Fx */
/*     s8Gcopy           s8getFeasLC  s8getR     s8getWeights  s8Gloc */
/*     s8Gprod           s8Gsize */
/*     s8Infs            s8InitH      s8InitPen */
/*     s8Merit */
/*     s8optimizeSlacks */
/*     s8HxLP            s8HxPP       s8HxQP     s8HxNull */
/*     s8rand            s8rc */
/*     s8scaleG          s8scaleJ     s8sInf     s8solveQP   s8stepLimits */


/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Subroutine */ int s8fdswitch_(integer *nncon0, integer *nncon, integer *
	nnobj, integer *itn, integer *cditns, logical *goodg, logical *
	needderivs, logical *usefd, doublereal *dualinf, doublereal *fcon, 
	doublereal *fobj, integer *iw, integer *leniw, doublereal *rw, 
	integer *lenrw)
{
    /* Format strings */
    static char fmt_1000[] = "(\002 Itn\002,i7,\002 -- Central differences i"
	    "nvoked.\002,\002  Small reduced gradient.\002)";

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static char str[80];
    static doublereal epsrf, cnorm, fdint1;
    extern doublereal dnrm1s_(integer *, doublereal *, integer *);
    static doublereal rgnorm, rgtest;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);
    static logical central;
    static doublereal objsize;

    /* Fortran I/O blocks */
    static icilist io___9 = { 0, str, 0, fmt_1000, 80, 1 };


/*     ================================================================== */
/*     s8fdSwitch   controls the switch between forward and central */
/*     differences. */

/*     If the forward-difference estimate of the reduced gradient of the */
/*     Lagrangian is small,  a switch is made to central differences. */
/*     In this case, the derivatives are recomputed and the QP is solved */
/*     again. */

/*     On the other hand, if central differences have produced a large */
/*     reduced-gradient norm, switch back to forward differences. */

/*     31 Mar 2000: First version of s8fdSwitch written for SNOPT 6.1. */
/*     03 Aug 2003: snPRNT adopted. */
/*     03 Aug 2003: Current version of s8fdSwitch. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/* forwd diffs or cntrl diffs */
/* infoTag(6) */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --fcon;
    --iw;
    --rw;

    /* Function Body */
    epsrf = rw[73];
/* relative function precision. */
    fdint1 = rw[76];
/* (1) forwrd diff. interval */
    central = iw[181] == 2;
    if (*nncon == 0) {
	cnorm = 0.;
    } else {
	cnorm = dnrm1s_(nncon, &fcon[1], &c__1);
    }
    if (*nnobj == 0) {
	objsize = 0.;
    } else {
	objsize = abs(*fobj);
    }
    *goodg = TRUE_;
    rgtest = (objsize + 1. + cnorm) * epsrf / fdint1;
    rgnorm = *dualinf;
    if (central) {
	if (rgnorm > rgtest * 10. && *cditns > 0) {
	    iw[181] = 1;
	    central = FALSE_;
	    if (*usefd) {
		iw[242] = 1;
/* Forward differences */
	    }
	}
    } else {
	if (rgnorm <= rgtest) {
	    *cditns = 0;
	    iw[181] = 2;
	    if (*usefd) {
		*goodg = FALSE_;
		*needderivs = TRUE_;
		iw[242] = 2;
/* Central differences */
		s_wsfi(&io___9);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
	    }
	}
    }
    return 0;
} /* s8fdswitch_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8fdSwitch */
/* Subroutine */ int s8fv_(logical *elastic, integer *n, integer *nncon, 
	doublereal *tolz, doublereal *wtinf, doublereal *bl, doublereal *bu, 
	doublereal *fv, doublereal *x, doublereal *ycon, doublereal *fx)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static integer i__, j;
    static doublereal xj, blj, buj, fvi, fxi, fvl, fvu, yconi, yconv, yconw;

/*     ================================================================== */
/*     s8Fv  computes the vector of nonlinear constraint violations: */
/*        Fv = fCon + A(linear)*x - (nonlinear slacks) */

/*     If the Lagrange multiplier is zero, the violation can be set to */
/*     any value without changing the merit function.  In this case we */
/*     try and set the slack so that Fv is zero (subject to the slack */
/*     being feasible). */

/*     In elastic mode we implicitly adjust the variables v and w such */
/*     that   c - s(feas) + v - w = 0,  with  v >= 0  and  w >= 0. */

/*     On entry, */
/*        x   =  the current x. */
/*        Fx  =  fCon + A(linear)*x,   defined in s8Fx. */

/*     On exit, */
/*        x   =  x containing the modified slacks. */
/*        Fv  =  fCon + A(linear)*x -  slacks. */
/*        Fx  =  unaltered. */

/*     19 Apr 2001: First version. */
/*     19 Apr 2001: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --fx;
    --ycon;
    --x;
    --fv;
    --bu;
    --bl;

    /* Function Body */
    i__1 = *nncon;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = *n + i__;
	xj = x[j];
	fxi = fx[i__];
	fvi = fxi - xj;
	yconi = ycon[i__];
	blj = bl[j];
	buj = bu[j];
	fvu = fxi - buj;
	fvl = fxi - blj;
	yconv = (d__1 = *wtinf - yconi, abs(d__1));
/* Multiplier for v in elastic mode */
	yconw = (d__1 = *wtinf + yconi, abs(d__1));
/* Multiplier for w in elastic mode */
	if (*elastic && xj <= blj && yconv <= *tolz) {
	    if (fvi > 0.) {
		fvi = max(0.,fvl);
	    } else {
		fvi = 0.;
	    }
	} else if (*elastic && xj >= buj && yconw <= *tolz) {
	    if (fvi < 0.) {
		fvi = min(0.,fvu);
	    } else {
		fvi = 0.;
	    }
	} else {
	    if (yconi <= *tolz && fvi > 0.) {
		fvi = max(0.,fvu);
	    } else if (yconi >= -(*tolz) && fvi < 0.) {
		fvi = min(0.,fvl);
	    }
	}
	xj = fxi - fvi;
	fv[i__] = fvi;
	x[j] = xj;
    }
    return 0;
} /* s8fv_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8Fv */
/* Subroutine */ int s8fx_(integer *n, integer *nncon, integer *nnjac, 
	doublereal *tolz, integer *nej, integer *nlocj, integer *locj, 
	integer *indj, doublereal *jcol, doublereal *fcon, doublereal *x, 
	doublereal *fx)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer nlin;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), s2aprod_(integer *, doublereal *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *);

/*     ================================================================== */
/*     s8Fx  defines the nonlinear constraint values */
/*       Fx  =  true nonlinear slack = fCon + A(linear)*x, */

/*     09 Jan 1992: First version based on Minos routine m8viol. */
/*     16 Nov 1998: Norm x changed to include only columns. */
/*     21 Oct 2000: Made compatible with SNOPT 6.1 */
/*     21 Oct 2000: Current version of s8Fx */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
/*     Compute the nonlinear constraint value. */
/*     Set  Fx  =  fCon + (linear A)*x,   excluding slacks. */
    /* Parameter adjustments */
    --fx;
    --x;
    --fcon;
    --jcol;
    --indj;
    --locj;

    /* Function Body */
    dcopy_(nncon, &fcon[1], &c__1, &fx[1], &c__1);
    nlin = *n - *nnjac;
    if (nlin > 0) {
	i__1 = nlin + 1;
	s2aprod_(&c__0, tolz, nej, &i__1, &locj[*nnjac + 1], &indj[1], &jcol[
		1], &c_b12, &x[*nnjac + 1], &nlin, &c_b12, &fx[1], nncon);
    }
    return 0;
} /* s8fx_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8Fx */
/* Subroutine */ int s8gcopy_(integer *nncon, integer *nnjac, integer *nej, 
	integer *nlocj, integer *locj, integer *indj, integer *neg1, integer *
	nlocg1, integer *locg1, doublereal *g1, integer *neg2, integer *
	nlocg2, integer *locg2, doublereal *g2)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer j, k, l1, l2, ir;

/*     ================================================================== */
/*     s8Gcopy  copies G1 into G2 when either  G1 or  G2 */
/*     is stored in the upper-left hand corner of J. */

/*     16 Sep 1993: First version. */
/*     26 Oct 2000: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --indj;
    --locj;
    --g1;
    --locg1;
    --g2;
    --locg2;

    /* Function Body */
    i__1 = *nnjac;
    for (j = 1; j <= i__1; ++j) {
	l1 = locg1[j];
	l2 = locg2[j];
	i__2 = locj[j + 1] - 1;
	for (k = locj[j]; k <= i__2; ++k) {
	    ir = indj[k];
	    if (ir > *nncon) {
		goto L100;
	    }
	    g2[l2] = g1[l1];
	    ++l1;
	    ++l2;
	}
L100:
	;
    }
    return 0;
} /* s8gcopy_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8Gcopy */
/* Subroutine */ int s8getfeaslc_(integer *iexit, integer *starttype, U_fp 
	mnrlog, integer *lenr, integer *m, integer *maxs, integer *mbs, 
	integer *n, integer *nb, integer *nncon0, integer *nncon, integer *
	nnh0, integer *nnh, integer *ndegen, integer *ns, integer *numlc, 
	integer *numliq, integer *itn, integer *itnlim, integer *itqp, 
	integer *mnrprtlvl, doublereal *scaleobj, doublereal *toloptqp, 
	doublereal *tolx, integer *ninf, doublereal *sinf, integer *ninfe, 
	doublereal *sinfe, doublereal *wtinf, doublereal *pinorm, doublereal *
	rgnorm, integer *nej, integer *nlocj, integer *locj, integer *indj, 
	doublereal *jcol, integer *neh, integer *nloch, integer *loch, 
	integer *indh, doublereal *hcol, integer *etype, integer *estate, 
	integer *feastype, integer *hs, integer *kbs, doublereal *bl, 
	doublereal *bu, doublereal *blqp, doublereal *buqp, doublereal *blbs, 
	doublereal *bubs, doublereal *blsave, doublereal *busave, doublereal *
	gbs, doublereal *gqp, doublereal *hdx, doublereal *pbs, doublereal *
	pi, doublereal *r__, doublereal *rc, doublereal *rg, doublereal *
	qprhs, doublereal *scales, doublereal *x0, doublereal *x, doublereal *
	xbs, integer *iy, integer *iy1, doublereal *y, doublereal *y1, 
	doublereal *y2, char *cw, integer *lencw, integer *iw, integer *leniw,
	 doublereal *rw, integer *lenrw, ftnlen cw_len)
{
    /* Format strings */
    static char fmt_8000[] = "(\002 Itn\002,i7,\002: Feasible linear rows"
	    "\002)";
    static char fmt_8100[] = "(\002 Itn\002,i7,\002: PP\002,i1,\002.  Minimi"
	    "zing  Norm(x-x0)\002)";
    static char fmt_8200[] = "(\002 Itn\002,i7,\002: PP\002,i1,\002.  Norm(x"
	    "-x0) approximately minimized  (\002,1p,e8.2,\002)\002)";
    static char fmt_8300[] = "(\002 Itn\002,i7,\002: PP1.  Making nonlinear "
	    "variables feasible\002)";
    static char fmt_8400[] = "(\002 Itn\002,i7,\002: PP1. \002,i7,\002 nonli"
	    "near variables made feasible\002)";

    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static doublereal hcondbnd;
    extern /* Subroutine */ int s8hxnull_();
    static integer elastics, minimize, probtype;
    static doublereal toloptfp, toloptpp;
    static integer j;
    static doublereal x0j;
    static integer itqptarget;
    static char str[80];
    static doublereal eps0, eps2;
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
	    s5qp_(integer *, integer *, char *, integer *, U_fp, U_fp, U_fp, 
	    integer *, integer *, logical *, logical *, logical *, integer *, 
	    logical *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, integer *, integer *, doublereal *, doublereal *, 
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
    static integer suboptimize;
    static doublereal obja;
    static integer eigh;
    static logical gotr;
    extern /* Subroutine */ int iload_(integer *, integer *, integer *, 
	    integer *);
    static integer emode;
    static logical needx;
    static doublereal objpp;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer nviol;
    extern /* Subroutine */ int s8hxpp_(), s8hxqp_();
    static doublereal zcndbd;
    static logical needlu;
    static integer nobjfp, iobjpp, nobjpp, inform__, lvlppm, typelu;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);
    static integer nobjfp0, nobjpp0;
    extern /* Subroutine */ int s2aprod_(integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *);
    static logical elastic;
    static integer hvcalls, lvlobje;
    static char probtag[20];
    static integer msbsave, itqpmax;

    /* Fortran I/O blocks */
    static icilist io___55 = { 0, str, 0, fmt_8000, 80, 1 };
    static icilist io___56 = { 0, str, 0, fmt_8100, 80, 1 };
    static icilist io___62 = { 0, str, 0, fmt_8200, 80, 1 };
    static icilist io___63 = { 0, str, 0, fmt_8300, 80, 1 };
    static icilist io___64 = { 0, str, 0, fmt_8400, 80, 1 };
    static icilist io___69 = { 0, str, 0, fmt_8200, 80, 1 };


/*     ================================================================== */
/*     s8getfeasLC finds a feasible point for the linear constraints. */

/*     On entry: */
/*     --------- */
/*     startType     is 0,1,2,3: an integer version of the solver's Start */
/*                   (character for sqopt, snoptb, snoptc, npopt */
/*                   integer    for snopta). */

/*     nInf          is the number of infeasible equality rows. */

/*     nInfE         is zero (elastic mode is turned off in s5getB). */

/*     kBS(m+1:m+nS) contains a basis specified by nS, hs(*) and x(*). */
/*                   In particular, there must be nS values hs(j) = 2, */
/*                   and the corresponding j's must be listed in */
/*                   kBS(m+1:m+nS). If nInf = 0, the basis is infeasible. */

/*     bl, bu        contain the true bounds (optionally scaled). */

/*     blSave, buSave are copies of the (optionally scaled) upper */
/*                   and bounds set in s5getB. Entries for the nonlinear */
/*                   rows are relaxed. */

/*     blQP, buQP    contain the working lower and upper bounds. */


/*      iExit       Result */
/*      -----       ------ */
/*       >0         Fatal error */
/*        0         Feasible point found */

/*     11 May 1994: First version of s8feasLC. */
/*     19 Aug 1996: First minsum version. */
/*     05 Feb 1998: Proximal point norm changed to one-norm. */
/*     23 Dec 1999: Optional Proximal Point methods 0 and  2 added. */
/*     03 Aug 2003: snPRNT and snEXIT adopted. */
/*     16 May 2006: Explicit target itQP added. */
/*     18 Jun 2008: Hdx, pBS and rg added as arguments. */
/*     22 Feb 2015: First FP call to s5QP instead of s5LP (for nS > 0). */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/* max # of new superbasics */
/* >0 => Minor heading for iPrint */
/* >0 => Minor heading for iSumm */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --r__;
    --pi;
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
    --x;
    --x0;
    --scales;
    --rc;
    --pbs;
    --busave;
    --blsave;
    --buqp;
    --blqp;
    --bu;
    --bl;
    --hs;
    --estate;
    --etype;
    --qprhs;
    --hdx;
    --gqp;
    --jcol;
    --indj;
    --locj;
    --hcol;
    --indh;
    --loch;
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    lvlppm = iw[79];
/* 1(2)-norm proximal point method for x0 */
    eps0 = rw[2];
/* eps**(4/5) */
    eps2 = rw[4];
/* eps**(1/2) */
    hcondbnd = rw[85];
/* bound on the condition of ZHZ */
    zcndbd = rw[86];
/* bound on the condition of Z */
    obja = 0.;
    iobjpp = 0;
    minimize = 1;
/* Local value */
    needlu = TRUE_;
    needx = needlu;
/*     Set the QP rhs to make x satisfy the (relaxed) nonlinear rows. */
/*     The array  QPrhs  contains the rhs. */
/*     Use a fairly tight optimality tolerance for phase 1. */
    if (*nncon > 0) {
	dcopy_(nncon, &x[*n + 1], &c__1, &qprhs[1], &c__1);
	s2aprod_(&c__0, &eps0, nej, nlocj, &locj[1], &indj[1], &jcol[1], &
		c_b12, &x[1], n, &c_b21, &qprhs[1], nncon);
    }
    if (*numliq > 0 || *ninf > 0) {
/*        --------------------------------------------------------------- */
/*        s5getB has found a basis for the linear constraints that may */
/*        or may not be feasible.  Find a feasible point for all the */
/*        linear constraints.  If none exists, minimize the sum of */
/*        infeasibilities of the linear rows, subject to the bounds. */
/*        --------------------------------------------------------------- */
/*        Make the bounds on the linear rows elastic. */
/*        The bounds on the linear rows are elastic */
	iload_(numlc, &c__3, &etype[*n + *nncon + 1], &c__1);
	if (*ninf == 0) {
/*           s5getB found feasible linear equality constraints. */
	    elastic = FALSE_;
/* Start in normal mode */
	    emode = 1;
/* Enter elastic mode if infeasible */
	    probtype = 0;
	} else {
/*           s5getB found infeasible linear equality constraints. */
/*           Minimize the sum of infeasibilities. */
	    elastic = TRUE_;
/* Start in elastic mode */
	    emode = 2;
/* elastic mode throughout */
	    elastics = 0;
	    probtype = 1;
	}
	lvlobje = 2;
/* W1 = 0, W2 = 1   W1*trueObj + W2*sInfE */
/* In elastic mode, use: */
	s_copy(probtag, "linear rows", (ftnlen)20, (ftnlen)11);
	suboptimize = -1;
	toloptfp = max(1e-6,eps2);
/* Default optimality tol. */
	eigh = 0;
	iw[223] = 1;
/* Switch to QP print   heading */
	iw[225] = 1;
/* Switch to QP summary heading */
	nobjfp = 0;
/* No explicit gradient in proximal point */
	nobjfp0 = 1;
	gotr = FALSE_;
	typelu = 3;
	hvcalls = 0;
	s5qp_(&inform__, &probtype, probtag, &suboptimize, (U_fp)mnrlog, (
		U_fp)s8hxnull_, (U_fp)s8hxnull_, &hvcalls, &eigh, &elastic, &
		gotr, &needlu, &typelu, &needx, lenr, m, maxs, mbs, n, nb, 
		ndegen, nnh0, nnh, &nobjfp0, &nobjfp, nnh0, nnh, ns, itqp, 
		itnlim, itnlim, itn, &emode, &lvlobje, mnrprtlvl, &minimize, &
		iobjpp, scaleobj, &obja, &objpp, &hcondbnd, &zcndbd, &
		toloptfp, toloptqp, tolx, ninf, sinf, &elastics, ninfe, sinfe,
		 wtinf, pinorm, rgnorm, nej, nlocj, &locj[1], &indj[1], &jcol[
		1], neh, nloch, &loch[1], &indh[1], &hcol[1], &etype[1], &
		estate[1], &feastype[1], &hs[1], &kbs[1], &bl[1], &bu[1], &
		blqp[1], &buqp[1], &blbs[1], &bubs[1], &gbs[1], &gqp[1], &gqp[
		1], &hdx[1], &pbs[1], &pi[1], &r__[1], &rc[1], &rg[1], nncon0,
		 nncon, &qprhs[1], &scales[1], nnh0, nnh, &x0[1], &x[1], &xbs[
		1], &x0[1], &iy[1], &iy1[1], &y[1], &y1[1], &y2[1], cw + 8, 
		lencw, &iw[1], leniw, &rw[1], lenrw, cw + 8, lencw, &iw[1], 
		leniw, &rw[1], lenrw, (ftnlen)20, (ftnlen)8, (ftnlen)8);
/*        Check for trouble in s5QP. */
/*        iExit        Status */
/*        -----        ------ */
/*         -3          Too many iterations */
/*         -2          Phase 1 is unbounded */
/*         -1          infeasible nonelastics */
/*          0          infeasibilities minimized */
/*         >0          Fatal error */

/*        Time to stop if the linear constraints are infeasible. */
/*        If inform = 0, the sum of infeasibilities will have been */
/*        minimized. */
	if (inform__ != 0 || *ninf > 0 || *ninfe > 0) {
	    if (inform__ > 0) {
		*iexit = inform__;
/* Fatal error */
	    } else if (inform__ == -3) {
		*iexit = 31;
/* iterations limit */
	    } else if (*ninf > 0) {
		*iexit = 11;
/* infeasible linear constraints */
	    } else if (*ninfe > 0) {
		*iexit = 14;
/* infeasible linear constraints */
	    }
	    if (*iexit != 0) {
		goto L800;
	    }
	}
/*        Print something brief if s5LP didn't already do so. */
	if (*mnrprtlvl >= 1) {
	    s_wsfi(&io___55);
	    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__31, str, &iw[1], leniw, (ftnlen)80);
	    snprnt_(&c__22, str, &iw[1], leniw, (ftnlen)80);
	}
	needlu = FALSE_;
    }
/*     ------------------------------------------------------------------ */
/*     x is feasible for the linear constraints. */
/*     The linear rows are not allowed to go infeasible again. */
/*     ------------------------------------------------------------------ */
    iload_(numlc, &c__0, &etype[*n + *nncon + 1], &c__1);
    if (lvlppm > 0 && *nnh > 0 && *starttype == 0) {
/*        =============================================================== */
/*        Find a feasible point closest to x0. */
/*        Minimize norm(x - x0). */
/*        =============================================================== */
	if (*mnrprtlvl >= 1) {
	    s_wsfi(&io___56);
	    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&lvlppm, (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
	}
	if (lvlppm == 1) {
/*           ------------------------------------------------------------ */
/*           Minimize the one-norm of (x-x0) by fixing the nonlinear */
/*           variables so that bl = x0 = bu.  Any bl or bu that is moved */
/*           to  x0  is made elastic. */
/*           ------------------------------------------------------------ */
	    i__1 = *nnh;
	    for (j = 1; j <= i__1; ++j) {
		if (bl[j] == bu[j]) {
/*                 Relax */
		} else {
		    x0j = x0[j];
		    bl[j] = x0j;
		    bu[j] = x0j;
		    blqp[j] = bl[j];
		    buqp[j] = bu[j];
		    etype[j] = 3;
		    if (hs[j] <= 1) {
			x[j] = x0j;
		    }
		}
	    }
	    elastic = FALSE_;
/* Start in normal mode */
	    emode = 1;
/* Enter elastic mode if infeasible */
	    lvlobje = 2;
/* 0*trueObj + 1*sInfE */
/* In elastic mode, use: */
	    s_copy(probtag, "norm(x-x0) problem  ", (ftnlen)20, (ftnlen)20);
	    iw[223] = 1;
/* New LP print   header */
	    iw[225] = 1;
/* New LP summary header */
	    needx = TRUE_;
	    suboptimize = -1;
	    toloptfp = .01;
/* Sloppy phase1 optimality tol for */
	    toloptpp = .01;
/* Sloppy phase 2 opt tol */
	    itqpmax = *itqp + 200;
/* Limit the minor iterations */
	    s5lp_(&inform__, &c__1, probtag, &elastic, &suboptimize, (U_fp)
		    mnrlog, &needlu, &needx, m, n, nb, ndegen, itqp, &itqpmax,
		     itn, &emode, &lvlobje, mnrprtlvl, &minimize, &iobjpp, 
		    scaleobj, &obja, &toloptfp, &toloptpp, tolx, ninf, sinf, &
		    elastics, ninfe, sinfe, wtinf, pinorm, rgnorm, nej, nlocj,
		     &locj[1], &indj[1], &jcol[1], &etype[1], &estate[1], &
		    feastype[1], &hs[1], &kbs[1], &bl[1], &bu[1], &blqp[1], &
		    buqp[1], &blbs[1], &bubs[1], &gbs[1], &pi[1], &rc[1], 
		    nncon0, nncon, &qprhs[1], &scales[1], &x[1], &xbs[1], &x0[
		    1], &iy[1], &iy1[1], &y[1], &y1[1], cw + 8, lencw, &iw[1],
		     leniw, &rw[1], lenrw, (ftnlen)20, (ftnlen)8);
/*           Some elastic variables may have moved outside their bounds. */
/*           Count them.  Reset the true bounds. */
/*           If necessary,  get feasible again with the normal tolOptQP. */
	    nviol = 0;
	    i__1 = *nnh;
	    for (j = 1; j <= i__1; ++j) {
		bl[j] = blsave[j];
		bu[j] = busave[j];
		blqp[j] = bl[j];
		buqp[j] = bu[j];
		if (x[j] < bl[j] - *tolx || x[j] > bu[j] + *tolx) {
		    ++nviol;
		}
	    }
/*           Check for errors in s5LP. */
/*           inform values are = -3,-2,-1, 0, >0 */
	    if (inform__ != 0) {
		if (inform__ > 0) {
		    *iexit = inform__;
/* Fatal error */
		} else if (inform__ == -3 && *itn >= *itnlim) {
		    *iexit = 31;
/* iterations limit */
		}
		if (*iexit != 0) {
		    goto L800;
		}
	    }
	    if (inform__ == 0 && *mnrprtlvl >= 1) {
		s_wsfi(&io___62);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&lvlppm, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&(*sinf), (ftnlen)sizeof(doublereal));
		e_wsfi();
		snprnt_(&c__33, str, &iw[1], leniw, (ftnlen)80);
		if (nviol > 0) {
		    s_wsfi(&io___63);
		    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		    e_wsfi();
		    snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
		}
	    }
	    if (nviol > 0) {
/*              Return the infeasible variables to feasibility. */
		s_copy(probtag, "linear rows again   ", (ftnlen)20, (ftnlen)
			20);
		elastic = FALSE_;
		needlu = FALSE_;
		needx = TRUE_;
		suboptimize = -1;
		toloptfp = eps2;
/* Revert to accurate phase 1 opt tol */
		if (inform__ != 0) {
		    needlu = TRUE_;
		}
		s5lp_(&inform__, &c__0, probtag, &elastic, &suboptimize, (
			U_fp)mnrlog, &needlu, &needx, m, n, nb, ndegen, itqp, 
			itnlim, itn, &emode, &lvlobje, mnrprtlvl, &minimize, &
			iobjpp, scaleobj, &obja, &toloptfp, toloptqp, tolx, 
			ninf, sinf, &elastics, ninfe, sinfe, wtinf, pinorm, 
			rgnorm, nej, nlocj, &locj[1], &indj[1], &jcol[1], &
			etype[1], &estate[1], &feastype[1], &hs[1], &kbs[1], &
			bl[1], &bu[1], &blqp[1], &buqp[1], &blbs[1], &bubs[1],
			 &gbs[1], &pi[1], &rc[1], nncon0, nncon, &qprhs[1], &
			scales[1], &x[1], &xbs[1], &x0[1], &iy[1], &iy1[1], &
			y[1], &y1[1], cw + 8, lencw, &iw[1], leniw, &rw[1], 
			lenrw, (ftnlen)20, (ftnlen)8);
/*              Possible inform values are = -3,-2,-1, 0, >0 */
		if (inform__ != 0) {
		    if (inform__ > 0) {
			*iexit = inform__;
/* Fatal error */
		    } else if (inform__ == -3) {
			*iexit = 31;
/* iterations limit */
		    } else if (*ninf > 0) {
			*iexit = 11;
/* infeasible (should not happen here) */
		    }
		    if (*iexit != 0) {
			goto L800;
		    }
		}
		if (inform__ == 0 && *mnrprtlvl >= 1) {
		    s_wsfi(&io___64);
		    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&nviol, (ftnlen)sizeof(integer));
		    e_wsfi();
		    snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
		}
	    }
	    *ninf = 0;
	    *sinf = 0.;
	    *ninfe = 0;
	    *sinfe = 0.;
/*           x(1:nnH) are feasible and must remain so. */
	    iload_(nnh, &c__0, &etype[1], &c__1);
	} else if (lvlppm == 2) {
/*           ------------------------------------------------------------ */
/*           Minimize the two-norm of (x-x0). */
/*           ------------------------------------------------------------ */
	    nobjpp = 0;
/* No explicit gradient in proximal point */
	    nobjpp0 = 1;
	    gotr = FALSE_;
	    typelu = 3;
	    hvcalls = 0;
	    s_copy(probtag, "norm(x-x0) problem  ", (ftnlen)20, (ftnlen)20);
	    if (*nnh < *n) {
		eigh = 0;
	    } else {
		eigh = 1;
	    }
	    iw[223] = 1;
/* Switch to QP print   heading */
	    iw[225] = 1;
/* Switch to QP summary heading */
	    elastic = FALSE_;
/* Start in normal mode */
	    emode = 0;
/* elastic mode not needed. */
	    lvlobje = 0;
/* No objective in elastic mode */
	    needx = FALSE_;
	    itqpmax = *itqp + 100;
/* Limit the minor iterations */
	    itqptarget = *itqp + 100;
	    msbsave = iw[95];
	    iw[95] = 100;
/* and the number of new superbasics */
	    suboptimize = -1;
	    toloptfp = eps2;
	    toloptpp = .01;
/* Sloppy phase 2 opt tol */
	    s5qp_(&inform__, &c__6, probtag, &suboptimize, (U_fp)mnrlog, (
		    U_fp)s8hxpp_, (U_fp)s8hxqp_, &hvcalls, &eigh, &elastic, &
		    gotr, &needlu, &typelu, &needx, lenr, m, maxs, mbs, n, nb,
		     ndegen, nnh0, nnh, &nobjpp0, &nobjpp, nnh0, nnh, ns, 
		    itqp, &itqpmax, &itqptarget, itn, &emode, &lvlobje, 
		    mnrprtlvl, &minimize, &iobjpp, scaleobj, &obja, &objpp, &
		    hcondbnd, &zcndbd, &toloptfp, &toloptpp, tolx, ninf, sinf,
		     &elastics, ninfe, sinfe, wtinf, pinorm, rgnorm, nej, 
		    nlocj, &locj[1], &indj[1], &jcol[1], neh, nloch, &loch[1],
		     &indh[1], &hcol[1], &etype[1], &estate[1], &feastype[1], 
		    &hs[1], &kbs[1], &bl[1], &bu[1], &blqp[1], &buqp[1], &
		    blbs[1], &bubs[1], &gbs[1], &gqp[1], &gqp[1], &hdx[1], &
		    pbs[1], &pi[1], &r__[1], &rc[1], &rg[1], nncon0, nncon, &
		    qprhs[1], &scales[1], nnh0, nnh, &x0[1], &x[1], &xbs[1], &
		    x0[1], &iy[1], &iy1[1], &y[1], &y1[1], &y2[1], cw + 8, 
		    lencw, &iw[1], leniw, &rw[1], lenrw, cw + 8, lencw, &iw[1]
		    , leniw, &rw[1], lenrw, (ftnlen)20, (ftnlen)8, (ftnlen)8);
	    iw[95] = msbsave;
/*           Check for trouble. */
/*           Possible inform values are = -9(-1)-1, 0, >0 */
	    if (inform__ != 0 || *ninf > 0) {
		if (inform__ > 0) {
		    *iexit = inform__;
/* Fatal LU error or time limit */
		} else if (inform__ == -3 && *itn >= *itnlim) {
		    *iexit = 31;
/* iterations limit */
		} else if (inform__ == -1 || *ninf > 0) {
		    *iexit = 11;
/* infeasible (should not happen here) */
		}
		if (*iexit != 0) {
		    goto L800;
		}
	    }
/*           Note: objQP is an updated quantity that may be slightly */
/*           negative. */
	    if (*mnrprtlvl >= 1) {
		s_wsfi(&io___69);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&lvlppm, (ftnlen)sizeof(integer));
		d__1 = abs(objpp);
		do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
		e_wsfi();
		snprnt_(&c__31, str, &iw[1], leniw, (ftnlen)80);
		snprnt_(&c__22, str, &iw[1], leniw, (ftnlen)80);
	    }
	}
/* Proximal Point method 1 */
    }
/* nnH > 0 */
L800:
    return 0;
} /* s8getfeaslc_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8getFeasLC */
/* Subroutine */ int s8getweights_(integer *job, logical *boosted, integer *
	itn, doublereal *g0norm, doublereal *gnorm, doublereal *wtinf0, 
	doublereal *wtinf, doublereal *wtmax, doublereal *weight, doublereal *
	wtfactor, doublereal *wtscale, integer *iw, integer *leniw)
{
    /* Format strings */
    static char fmt_1000[] = "(\002 Itn\002,i7,\002: Elastic weight increase"
	    "d to \002,1p,e11.3)";

    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static doublereal newweight;
    static char str[80];
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___72 = { 0, str, 0, fmt_1000, 80, 1 };


/*     ================================================================== */
/*     s8getWeights  initializes or updates the elastic weight  wtInf. */
/*     The elastic weight is given by  wtInf = wtScale*weight, */
/*     where wtScale is some scale-dependent quantity (fObj here). */
/*     wtInf is increased by redefining weight as weight*wtFactor, where */
/*     wtFactor is a constant factor. */

/*     weight, wtFactor and wtScale are 'saved' local variables. */

/*     20 Feb 1997: First version of s8getWeights. */
/*     27 Apr 2001: wtMax introduced as parameter instead of local. */
/*     03 Aug 2003: snPRNT adopted. */
/*     03 Aug 2003: Current version of s8getWeights. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --iw;

    /* Function Body */
    if (*job == 0) {
/*        Set the weight. */
/*        weight is the ``unscaled'' weight on the infeasibilities. */
/*        wtScale is a scale factor based on the current gradient. */
/* Computing MAX */
	d__1 = 1., d__2 = *g0norm + *gnorm;
	*wtscale = max(d__1,d__2);
/*        wtScale  = max( 1.0d+2, g0Norm + gNorm) */
/*        wtScale  = max( 1.0d+0, g0Norm + gNorm) */
	*wtfactor = 10.;
	*weight = *wtinf0;
	*wtinf = *wtscale * *weight;
    } else if (*job == 1) {
/*        If possible, boost the weight. */
/* Computing MIN */
	d__1 = *wtfactor * *weight;
	newweight = min(d__1,*wtmax);
	*boosted = newweight > *weight;
	if (*boosted) {
	    *weight = newweight;
	    *wtinf = *weight * *wtscale;
	    *wtfactor *= 10.;
	    s_wsfi(&io___72);
	    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*wtinf), (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
	}
    }
    return 0;
} /* s8getweights_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8getWeights */
/* Subroutine */ int s8gloc_(integer *nncon, integer *nnjac, integer *nej, 
	integer *nlocj, integer *locj, integer *indj, integer *negcon, 
	integer *nlocg, integer *locg)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer j, k, ir, neg;

/*     ================================================================== */
/*     s8Gloc  counts the number of nonlinear Jacobian elements and */
/*     assembles their column pointers in locG. */

/*     29 Oct 2000: First version of s8Gloc. */
/*     31 Aug 2008: Local variable used for negCon. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --indj;
    --locj;
    --locg;

    /* Function Body */
    neg = 0;
    locg[1] = 1;
    i__1 = *nnjac;
    for (j = 1; j <= i__1; ++j) {
	i__2 = locj[j + 1] - 1;
	for (k = locj[j]; k <= i__2; ++k) {
	    ir = indj[k];
	    if (ir > *nncon) {
		goto L100;
	    }
	    ++neg;
	}
L100:
	locg[j + 1] = neg + 1;
    }
    return 0;
} /* s8gloc_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8Gloc */
/* Subroutine */ int s8gprod_(integer *task, doublereal *tolz, integer *nej, 
	integer *nlocj, integer *locj, integer *indj, integer *negcon, 
	integer *nlocg, integer *locg, doublereal *gcon, doublereal *alpha, 
	doublereal *x, integer *lenx, doublereal *beta, doublereal *y, 
	integer *leny)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, ig, ij, ir;
    static doublereal xj, sum, alphxj;

/*     ================================================================== */
/*     s8Gprod computes matrix-vector products involving J and x.  The */
/*     variable task specifies the operation to be performed as follows: */
/*       task = 'N' (normal)          y := alpha*J *x + beta*y, */
/*       task = 'T' (transpose)       y := alpha*J'*x + beta*y, */
/*     where alpha and beta are scalars, x and y are vectors and J is a */
/*     sparse matrix whose columns are in natural order. */

/*     26 Oct 2000: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --indj;
    --locj;
    --gcon;
    --locg;
    --x;
    --y;

    /* Function Body */
    if (*alpha == 0. && *beta == 1.) {
	return 0;
    }
/*     First form  y := beta*y. */
    if (*beta != 1.) {
	if (*beta == 0.) {
	    i__1 = *leny;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		y[i__] = 0.;
	    }
	} else {
	    i__1 = *leny;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		y[i__] = *beta * y[i__];
	    }
	}
    }
    if (*alpha == 0.) {
/*        Relax */
    } else if (*alpha == -1.) {
	if (*task == 0) {
	    i__1 = *lenx;
	    for (j = 1; j <= i__1; ++j) {
		xj = x[j];
		if (abs(xj) > *tolz) {
		    ig = locg[j];
		    i__2 = locj[j + 1] - 1;
		    for (ij = locj[j]; ij <= i__2; ++ij) {
			ir = indj[ij];
			if (ir > *leny) {
			    goto L100;
			}
			y[ir] -= gcon[ig] * xj;
			++ig;
		    }
		}
L100:
		;
	    }
	} else if (*task == 1) {
	    i__1 = *leny;
	    for (j = 1; j <= i__1; ++j) {
		sum = y[j];
		ig = locg[j];
		i__2 = locj[j + 1] - 1;
		for (ij = locj[j]; ij <= i__2; ++ij) {
		    ir = indj[ij];
		    if (ir > *lenx) {
			goto L200;
		    }
		    sum -= gcon[ig] * x[ir];
		    ++ig;
		}
L200:
		y[j] = sum;
	    }
	}
    } else {
/* General alpha */
	if (*task == 0) {
	    i__1 = *lenx;
	    for (j = 1; j <= i__1; ++j) {
		alphxj = *alpha * x[j];
		if (abs(alphxj) > *tolz) {
		    ig = locg[j];
		    i__2 = locj[j + 1] - 1;
		    for (ij = locj[j]; ij <= i__2; ++ij) {
			ir = indj[ij];
			if (ir > *leny) {
			    goto L300;
			}
			y[ir] += gcon[ig] * alphxj;
			++ig;
		    }
		}
L300:
		;
	    }
	} else if (*task == 1) {
	    i__1 = *leny;
	    for (j = 1; j <= i__1; ++j) {
		sum = 0.;
		ig = locg[j];
		i__2 = locj[j + 1] - 1;
		for (ij = locj[j]; ij <= i__2; ++ij) {
		    ir = indj[ij];
		    if (ir > *lenx) {
			goto L400;
		    }
		    sum += gcon[ig] * x[ir];
		    ++ig;
		}
L400:
		y[j] += *alpha * sum;
	    }
	}
/* task .eq. Normal */
    }
/* general alpha */
    return 0;
} /* s8gprod_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8Gprod */
/* Subroutine */ int s8gsize_(integer *m, integer *nncon, integer *nnjac, 
	integer *nej, integer *nlocj, integer *locj, integer *indj, integer *
	negcon)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer k, ir, neg, last, nlocg;

/* ================================================================= */
/* s8Gsize  counts the number of nonlinear Jacobian elements. */

/* 04 Nov 2000: First version of s8Gsize */
/* 31 Aug 2008: Local variable used for negCon. */
/* ================================================================= */
/* ----------------------------------------------------------------- */
    /* Parameter adjustments */
    --indj;
    --locj;

    /* Function Body */
    neg = 0;
    nlocg = *nnjac + 1;
    if (*nncon > 0) {
	last = locj[nlocg] - 1;
	if (*nncon == *m) {
	    neg = last;
	} else {
	    i__1 = last;
	    for (k = 1; k <= i__1; ++k) {
		ir = indj[k];
		if (ir <= *nncon) {
		    ++neg;
		}
	    }
	}
    }
    *negcon = max(1,neg);
    return 0;
} /* s8gsize_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8Gsize */
/* Subroutine */ int s8infs_(logical *elastic, integer *n, integer *nb, 
	integer *nncon0, integer *nncon, doublereal *tolx, doublereal *wtinf, 
	doublereal *prinf, doublereal *dualinf, integer *jprinf, integer *
	jdualinf, doublereal *bl, doublereal *bu, doublereal *fx, doublereal *
	rc, doublereal *x)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j;
    static doublereal v, w, dj, xj, tol, viol, slack;

/*     ================================================================== */
/*     s8Infs computes the maximum primal and dual infeasibilities, */
/*     using bl, bu, rc, x and the true nonlinear slacks Fxslk. */
/*     The linear constraints and bounds are assumed to be satisfied. */
/*     The primal infeasibility is therefore the maximum violation of */
/*     the nonlinear constraints. */
/*     The dual infeasibility is the maximum complementarity gap */
/*     for the bound constraints (with bounds assumed to be no further */
/*     than 1.0 from each x(j)). */

/*     prInf, dualInf   return the max primal and dual infeas. */

/*     20 Feb 1994: First version based on Minos 5.5 routine m8infs. */
/*     25 Oct 1996: Elastic mode added. */
/*     11 Sep 2014: Nonlinear Constraint violations replace slack values. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --x;
    --rc;
    --bu;
    --bl;
    --fx;

    /* Function Body */
    *jprinf = 0;
    *prinf = 0.;
    tol = *tolx;
/*     See how much  Fx  violates the bounds on the nonlinear slacks. */
/*     prInf is the maximum violation. */
    i__1 = *nncon;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = *n + i__;
	slack = fx[i__];
/* Computing MAX */
	d__1 = 0., d__2 = bl[j] - slack, d__1 = max(d__1,d__2), d__2 = slack 
		- bu[j];
	viol = max(d__1,d__2);
	if (*prinf < viol) {
	    *prinf = viol;
	    *jprinf = j;
	}
    }
/*     ------------------------------------------------------------------ */
/*     + rc(j)  is the multiplier for lower bound constraints. */
/*     - rc(j)  is the multiplier for upper bound constraints. */
/*     dualInf is the maximum component-wise complementarity measure. */
/*     ------------------------------------------------------------------ */
    *jdualinf = 0;
    *dualinf = 0.;
    i__1 = *nb;
    for (j = 1; j <= i__1; ++j) {
	dj = rc[j];
	if (dj != 0.) {
	    if (j <= *n || j > *n + *nncon) {
		xj = x[j];
	    } else {
		xj = fx[j - *n];
	    }
	    if (dj > 0.) {
/* Computing MIN */
		d__1 = xj - bl[j];
		dj *= min(d__1,1.);
	    } else if (dj < 0.) {
/* Computing MIN */
		d__1 = bu[j] - xj;
		dj = -dj * min(d__1,1.);
	    }
	    if (*dualinf < dj) {
		*dualinf = dj;
		*jdualinf = j;
	    }
	}
/* dj nonzero */
    }
/*     ------------------------------------------------------------------ */
/*     Include contributions from the elastic variables. */
/*     ------------------------------------------------------------------ */
    if (*elastic) {
	i__1 = *n + *nncon;
	for (j = *n + 1; j <= i__1; ++j) {
	    dj = rc[j];
	    v = bl[j] - x[j];
	    w = x[j] - bu[j];
	    if (v > tol) {
		dj = (d__1 = *wtinf - dj, abs(d__1)) * min(v,1.);
	    } else if (w > tol) {
		dj = (d__1 = *wtinf + dj, abs(d__1)) * min(w,1.);
	    } else {
		dj = 0.;
	    }
	    if (*dualinf < dj) {
		*dualinf = dj;
		*jdualinf = j;
	    }
	}
    }
    return 0;
} /* s8infs_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8Infs */
/* Subroutine */ int s8inith_(integer *nnh0, integer *nnh, doublereal *gnorm0,
	 doublereal *gnorm, doublereal *ud0, doublereal *hd, integer *iw, 
	integer *leniw, doublereal *rw, integer *lenrw)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    extern /* Subroutine */ int s8reseth_(integer *, doublereal *, doublereal 
	    *, integer *, integer *, doublereal *, integer *);

/*     ================================================================== */
/*     s8InitH  defines the initial approximate Hessian H. */

/*     On entry */
/*     -------- */
/*     gNorm0  is the term from the scaled constant objective row. */
/*     gNorm   approximates the two-norm of the nonlinear obj. gradient. */

/*     On exit */
/*     -------- */
/*     UD0     is a positive scalar such that H = UD0*UD0*I */
/*     HD      is the diagonal marix such that diag(HD) = H. */

/*     19 Aug 2013: First version based on dnInit. */
/*     20 Aug 2013: Last update. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
/* $$$      if (nnH .eq. 0) then */
/* $$$         UD0  = one */
/* $$$      else */
/* $$$         rnnH = nnH */
/* $$$         if (gNorm .gt. zero) then */
/* $$$            rnnH = nnH */
/* $$$            UD0  = sqrt(gNorm/sqrt(rnnH)) */
/* $$$         else */
/* $$$            UD0  = one */
/* $$$         end if */
/* $$$      end if */
/* Approximate Hessian type */
    /* Parameter adjustments */
    --hd;
    --iw;
    --rw;

    /* Function Body */
    if (*nnh == 0) {
	*ud0 = 1.;
/* UD0 is not used */
    } else {
/* Computing MAX */
	d__1 = 1., d__2 = *gnorm0 + *gnorm;
	*ud0 = sqrt(max(d__1,d__2) * 1.);
    }
/* Computing MIN */
    d__1 = max(*ud0,.01);
    *ud0 = min(d__1,10.);
    if (*nnh > 0) {
	iw[243] = -1;
	s8reseth_(nnh, ud0, &hd[1], &iw[1], leniw, &rw[1], lenrw);
    }
    return 0;
} /* s8inith_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8InitH */
/* Subroutine */ int s8initpen_(integer *nncon, doublereal *penparm, 
	doublereal *xpen0, doublereal *xpen, doublereal *rw, integer *lenrw)
{
    static doublereal eps;
    extern /* Subroutine */ int dload_(integer *, doublereal *, doublereal *, 
	    integer *);

/*     ================================================================== */
/*     s8InitPen  defines initial values for the penalty parameters. */

/*     19 Aug 2013: First version based on dnInit. */
/*     20 Aug 2013: Last update. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --xpen;
    --penparm;
    --rw;

    /* Function Body */
    eps = rw[1];
/*     --------------------------------------------- */
/*     Initialize the penalty parameters. */
/*     --------------------------------------------- */
/* machine precision.  IEEE DP  2.22e-16 */
    penparm[4] = 1.;
    penparm[1] = 1.;
    penparm[2] = 1. / eps;
    penparm[3] = *xpen0;
    dload_(nncon, xpen0, &xpen[1], &c__1);
    return 0;
} /* s8initpen_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8InitPen */
/* Subroutine */ int s8merit_(integer *nncon, doublereal *fmerit, doublereal *
	gmerit, doublereal *hmerit, doublereal *penparm, doublereal *fv, 
	doublereal *xpen, doublereal *w, doublereal *rw, integer *lenrw)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static logical overflow;
    static integer i__;
    static doublereal eps0;
    extern doublereal ddiv_(doublereal *, doublereal *, logical *), ddot_(
	    integer *, doublereal *, integer *, doublereal *, integer *), 
	    dnrm2_(integer *, doublereal *, integer *);
    static doublereal xpen0;
    extern /* Subroutine */ int ddscl_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dcopy_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static doublereal ppscl, xpeni;
    static logical boost;
    static doublereal wnorm, penold, penmax, penmin, pennew, incrun, rtundf, 
	    penlty, pendamp, pennorm;

/*     ================================================================== */
/*     s8Merit  computes the contributions to the merit function and its */
/*     directional derivative from the nonlinear constraints. */
/*     The penalty parameters  xPen(j)  are increased if */
/*     the directional derivative is not sufficiently negative. */

/*     On entry: */
/*         Fv     is the violation c(x) + A(linear)x - s,  where */
/*                s  minimizes the merit function with respect to the */
/*                nonlinear slacks only. */

/*     30 Dec 1991: First version based on Npsol 4.0 routine npmrt. */
/*     02 Nov 1996: Multipliers no longer updated here. */
/*     19 Jul 1997: Thread-safe version. */
/*     21 Oct 2000: Made compatible with SNOPT 6.1 */
/*     21 Oct 2000: Current version of s8Merit. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --w;
    --xpen;
    --fv;
    --penparm;
    --rw;

    /* Function Body */
    eps0 = rw[2];
    rtundf = rw[10];
    xpen0 = rw[89];
    pendamp = penparm[1];
    penmax = penparm[2];
    pennorm = penparm[3];
    incrun = penparm[4];
    overflow = FALSE_;
/*     Find the quantities that define  penMin, the vector of minimum */
/*     two-norm such that the directional derivative is one half of */
/*     approximate curvature   - (p)'H(p). */
/*     The factor  rtUndf  tends to keep  xPen  sparse. */
    i__1 = *nncon;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if ((d__1 = fv[i__], abs(d__1)) <= rtundf) {
	    w[i__] = 0.;
	} else {
/* Computing 2nd power */
	    d__1 = fv[i__];
	    w[i__] = d__1 * d__1;
	}
    }
    wnorm = dnrm2_(nncon, &w[1], &c__1);
    d__1 = *gmerit + *hmerit * .5;
    ppscl = ddiv_(&d__1, &wnorm, &overflow);
    if (abs(ppscl) <= penmax && ! overflow) {
/*        --------------------------------------------------------------- */
/*        Bounded  penMin  found.  The final value of  xPen(i)  will */
/*        never be less than  penMin(i).  A trial value  penNew  is */
/*        computed that is equal to the geometric mean of the previous */
/*        xPen  and a damped value of penMin.  The new  xPen  is defined */
/*        as  penNew  if it is less than half the previous  xPen  and */
/*        greater than  penMin. */
/*        --------------------------------------------------------------- */
	i__1 = *nncon;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	    d__1 = w[i__] / wnorm * ppscl;
	    penmin = max(d__1,0.);
	    xpeni = xpen[i__];
	    pennew = sqrt(xpeni * (pendamp + penmin));
	    if (pennew < xpeni * .5) {
		xpeni = pennew;
	    }
	    xpeni = max(xpeni,penmin);
	    xpen[i__] = max(xpeni,xpen0);
	}
	penold = pennorm;
	pennorm = dnrm2_(nncon, &xpen[1], &c__1);
/*        --------------------------------------------------------------- */
/*        If  IncRun = true,  there has been a run of iterations in */
/*        which the norm of  xPen  has not decreased.  Conversely, */
/*        IncRun = false  implies that there has been a run of */
/*        iterations in which the norm of xPen has not increased.  If */
/*        IncRun changes during this iteration the damping parameter */
/*        penDamp is increased by a factor of two.  This ensures that */
/*        xPen(j) will oscillate only a finite number of times. */
/*        --------------------------------------------------------------- */
	boost = FALSE_;
	if (incrun > 0. && pennorm < penold) {
	    boost = TRUE_;
	}
	if (incrun < 0. && pennorm > penold) {
	    boost = TRUE_;
	}
	if (boost) {
/* Computing MIN */
	    d__1 = 1 / eps0, d__2 = pendamp * 2.;
	    pendamp = min(d__1,d__2);
	    incrun = -incrun;
	}
    }
/*     ------------------------------------------------------------------ */
/*     Compute the new value and directional derivative of the */
/*     merit function. */
/*     ------------------------------------------------------------------ */
    dcopy_(nncon, &fv[1], &c__1, &w[1], &c__1);
    ddscl_(nncon, &xpen[1], &c__1, &w[1], &c__1);
    penlty = ddot_(nncon, &w[1], &c__1, &fv[1], &c__1);
    *fmerit += penlty * .5;
    *gmerit -= penlty;
    penparm[1] = pendamp;
    penparm[2] = penmax;
    penparm[3] = pennorm;
    penparm[4] = incrun;
    return 0;
} /* s8merit_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine  s8Merit */
/* Subroutine */ int s8optimizeslacks_(logical *elastic, integer *n, integer *
	nb, integer *nncon, integer *ninfe, doublereal *sinfe, doublereal *
	featol, doublereal *wtinf, doublereal *bl, doublereal *bu, doublereal 
	*fv, doublereal *x, doublereal *ycon, doublereal *xpen, doublereal *
	fx)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static integer i__, j;
    static doublereal dm, vl, xj, vu, blj, con, buj, yconi, xpeni;

/*     ================================================================== */
/*     s8optimizeSlacks computes the nonlinear constraint violations: */
/*        Fv = fCon + A(linear)*x - (optimal nonlinear slacks) */
/*           =        Fx          - (optimal nonlinear slacks) */

/*     The optimal nonlinear slacks are computed as follows: */
/*     (1) Feasible  nonlinear slacks are adjusted so that they minimize */
/*         the merit function subject to  x  and  yCon  being held */
/*         constant. */
/*     (2) Infeasible slacks are compared with the true nonlinear slacks, */
/*         and, if necessary, they are adjusted so that the sum of */
/*         infeasibilities is reduced. */

/*     If yCon is zero, the violation can be set to any value without */
/*     changing the merit function.  In this case we choose the slack to */
/*     so that  the violation is zero (subject to the constraints above). */

/*     On entry, */
/*        x     =  the current x. */
/*        Fx    =  fCon + A(linear)*x,   defined in s8Fx. */

/*     On exit, */
/*        x     =  contains the optimal slacks in x(n+1:n+nnCon). */
/*                 Other elements are left unchanged. */
/*        Fv    =  fCon + A(linear)*x - optimal slacks. */
/*        Fx    =  unaltered. */
/*        nInfE =  number of optimal slack infeasibilities. */
/*        sInfE =  sum of optimal slack infeasibilities. */

/*     19 Nov 2012: First version based on SNOPT routine s8sOpt. */
/*     21 Sep 2013: Added the sum and number of slack infeasibilities. */
/*     01 Jan 2014: Reworked to handle elastic mode better. */
/*     20 Jan 2014: Reworked and simplified to compute the correct sInfE. */
/*     20 Jan 2014: Last update. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --x;
    --bu;
    --bl;
    --fx;
    --xpen;
    --ycon;
    --fv;

    /* Function Body */
    *ninfe = 0;
    *sinfe = 0.;
    i__1 = *nncon;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = *n + i__;
	con = fx[i__];
	xj = x[j];
	xpeni = xpen[i__];
	yconi = ycon[i__];
	blj = bl[j];
	buj = bu[j];
/* -------------------------------------------------------------- */
/* Redefine  xj  so that it minimizes the merit function */
/* subject to upper and lower bounds determined by the current */
/* multipliers. */
/* -------------------------------------------------------------- */
	if (*elastic && xj < blj - *featol) {
/* ----------------------------------------------------------- */
/* This slack is below its lower bound. */
/* ----------------------------------------------------------- */
	    yconi -= *wtinf;
	    dm = yconi - xpeni * (con - xj);
	    if (dm < 0.) {
/* xj must increase to reduce the merit function. */

/* s = max( s0, min(c,c - Mul/Pen) ) */

/* If Mul <= 0             then c < c - Mul/Pen */
/* If Mul >  0 and  dM < 0 then Pen is nonzero. */
		if (yconi <= 0.) {
		    xj = max(xj,con);
		} else {
/* Computing MAX */
		    d__1 = xj, d__2 = con - yconi / xpeni;
		    xj = max(d__1,d__2);
		}
	    } else if (dm > 0.) {
/* xj must decrease to reduce the merit function. */

/* s = min( s0, max(c,c - Mul/Pen) ) */

/* If Mul >= 0             then c > c - Mul/Pen */
/* If Mul <  0 and  dM > 0 then Pen is nonzero. */
		if (yconi >= 0.) {
		    xj = min(xj,con);
		} else {
/* Computing MIN */
		    d__1 = xj, d__2 = con - yconi / xpeni;
		    xj = min(d__1,d__2);
		}
	    }
	} else if (*elastic && xj > buj + *featol) {
/* ----------------------------------------------------------- */
/* This slack is at or above its upper bound in elastic mode. */
/* ----------------------------------------------------------- */
	    yconi += *wtinf;
	    dm = yconi - xpeni * (con - xj);
	    if (dm < 0.) {
/* xj must increase to reduce the merit function. */

/* s = max( s0, min(c,c - Mul/Pen) ) */

/* If Mul <= 0             then c < c - Mul/Pen */
/* If Mul >  0 and  dM < 0 then Pen is nonzero. */
		if (yconi <= 0.) {
		    xj = max(xj,con);
		} else {
/* Computing MAX */
		    d__1 = xj, d__2 = con - yconi / xpeni;
		    xj = max(d__1,d__2);
		}
	    } else if (dm > 0.) {
/* xj must decrease to reduce the merit function. */

/* s = min( s0, max(c,c - Mul/Pen) ) */

/* If Mul >= 0             then c > c - Mul/Pen */
/* If Mul <  0 and  dM > 0 then Pen is nonzero. */
		if (yconi >= 0.) {
		    xj = min(xj,con);
		} else {
/* Computing MIN */
		    d__1 = xj, d__2 = con - yconi / xpeni;
		    xj = min(d__1,d__2);
		}
	    }
	} else {
/* ----------------------------------------------------------- */
/* Feasible slack. */
/* ----------------------------------------------------------- */
/* Define  dM, the derivative of the merit function. */
/* Require      max( bl, tbl ) <=  xj <= min( bu,tbu ). */
	    d__1 = con - xj;
	    dm = yconi - pow_dd(&xpeni, &d__1);
	    if (dm < 0.) {
/* xj must increase to reduce the merit function. */

/* s = max( s0, min(c,c - Mul/Pen, bu) ) */

/* If Mul <= 0             then c < c - Mul/Pen */
/* If Mul >  0 and  dM < 0 then Pen is nonzero. */
		if (yconi <= 0.) {
/* Computing MAX */
		    d__1 = xj, d__2 = min(buj,con);
		    xj = max(d__1,d__2);
		} else {
/* Computing MAX */
/* Computing MIN */
		    d__3 = buj, d__4 = con - yconi / xpeni;
		    d__1 = xj, d__2 = min(d__3,d__4);
		    xj = max(d__1,d__2);
		}
	    } else if (dm > 0.) {
/* xj must decrease to reduce the merit function. */

/* s = min( s0, max(c,c - Mul/Pen, bl) ) */

/* If Mul >= 0             then c > c - Mul/Pen */
/* If Mul <  0 and  dM > 0 then Pen is nonzero. */
		if (yconi >= 0.) {
/* Computing MIN */
		    d__1 = xj, d__2 = max(blj,con);
		    xj = min(d__1,d__2);
		} else {
/* Computing MIN */
/* Computing MAX */
		    d__3 = blj, d__4 = con - yconi / xpeni;
		    d__1 = xj, d__2 = max(d__3,d__4);
		    xj = min(d__1,d__2);
		}
	    }
	}
	vl = blj - xj;
	vu = xj - buj;
	if (vl > *featol || vu > *featol) {
	    ++(*ninfe);
	    *sinfe += max(vl,vu);
	}
	fv[i__] = con - xj;
	x[j] = xj;
    }
    return 0;
} /* s8optimizeslacks_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8optimizeSlacks */
/* Subroutine */ int s8hxlp_(integer *nnh, doublereal *x, doublereal *hx, 
	integer *status, char *cu, integer *lencu, integer *iu, integer *
	leniu, doublereal *ru, integer *lenru, ftnlen cu_len)
{
/*     ================================================================== */
/*     s8HxLP is the argument qpHx for s5solve when s5solve is called */
/*     from one of the snOpt wrappers. */

/*     04 Dec 2004: First version of s8HxLP. */
/*     04 Dec 2004: Current version of s8HxLP. */
/*     ================================================================== */
/* Relax */
    /* Parameter adjustments */
    --hx;
    --x;
    cu -= 8;
    --iu;
    --ru;

    /* Function Body */
    return 0;
} /* s8hxlp_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8HxLP */
/* Subroutine */ int s8hxnp_(S_fp userhv, integer *nnh, integer *neh, integer 
	*nloch, integer *loch, integer *indh, doublereal *hcol, doublereal *v,
	 doublereal *hv, integer *status, char *cu, integer *lencu, integer *
	iu, integer *leniu, doublereal *ru, integer *lenru, char *cw, integer 
	*lencw, integer *iw, integer *leniw, doublereal *rw, integer *lenrw, 
	ftnlen cu_len, ftnlen cw_len)
{
    static integer lxscaled, lvlscale, lx;
    extern /* Subroutine */ int s8callstatus_(integer *, integer *, integer *)
	    , ddscl_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer nncon;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer lycon;
    static logical scaled;
    static integer hvmode, lscales;

/* ================================================================= */
/* s8HxNP (aka Hprod)  is called by the QP solver. */
/* It is a wrapper that calls userHv  (aka Hprod1) to compute Hv, */
/* the product of the Hessian and v and scales it. */

/* 30 Apr 2011: First version of s8npHx */
/* 07 May 2011: HvMode added to the argument list of userHv. */
/* ================================================================= */
/* ----------------------------------------------------------------- */
    /* Parameter adjustments */
    --hv;
    --v;
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
    nncon = iw[23];
/* # of nonlinear constraints */
    lvlscale = iw[75];
/* scale option */
    lx = iw[300];
/* x(nb) at which H is defined */
    lscales = iw[296];
/* scales(nb)  = row and column scales */
    lxscaled = iw[302];
/* xScaled(n)  = copy of the scaled x */
    lycon = iw[348];
/* yCon (nnCon) = multipliers for fCon */
    scaled = lvlscale > 0;
/* Determine the call-status for the NP subproblem. */
    s8callstatus_(status, &iw[1], leniw);
    if (scaled) {
	dcopy_(nnh, &v[1], &c__1, &rw[lxscaled], &c__1);
	ddscl_(nnh, &rw[lscales], &c__1, &v[1], &c__1);
/* Scale the base x and the multipliers */
    }
    hvmode = 1;
/* Do NOT compute a new Hessian. */
    (*userhv)(&hvmode, &nncon, nnh, neh, nloch, &loch[1], &indh[1], &hcol[1], 
	    &rw[lycon], &rw[lx], &v[1], &hv[1], status, cu + 8, lencu, &iu[1],
	     leniu, &ru[1], lenru, (ftnlen)8);
    if (scaled) {
	dcopy_(nnh, &rw[lxscaled], &c__1, &v[1], &c__1);
	ddscl_(nnh, &rw[lscales], &c__1, &hv[1], &c__1);
/* Unscale the base x and the multipliers! */
    }
    return 0;
} /* s8hxnp_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8HxNP */
/* Subroutine */ int s8hxpp_(U_fp hprod, integer *nnh, integer *neh, integer *
	nloch, integer *loch, integer *indh, doublereal *hcol, doublereal *x, 
	doublereal *hx, integer *status, char *cu, integer *lencu, integer *
	iu, integer *leniu, doublereal *ru, integer *lenru, char *cw, integer 
	*lencw, integer *iw, integer *leniw, doublereal *rw, integer *lenrw, 
	ftnlen cu_len, ftnlen cw_len)
{
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);

/*     ================================================================== */
/*     s8HxPP  defines the product  H*x  for the proximal-point QP */
/*     subproblem of snopt. */

/*     On exit,    Hx   = x. */

/*     23 Oct 1993: First version of s8HxPP. */
/*     02 Aug 2000: Current version. */
/*     ================================================================== */
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
    dcopy_(nnh, &x[1], &c__1, &hx[1], &c__1);
    return 0;
} /* s8hxpp_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8HxPP */
/* Subroutine */ int s8hxqp_(S_fp usrhx, integer *nnh, integer *neh, integer *
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
/*     s8HxQP  computes the user-defined product  Hx  and scales it. */

/*     07 Nov 2014: First   version of s8Hxqp. */
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
} /* s8hxqp_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8HxQP */
/* Subroutine */ int s8hxnull_(integer *nnh, doublereal *x, doublereal *hx, 
	integer *status, char *cu, integer *lencu, integer *iu, integer *
	leniu, doublereal *ru, integer *lenru, ftnlen cu_len)
{
    /* Format strings */
    static char fmt_1000[] = "(//\002 XXX  The default (dummy) version of su"
	    "broutine Hx\002,\002     has been called from SNOPT. \002)";

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void);

    /* Local variables */
    static integer nout;

    /* Fortran I/O blocks */
    static cilist io___142 = { 0, 0, 0, fmt_1000, 0 };


/*     ================================================================== */
/*     This is the dummy (empty) version of the routine qpHx. */
/*     It should never be called. */

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
	    io___142.ciunit = nout;
	    s_wsfe(&io___142);
	    e_wsfe();
	}
    }
    return 0;
} /* s8hxnull_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8HxNull */
/* Subroutine */ int s8rand_(integer *leng, integer *neg, doublereal *g)
{
    static integer seeds[3];
    extern /* Subroutine */ int ddrand_(integer *, doublereal *, integer *, 
	    integer *);

/*     ================================================================== */
/*     s8rand  fills the array g with random numbers. */

/*     15 Nov 1991: First version of s8rand in s8aux. */
/*     30 Jun 1999: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --g;

    /* Function Body */
    if (*neg <= 0) {
	return 0;
    }
    seeds[0] = 1547;
    seeds[1] = 2671;
    seeds[2] = 3770;
    ddrand_(neg, &g[1], &c__1, seeds);
    return 0;
} /* s8rand_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8rand */
/* Subroutine */ int s8rc_(doublereal *scaleobj, integer *minimize, integer *
	iobj, integer *m, integer *n, integer *nb, integer *nnobj0, integer *
	nnobj, integer *nncon, integer *nnjac, integer *negcon, integer *nej, 
	integer *nlocj, integer *locj, integer *indj, doublereal *jcol, 
	doublereal *gobj, doublereal *gcon, doublereal *pi, doublereal *rc)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer j, k, l;
    static doublereal dj;
    static integer ir;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    static doublereal signobj;

/*     ================================================================== */
/*     s8rc   computes reduced costs rc = gObj - ( A  -I )'*pi, */
/*     using  gCon  as the top left-hand corner of A. */
/*     gCon, gObj and pi are assumed to exist. */

/*     s8rc   is called by s8SQP. */

/*     28 Sep 1993: First version, derived from m4rc. */
/*     31 Oct 1996: Min sum option added. */
/*     30 Oct 2000: Current version of s8rc. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --pi;
    --rc;
    --gobj;
    --gcon;
    --jcol;
    --indj;
    --locj;

    /* Function Body */
    l = 0;
    i__1 = *nnjac;
    for (j = 1; j <= i__1; ++j) {
	dj = 0.;
	i__2 = locj[j + 1] - 1;
	for (k = locj[j]; k <= i__2; ++k) {
	    ir = indj[k];
	    if (ir <= *nncon) {
		++l;
		dj += pi[ir] * gcon[l];
	    } else {
		dj += pi[ir] * jcol[k];
	    }
	}
	rc[j] = -dj;
    }
    i__1 = *n;
    for (j = *nnjac + 1; j <= i__1; ++j) {
	dj = 0.;
	i__2 = locj[j + 1] - 1;
	for (k = locj[j]; k <= i__2; ++k) {
	    ir = indj[k];
	    dj += pi[ir] * jcol[k];
	}
	rc[j] = -dj;
    }
    dcopy_(m, &pi[1], &c__1, &rc[*n + 1], &c__1);
/*     Include the nonlinear objective gradient. */
    signobj = (doublereal) (*minimize);
    if (*nnobj > 0) {
	daxpy_(nnobj, &signobj, &gobj[1], &c__1, &rc[1], &c__1);
    }
    if (*iobj > 0) {
	rc[*n + *iobj] += signobj * *scaleobj;
    }
    return 0;
} /* s8rc_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8rc */
/* Subroutine */ int s8scaleg_(integer *nnobj, doublereal *scales, doublereal 
	*gobj, doublereal *rw, integer *lenrw)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j;
    static doublereal grad, gdummy;

/*     ================================================================== */
/*     s8scaleG  scales the objective gradient. */
/*     s8scaleG is called by funwrapper only if modefg = 2. */
/*     Hence, it is used to scale known gradient elements (if any), */
/*     but is not called when missing gradients are being estimated */
/*     by s6dobj. */

/*     17 Feb 1992: First version. */
/*     16 Jul 1997: Thread-safe version. */
/*     02 Jan 2001: Current version of s8scaleG. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --gobj;
    --scales;
    --rw;

    /* Function Body */
    gdummy = rw[69];
/* definition of 'unset' value */
    i__1 = *nnobj;
    for (j = 1; j <= i__1; ++j) {
	grad = gobj[j];
	if (grad != gdummy) {
	    gobj[j] = grad * scales[j];
	}
    }
    return 0;
} /* s8scaleg_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8scaleG */
/* Subroutine */ int s8scalej_(integer *nncon, integer *nnjac, integer *
	negcon, integer *n, doublereal *scales, integer *nej, integer *nlocj, 
	integer *locj, integer *indj, doublereal *gcon, doublereal *rw, 
	integer *lenrw)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer j, k, l, ir;
    static doublereal grad, cscale, gdummy;

/*     ================================================================== */
/*     s8scaleJ  scales the Jacobian. */
/*     s8scaleJ is called by funwrapper only if modefg = 2. */
/*     Hence, it is used to scale known gradient elements (if any), */
/*     but is not called when missing gradients are being estimated */
/*     by s6dcon. */

/*     17 Feb 1992: First version based on Minos routine m8sclj. */
/*     16 Jul 1997: Thread-safe version. */
/*     02 Dec 2001: Current version of s8scaleJ. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --gcon;
    --scales;
    --indj;
    --locj;
    --rw;

    /* Function Body */
    gdummy = rw[69];
/* definition of 'unset' value */
    l = 0;
    i__1 = *nnjac;
    for (j = 1; j <= i__1; ++j) {
	cscale = scales[j];
	i__2 = locj[j + 1] - 1;
	for (k = locj[j]; k <= i__2; ++k) {
	    ir = indj[k];
	    if (ir > *nncon) {
		goto L300;
	    }
	    ++l;
	    grad = gcon[l];
	    if (grad != gdummy) {
		gcon[l] = grad * cscale / scales[*n + ir];
	    }
	}
L300:
	;
    }
    return 0;
} /* s8scalej_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8scaleJ */
/* Subroutine */ int s8sinf_(integer *n, integer *nb, integer *nncon, 
	doublereal *tolx, integer *ninf, doublereal *sinf, doublereal *bl, 
	doublereal *bu, doublereal *x)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j;
    static doublereal tol, slack, violl, violu;

/*     ================================================================== */
/*     s8sInf computes the sum of infeasibilities of the nonlinear slacks */
/*     using bl, bu and x. */

/*     10 Jan 1997: First version of s8sInf. */
/*     30 Oct 2000: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --x;
    --bu;
    --bl;

    /* Function Body */
    *ninf = 0;
    *sinf = 0.;
    tol = *tolx;
/*     See how much  x(n+1:n+nnCon) violates its bounds. */
    i__1 = *nncon;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = *n + i__;
	slack = x[j];
	violl = bl[j] - slack;
	violu = slack - bu[j];
	if (violl > tol || violu > tol) {
	    ++(*ninf);
	    *sinf += max(violl,violu);
	}
    }
    return 0;
} /* s8sinf_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8sInf */
/* Subroutine */ int s8solveqp_(integer *iexit, U_fp mnrlog, U_fp hprod, U_fp 
	hprod1, integer *hvcalls, logical *elastic, logical *gotr, integer *
	itn, integer *itqp, integer *lenr, integer *m, integer *maxs, integer 
	*mbs, integer *n, integer *nb, integer *nncon0, integer *nncon, 
	integer *nnobj0, integer *nnobj, integer *nnh0, integer *nnh, integer 
	*ns, integer *ndegen, integer *mjrprtlvl, integer *mnrprtlvl, integer 
	*minimize, integer *iobj, doublereal *scaleobj, doublereal *objadd, 
	doublereal *objqp, doublereal *condzhz, doublereal *toloptfp, 
	doublereal *toloptqpk, doublereal *tolx, integer *ninf, doublereal *
	sinf, integer *elastics, integer *ninfe, doublereal *sinfe, 
	doublereal *wtinf, doublereal *ud0, doublereal *pinorm, integer *nej, 
	integer *nlocj, integer *locj, integer *indj, doublereal *jcol, 
	integer *neh, integer *nloch, integer *loch, integer *indh, 
	doublereal *hcol, integer *etype, integer *estate, integer *feastype, 
	integer *hs, integer *kbs, doublereal *bl, doublereal *bu, doublereal 
	*blqp, doublereal *buqp, doublereal *blbs, doublereal *bubs, 
	doublereal *gbs, doublereal *gqp, doublereal *gobj, doublereal *hd, 
	doublereal *hdx, doublereal *pbs, doublereal *pi, doublereal *r__, 
	doublereal *rc, doublereal *rg, doublereal *rg2, doublereal *qprhs, 
	doublereal *scales, doublereal *x, doublereal *xbs, doublereal *xqp0, 
	doublereal *xqp, integer *iy, integer *iy1, doublereal *y, doublereal 
	*y1, doublereal *y2, char *cu, integer *lencu, integer *iu, integer *
	leniu, doublereal *ru, integer *lenru, char *cw, integer *lencw, 
	integer *iw, integer *leniw, doublereal *rw, integer *lenrw, ftnlen 
	cu_len, ftnlen cw_len)
{
    /* Format strings */
    static char fmt_1500[] = "(\002 Itn\002,i7,\002: Infeasible nonelastics "
	    "in QP feasibility\002,\002 phase\002)";
    static char fmt_1100[] = "(\002 Itn\002,i7,\002: Infeasible subproblem"
	    ".\002,\002 Elastic mode started with weight = \002,1p,e8.1)";
    static char fmt_1200[] = "(\002 Itn\002,i7,\002: Feasible QP non-elast"
	    "ics\002)";
    static char fmt_1300[] = "(\002 Itn\002,i7,\002: Feasible QP subproblem"
	    " \002)";
    static char fmt_1400[] = "(\002 Itn\002,i7,\002: Large multipliers.\002"
	    ",\002 Elastic mode started with weight = \002,1p,e8.1)";
    static char fmt_2001[] = "(\002 Itn\002,i7,\002: Infeasible nonelastics "
	    "in QP optimality\002,\002 phase\002)";
    static char fmt_2002[] = "(\002 Itn\002,i7,\002: Unbounded QP subproble"
	    "m\002)";
    static char fmt_3000[] = "(\002 Itn\002,i7,\002: Hessian reset\002)";
    static char fmt_2008[] = "(\002 Itn\002,i7,\002: Ill-conditioned QP null"
	    "-space basis.\002,\002 Cond = \002,1p,e8.1)";
    static char fmt_2006[] = "(\002 Itn\002,i7,\002: Indefinite QP reduced H"
	    "essian\002)";
    static char fmt_2007[] = "(\002 Itn\002,i7,\002: Large QP reduced gradie"
	    "nt\002)";
    static char fmt_2010[] = "(\002 Itn\002,i7,\002: Too many CG subspace it"
	    "erations.\002)";

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    extern /* Subroutine */ int s5zhzfac_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *), 
	    s8reseth_(integer *, doublereal *, doublereal *, integer *, 
	    integer *, doublereal *, integer *);
    static logical feasible;
    static doublereal hcondbnd;
    static integer lvlobjfp;
    static logical normalin;
    static integer qpsolver;
    static logical terminate;
    static integer lurequest, itqptarget;
    static char str[80];
    extern /* Subroutine */ int s5qn_(integer *, integer *, char *, integer *,
	     U_fp, U_fp, U_fp, integer *, logical *, logical *, logical *, 
	    integer *, logical *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, integer *
	    , doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, doublereal *, doublereal *, doublereal *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, ftnlen, 
	    ftnlen, ftnlen), s5qp_(integer *, integer *, char *, integer *, 
	    U_fp, U_fp, U_fp, integer *, integer *, logical *, logical *, 
	    logical *, integer *, logical *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, integer *, integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, doublereal *, integer *, integer 
	    *, integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, char *, integer *, integer *, integer *, doublereal 
	    *, integer *, char *, integer *, integer *, integer *, doublereal 
	    *, integer *, ftnlen, ftnlen, ftnlen);
    static integer suboptimize, eigh;
    static logical newb;
    static integer ngqp, maxr;
    static logical luok, nonlinearcon;
    static integer ngqp0;
    extern /* Subroutine */ int s5zhz_(integer *, U_fp, U_fp, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, char *, integer *, integer *, integer *, 
	    doublereal *, integer *, char *, integer *, integer *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen);
    static integer emode;
    static doublereal objfp;
    static logical needx;
    static integer maxsb;
    static doublereal flmax;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer nswap;
    static logical newlu;
    extern /* Subroutine */ int s2bfac_(integer *, integer *, logical *, 
	    logical *, logical *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, integer *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, integer *);
    static doublereal normz;
    static logical rcheck;
    static doublereal zcndbd;
    static logical needlu;
    static integer inform__;
    static logical solved;
    static integer itnlim, mminor, lvlpre;
    static doublereal plinfy, rgnorm, zhzmin, zhzmax;
    static integer typelu;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen), s2trylu_(integer *, integer *, integer *, integer *, 
	    logical *, integer *, integer *, integer *, doublereal *, integer 
	    *);
    static integer lvlobje;
    static char probtag[20];
    static doublereal targeth;
    static logical resolve;
    static doublereal targetz;
    static integer itqpmax, rankzhz;
    extern /* Subroutine */ int s5rcheck_(integer *, U_fp, U_fp, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, integer *, integer 
	    *, integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___199 = { 0, str, 0, fmt_1500, 80, 1 };
    static icilist io___200 = { 0, str, 0, fmt_1100, 80, 1 };
    static icilist io___201 = { 0, str, 0, fmt_1200, 80, 1 };
    static icilist io___202 = { 0, str, 0, fmt_1300, 80, 1 };
    static icilist io___215 = { 0, str, 0, fmt_1400, 80, 1 };
    static icilist io___216 = { 0, str, 0, fmt_2001, 80, 1 };
    static icilist io___217 = { 0, str, 0, fmt_2002, 80, 1 };
    static icilist io___218 = { 0, str, 0, fmt_3000, 80, 1 };
    static icilist io___219 = { 0, str, 0, fmt_2008, 80, 1 };
    static icilist io___220 = { 0, str, 0, fmt_2006, 80, 1 };
    static icilist io___221 = { 0, str, 0, fmt_2007, 80, 1 };
    static icilist io___222 = { 0, str, 0, fmt_2010, 80, 1 };
    static icilist io___223 = { 0, str, 0, fmt_3000, 80, 1 };


/*     ================================================================== */
/*     s8solveQP  computes  xQP, the solution of the QP subproblem. */
/*     By construction, the problem has  nnH  nonlinear variables, */

/*     The SQP base point  x  is not altered. */

/*     On entry: */
/*     --------- */
/*     eType    contains the allowed type of elastic violation. */
/*     eState   contains the elastic states. */
/*     The LU factorization is assumed to be known. */
/*     The arrays  xBS, blBS and buBS are defined. */

/*     iExit     Status */
/*     -----     ------ */
/*      >0         Fatal error */
/*       0         QP solution found */
/*      -1         Too many iterations */
/*      -2         Too many superbasics */

/*     LUrequest =  1  Frequency */
/*     LUrequest =  2  LU nonzeros increased */
/*     LUrequest =  3 */
/*     LUrequest =  4 */
/*     LUrequest =  5  Singular after LU mod */
/*     LUrequest =  6  Unstable LU mod (growth in new column of U) */
/*     LUrequest =  7  Not enough memory */
/*     LUrequest =  8 */
/*     LUrequest =  9 */
/*     LUrequest = 10  Row error in setx */
/*     LUrequest = 11  Big  dx   in setx */

/*     LUrequest = 20  Infeasible nonelastics in phase 2. */
/*     LUrequest = 21  Iterative refinement failed in QP */
/*     LUrequest = 22  Unbounded QP */
/*     LUrequest = 23 */
/*     LUrequest = 24  Small directional derivative in QP */
/*     LUrequest = 25  Ill-conditioned Z */
/*     LUrequest = 26  Indefinite Z'HZ in QP */
/*     LUrequest = 27  R singular after bound swap in QP */

/*     On output, */
/*     QPerr points to ' ', 't', 'u' or 'w'. */
/*     QPfea points to ' '  or 'i'. */

/*     30 Dec 1991: First version of s8solveQP (as s8iQP). */
/*     19 Jul 1997: Thread-safe version. */
/*     31 Jul 2003: dnormj used for norm of the nonlinear pis. */
/*     03 Aug 2003: snEXIT and snPRNT  adopted. */
/*     19 Jun 2008: Hprod, Hprod1 added as arguments. */
/*     07 Oct 2014: infoTags added. */
/*     18 Oct 2014: Switch to elastic mode based on piNorm. */
/*     03 Nov 2014: neH, indH, locH, Hcol added as arguments. */
/*     29 Dec 2014: Nonzero nnH argument for FP phase. */
/*     29 Dec 2014: Merged with s8iQN. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/* condition estimate of Z */
/* Current QP solver */
/* Current precon mode */
/* itns since last factorize */
/* # lines in log     file */
/* # lines in summary file */
/* >0 => Minor heading for iPrint */
/* >0 => Minor heading for iSumm */
/* TagInfo(1): QN update type */
/* TagInfo(2): */
/* TagInfo(3): Line search outcome */
/* TagInfo(4): QP Feasibility */
/* TagInfo(5)  QP Optimality */
/* TagInfo(6) */
/* TagInfo(7): Approx Hessian type */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --r__;
    --pi;
    --rg2;
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
    --xqp;
    --xqp0;
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
    --qprhs;
    --gobj;
    --hdx;
    --hd;
    --gqp;
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
    maxr = iw[52];
/* max columns of R */
    qpsolver = iw[55];
/* = 0:1:2   => QPChol:CG:QN QP solver */
    lvlpre = iw[77];
/* >0     => QN preconditioned CG */
    itnlim = iw[89];
/* limit on total iterations */
    mminor = iw[91];
/* limit on minor iterations */
    flmax = rw[8];
/* est. of the largest pos. real */
    hcondbnd = rw[85];
/* bound on the condition of ZHZ */
    zcndbd = rw[86];
/* bound on the condition of Z */
    *iexit = 0;
    *itqp = 0;
    s_copy(probtag, "QP subproblem", (ftnlen)20, (ftnlen)13);
    plinfy = flmax;
    targeth = hcondbnd;
/* condH > targetH => H is ill-conditioned */
    targetz = zcndbd;
/* condZ > targetZ => Z is ill-conditioned */
    *condzhz = 0.;
    normalin = ! (*elastic);
/* eigH    encodes the inertia of the Hessian H */
/*         eigH        Hessian */
/*         -----   --------------------- */
/*          -1                indefinite */
/*           0     positive semidefinite */
/*           1     positive semidefinite */
    if (*nnh < *n) {
	eigh = 0;
    } else {
	eigh = 1;
    }
    iw[220] = 0;
    iw[221] = 0;
    iw[223] = 1;
    iw[225] = 1;
    iw[241] = 0;
/* Optimal  QP subproblem */
    iw[240] = 0;
/*     If nS dips below maxR, go back to using the default solver. */
/* Feasible QP subproblem */
    if (*ns <= maxr) {
	if (qpsolver == 0 || qpsolver == 2) {
	    iw[208] = qpsolver;
	}
    }
    iw[209] = lvlpre;
/* Current precon mode */
    ngqp = *nnh;
    ngqp0 = max(ngqp,1);
    typelu = 3;
    lurequest = 0;
/* The first LU has been computed. */
    nonlinearcon = *nncon > 0;
    feasible = *nncon == 0;
/* ----------------------------------------------------------------- */
/* Find a feasible point. */
/* If the constraints are linear, x is already feasible. */
/* Otherwise, find a feasible x for this linearization. */
/* Minimize the sum of the elastic variables */
/* subject to keeping the non-elastic variables feasible. */
/* Elastic variables can move outside their bounds. */
/* ----------------------------------------------------------------- */
/* Loop back here in the unlikely event that the nonelastic bounds */
/* become infeasible while solving the QP subproblem. */
/* Status of the non-elastic variables */
L100:
    if (! feasible) {
	itqpmax = itnlim;
	itqptarget = itnlim;
	suboptimize = -1;
	if (normalin) {
/* Set eMode to switch to Elastic mode on infeasibility. */
	    emode = 1;
	} else {
/* Already in elastic mode. */
	    emode = 2;
	}
/* When getting feasible in elastic mode, set lvlObjFP to use */
/* the composite objective: */
/*  w1*Obj + w2*sInf,  with w1 = 0, w2 = wtInf. */
/* This minimizes the sum of the infeasibilities of the */
/* elastic variables subject to the nonelastic constraints. */
	lvlobjfp = 2;
	if (iw[208] == 0) {
	    *gotr = FALSE_;
	}
	luok = TRUE_;
	terminate = FALSE_;
/*        =============================================================== */
/* while (.not. Terminate  .and.  LUok) do */
L500:
	if (! terminate && luok) {
	    needlu = lurequest > 0;
	    needx = needlu;
	    maxsb = maxr;
	    s5qn_(&inform__, &c__4, probtag, &suboptimize, (U_fp)mnrlog, (
		    U_fp)hprod, (U_fp)hprod1, hvcalls, elastic, gotr, &needlu,
		     &typelu, &needx, lenr, m, &maxsb, mbs, n, nb, ndegen, &
		    ngqp0, &ngqp, nnobj0, nnobj, nnh0, nnh, ns, itqp, &
		    itqpmax, &itqptarget, itn, &emode, &lvlobjfp, mnrprtlvl, 
		    minimize, iobj, scaleobj, objadd, &objfp, condzhz, &
		    targetz, toloptfp, toloptqpk, tolx, ninf, sinf, elastics, 
		    ninfe, sinfe, wtinf, pinorm, nej, nlocj, &locj[1], &indj[
		    1], &jcol[1], neh, nloch, &loch[1], &indh[1], &hcol[1], &
		    etype[1], &estate[1], &feastype[1], &hs[1], &kbs[1], &bl[
		    1], &bu[1], &blqp[1], &buqp[1], &blbs[1], &bubs[1], &gbs[
		    1], &gobj[1], &gqp[1], &hdx[1], &pbs[1], &pi[1], &r__[1], 
		    &rc[1], &rg[1], &rg2[1], nncon0, nncon, &qprhs[1], &
		    scales[1], nnh0, nnh, &x[1], &xqp[1], &xbs[1], &x[1], &iy[
		    1], &iy1[1], &y[1], &y1[1], &y2[1], cu + 8, lencu, &iu[1],
		     leniu, &ru[1], lenru, cw + 8, lencw, &iw[1], leniw, &rw[
		    1], lenrw, (ftnlen)20, (ftnlen)8, (ftnlen)8);
/* Check for trouble.  Here are the possibilities: */

/* inform      Result */
/* ------      ------ */
/*  >0         Fatal LU error */
/*   0         Found a feasible point for the nonelastics */
/*  -1         The nonelastics are infeasible */
/*  -2         Phase 1 is unbounded */
/*  -3         Too many iterations */
/*  -4         Void */
/*  -5         Superbasic limit exceeded */
	    if (inform__ > 0) {
		*iexit = inform__;
/* Fatal LU error or time limit */
		goto L900;
	    }
	    terminate = inform__ == 0 || inform__ == -3 || inform__ == -5;
	    if (! terminate) {
/* ======================================================== */
/* Trouble. */
/* inform = -2 implies that phase 1 was unbounded, */
/*             which can only occur if a bad basis gives */
/*             a large search direction */
/* inform = -1 implies the nonelastics are infeasible, */
/*             which should not happen since we already */
/*             know a feasible point for the nonelastics. */
/* ======================================================== */
/* Treat both cases as infeasible. Repeatedly refactorize */
/* with tighter tols before declaring LC infeasibility. */
		inform__ = -1;
		s_wsfi(&io___199);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
		s2trylu_(itn, &c__22, ns, &lurequest, &luok, &typelu, &iw[1], 
			leniw, &rw[1], lenrw);
		if (! luok) {
		    *iexit = 15;
/* infeasible linear constraints in QP. */
		    goto L900;
		}
	    }
	    goto L500;
	}
/* end while */
/*        --------------------------------------------------------------- */
	if (inform__ < 0) {
	    goto L800;
	}
/* Itns or time limit exceeded */
	if (*elastic && normalin) {
/* The QP switched to elastic mode. */
/* The linearized constraints are infeasible. */
	    if (*mjrprtlvl >= 1 || *mnrprtlvl >= 10) {
		s_wsfi(&io___200);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&(*wtinf), (ftnlen)sizeof(doublereal));
		e_wsfi();
		snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
		iw[223] = 1;
		iw[225] = 1;
	    }
	} else if (*mjrprtlvl > 10 && *mnrprtlvl > 10) {
/* No change in mode. */
	    if (*elastic) {
		s_wsfi(&io___201);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__21, str, &iw[1], leniw, (ftnlen)80);
	    } else {
		s_wsfi(&io___202);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__21, str, &iw[1], leniw, (ftnlen)80);
	    }
	}
    }
/*     ------------------------------------------------------------------ */
/*     The inelastic variables (x's and linear slacks) are now feasible. */
/*     Save them in xQP0 for use with the BFGS update. */

/*     Solve the QP subproblem. */
/*     Loop back sometimes if we need a BS factorize. */
/*     ------------------------------------------------------------------ */
/* .not. Feasible */
    dcopy_(nb, &xqp[1], &c__1, &xqp0[1], &c__1);
    feasible = TRUE_;
/* For nonlinear constraints. */
/* Set lvlObjE to use the composite objective  Obj + wtInf*sInf */
/* after any switch to elastic mode. */
/* Linear constraints.  In theory, the subproblem should be */
/* feasible. If it is not, do not switch to Elastic mode. */
/* the nonelastics are feasible */
    if (nonlinearcon) {
	lvlobje = 1;
    } else {
	emode = 0;
	lvlobje = 0;
    }
    if (*nnh > 0) {
	suboptimize = 0;
    } else {
	suboptimize = -1;
    }
    itqpmax = itnlim;
    itqptarget = *itqp + mminor;
    lurequest = 0;
    typelu = 3;
    luok = TRUE_;
    terminate = FALSE_;
    solved = FALSE_;
/*     ================================================================== */
/*     while (.not. (Solved  .or.  Terminate)  .and.  LUok) do */
L600:
    if (! (solved || terminate) && luok) {
	inform__ = 0;
	needlu = lurequest > 0;
	needx = needlu;
	if (needlu) {
	    s2bfac_(iexit, &typelu, &needlu, &newlu, &newb, iobj, itn, 
		    mjrprtlvl, &lurequest, m, mbs, n, nb, nnh, ns, &nswap, 
		    nej, nlocj, &locj[1], &indj[1], &jcol[1], &kbs[1], &hs[1],
		     &blqp[1], &buqp[1], &blbs[1], &bubs[1], nncon0, nncon, &
		    qprhs[1], &xqp[1], &xbs[1], &iy[1], &iy1[1], &y[1], &y1[1]
		    , &iw[1], leniw, &rw[1], lenrw);
	    if (*iexit != 0) {
		goto L900;
	    }
	    if (nswap > 0) {
		*gotr = FALSE_;
	    }
/* Reset R. */
	    lurequest = 0;
	}
	if (*mnrprtlvl >= 1) {
	    iw[223] = 1;
/* QP print   header */
	    iw[225] = 1;
/* QP summary header */
	}
	if (iw[208] == 0) {
	    if (! (*gotr)) {
		if (*ns > 0) {
/* ----------------------------------------------------- */
/* Compute and factorize  Z'HZ using Cholesky. */
/* ----------------------------------------------------- */
		    s5zhz_(&inform__, (U_fp)hprod, (U_fp)hprod1, hvcalls, &
			    maxr, lenr, minimize, m, mbs, n, nb, nnh, ns, nej,
			     nlocj, &locj[1], &indj[1], &jcol[1], neh, nloch, 
			    &loch[1], &indh[1], &hcol[1], &zhzmin, &zhzmax, &
			    rw[192], &normz, &kbs[1], &r__[1], &y[1], &y1[1], 
			    &y2[1], cu + 8, lencu, &iu[1], leniu, &ru[1], 
			    lenru, cw + 8, lencw, &iw[1], leniw, &rw[1], 
			    lenrw, (ftnlen)8, (ftnlen)8);
		    if (inform__ != 0) {
			*iexit = inform__;
/* Fatal error when solving with LU. */
			goto L900;
		    }
		    s5zhzfac_(&inform__, &eigh, &c__0, itn, lenr, m, &maxr, 
			    mbs, nb, ns, &targeth, &zhzmin, &zhzmax, &normz, &
			    rankzhz, &hs[1], &kbs[1], &iy[1], &blqp[1], &buqp[
			    1], &blbs[1], &bubs[1], &xqp[1], &xbs[1], &r__[1],
			     &y1[1], &iw[1], leniw);
/* Possible values of inform are: */
/* inform  Status */
/* ------  ------ */
/*  -2     H singular (but should be positive definite) */
/*  -1     H indefinite */
/*   0     normal exit */
		    if (inform__ < 0) {
			inform__ = -9;
		    }
		    rcheck = FALSE_;
		    if (rcheck) {
			s5rcheck_(iexit, (U_fp)hprod, (U_fp)hprod1, hvcalls, &
				eigh, itn, minimize, &maxr, lenr, m, mbs, n, 
				nb, nnh, ns, nej, nlocj, &locj[1], &indj[1], &
				jcol[1], neh, nloch, &loch[1], &indh[1], &
				hcol[1], &kbs[1], &r__[1], &y[1], &y1[1], &y2[
				1], cu + 8, lencu, &iu[1], leniu, &ru[1], 
				lenru, cw + 8, lencw, &iw[1], leniw, &rw[1], 
				lenrw, (ftnlen)8, (ftnlen)8);
			if (*iexit != 0) {
			    goto L900;
			}
		    }
		}
/* nS > 0 */
		*gotr = inform__ == 0;
	    }
/* not gotR */
	    if (inform__ == 0) {
/* -------------------------------------------------------- */
/* Solve the QP subproblem using the Cholseky QP solver.. */
/* maxS = maxR to force termination at SB limit. */
/* -------------------------------------------------------- */
		needlu = FALSE_;
		needx = needlu;
		maxsb = maxr;
		s5qp_(&inform__, &c__5, probtag, &suboptimize, (U_fp)mnrlog, (
			U_fp)hprod, (U_fp)hprod1, hvcalls, &eigh, elastic, 
			gotr, &needlu, &typelu, &needx, lenr, m, &maxsb, mbs, 
			n, nb, ndegen, &ngqp0, &ngqp, nnobj0, nnobj, nnh0, 
			nnh, ns, itqp, &itqpmax, &itqptarget, itn, &emode, &
			lvlobje, mnrprtlvl, minimize, iobj, scaleobj, objadd, 
			objqp, &targeth, &targetz, toloptfp, toloptqpk, tolx, 
			ninf, sinf, elastics, ninfe, sinfe, wtinf, pinorm, &
			rgnorm, nej, nlocj, &locj[1], &indj[1], &jcol[1], neh,
			 nloch, &loch[1], &indh[1], &hcol[1], &etype[1], &
			estate[1], &feastype[1], &hs[1], &kbs[1], &bl[1], &bu[
			1], &blqp[1], &buqp[1], &blbs[1], &bubs[1], &gbs[1], &
			gobj[1], &gqp[1], &hdx[1], &pbs[1], &pi[1], &r__[1], &
			rc[1], &rg[1], nncon0, nncon, &qprhs[1], &scales[1], 
			nnh0, nnh, &x[1], &xqp[1], &xbs[1], &x[1], &iy[1], &
			iy1[1], &y[1], &y1[1], &y2[1], cu + 8, lencu, &iu[1], 
			leniu, &ru[1], lenru, cw + 8, lencw, &iw[1], leniw, &
			rw[1], lenrw, (ftnlen)20, (ftnlen)8, (ftnlen)8);
	    }
	} else if (iw[208] == 2 || iw[208] == 1) {
/*           ------------------------------------------------------------ */
/*           Solve the QP using the quasi-Newton or CG solver.. */
/*           maxS = maxR to avoid automatic switch to CG mode */
/*           ------------------------------------------------------------ */
	    if (iw[208] == 2) {
		maxsb = maxr;
	    } else {
		maxsb = *maxs;
	    }
	    s5qn_(&inform__, &c__5, probtag, &suboptimize, (U_fp)mnrlog, (
		    U_fp)hprod, (U_fp)hprod1, hvcalls, elastic, gotr, &needlu,
		     &typelu, &needx, lenr, m, &maxsb, mbs, n, nb, ndegen, &
		    ngqp0, &ngqp, nnobj0, nnobj, nnh0, nnh, ns, itqp, &
		    itqpmax, &itqptarget, itn, &emode, &lvlobje, mnrprtlvl, 
		    minimize, iobj, scaleobj, objadd, objqp, condzhz, &
		    targetz, toloptfp, toloptqpk, tolx, ninf, sinf, elastics, 
		    ninfe, sinfe, wtinf, pinorm, nej, nlocj, &locj[1], &indj[
		    1], &jcol[1], neh, nloch, &loch[1], &indh[1], &hcol[1], &
		    etype[1], &estate[1], &feastype[1], &hs[1], &kbs[1], &bl[
		    1], &bu[1], &blqp[1], &buqp[1], &blbs[1], &bubs[1], &gbs[
		    1], &gobj[1], &gqp[1], &hdx[1], &pbs[1], &pi[1], &r__[1], 
		    &rc[1], &rg[1], &rg2[1], nncon0, nncon, &qprhs[1], &
		    scales[1], nnh0, nnh, &x[1], &xqp[1], &xbs[1], &x[1], &iy[
		    1], &iy1[1], &y[1], &y1[1], &y2[1], cu + 8, lencu, &iu[1],
		     leniu, &ru[1], lenru, cw + 8, lencw, &iw[1], leniw, &rw[
		    1], lenrw, (ftnlen)20, (ftnlen)8, (ftnlen)8);
	}
/* inform   Result */
/* ------   ------ */
/*   >0     Fatal LU error or time limit */
/*    0     QP solution found */
/*   -1     The nonelastics are infeasible */
/*   -2     The QP subproblem is unbounded */
/*   -3     Too many iterations */
/*   -4     The QP subproblem has a weak minimizer   (s5QP) */
/*   -5     Too many superbasics */
/*   -6     Reduced Hessian not psd with new column  (s5QP) */
/*   -7     Z'g could not be made sufficiently small (s5QP) */
/*   -8     Ill-conditioned Z */
/*   -9     Reduced Hessian not psd when factored    (s5QP) */
/*  -10     Too many CG subspace iterations          (s5QN) */
/* QPSolver = QN or CG */
	if (inform__ > 0) {
	    *iexit = inform__;
/* Fatal LU error or time limit */
	    goto L900;
	}
	solved = inform__ == 0;
	terminate = inform__ == -3;
	resolve = inform__ == -1 || inform__ == -2 || inform__ == -4 || 
		inform__ == -5 || inform__ == -6 || inform__ == -7 || 
		inform__ == -8 || inform__ == -9 || inform__ == -10;
	if (solved) {
/* Finish if there are no large multipliers. */
/* Otherwise, set elastic mode and solve the QP again, */
	    if (! (*elastic) && nonlinearcon) {
		if (*pinorm > *wtinf) {
		    *elastic = TRUE_;
		    s_wsfi(&io___215);
		    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&(*wtinf), (ftnlen)sizeof(
			    doublereal));
		    e_wsfi();
		    snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
		    iw[223] = 1;
/* New Minor print   header */
		    iw[225] = 1;
/* New MinorLP summary header */
		    solved = FALSE_;
		}
	    }
	} else if (terminate) {
/* Relax */
	} else if (resolve) {
/* =========================================================== */
/* Neither solved nor terminated. */
/* inform = -1, -2, -4, -5, -6, -7, -8, -9, -10 */
/* Try to solve again with one or more of the following: */
/*  (1) with a different H */
/*  (2) a different QP solver */
/*  (3) or more accurate initial LU. */
/* =========================================================== */
	    if (inform__ == -1) {
/* -------------------------------------------------------- */
/* The nonelastics are infeasible. This should not happen */
/* because Phase 1 has already found a feasible point for */
/* the nonelastics. The basis must be ill-conditioned. */
/* Refactorize with tighter tols and restart at the known */
/* feasible point. */
/* -------------------------------------------------------- */
		s_wsfi(&io___216);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
		s2trylu_(itn, &c__20, ns, &lurequest, &luok, &typelu, &iw[1], 
			leniw, &rw[1], lenrw);
		if (! luok) {
		    *iexit = 15;
/* infeasible linear constraints in QP. */
		    goto L900;
		}
		*elastic = FALSE_;
		feasible = FALSE_;
		goto L100;
	    } else if (inform__ == -2) {
/*              --------------------------------------------------------- */
/*              The QP is unbounded. */
/*              ------------------------------------------------------ */
		s_wsfi(&io___217);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
		if (iw[208] == 0) {
/*                 The Reduced Hessian is positive semidefinite. */
/*                 Refactor for safety. */
		    if (iw[215] > 0) {
			s2trylu_(itn, &lurequest, ns, &lurequest, &luok, &
				typelu, &iw[1], leniw, &rw[1], lenrw);
			if (! luok) {
			    terminate = TRUE_;
			}
		    } else {
			terminate = TRUE_;
		    }
		} else if (iw[208] == 2) {
/*                 The Hessian is positive definite. */
/*                 Reset both the full and reduced Hessian. */
		    if (iw[243] != 2) {
			if (*nnh > 0) {
			    iw[243] = 2;
			    s8reseth_(nnh, ud0, &hd[1], &iw[1], leniw, &rw[1],
				     lenrw);
			}
			s_wsfi(&io___218);
			do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer)
				);
			e_wsfi();
			snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
			*gotr = FALSE_;
		    }
		    s2trylu_(itn, &c__25, ns, &lurequest, &luok, &typelu, &iw[
			    1], leniw, &rw[1], lenrw);
		    if (! luok) {
			terminate = TRUE_;
/* Done all we can. Terminate. */
		    }
		}
	    } else if (inform__ == -4) {
/*              --------------------------------------------------------- */
/*              Weak QP minimizer. */
/*              ------------------------------------------------------ */
		if (iw[215] > 0) {
		    s2trylu_(itn, &lurequest, ns, &lurequest, &luok, &typelu, 
			    &iw[1], leniw, &rw[1], lenrw);
		    if (! luok) {
			terminate = TRUE_;
		    }
		} else {
		    terminate = TRUE_;
		}
	    } else if (inform__ == -5) {
/*              --------------------------------------------------------- */
/*              Too many superbasics. */
/*              Switch to CG mode if possible. */
/*              ------------------------------------------------------ */
		if (maxr < *maxs) {
		    iw[208] = 1;
/* Switch to CG */
		    iw[209] = 0;
/* with no preconditioning */
		    *gotr = FALSE_;
		    iw[241] = 0;
		} else {
		    terminate = TRUE_;
		}
	    } else if (inform__ == -8) {
/*              --------------------------------------------------------- */
/*              condZ > targetZ  while computing the search direction. */
/*              Refactorize B, possibly with a reduced factor tol. If */
/*              the factor tol is already tight, accept Z, however bad. */
/*              --------------------------------------------------------- */
		s_wsfi(&io___219);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&rw[192], (ftnlen)sizeof(doublereal));
		e_wsfi();
		snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
		s2trylu_(itn, &c__25, ns, &lurequest, &luok, &typelu, &iw[1], 
			leniw, &rw[1], lenrw);
		if (! luok) {
		    targetz = plinfy;
		    luok = TRUE_;
		}
	    } else if (inform__ == -6 || inform__ == -7 || inform__ == -9 || 
		    inform__ == -10) {
/*              --------------------------------------------------------- */
/*              Either Z'HZ is not positive definite or Z'g is large. */
/*              Most likely, Z'HZ is ill-conditioned. */

/*              If the size of norm of Z (relative to J) is large, try */
/*              to get a better Z by refactorizing B. This is tried once. */
/*              Otherwise, the off-diagonals of H are discarded, or H is */
/*              reset to a diagonal if it is already diagonal. */
/*              Refactorize B, possibly with a reduced factor tol. If */
/*              the factor tol is already tight, accept Z, however bad. */
/*              ---------------------------------------------------------- */
		if (inform__ == -6 || inform__ == -9) {
		    s_wsfi(&io___220);
		    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		    e_wsfi();
		} else if (inform__ == -7) {
		    s_wsfi(&io___221);
		    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		    e_wsfi();
		} else if (inform__ == -10) {
		    s_wsfi(&io___222);
		    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		    e_wsfi();
		}
		snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
		if (iw[243] != 2) {
/*                 Discard the off-diagonals of H. */
/*                 If H is already diagonal, H is set to the identity. */
		    s8reseth_(nnh, ud0, &hd[1], &iw[1], leniw, &rw[1], lenrw);
		    s_wsfi(&io___223);
		    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		    e_wsfi();
		    snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
		    *gotr = FALSE_;
		} else if (iw[243] == 2) {
/*                 H = I, so Z'Z is not positive definite. */
		    if (inform__ == -6 || inform__ == -9) {
			lurequest = 25;
		    } else if (inform__ == -7) {
			lurequest = 21;
		    }
		    s2trylu_(itn, &lurequest, ns, &lurequest, &luok, &typelu, 
			    &iw[1], leniw, &rw[1], lenrw);
		    if (! luok) {
			*iexit = 44;
/* ill-conditioned null-space basis */
			goto L900;
		    }
		}
	    }
/* inform ne 0 */
	}
	goto L600;
    }
/* end while not solved or terminated */
/*     ------------------------------------------------------------------ */
L800:
    if (*ninfe > 0) {
	iw[240] = 1;
    }
    if (inform__ == 0) {
	iw[241] = max(suboptimize,0);
/* = 0,1 or 2 */
    } else if (inform__ == -1) {
	*iexit = 15;
/* infeasible nonelastics */
    } else if (inform__ == -2) {
	iw[241] = 3;
/* unbounded subproblem */
    } else if (inform__ == -3) {
	iw[241] = 2;
	*iexit = -1;
/* too many iterations */
    } else if (inform__ == -4) {
	iw[241] = 4;
/* weak QP solution */
    } else if (inform__ == -5 && feasible) {
	iw[241] = 5;
	*iexit = -2;
/* superbasic limit */
    } else if (inform__ == -5) {
	*iexit = 33;
/* superbasic limit */
    }
L900:
    return 0;
} /* s8solveqp_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8solveQP */
/* Subroutine */ int s8steplimits_(logical *fonlyls, integer *nb, integer *
	nncon, integer *nnobj, integer *majors, integer *qnskips, doublereal *
	step, doublereal *stepmin, doublereal *steplimit, doublereal *stepmax,
	 doublereal *tolz, doublereal *dxnorm, doublereal *xnorm, doublereal *
	bl, doublereal *bu, doublereal *x, doublereal *dx, integer *iw, 
	integer *leniw, doublereal *rw, integer *lenrw)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Local variables */
    static logical overflow;
    static integer j;
    static doublereal res;
    extern doublereal ddiv_(doublereal *, doublereal *, logical *);
    static doublereal tolp, bigdx;
    static integer gotfd;
    static doublereal xdlim, pivot, fdint1;
    static integer hdinfo, lvldif;
    static doublereal pivabs;
    static logical switch__;
    static doublereal stepqp, tolpiv;

/*     ================================================================== */
/*     s8step  finds the maximum, minimum and initial value for the */
/*     linesearch step. */

/*     For problems with nonlinear constraints, the maximum step stepMax */
/*     is one.  If there are only linear constraints the maximum step is */
/*     the largest step such that x + step*dx  reaches one of its bounds. */

/*     All step sizes are subject to the user-specified limit  stepLimit. */

/*     04 Dec 1992: First version of s8step based on npsol routine npalf. */
/*     31 Mar 2000: Updated for SNOPT 6.1. */
/*     19 Mar 2001: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --dx;
    --x;
    --bu;
    --bl;
    --iw;
    --rw;

    /* Function Body */
    tolpiv = rw[60];
/* excludes small elements of y. */
    bigdx = rw[72];
/* unbounded step. */
    fdint1 = rw[76];
/* (1) forwrd diff. interval */
    xdlim = rw[80];
/* Step limit */
    lvldif = iw[181];
/* 1(2) for forwd (cntrl) diffs */
    gotfd = iw[183];
/* > 0 => some differences needed */
    hdinfo = iw[243];
/* infoTag(7): Approximate Hessian type */
    overflow = FALSE_;
/*     ================================================================== */
/*     Switch   indicates if there is an option to switch to */
/*              central differences to get a better search direction. */
/*     stepQP   is the step predicted by the QP subproblem (usually 1). */
/*     stepMax  is the largest feasible steplength subject to a */
/*              user-defined limit, bigdx, on the change in  x. */
/*     step     is initialized subject to a user-defined limit, xdlim. */
/*     ================================================================== */
    if (*nncon == 0 && *nnobj == 0) {
/*        Linear program */
/* LP !! */
	*step = 1.;
	*stepmin = 1.;
	*steplimit = 1.;
	*stepmax = 1.;
    } else {
/*        Nolinear program */
	switch__ = gotfd > 0 && lvldif != 2;
	*stepmin = 0.;
	if (*fonlyls && switch__) {
	    *stepmin = fdint1 * (*xnorm + 1.) / *dxnorm;
	}
	stepqp = 1.;
	if (*nncon > 0 && (*qnskips == 0 || hdinfo != 2)) {
	    *stepmax = 1.;
	} else {
	    tolp = tolpiv * *dxnorm;
	    *stepmax = ddiv_(&bigdx, dxnorm, &overflow);
	    *step = *stepmax;
	    j = 1;
/* +          while (j .le. nb  .and.  step .gt. stepQP) do */
L100:
	    if (j <= *nb && *step > stepqp) {
		pivot = dx[j];
		pivabs = abs(pivot);
		if (pivabs > tolp) {
		    if (pivot <= 0.) {
			res = x[j] - bl[j];
			if (*step * pivabs > res) {
			    *step = res / pivabs;
			}
		    } else {
			res = bu[j] - x[j];
			if (*step * pivabs > res) {
			    *step = res / pivabs;
			}
		    }
		}
		++j;
		goto L100;
/* +          end while */
	    }
	    *step = max(*step,stepqp);
	    if (*step < stepqp + *tolz) {
		*step = stepqp;
	    }
	    *stepmax = *step;
	}
	d__1 = (*xnorm + 1.) * xdlim;
	*steplimit = ddiv_(&d__1, dxnorm, &overflow);
	if (*majors <= 1) {
/* Computing MIN */
	    d__1 = *steplimit, d__2 = ddiv_(&c_b12, dxnorm, &overflow);
	    *steplimit = min(d__1,d__2);
	}
	*stepmax = min(*steplimit,*stepmax);
	*step = min(*steplimit,1.);
    }
    return 0;
} /* s8steplimits_ */

