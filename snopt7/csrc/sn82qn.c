/* ./src/sn82qn.f -- translated by f2c (version 20100827).
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

static doublereal c_b2 = 0.;
static integer c__1 = 1;
static integer c__0 = 0;
static doublereal c_b13 = -1.;
static doublereal c_b25 = 1.;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*     File  sn82qn.f */

/*     s8getH     s8ResetH   s8HQN    s8Hwrapper   s8Hx   s8xHx */
/*     s8Hupdate  s8HmodA    s8HmodB */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Subroutine */ int s8geth_(integer *nnh, integer *lenh, doublereal *u, 
	doublereal *h__, doublereal *y, doublereal *y1)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j;
    extern /* Subroutine */ int dload_(integer *, doublereal *, doublereal *, 
	    integer *), dcopy_(integer *, doublereal *, integer *, doublereal 
	    *, integer *);
    static integer jthcol;
    extern /* Subroutine */ int s6rprod_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *);

/*     ================================================================== */
/*     s8getH  computes the product H = U'U, where  U is the Cholesky */
/*     factor of the approximate Hessian of the Lagrangian.  The matrix */
/*     U is stored by rows in the one-dimensional array  U. */
/*     lenH defines the length of U.  lenH must be at least */
/*     nnH*(nnH + 1)/2.  The result is stored by columns in the upper */
/*     triangular array H. */

/*     03 Sep 2006: First version of s8getH. */
/*     03 Sep 2006: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
/*     Compute the product y1 = U'Uy, where  U is an upper- */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --y1;
    --y;
    --h__;
    --u;

    /* Function Body */
    jthcol = 1;
    i__1 = *nnh;
    for (j = 1; j <= i__1; ++j) {
	jthcol = jthcol + j - 1;
	dload_(nnh, &c_b2, &y[1], &c__1);
	y[j] = 1.;
	s6rprod_(&c__0, nnh, nnh, lenh, &u[1], &y[1], &y1[1]);
	s6rprod_(&c__1, nnh, nnh, lenh, &u[1], &y1[1], &y[1]);
	dcopy_(&j, &y[1], &c__1, &h__[jthcol], &c__1);
    }
    return 0;
} /* s8geth_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8getH */
/* Subroutine */ int s8reseth_(integer *nnh, doublereal *ud0, doublereal *hd, 
	integer *iw, integer *leniw, doublereal *rw, integer *lenrw)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static doublereal condumax;
    static logical overflow;
    static integer j, lu, lu0;
    extern doublereal ddiv_(doublereal *, doublereal *, logical *);
    static integer lenu;
    extern /* Subroutine */ int dload_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal hdmin, hdmax;
    extern /* Subroutine */ int s8fmh0_(integer *, doublereal *, integer *, 
	    doublereal *), s8lmh0_(integer *, doublereal *, doublereal *);
    static doublereal condhd;
    static integer hdtype, lvlhess;

/*     ================================================================== */
/*     s8ResetH  initializes or resets  the BFGS approximate Hessian. */

/*     The BFGS Hessian is  H = U'U,  where the form of U is defined */
/*     according to the definition of lvlHess. Copies of the diagonals of */
/*     H are held in HD. */

/*     At a reset,  H is replaced by a diagonal matrix. */

/*     If HD is well-conditioned, H is reset to HD (i.e., the off */
/*     diagonals of H are discarded). This strategy implies that HD is */
/*     always the diagonal of H,  regardless of resets. */

/*     If HD is ill-conditioned, H and HD are reset to the scaled */
/*     identity matrix. */

/*     On entry, */
/*     --------- */
/*      UD0    is the new diagonal of U if U is reset to a diagonal. */
/*             Its value reflects the scale of the gradient and is set */
/*             in s8SQP. */

/*      HD     are the diagonal elements of the BFGS QN Hessian. */
/*             i.e.,  HD = diag(U*U) */

/*      HDInfo defines the input status of H. */
/*             UnSet  (-1)  => U not set */
/*             Normal ( 0)  => U'U is a non-diagonal BFGS approx */
/*             Diag   ( 1)  => U'U is general diagonal */
/*             Unit   ( 2)  => U'U is a multiple of the identity */

/*     On exit */
/*     ------- */
/*      HD     is defined according to the input value of HDInfo */

/*             HDInfo */
/*             ----- */
/*             UnSet   =>    U =     U0*I,   HD = (U0*U0)*I */
/*             Normal  =>    U = sqrt(HD),   HD   unchanged if HD is good */
/*                        or U =     U0*I,   HD = (U0*U0)*I if HD is bad */
/*             Diag    =>    U =     U0*I,   HD = (U0*U0)*I */
/*             Unit    =>    U =     U0*I,   HD = (U0*U0)*I */

/*     HDInfo  indicates what happened. */
/*             Diag    =>    normal reset,   HD not changed. */
/*             Unit    =>    Unit   reset or HD was bad. */

/*             In all cases, diag(U'U) = HD */

/*     iw(QNmods), the number of BFGS updates since the last reset, */
/*                 is set to 0. */

/*     19 Nov 2012: First version based on SNOPT routine s8H0. */
/*     15 Aug 2013: Last update. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/* Approximate Hessian type */
/*     ------------------------------------------------------------------ */
/* BFGS updates since last reset */
    /* Parameter adjustments */
    --hd;
    --iw;
    --rw;

    /* Function Body */
    condumax = rw[87];
/* max cond estimator for U with H = U'U */
    lvlhess = iw[72];
/* LM, FM or Exact Hessian */
    hdtype = iw[243];
    if (hdtype == 0) {
/*        --------------------------------------------------------------- */
/*        Try and set U so that U'U = HD, where HD is the diagonal of the */
/*        Hessian.  First, check the condition of HD and reset it to */
/*        UD0*UD0*I  if its ill-conditioned or not positive definite. */
/*        --------------------------------------------------------------- */
	overflow = FALSE_;
	hdmin = hd[1];
/* strictly positive in exact arithmetic */
	hdmax = hdmin;
	i__1 = *nnh;
	for (j = 2; j <= i__1; ++j) {
/* Computing MIN */
	    d__1 = hd[j];
	    hdmin = min(d__1,hdmin);
/* Computing MAX */
	    d__1 = hd[j];
	    hdmax = max(d__1,hdmax);
	}
	condhd = ddiv_(&hdmax, &hdmin, &overflow);
	if (hdmin <= 0. || condhd >= condumax * condumax) {
	    hdtype = 2;
	}
    }
    if (hdtype == 0) {
	iw[243] = 1;
/* Set U to sqrt(HD) */
    } else {
/* HDInfo == Unset, Diag or Unit */
	iw[243] = 2;
/* Set diag(U) to UD0*I */
	d__1 = *ud0 * *ud0;
	dload_(nnh, &d__1, &hd[1], &c__1);
    }
/*     ------------------------------------------------------------------ */
/*     Zero the off-diagonal elements of U. */
/*     How this is done depends on the way that U is stored. */
/*     ------------------------------------------------------------------ */
    if (lvlhess == 0) {
/*        ----------------------- */
/*        Limited-memory Hessian. */
/*        ----------------------- */
	lu0 = iw[346];
/* Square root of initial BFGS diagonal */
	s8lmh0_(nnh, &hd[1], &rw[lu0]);
    } else if (lvlhess == 1) {
/*        ----------------------- */
/*        Full-memory Hessian. */
/*        ----------------------- */
	lu = iw[391];
/* U(lenU), full-memory BFGS H = U'U */
	lenu = iw[392];

	s8fmh0_(nnh, &hd[1], &lenu, &rw[lu]);
    }
    iw[381] = 0;
    return 0;
} /* s8reseth_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8ResetH */
/* Subroutine */ int s8hqn_(integer *iexit, U_fp funwrapper, U_fp funcon, 
	U_fp funobj, U_fp userhv, logical *usefd, integer *qptype, integer *
	starttype, integer *lenr, integer *m, integer *mbs, integer *n, 
	integer *nb, integer *nncon0, integer *nncon, integer *nnjac, integer 
	*nnh, integer *nnobj0, integer *nnobj, integer *ns, integer *nmajor, 
	integer *qnskips, doublereal *ud0, doublereal *step, integer *
	minimize, doublereal *dxhdx, integer *rtrmods, logical *gotr, 
	doublereal *penparm, doublereal *fobj, doublereal *fcon, doublereal *
	gcon, doublereal *gobj, doublereal *fcon1, doublereal *gcon1, 
	doublereal *gobj1, integer *nej, integer *nlocj, integer *locj, 
	integer *indj, doublereal *jcol, integer *neh, integer *nloch, 
	integer *loch, integer *indh, doublereal *hcol, integer *negcon, 
	integer *nlocg, integer *locg, integer *kbs, doublereal *bl, 
	doublereal *bu, doublereal *dx, doublereal *dg, doublereal *udx, 
	doublereal *hdx, doublereal *hd, doublereal *ycon1, doublereal *r__, 
	doublereal *x, doublereal *x1, doublereal *xqp0, doublereal *xpen, 
	doublereal *y, doublereal *y1, doublereal *y2, char *cu, integer *
	lencu, integer *iu, integer *leniu, doublereal *ru, integer *lenru, 
	char *cw, integer *lencw, integer *iw, integer *leniw, doublereal *rw,
	 integer *lenrw, ftnlen cu_len, ftnlen cw_len)
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    extern /* Subroutine */ int s8reseth_(integer *, doublereal *, doublereal 
	    *, integer *, integer *, doublereal *, integer *);
    static logical overflow;
    static integer mqnskips;
    extern /* Subroutine */ int s8hupdate_(integer *, integer *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, integer *), s6rupdate_(integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, doublereal *, integer *), s8initpen_(integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *)
	    ;
    static integer nbs;
    static doublereal ydx, eps0, eps1;
    static integer kfac;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *), ddiv_(doublereal *, doublereal *, logical *);
    static doublereal rnnh;
    static integer maxr;
    static doublereal rydx;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static logical nonlinearcon;
    static doublereal xpen0;
    extern /* Subroutine */ int dload_(integer *, doublereal *, doublereal *, 
	    integer *), dscal_(integer *, doublereal *, doublereal *, integer 
	    *), dcopy_(integer *, doublereal *, integer *, doublereal *, 
	    integer *), daxpy_(integer *, doublereal *, doublereal *, integer 
	    *, doublereal *, integer *);
    static integer nzero;
    static doublereal snorm;
    static integer inform__;
    static doublereal glnorm, rdxhdx, penunm, ydxmin, u0scale;
    extern /* Subroutine */ int s8hmoda_(integer *, U_fp, U_fp, U_fp, U_fp, 
	    logical *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     integer *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, ftnlen, ftnlen), 
	    s8hmodb_(integer *, integer *, doublereal *, integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), s8gprod_(
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *);
    static logical updated;
    static doublereal signobj;

/*     ================================================================== */
/*     s8HQN  does the quasi-Newton update with vectors */
/*        dx = x1 - x   and   dg = gL(x1) - gL(x). */

/*     On entry: */
/*      xQP is the QP solution. */

/*     23 Apr 1999: First version of s8HQN, */
/*     18 Feb 2001: LM H stored in product form. */
/*     12 Oct 2003: snEXIT and SNPRNT adopted */
/*     10 Jan 2005: FM H stored in product form. */
/*     23 Jun 2008: y3 no longer an argument. */
/*     19 Oct 2014: Added user-defined Hessian. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/* Tags(1) */
/* Hessian mod type */
/*     ------------------------------------------------------------------ */
/* Approximate Hessian type */
    /* Parameter adjustments */
    --r__;
    --kbs;
    --x;
    --y2;
    --y1;
    --y;
    --xqp0;
    --bu;
    --bl;
    --xpen;
    --ycon1;
    --fcon1;
    --fcon;
    --x1;
    --hd;
    --hdx;
    --udx;
    --dg;
    --dx;
    --gobj1;
    --gobj;
    --penparm;
    --jcol;
    --indj;
    --locj;
    --hcol;
    --indh;
    --loch;
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
    maxr = iw[52];
/* max columns of R */
    kfac = iw[59];
/* factorization frequency */
    mqnskips = iw[67];
/* # largest allowable  QNskips */
    eps0 = rw[2];
/* eps**(4/5) */
    eps1 = rw[3];
/* eps**(2/3) */
    xpen0 = rw[89];
/* initial penalty parameter. */
    *iexit = 0;
    overflow = FALSE_;
    nonlinearcon = *nncon > 0;
    nbs = *m + *ns;
    signobj = (doublereal) (*minimize);
    iw[237] = 0;
    iw[238] = 0;
    ydx = 0.;
/*     --------------------------------------------------------------- */
/*     Compute  dx = x1 - x  and  dg = gL1 - gL. */
/*     Compute the approx. curvature ydx and new scale factor U0. */
/*     --------------------------------------------------------------- */
    dcopy_(nnh, &x1[1], &c__1, &dx[1], &c__1);
    daxpy_(nnh, &c_b13, &x[1], &c__1, &dx[1], &c__1);
    dscal_(nnh, step, &hdx[1], &c__1);
    dscal_(nnh, step, &udx[1], &c__1);
    *dxhdx = *dxhdx * *step * *step;
    if (*nnobj > 0) {
	dcopy_(nnobj, &gobj1[1], &c__1, &dg[1], &c__1);
	if (*minimize < 0) {
	    dscal_(nnobj, &signobj, &dg[1], &c__1);
	}
    }
    nzero = *nnh - *nnobj;
    if (nzero > 0) {
	dload_(&nzero, &c_b2, &dg[*nnobj + 1], &c__1);
    }
    if (*nncon > 0) {
	s8gprod_(&c__1, &eps0, nej, nlocj, &locj[1], &indj[1], negcon, nlocg, 
		&locg[1], &gcon1[1], &c_b13, &ycon1[1], nncon, &c_b25, &dg[1],
		 nnjac);
    }
/*     gLnorm = dnormi( nnH, dg, 1 ) */
    glnorm = dnrm2_(nnh, &dg[1], &c__1);
    if (*nnobj > 0) {
	d__1 = -signobj;
	daxpy_(nnobj, &d__1, &gobj[1], &c__1, &dg[1], &c__1);
    }
    if (*nncon > 0) {
	s8gprod_(&c__1, &eps0, nej, nlocj, &locj[1], &indj[1], negcon, nlocg, 
		&locg[1], &gcon[1], &c_b25, &ycon1[1], nncon, &c_b25, &dg[1], 
		nnjac);
    }
    ydx = ddot_(nnh, &dg[1], &c__1, &dx[1], &c__1);
    if (*nmajor == 1 && *starttype != 3) {
/*        =============================================================== */
/*        The BFGS update is not applied on the first iteration, but the */
/*        latest curvature information is used to get a better scaled H. */
/*        =============================================================== */
	if (glnorm > 0.) {
	    rnnh = (doublereal) (*nnh);
	    *ud0 = sqrt(glnorm / sqrt(rnnh));
	} else {
	    *ud0 = 1.;
	}
	*gotr = FALSE_;
/* Computing MIN */
	d__1 = max(*ud0,.01);
	*ud0 = min(d__1,10.);
	s8reseth_(nnh, ud0, &hd[1], &iw[1], leniw, &rw[1], lenrw);
    } else {
/*        =============================================================== */
/*        Except on the first iteration, attempt a BFGS update. */
/*        Compute the smallest allowable curvature. */
/*        If the update cannot be done, s8HmodA attempts to find a */
/*        modified update using  dx = x1 - x defined with a new x. */
/*        Arrays fCon, gCon and gObj must be redefined at the new x. */
/*        =============================================================== */
	snorm = dnrm2_(nnh, &dx[1], &c__1);
	d__1 = sqrt((abs(ydx)));
	*ud0 = ddiv_(&d__1, &snorm, &overflow);
/* Computing MIN */
	d__1 = max(*ud0,.01);
	*ud0 = min(d__1,10.);
	penunm = 0.;
	ydxmin = *dxhdx * .001;
	updated = *dxhdx > 0. && (ydx >= ydxmin || ydx >= eps1);
	if (nonlinearcon && ! updated) {
/*           ------------------------------------------------------------ */
/*           Redefine  x, dx, Hdx and dg. */
/*           The problem functions are recomputed at x. */
/*           ------------------------------------------------------------ */
	    s8hmoda_(&inform__, (U_fp)funwrapper, (U_fp)funcon, (U_fp)funobj, 
		    (U_fp)userhv, usefd, n, nb, nncon0, nncon, nnjac, nnobj0, 
		    nnobj, nnh, minimize, step, dxhdx, &ydx, fobj, &fcon[1], &
		    gcon[1], &gobj[1], &gcon1[1], &gobj1[1], nej, nlocj, &
		    locj[1], &indj[1], neh, nloch, &loch[1], &indh[1], &hcol[
		    1], negcon, nlocg, &locg[1], &bl[1], &bu[1], &dx[1], &dg[
		    1], &udx[1], &hdx[1], &ycon1[1], &y1[1], &x[1], &xqp0[1], 
		    &y[1], cu + 8, lencu, &iu[1], leniu, &ru[1], lenru, cw + 
		    8, lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)8, (
		    ftnlen)8);
	    if (inform__ > 0) {
		*iexit = inform__;
/* User wants to stop */
		goto L999;
	    }
	    ydxmin = *dxhdx * .001;
	    updated = *dxhdx > 0. && (ydx >= ydxmin || ydx >= eps1);
	    if (updated) {
		iw[238] = 1;
/* Mod A  succeeded */
	    }
	    if (! updated && *dxhdx > 0.) {
/*              --------------------------------------------------------- */
/*              If all else fails, attempt to update the Hessian of */
/*              the augmented Lagrangian. */
/*              The target ydx is defined via tolg2. */
/*              --------------------------------------------------------- */
		ydxmin = *dxhdx * .1;
		s8hmodb_(nncon, nnjac, &eps0, nej, nlocj, &locj[1], &indj[1], 
			negcon, nlocg, &locg[1], &ydx, &ydxmin, &penunm, &
			fcon[1], &fcon1[1], &gcon[1], &gcon1[1], &dx[1], &dg[
			1], &y[1], &y1[1], &y2[1]);
		updated = ydx >= ydxmin;
		if (updated) {
		    iw[238] = 2;
/* Mod A + B succeeded */
		}
	    }
	}
/* NonlinearCon */
	if (updated) {
/*           ------------------------------------------------------------ */
/*           Update the approximate Hessian using (dg,Hdx). */
/*           If there are no nonlinear constraints,  apply the update */
/*           to the reduced Hessian. */
/*           ------------------------------------------------------------ */
	    *qnskips = 0;
	    if (ydx >= ydxmin && iw[243] == 0) {
		iw[237] = 0;
/* conventional BFGS */
	    } else {
		iw[237] = 1;
/* self-scaled  BFGS */
	    }
	    rydx = sqrt(ydx);
	    rdxhdx = sqrt(*dxhdx);
	    s8hupdate_(&iw[237], nnh, ud0, &u0scale, &rydx, &rdxhdx, &hd[1], &
		    dx[1], &hdx[1], &dg[1], &iw[1], leniw, &rw[1], lenrw);
	    if (*qptype == 0) {
		*gotr = *gotr && *nncon == 0 && *ns > 0 && iw[243] == 0 && *
			rtrmods < kfac;
	    } else {
		*gotr = *gotr && *ns > 0 && iw[243] == 0;
	    }
	    if (*gotr) {
		s6rupdate_(&iw[237], &maxr, lenr, m, n, &nbs, nnh, ns, &
			u0scale, &rdxhdx, nej, nlocj, &locj[1], &indj[1], &
			jcol[1], &kbs[1], &dg[1], &hdx[1], &r__[1], &y[1], &
			y1[1], &y2[1], &iw[1], leniw, &rw[1], lenrw);
		++(*rtrmods);
	    } else {
		*rtrmods = 0;
	    }
	} else {
/*           ------------------------------------------------------------ */
/*           No suitable update pair (dg,Hdx) could be found. */
/*           Skip the update.  Too many skips and we reset. */
/*           ------------------------------------------------------------ */
	    ++(*qnskips);
/*           Apply all updates to H and discard the off-diagonals. */
	    if (*qnskips % mqnskips == 0) {
		s8reseth_(nnh, ud0, &hd[1], &iw[1], leniw, &rw[1], lenrw);
		if (*qnskips % (mqnskips << 1) == 0) {
/*                 ------------------------------------------------------ */
/*                 Reset the multipliers and penalty parameters */
/*                 ------------------------------------------------------ */
		    s8initpen_(nncon, &penparm[1], &xpen0, &xpen[1], &rw[1], 
			    lenrw);
		    dload_(nncon, &c_b2, &ycon1[1], &c__1);
		}
	    }
	}
    }
/* nMajor > 1 */
L999:
    return 0;
} /* s8hqn_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8HQN */
/* Subroutine */ int s8hwrapper_(S_fp hprod, integer *nnh, integer *neh, 
	integer *nloch, integer *loch, integer *indh, doublereal *hcol, 
	doublereal *x, doublereal *hx, integer *status, char *cu, integer *
	lencu, integer *iu, integer *leniu, doublereal *ru, integer *lenru, 
	char *cw, integer *lencw, integer *iw, integer *leniw, doublereal *rw,
	 integer *lenrw, ftnlen cu_len, ftnlen cw_len)
{
/*     ================================================================== */
/*     s8Hwrapper wraps Hprod, which multiplies the QP Hessian H by the */
/*     vector  x.   It is called by the QP solver. */

/*     On entry: */
/*        Status  = 0  => a normal call for H*x. */
/*        Status  = 1  => the first entry for a given QP. */
/*        Status ge 2  => last call for a given QP. Status = 2+iExit. */

/*     On exit: */
/*        Status lt 0   the user wants to stop. */

/*     03 Nov 2000: First version of s8Hwrapper. */
/*     03 Nov 2014: neH, indH, locH, Hcol added as arguments. */
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
    (*hprod)(nnh, &x[1], &hx[1], cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, 
	    (ftnlen)8);
    return 0;
} /* s8hwrapper_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8Hwrapper */
/* Subroutine */ int s8hx_(integer *nnh, doublereal *x, doublereal *hx, char *
	cw, integer *lencw, integer *iw, integer *leniw, doublereal *rw, 
	integer *lenrw, ftnlen cw_len)
{
    static integer minimize, ls, lu, lv, lu0, lux, lenu;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), s8fmhx_(integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *), s8lmhx_(integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *);
    static doublereal signobj;
    static integer mqnmods, lvlhess;

/*     ================================================================== */
/*     s8Hx  multiplies the QP Hessian  H by the vector  x. */
/*     It is used to define Hx for the QP subproblem. */

/*     This routine is called by a general QP solver, which will rescale */
/*     Hx by  signObj  when maximizing. */

/*     s8Hx calls one of the Hessian routines s8LMH, s8FMH, s8SDH, ... */
/*     according to the value of the options lvlDer and lvlHess. */
/*     Each of these routines defines a particular form of the Hessian. */
/*     At the moment the options are: */

/*        lvlHess = LM      Limited-Memory (LM) BFGS  (the default). */
/*        lvlHess = FM      Full-Memory    (FM) BFGS */
/*        lvlHess = Exact   FD or exact Hessian */

/*     30 Dec 1991: First version of s8Hx. */
/*     12 Jan 1996: Full memory Hessian option added. */
/*     04 Apr 1999: Exact and FD Hessian option added. */
/*     18 Feb 2001: LM H stored in product form. */
/*     10 Jan 2005: FM H stored in product form. */
/*     23 Jun 2008: Exact option added. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --hx;
    --x;
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    lvlhess = iw[72];
/* LM, FM or Exact Hessian */
    minimize = iw[199];
/* (-1)(+1)    => (max)(min) */
    lux = iw[345];
/* Ux(nnH)     = product of U with x */
    signobj = (doublereal) minimize;
    if (lvlhess == 0) {
/*        ----------------------- */
/*        Limited-memory Hessian. */
/*        ----------------------- */
	mqnmods = iw[54];
/* (ge 0) max # of BFGS updates */
	lu0 = iw[346];
/* Square root of initial BFGS diagonal */
	ls = iw[401];
/* sk's for BFGS products: (I + sk*vk') */
	lv = iw[402];
/* vk's for BFGS products: (I + sk*vk') */
	s8lmhx_(nnh, &x[1], &rw[lux], &hx[1], &mqnmods, &iw[381], &rw[lu0], &
		rw[ls], &rw[lv]);
    } else if (lvlhess == 1) {
/*        ----------------------- */
/*        Full-memory Hessian. */
/*        ----------------------- */
	lu = iw[391];

	lenu = iw[392];

	s8fmhx_(nnh, &x[1], &rw[lux], &hx[1], &lenu, &rw[lu]);
    }
    if (minimize < 0) {
	dscal_(nnh, &signobj, &hx[1], &c__1);
    }
    return 0;
} /* s8hx_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8Hx */
/* Subroutine */ int s8xhx_(integer *nnh, doublereal *x, doublereal *ux, 
	doublereal *hx, doublereal *xhx, integer *iw, integer *leniw, 
	doublereal *rw, integer *lenrw)
{
    static integer ls, lu, lv, lu0;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer lenu;
    extern /* Subroutine */ int s8fmhx_(integer *, doublereal *, doublereal *,
	     doublereal *, integer *, doublereal *), s8lmhx_(integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *);
    static integer mqnmods, lvlhess;

/*     ================================================================== */
/*     s8xHx  computes x'Hx and Hx, where H = U'U. */

/*     s8xHx calls one of the Hessian routines s8LMH, s8FMH, s8SDH, ... */
/*     according to the value of the options lvlDer and lvlHess. */
/*     Each of these routines defines a particular form of the Hessian. */
/*     At the moment the options are: */

/*        lvlHess = LM      Limited-Memory (LM) BFGS  (the default). */
/*        lvlHess = FM      Full-Memory    (FM) BFGS */
/*        lvlHess = Exact   FD or exact Hessian */

/*     10 Jan 2005: First version of s8Hx based on s8Hx */
/*     18 Feb 2001: LM H stored in product form. */
/*     10 Jan 2005: FM H stored in product form. */
/*     23 Jun 2008: Exact option added. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --hx;
    --ux;
    --x;
    --iw;
    --rw;

    /* Function Body */
    lvlhess = iw[72];
/* LM, FM or Exact Hessian */
    if (lvlhess == 0) {
/*        ----------------------- */
/*        Limited memory Hessian. */
/*        ----------------------- */
	mqnmods = iw[54];
/* (ge 0) max # of BFGS updates */
	lu0 = iw[346];
/* Square root of initial BFGS diagonal */
	ls = iw[401];
/* sk's for BFGS products: (I + sk*vk') */
	lv = iw[402];
/* vk's for BFGS products: (I + sk*vk') */
	s8lmhx_(nnh, &x[1], &ux[1], &hx[1], &mqnmods, &iw[381], &rw[lu0], &rw[
		ls], &rw[lv]);
    } else if (lvlhess == 1) {
/*        ----------------------- */
/*        Full memory Hessian. */
/*        ----------------------- */
	lu = iw[391];
/* U(lenU), dense Hessian factor */
	lenu = iw[392];

	s8fmhx_(nnh, &x[1], &ux[1], &hx[1], &lenu, &rw[lu]);
    }
    *xhx = ddot_(nnh, &ux[1], &c__1, &ux[1], &c__1);
    return 0;
} /* s8xhx_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8xHx */
/* Subroutine */ int s8hupdate_(integer *updatetype, integer *nnh, doublereal 
	*ud0, doublereal *u0scale, doublereal *rydx, doublereal *rdxhdx, 
	doublereal *hd, doublereal *dx, doublereal *hdx, doublereal *y, 
	integer *iw, integer *leniw, doublereal *rw, integer *lenrw)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    extern /* Subroutine */ int s8reseth_(integer *, doublereal *, doublereal 
	    *, integer *, integer *, doublereal *, integer *);
    static integer i__, ls, lu, lv;
    static doublereal yi;
    extern /* Subroutine */ int s8fmupdate_(integer *, integer *, doublereal *
	    , doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *), s8lmupdate_(integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    static integer lu0;
    static doublereal hdxi;
    static integer lenu, ludx;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal h0scale;
    static integer mqnmods, lvlhess;

/*     ================================================================== */
/*     s8Hupdate  applies the pair of vectors that define the BFGS update */
/*     or self-scaled BFGS update. */

/*     On entry: */
/*     --------- */
/*     Hdx   contains  H  times the difference x1 - x. */
/*     y     contains the gradient  difference g1 - g. */

/*     On exit: */
/*     --------- */
/*     Hdx   contains  H  times the difference x1 - x(new). */
/*     y     contains the gradient  difference g1 - g(new). */

/*     s8Hupdate calls one of the Hessian routines s8LMH, s8FMH, s8SDH, ... */
/*     according to the value of the option lvlHess. */
/*     At the moment the options are: */

/*        lvlHess = LM      Limited-Memory (LM) BFGS  (the default). */
/*        lvlHess = FM      Full-Memory    (FM) BFGS */
/*        lvlHess = Exact   FD or exact Hessian */

/*     19 Jul 1995: First version (s8Hupd) of s8Hupdate. */
/*     12 Jan 1996: Full-memory Hessian option added. */
/*     18 Feb 2001: LM H stored in product form. */
/*     12 Jan 2005: FM H stored in product form. */
/*     16 Jan 2005: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/* # of updates since last reset */
/*     ------------------------------------------------------------------ */
/* Tags(7): Approx Hessian type */
    /* Parameter adjustments */
    --y;
    --hdx;
    --dx;
    --hd;
    --iw;
    --rw;

    /* Function Body */
    mqnmods = iw[54];
/* (ge 0) max # of BFGS updates */
    lvlhess = iw[72];
/* LM, FM or Exact Hessian */
    if (*updatetype == 1) {
	*u0scale = *rydx / *rdxhdx;
	h0scale = *u0scale * *u0scale;
	*rdxhdx *= *u0scale;
/* Used later for the LC update. */
	dscal_(nnh, &h0scale, &hd[1], &c__1);
	dscal_(nnh, &h0scale, &hdx[1], &c__1);
    }
/*     Include the latest update in the Hessian diagonal. */
    i__1 = *nnh;
    for (i__ = 1; i__ <= i__1; ++i__) {
	hdxi = hdx[i__] / *rdxhdx;
	yi = y[i__] / *rydx;
/* Computing 2nd power */
	d__1 = hdxi;
/* Computing 2nd power */
	d__2 = yi;
	hd[i__] = hd[i__] - d__1 * d__1 + d__2 * d__2;
    }
    if (iw[381] >= mqnmods) {
/*        --------------------------------------------------------------- */
/*        Too many modifications. */
/*        Reset U to be the root of the diagonal of the current H. */
/*        --------------------------------------------------------------- */
	s8reseth_(nnh, ud0, &hd[1], &iw[1], leniw, &rw[1], lenrw);
    } else {
/*        --------------------------------------------------------------- */
/*        Apply the update to H = U'U. */
/*        The form of U depends on the selected implementation type. */
/*        --------------------------------------------------------------- */
	++iw[381];
	iw[243] = 0;
	if (lvlhess == 0) {
	    lu0 = iw[346];
/* Square root of initial BFGS diagonal */
	    ls = iw[401];
/* sk's for BFGS products: (I + sk*vk') */
	    lv = iw[402];
/* vk's for BFGS products: (I + sk*vk') */
	    s8lmupdate_(updatetype, &iw[381], &mqnmods, nnh, u0scale, rydx, 
		    rdxhdx, &hdx[1], &y[1], &dx[1], &rw[lu0], &rw[ls], &rw[lv]
		    );
	} else if (lvlhess == 1) {
	    ludx = iw[345];
/* Ux(nnH)      = product of U with x */
	    lu = iw[391];
/* U(lenU), full-memory BFGS Hessian H = */
	    lenu = iw[392];

	    s8fmupdate_(updatetype, nnh, u0scale, rydx, rdxhdx, &hdx[1], &y[1]
		    , &rw[ludx], &lenu, &rw[lu]);
	}
    }
    return 0;
} /* s8hupdate_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8HUpdate */
/* Subroutine */ int s8hmoda_(integer *iexit, S_fp funwrapper, U_fp funcon, 
	U_fp funobj, U_fp userhv, logical *usefd, integer *n, integer *nb, 
	integer *nncon0, integer *nncon, integer *nnjac, integer *nnobj0, 
	integer *nnobj, integer *nnh, integer *minimize, doublereal *step, 
	doublereal *dxhdx, doublereal *ydx, doublereal *fobj, doublereal *
	fcon, doublereal *gcon, doublereal *gobj, doublereal *gcon1, 
	doublereal *gobj1, integer *nej, integer *nlocj, integer *locj, 
	integer *indj, integer *neh, integer *nloch, integer *loch, integer *
	indh, doublereal *hcol, integer *negcon, integer *nlocg, integer *
	locg, doublereal *bl, doublereal *bu, doublereal *dx, doublereal *dg, 
	doublereal *udx, doublereal *hdx, doublereal *ycon1, doublereal *tdx, 
	doublereal *x, doublereal *xqp0, doublereal *y, char *cu, integer *
	lencu, integer *iu, integer *leniu, doublereal *ru, integer *lenru, 
	char *cw, integer *lencw, integer *iw, integer *leniw, doublereal *rw,
	 integer *lenrw, ftnlen cu_len, ftnlen cw_len)
{
    /* System generated locals */
    doublereal d__1;

    /* Local variables */
    static doublereal eps, eps0;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern /* Subroutine */ int s6getmissing_(integer *, integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     S_fp, U_fp, U_fp, U_fp, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, char *, integer *, integer *, integer *, doublereal 
	    *, integer *, char *, integer *, integer *, integer *, doublereal 
	    *, integer *, ftnlen, ftnlen);
    static logical nonlinearobj, nonlinearcon;
    extern /* Subroutine */ int s8xhx_(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *), dload_(integer *, doublereal *, doublereal *, integer 
	    *), dscal_(integer *, doublereal *, doublereal *, integer *), 
	    dcopy_(integer *, doublereal *, integer *, doublereal *, integer *
	    ), daxpy_(integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer nzero, modefg;
    extern /* Subroutine */ int s8gprod_(integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, integer *);
    static doublereal signobj;

/*     ================================================================== */
/*     s8HmodA   redefines the quantities x and dx, Hdx and dg  used for the */
/*     quasi-Newton update.  The problem functions are recomputed at x. */

/*     The new  x1  is  x1 + step*(xQP0 - x1),  where xQP0 is a */
/*     (nonelastic) feasible point from the QP subproblem. */

/*     s8HmodA is always called with nnH > 0. */

/*     02 Dec 1994: First version (s8x1) of s8HmodA. */
/*     20 Jul 1998: s8HmodA made self-contained */
/*     24 Aug 1998: Fixed bug found by Alan Brown at Nag. */
/*                  FD derivatives now computed correctly. */
/*                  Parameter useFD added. */
/*     11 Oct 1998: Facility to combine funobj and funcon added. */
/*     12 Oct 2003: snEXIT and snPRNT adopted. */
/*     14 Jan 2005: Argument Udx added for call to s8xHx. */
/*     16 Jun 2008: Call-status implemented correctly */
/*     19 Oct 2014: Added user-defined Hessian. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --x;
    --bu;
    --bl;
    --y;
    --ycon1;
    --fcon;
    --gobj1;
    --gobj;
    --xqp0;
    --tdx;
    --hdx;
    --udx;
    --dg;
    --dx;
    --indj;
    --locj;
    --hcol;
    --indh;
    --loch;
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
    eps = rw[1];
/* machine precision.  IEEE DP  2.22e-16 */
    eps0 = rw[2];
/* eps**(4/5) */
    *iexit = 0;
    modefg = 2;
    signobj = (doublereal) (*minimize);
    nonlinearcon = *nncon > 0;
    nonlinearobj = *nnobj > 0;
/*     Save dx in case a better update pair cannot be found. */
    dcopy_(nnh, &dx[1], &c__1, &tdx[1], &c__1);
/*     ------------------------------------------------- */
/*     dx =  dx - step*y,  with  y = xQP0 - x */
/*     ------------------------------------------------- */
    daxpy_(nnh, &c_b13, &x[1], &c__1, &xqp0[1], &c__1);
    d__1 = -(*step);
    daxpy_(nnh, &d__1, &xqp0[1], &c__1, &dx[1], &c__1);
/*     ------------------------------------------------- */
/*     Compute the minimum curvature. */
/*     If nnH < n, dxHdx may be zero (or negative fuzz). */
/*     ------------------------------------------------- */
    s8xhx_(nnh, &dx[1], &udx[1], &hdx[1], dxhdx, &iw[1], leniw, &rw[1], lenrw)
	    ;
    if (*dxhdx >= eps) {
/*        ----------------------------------------------- */
/*        Redefine  x  as   x + step*y  (y held in xQP0.) */
/*        Evaluate the functions at the new x. */
/*        ----------------------------------------------- */
	daxpy_(nnh, step, &xqp0[1], &c__1, &x[1], &c__1);
	(*funwrapper)(iexit, &modefg, &nonlinearcon, &nonlinearobj, n, negcon,
		 nncon0, nncon, nnjac, nnh, nnobj0, nnobj, (U_fp)funcon, (
		U_fp)funobj, (U_fp)userhv, &x[1], &ycon1[1], nej, nlocj, &
		locj[1], &indj[1], neh, nloch, &loch[1], &indh[1], &hcol[1], &
		fcon[1], fobj, &gcon[1], &gobj[1], cu + 8, lencu, &iu[1], 
		leniu, &ru[1], lenru, cw + 8, lencw, &iw[1], leniw, &rw[1], 
		lenrw, (ftnlen)8, (ftnlen)8);
	if (*iexit > 0) {
	    return 0;
	}
	if (*iexit == 0 && *usefd) {
	    s6getmissing_(iexit, n, negcon, nncon0, nncon, nnjac, nnh, nnobj0,
		     nnobj, (S_fp)funwrapper, (U_fp)funcon, (U_fp)funobj, (
		    U_fp)userhv, &bl[1], &bu[1], &x[1], &ycon1[1], nej, nlocj,
		     &locj[1], &indj[1], neh, nloch, &loch[1], &indh[1], &
		    hcol[1], &fcon[1], fobj, &gcon[1], &gobj[1], &y[1], cu + 
		    8, lencu, &iu[1], leniu, &ru[1], lenru, cw + 8, lencw, &
		    iw[1], leniw, &rw[1], lenrw, (ftnlen)8, (ftnlen)8);
	    if (*iexit > 0) {
		return 0;
	    }
	}
	if (*iexit == 0) {
/*           ------------------------------------------------------------ */
/*           The functions have been computed at x. */
/*           ------------------------------------------------------------ */
	    if (*nnobj > 0) {
		dcopy_(nnobj, &gobj1[1], &c__1, &dg[1], &c__1);
		daxpy_(nnobj, &c_b13, &gobj[1], &c__1, &dg[1], &c__1);
		if (*minimize < 0) {
		    dscal_(nnobj, &signobj, &dg[1], &c__1);
		}
	    }
	    nzero = *nnh - *nnobj;
	    if (nzero > 0) {
		dload_(&nzero, &c_b2, &dg[*nnobj + 1], &c__1);
	    }
	    if (*nncon > 0) {
		s8gprod_(&c__1, &eps0, nej, nlocj, &locj[1], &indj[1], negcon,
			 nlocg, &locg[1], &gcon1[1], &c_b13, &ycon1[1], nncon,
			 &c_b25, &dg[1], nnjac);
		s8gprod_(&c__1, &eps0, nej, nlocj, &locj[1], &indj[1], negcon,
			 nlocg, &locg[1], &gcon[1], &c_b25, &ycon1[1], nncon, 
			&c_b25, &dg[1], nnjac);
	    }
	    *ydx = ddot_(nnh, &dg[1], &c__1, &dx[1], &c__1);
	}
    }
    if (*dxhdx < eps || *iexit != 0) {
	dcopy_(nnh, &tdx[1], &c__1, &dx[1], &c__1);
	s8xhx_(nnh, &dx[1], &udx[1], &hdx[1], dxhdx, &iw[1], leniw, &rw[1], 
		lenrw);
    }
    return 0;
} /* s8hmoda_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8HmodA */
/* Subroutine */ int s8hmodb_(integer *nncon, integer *nnjac, doublereal *
	tolz, integer *nej, integer *nlocj, integer *locj, integer *indj, 
	integer *negcon, integer *nlocg, integer *locg, doublereal *ydx, 
	doublereal *ydxmin, doublereal *penunm, doublereal *fcon, doublereal *
	fcon1, doublereal *gcon, doublereal *gcon1, doublereal *dx, 
	doublereal *gd, doublereal *penu, doublereal *v, doublereal *w)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static logical overflow;
    static integer i__;
    static doublereal wi, diff, beta;
    extern doublereal ddiv_(doublereal *, doublereal *, logical *);
    static doublereal peni, wmax;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), daxpy_(integer *, doublereal *, doublereal *, integer 
	    *, doublereal *, integer *);
    static doublereal wnorm;
    static logical gotpen;
    extern /* Subroutine */ int s8gprod_(integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, integer *);

/*     ================================================================== */
/*     s8HmodB  attempts to find the a vector xPen  of minimum two-norm */
/*     such that there exists a BFGS update for the modified Lagrangian */
/*       La   = f(x) - lambda'(fCon1 - LfCon) */
/*                   + 1/2 (fCon1 - LfCon)'*diag(PenU)*(fCon1 - LfCon), */

/*     where  LfCon = fCon + J(x1)*dx. */

/*     On entry: */
/*     --------- */
/*     dx      is the nonlinear part of the search direction x2 - x1. */
/*     gd      is the Lagrangian gradient difference. */
/*     gCon    is the Jacobian at the old x. */
/*     gCon1   is the Jacobian at the new x. */
/*     ydx     is the approximate curvature of the Lagrangian. */
/*     ydxmin  (ydx < ydxmin) is the smallest acceptable approximate */
/*             curvature. */

/*     On exit, */
/*     --------- */
/*     gd      is the augmented Lagrangian gradient difference. */
/*     PenU    are the penalty parameters. */
/*     ydx     is unchanged unless gotPen is true, in which case */
/*              ydx = ydxmin. */

/*     08 Dec 1991: First version based on  Npsol  routine npupdt. */
/*     26 Oct 2000: Current version of s8HmodB. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/* -->  parameter         (PenMax = 1.0d+5) */
/* -->  parameter         (PenMax = 1.0d+16) */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --w;
    --v;
    --penu;
    --fcon1;
    --fcon;
    --gd;
    --dx;
    --indj;
    --locj;
    --gcon1;
    --gcon;
    --locg;

    /* Function Body */
    overflow = FALSE_;
/*     Try an augmented Lagrangian term to increase ydx. */
    *penunm = 0.;
/*     Compute  v = J1*dx and w = (J2 - J1)*dx = J2*dx - v. */
    s8gprod_(&c__0, tolz, nej, nlocj, &locj[1], &indj[1], negcon, nlocg, &
	    locg[1], &gcon[1], &c_b25, &dx[1], nnjac, &c_b2, &v[1], nncon);
    s8gprod_(&c__0, tolz, nej, nlocj, &locj[1], &indj[1], negcon, nlocg, &
	    locg[1], &gcon1[1], &c_b25, &dx[1], nnjac, &c_b2, &w[1], nncon);
    daxpy_(nncon, &c_b13, &v[1], &c__1, &w[1], &c__1);
/*     Compute the difference between c and its linearization. */
/*     v  =  c - cL = fCon1 - (fCon + J1*s) = fCon1 - fCon - J1*s. */
    daxpy_(nncon, &c_b13, &fcon1[1], &c__1, &v[1], &c__1);
    daxpy_(nncon, &c_b25, &fcon[1], &c__1, &v[1], &c__1);
    dscal_(nncon, &c_b13, &v[1], &c__1);
/*     --------------------------------------------------------- */
/*     Compute the minimum-length vector of penalty parameters */
/*     that makes the approximate curvature equal to  ydxmin. */
/*     --------------------------------------------------------- */
/*     Use w to hold the constraint on PenU. */
/*     Minimize            norm(PenU) */
/*     subject to   ( Sum( w(i)*PenU(i) )  =   const, */
/*                  (           PenU(i)   .ge. 0. */
    wmax = 0.;
    i__1 = *nncon;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wi = w[i__] * v[i__];
	wmax = max(wmax,wi);
	w[i__] = max(0.,wi);
    }
    wnorm = dnrm2_(nncon, &w[1], &c__1);
    diff = *ydxmin - *ydx;
    d__1 = wmax * diff;
/* Computing 2nd power */
    d__3 = wnorm;
    d__2 = d__3 * d__3;
    beta = ddiv_(&d__1, &d__2, &overflow);
    gotpen = ! overflow && wmax > 0. && beta < 1e5;
    if (gotpen) {
/* Computing 2nd power */
	d__1 = wnorm;
	beta = diff / (d__1 * d__1);
	i__1 = *nncon;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    wi = w[i__];
	    peni = beta * wi;
	    v[i__] = peni * v[i__];
	    *ydx += peni * wi;
	    penu[i__] = peni;
	}
	*ydx = max(*ydx,*ydxmin);
	*penunm = dnrm2_(nncon, &penu[1], &c__1);
/*        Update  gd  by the term  (J2' - J1')*v, */
/*        with v = diag(PenU)*(fCon1 - fCon - J1*s) from above. */
	s8gprod_(&c__1, tolz, nej, nlocj, &locj[1], &indj[1], negcon, nlocg, &
		locg[1], &gcon1[1], &c_b25, &v[1], nncon, &c_b25, &gd[1], 
		nnjac);
	s8gprod_(&c__1, tolz, nej, nlocj, &locj[1], &indj[1], negcon, nlocg, &
		locg[1], &gcon[1], &c_b13, &v[1], nncon, &c_b25, &gd[1], 
		nnjac);
    }
/* gotPen */
    return 0;
} /* s8hmodb_ */

