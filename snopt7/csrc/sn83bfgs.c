/* ./src/sn83bfgs.f -- translated by f2c (version 20100827).
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

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*     File  sn83bfgs.f.   Full and limited memory BFGS routines. */

/*     s8FMH0   s8FMupdate   s8FMHx */
/*     s8LMH0   s8LMupdate   s8LMHx */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Subroutine */ int s8fmh0_(integer *nnh, doublereal *hd, integer *lenu, 
	doublereal *u)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer j, l, incr;
    extern /* Subroutine */ int dload_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer nzeros;

/*     ================================================================== */
/*     s8FMH0 zeros the off-diagonal elements of H such that H = U'U. */
/*     The diagonals of U are set to the square roots of the diagonal HD. */

/*     19 Jul 1995: First version of s8FMH0. */
/*     06 Sep 1998: Pre- and post-QP diagonal Hessian scaling added. */
/*     13 Jan 2005: HD always positive semidefinite. */
/*     05 Oct 2014: Reorganized to reflect qnInit in dnopt. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
/*     Zero the off-diagonal elements of U. */
/*     ------------------------------------------------------------ */
    /* Parameter adjustments */
    --hd;
    --u;

    /* Function Body */
    incr = *nnh;
    nzeros = *nnh - 1;
    l = 1;
    i__1 = *nnh;
    for (j = 1; j <= i__1; ++j) {
	u[l] = sqrt(hd[j]);
	if (j < *nnh) {
	    dload_(&nzeros, &c_b2, &u[l + 1], &c__1);
	    l += incr;
	    --incr;
	    --nzeros;
	}
    }
    return 0;
} /* s8fmh0_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8FMH0 */
/* Subroutine */ int s8fmupdate_(integer *updatetype, integer *nnh, 
	doublereal *u0scale, doublereal *rydx, doublereal *rdxhdx, doublereal 
	*hdx, doublereal *y, doublereal *udx, integer *lenu, doublereal *u)
{
    /* System generated locals */
    doublereal d__1;

    /* Local variables */
    static doublereal t;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal told;
    static integer numu;
    static doublereal tolz;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer iexit;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static doublereal ulast;
    extern /* Subroutine */ int s6rmod_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *);
    static integer lastnz;

/*     ================================================================== */
/*     s8FMupdate applies the full-memory BFGS update to H = U'U. */
/*     If defined, the self-scaling BFGS update parameter is saved. */
/*     It is needed to update the reduced Hessian when there are only */
/*     linear constraints. */

/*     19 Jul 1995: First version of s8FMup. */
/*     06 Sep 1998: Pre- and post-QP diagonal Hessian scaling added. */
/*     18 Feb 2001: LM H stored in product form. */
/*     13 Jan 2005: FM H stored in product form. */
/*     15 Jan 2005: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --udx;
    --y;
    --hdx;
    --u;

    /* Function Body */
    told = 0.;
    tolz = 0.;
    if (*updatetype == 1) {
	numu = *nnh * (*nnh + 1) / 2;
	dscal_(&numu, u0scale, &u[1], &c__1);
/* multiplies U by U0scale. */
	dscal_(nnh, u0scale, &udx[1], &c__1);
    }
/*     ------------------------------------------------------------------ */
/*     Overwrite (Udx,y) with the vectors (Us,v)  such that */
/*       Us  = Udx / rdxHdx,    v = (1/rydx) gdif - (1/rdxHdx) Hdx. */

/*     Then, U(new) = U + Us v',  with H = U'U. */

/*     Hdx and v  are saved to update R for LC problems. */
/*     ------------------------------------------------------------------ */
    t = ddot_(nnh, &y[1], &c__1, &hdx[1], &c__1);
    if (t >= 0.) {
	d__1 = 1. / *rydx;
	dscal_(nnh, &d__1, &y[1], &c__1);
    } else {
	d__1 = -1. / *rydx;
	dscal_(nnh, &d__1, &y[1], &c__1);
    }
    d__1 = -1. / *rdxhdx;
    daxpy_(nnh, &d__1, &hdx[1], &c__1, &y[1], &c__1);
    d__1 = 1. / *rdxhdx;
    dscal_(nnh, &d__1, &udx[1], &c__1);
/*     ------------------------------------------------------------------ */
/*     Restore  U + Us y' to triangular form  (overwriting Udx). */
/*     ------------------------------------------------------------------ */
    ulast = 0.;
    lastnz = *nnh;
    s6rmod_(&iexit, nnh, nnh, lenu, &u[1], &udx[1], &y[1], &lastnz, &ulast, &
	    told, &tolz);
    return 0;
} /* s8fmupdate_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8FMupdate */
/* Subroutine */ int s8fmhx_(integer *nnh, doublereal *x, doublereal *ux, 
	doublereal *hx, integer *lenu, doublereal *u)
{
    extern /* Subroutine */ int s6rprod_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *);

/*     ================================================================== */
/*     s8FMHx  computes the product Hx = U'Ux, where  U is an upper- */
/*     triangular matrix stored by rows in the one-dimensional array  U. */
/*     lenU defines the length of U.  lenU must be at least */
/*     nnH*(nnH + 1)/2. */

/*     12 Jan 1996: First version of s8FMHx */
/*     12 Jan 2005: H held as U'U. */
/*     15 Jan 2005: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --hx;
    --ux;
    --x;
    --u;

    /* Function Body */
    s6rprod_(&c__0, nnh, nnh, lenu, &u[1], &x[1], &ux[1]);
    s6rprod_(&c__1, nnh, nnh, lenu, &u[1], &ux[1], &hx[1]);
    return 0;
} /* s8fmhx_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8FMHx */
/* Subroutine */ int s8lmh0_(integer *nnh, doublereal *hd, doublereal *u0)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__;

/*     ================================================================== */
/*     s8LMH0  zeros the off-diagonal elements of H such that H = U'U. */
/*     The diagonals of U are set to the square roots of the diagonal HD. */

/*     19 Jul 1995: First version of s8LMH0. */
/*     06 Sep 1998: Pre- and post-QP diagonal Hessian scaling added. */
/*     18 Feb 2001: H stored in product form. */
/*     13 Jan 2005: Hd always positive semidefinite. */
/*     05 Oct 2014: Reorganized to reflect qnInit in dnopt. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --u0;
    --hd;

    /* Function Body */
    i__1 = *nnh;
    for (i__ = 1; i__ <= i__1; ++i__) {
	u0[i__] = sqrt(hd[i__]);
    }
    return 0;
} /* s8lmh0_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8LMH0 */
/* Subroutine */ int s8lmupdate_(integer *updatetype, integer *qnmods, 
	integer *mqnmods, integer *nnh, doublereal *u0scale, doublereal *rydx,
	 doublereal *rdxhdx, doublereal *hdx, doublereal *y, doublereal *dx, 
	doublereal *u0, doublereal *s, doublereal *v)
{
    /* System generated locals */
    integer s_dim1, s_offset, v_dim1, v_offset;
    doublereal d__1;

    /* Local variables */
    static doublereal t;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dcopy_(integer *, doublereal *, integer *, doublereal 
	    *, integer *), daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);

/*     ================================================================== */
/*     s8LMupdate computes the limited-memory BFGS update. */

/*     If defined, the self-scaling BFGS parameter U0scalee is saved. */
/*     It is needed to update the reduced Hessian when there are only */
/*     linear constraints. */

/*     19 Jul 1995: First version of s8LMup. */
/*     06 Sep 1998: Pre- and post-QP diagonal Hessian scaling added. */
/*     18 Feb 2001: LM H stored in product form. */
/*     13 Jan 2005: FM H stored in product form. */
/*     16 Jan 2005: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    v_dim1 = *nnh;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    s_dim1 = *nnh;
    s_offset = 1 + s_dim1;
    s -= s_offset;
    --u0;
    --dx;
    --y;
    --hdx;

    /* Function Body */
    if (*updatetype == 1) {
	dscal_(nnh, u0scale, &u0[1], &c__1);
/* multiplies U0 by U0scale. */
    }
/*     ------------------------------------------------------------------ */
/*     Space remains. Store s and v, where */
/*     U(new) = U(I + sv'), with H = U'U. */
/*       S  =  dx / rdxHdx,    V =  (1/rydx) gdif - (1/rdxHdx) Hdx. */

/*     Hdx and the modified y (= v) are used to update the reduced */
/*     Hessian for LC problems. */
/*     ------------------------------------------------------------------ */
    dcopy_(nnh, &dx[1], &c__1, &s[*qnmods * s_dim1 + 1], &c__1);
    d__1 = 1. / *rdxhdx;
    dscal_(nnh, &d__1, &s[*qnmods * s_dim1 + 1], &c__1);
    t = ddot_(nnh, &y[1], &c__1, &hdx[1], &c__1);
    if (t >= 0.) {
	d__1 = 1. / *rydx;
	dscal_(nnh, &d__1, &y[1], &c__1);
    } else {
	d__1 = -1. / *rydx;
	dscal_(nnh, &d__1, &y[1], &c__1);
    }
    d__1 = -1. / *rdxhdx;
    daxpy_(nnh, &d__1, &hdx[1], &c__1, &y[1], &c__1);
    dcopy_(nnh, &y[1], &c__1, &v[*qnmods * v_dim1 + 1], &c__1);
    return 0;
} /* s8lmupdate_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8LMupdate */
/* Subroutine */ int s8lmhx_(integer *nnh, doublereal *x, doublereal *ux, 
	doublereal *hx, integer *mqnmods, integer *qnmods, doublereal *u0, 
	doublereal *s, doublereal *v)
{
    /* System generated locals */
    integer s_dim1, s_offset, v_dim1, v_offset, i__1;

    /* Local variables */
    static doublereal c__;
    static integer k;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern /* Subroutine */ int ddscl_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dcopy_(integer *, doublereal *, integer 
	    *, doublereal *, integer *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);

/*     ================================================================== */
/*     s8LMHx forms the product  Hx  for the limited-memory s8Hx. */
/*     H = U'U, where U = U0*(I + s1*v1')*(I + s2*v2')...(I + sk*vk'). */
/*     with  k = QNmods */

/*     19 Jul 1995: First version of s8LMHx */
/*     18 Feb 2001: H stored in product form. */
/*     12 Jan 2005: Ux added as argument. */
/*     16 Jan 2005: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --u0;
    --hx;
    --ux;
    --x;
    v_dim1 = *nnh;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    s_dim1 = *nnh;
    s_offset = 1 + s_dim1;
    s -= s_offset;

    /* Function Body */
    dcopy_(nnh, &x[1], &c__1, &ux[1], &c__1);
/*     Multiply by U. */
    for (k = *qnmods; k >= 1; --k) {
	c__ = ddot_(nnh, &v[k * v_dim1 + 1], &c__1, &ux[1], &c__1);
	daxpy_(nnh, &c__, &s[k * s_dim1 + 1], &c__1, &ux[1], &c__1);
    }
    ddscl_(nnh, &u0[1], &c__1, &ux[1], &c__1);
/*     Multiply by U'. */
    dcopy_(nnh, &ux[1], &c__1, &hx[1], &c__1);
    ddscl_(nnh, &u0[1], &c__1, &hx[1], &c__1);
    i__1 = *qnmods;
    for (k = 1; k <= i__1; ++k) {
	c__ = ddot_(nnh, &s[k * s_dim1 + 1], &c__1, &hx[1], &c__1);
	daxpy_(nnh, &c__, &v[k * v_dim1 + 1], &c__1, &hx[1], &c__1);
    }
    return 0;
} /* s8lmhx_ */

