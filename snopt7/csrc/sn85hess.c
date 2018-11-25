/* ./src/sn85hess.f -- translated by f2c (version 20100827).
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

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*     File  sn85Hess.f */

/*     s8SDH0   s8SDIx */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Subroutine */ int s8sdh0_(integer *hqntype, integer *nnh, doublereal *
	u0pre, doublereal *u0)
{
    extern /* Subroutine */ int dload_(integer *, doublereal *, doublereal *, 
	    integer *);

/*     ================================================================== */
/*     s8SDH0 resets the approximate Hessian H to a diagonal matrix. */
/*     On entry, the value of HQNType is as follows: */

/*       HQNType */
/*       ------- */
/*       HUnset (-1)      H not set. */
/*       HNorml ( 0)      H is a Hessian of the form defined by  lvlHes. */
/*       HDiag  ( 1)      H is a diagonal matrix. */
/*       HUnit  ( 2)      H is an identity matrix. */

/*     19 Jul 1995: First version of s8SDH0. */
/*     06 Sep 1998: Pre- and post-QP diagonal Hessian scaling added. */
/*     25 Mar 2005: Resurrected for snopt8. */
/*     25 Mar 2005: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     Set H0 to a multiple of the identity. */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --u0;

    /* Function Body */
    *hqntype = 2;
    dload_(nnh, u0pre, &u0[1], &c__1);
    return 0;
} /* s8sdh0_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8SDH0 */
/* Subroutine */ int s8sdix_(integer *nnh, doublereal *u0, doublereal *x, 
	doublereal *hx)
{
    extern /* Subroutine */ int ddscl_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dcopy_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);

/*     ================================================================== */
/*     s8SDIx  multiplies the QP Hessian  H by the vector  x. */
/*     It is used to define Hx for the first QP subproblem. */

/*     19 Jul 1995: First version of s8SDIx. */
/*     13 Apr 2005: Current version. */
/*     ================================================================== */
    /* Parameter adjustments */
    --hx;
    --x;
    --u0;

    /* Function Body */
    dcopy_(nnh, &x[1], &c__1, &hx[1], &c__1);
    ddscl_(nnh, &u0[1], &c__1, &hx[1], &c__1);
    ddscl_(nnh, &u0[1], &c__1, &hx[1], &c__1);
    return 0;
} /* s8sdix_ */

