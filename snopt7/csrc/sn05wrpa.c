/* ./src/sn05wrpa.f -- translated by f2c (version 20100827).
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

static integer c__13 = 13;
static integer c__1 = 1;
static integer c__3 = 3;
static doublereal c_b26 = 0.;
static integer c__4 = 4;
static integer c__0 = 0;
static integer c_n4 = -4;
static doublereal c_b34 = -1.;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*     File  sn05wrpa.f  --- user-function wrapper for  snOptA. */

/*     s0fgA    s0fgA1 */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Subroutine */ int s0fga_(integer *iexit, integer *modefg, logical *
	nonlinearcon, logical *nonlinearobj, integer *n, integer *negcon, 
	integer *nncon0, integer *nncon, integer *nnjac, integer *nnh, 
	integer *nnobj0, integer *nnobj, U_fp userfg, U_fp dummyf, U_fp 
	dummyh, doublereal *x, doublereal *ycon, integer *nej, integer *nlocj,
	 integer *locj, integer *indj, integer *neh, integer *nloch, integer *
	loch, integer *indh, doublereal *hcol, doublereal *fcon, doublereal *
	fobj, doublereal *gcon, doublereal *gobj, char *cu, integer *lencu, 
	integer *iu, integer *leniu, doublereal *ru, integer *lenru, char *cw,
	 integer *lencw, integer *iw, integer *leniw, doublereal *rw, integer 
	*lenrw, ftnlen cu_len, ftnlen cw_len)
{
    /* Format strings */
    static char fmt_1100[] = "(\002 The user has defined\002,i8,\002   out"
	    " of\002,i8,\002   first  derivatives\002)";
    static char fmt_2100[] = "(\002 XXX  Some first  derivatives are missing"
	    " ---\002,\002 derivative level reduced to 0\002)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static integer lxscaled, lvlscale, k, kg, lf, lg, nf, neg, nkx, lxn;
    static char str[80];
    static integer leng, mode, lkxn;
    extern /* Subroutine */ int s8callstatus_(integer *, integer *, integer *)
	    , s0fga1_(integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, U_fp, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, ftnlen, ftnlen), dload_(
	    integer *, doublereal *, doublereal *, integer *), ddscl_(integer 
	    *, doublereal *, integer *, doublereal *, integer *), dddiv_(
	    integer *, doublereal *, integer *, doublereal *, integer *);
    static integer lindg, llocg, nlocg, lfmul;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer lycon;
    static logical scaled;
    static integer lgobju, lnglin, ligfun, ljgvar, lgconu;
    static doublereal gdummy;
    static integer knowng, status;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);
    static integer lscales;
    extern /* Subroutine */ int s8scaleg_(integer *, doublereal *, doublereal 
	    *, doublereal *, integer *), s8scalej_(integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, integer *, integer 
	    *, integer *, doublereal *, doublereal *, integer *);

    /* Fortran I/O blocks */
    static icilist io___30 = { 0, str, 0, fmt_1100, 80, 1 };
    static icilist io___31 = { 0, str, 0, fmt_2100, 80, 1 };


/*     ================================================================== */
/*     s0fgA  is an instance of funwrapper that calls the user-written */
/*     routine  userfg  to evaluate the problem functions and possibly */
/*     their gradients. */

/*     Subroutine  userfg  is called using modefg to control */
/*     the gradients as follows: */

/*     modefg        Task */
/*     ------        ---- */
/*       2     Assign fCon, fObj and all known elements of gCon and gObj. */
/*       1     Assign all known elements of gCon and gObj. */
/*             (fObj and fCon are ignored). */
/*       0     Assign fObj, fCon.  (gCon and gObj are ignored). */

/*     Since objective and constraints are computed simultaneously, */
/*     the input variables  NonlinearCon  and  NonlinearObj are ignored. */

/*     31 Oct 1998: First version based on s0fg in SNOPT 5.3-4. */
/*     03 Aug 2003: snEXIT and snPRNT adopted. */
/*     16 Jun 2008: Call-status implemented correctly. */
/*     15 Nov 2010: Call-status removed from the argument list. */
/*     07 Nov 2014: yCon added as argument. */
/*     13 Nov 2014: Scales applied for lvlScale > 0. */
/*     14 Feb 2015: proxWeight added. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/* = 0,1,2 or 3, deriv. level */
/* > 0 => some differences needed */
/* > 0 => some exact derivs */
/* calls to fCon: mode = 0 */
/* calls to fCon  mode > 0 */
/* calls to fObj: mode = 0 */
/*     ------------------------------------------------------------------ */
/* calls to fObj: mode > 0 */
    /* Parameter adjustments */
    --x;
    --gcon;
    --fcon;
    --ycon;
    --gobj;
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
    lvlscale = iw[75];
/* scale option */
    nf = iw[248];
/* # of components of user-defined F */
    neg = iw[249];
/* # of components of user-defined G */
    lkxn = iw[252];
/* jN = kxN(j ) => col j of Jcol is variable jN */
    llocg = iw[260];
/* locG(nlocG) = column pointers for indG */
    lindg = iw[261];
/* indG(neG) holds the row indices for gij */
    lnglin = iw[262];
/* nGlin(j) = # linear elems in col j of gCon */
    ligfun = iw[266];
/* iGfun(neG) row list of reordered G nonzeros */
    ljgvar = iw[267];
/* iGvar(neG) col list of reordered G nonzeros */
    lscales = iw[296];
/* scales(nb)  = row and column scales */
    lxscaled = iw[302];
/* xScaled(n)  = copy of scaled  x */
    lgconu = iw[319];
/* record of unknown derivatives and constants */
    lgobju = iw[323];
/* record of unknown derivatives */
    lxn = iw[327];
/* xN(n)       = variables in natural order */
    lf = iw[328];
/* F(nF)       = user-defined F */
    lfmul = iw[329];
/* Fmul(nF)    = user-defined multipliers */
    lg = iw[330];
/* G (lenG)    = problem derivatives */
    lycon = iw[348];
/* yCon(nnCon) = multipliers for F */
    gdummy = rw[69];
/* definition of an 'unset' value */
    *iexit = 0;
    scaled = lvlscale > 0;
    nlocg = *nnjac + 1;
    leng = max(neg,1);
    nkx = *n + nf;
/* Determine the user-function call-status. */
    s8callstatus_(&status, &iw[1], leniw);
    if (status == 1) {
/* -------------------------------------------------------------- */
/* First evaluation of the problem functions in snOptA */
/* On entry, lvlScale = 0. */
/* -------------------------------------------------------------- */
	iw[183] = 0;
/* > 0 => some differences needed */
	iw[184] = 0;
/* > 0 => some exact derivatives provided */
	snprnt_(&c__13, " ", &iw[1], leniw, (ftnlen)1);
	dload_(&neg, &gdummy, &rw[lg], &c__1);
    }
/* ----------------------------------------------------------------- */
/* Unscale x. */
/* ----------------------------------------------------------------- */
    if (scaled) {
	dcopy_(n, &x[1], &c__1, &rw[lxscaled], &c__1);
	ddscl_(n, &rw[lscales], &c__1, &x[1], &c__1);
    }
/* ----------------------------------------------------------------- */
/* Compute the user-defined functions and derivatives. */
/* ----------------------------------------------------------------- */
    mode = status;
    s0fga1_(modefg, &mode, &nf, n, &nkx, nncon0, nncon, nnjac, nnh, nnobj0, 
	    nnobj, &iw[lkxn], (U_fp)userfg, &rw[lxn], &x[1], &rw[lf], &rw[
	    lfmul], &fcon[1], &rw[lycon], fobj, &gobj[1], negcon, &nlocg, &iw[
	    llocg], &iw[lindg], &gcon[1], &rw[lgconu], &iw[lnglin], &iw[
	    ligfun], &iw[ljgvar], &leng, &neg, &rw[lg], cu + 8, lencu, &iu[1],
	     leniu, &ru[1], lenru, cw + 8, lencw, &iw[1], leniw, &rw[1], 
	    lenrw, (ftnlen)8, (ftnlen)8);
    ++iw[189];
    ++iw[194];
    if (*modefg > 0) {
	++iw[190];
	++iw[195];
    }
/* ----------------------------------------------------------------- */
/* Scale  x and the derivatives. */
/* ----------------------------------------------------------------- */
    if (scaled) {
	dcopy_(n, &rw[lxscaled], &c__1, &x[1], &c__1);
	dddiv_(nncon, &rw[lscales + *n], &c__1, &fcon[1], &c__1);
	if (*modefg > 0 && iw[184] > 0) {
	    s8scaleg_(nnobj, &rw[lscales], &gobj[1], &rw[1], lenrw);
	    s8scalej_(nncon, nnjac, negcon, n, &rw[lscales], nej, nlocj, &
		    locj[1], &indj[1], &gcon[1], &rw[1], lenrw);
	}
    }
    if (mode < 0) {
/* -------------------------------------------------------------- */
/* The user may be saying the function is undefined (mode = -1) */
/* or may just want to stop                         (mode < -1). */
/* -------------------------------------------------------------- */
	if (mode == -1) {
	    *iexit = -1;
	} else {
	    *iexit = 71;
	}
    }
/* ================================================================= */
/* Do some housekeeping after the first call. */
/* ================================================================= */
    if (status == 1 && *iexit == 0) {
/* -------------------------------------------------------------- */
/* Count how many Jacobian elements are provided. */
/* -------------------------------------------------------------- */
	knowng = 0;
	kg = lg;
	i__1 = neg;
	for (k = 1; k <= i__1; ++k) {
	    if (rw[kg] != gdummy) {
		++knowng;
	    }
	    ++kg;
	}
	s_wsfi(&io___30);
	do_fio(&c__1, (char *)&knowng, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&neg, (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)80);
	if (knowng < neg) {
	    iw[183] = 1;
	}
	if (knowng > 0) {
	    iw[184] = 1;
	}
	if (knowng < neg) {
/* ----------------------------------------------------------- */
/* Missing derivatives. */
/* Keep a record of them in gObjU and gConU. */
/* Reduce the derivative level if necessary. */
/* ----------------------------------------------------------- */
	    dcopy_(nnobj, &gobj[1], &c__1, &rw[lgobju], &c__1);
	    dcopy_(negcon, &gcon[1], &c__1, &rw[lgconu], &c__1);
	    if (iw[70] == 3) {
		iw[70] = 0;
		s_wsfi(&io___31);
		e_wsfi();
		snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)80);
	    }
	}
    }
    return 0;
} /* s0fga_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s0fgA */
/* Subroutine */ int s0fga1_(integer *modefg, integer *status, integer *nf, 
	integer *n, integer *nkx, integer *nncon0, integer *nncon, integer *
	nnjac, integer *nnh, integer *nnobj0, integer *nnobj, integer *kxn, 
	S_fp userfg, doublereal *xn, doublereal *x, doublereal *f, doublereal 
	*fmul, doublereal *fcon, doublereal *ycon, doublereal *fobj, 
	doublereal *gobj, integer *negcon, integer *nlocg, integer *locg, 
	integer *indg, doublereal *gcon, doublereal *gconu, integer *nglin, 
	integer *igfun, integer *jgvar, integer *leng, integer *neg, 
	doublereal *g, char *cu, integer *lencu, integer *iu, integer *leniu, 
	doublereal *ru, integer *lenru, char *cw, integer *lencw, integer *iw,
	 integer *leniw, doublereal *rw, integer *lenrw, ftnlen cu_len, 
	ftnlen cw_len)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer minimize, i__, j, k, in, jn, lx0;
    static doublereal proxweight;
    static integer iobj;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer needf, needg, needh;
    extern /* Subroutine */ int dload_(integer *, doublereal *, doublereal *, 
	    integer *), dscal_(integer *, doublereal *, doublereal *, integer 
	    *), dcopy_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer nextg;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *), s1time_(integer *, integer *,
	     integer *, integer *, doublereal *, integer *);
    static integer objrow, lvltim;
    static doublereal signobj;

/*     ================================================================== */
/*     s0fgA1   does the work for s0fgA. */

/*     08 Nov 1998: First version of s0fgA1. */
/*     09 Apr 1999: Updated for SnoptA. */
/*     01 Apr 2005: Cosmetic changes. */
/*     14 Feb 2015: proxWeight added. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --fmul;
    --f;
    --x;
    --xn;
    --kxn;
    --ycon;
    --fcon;
    --gobj;
    --gconu;
    --gcon;
    --indg;
    --nglin;
    --locg;
    --g;
    --jgvar;
    --igfun;
    cu -= 8;
    --iu;
    --ru;
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    proxweight = rw[91];
/* Proximal-point weight */
    minimize = iw[87];
/* 1, -1  => MIN, MAX */
    lvltim = iw[182];
/* Timing level */
    iobj = iw[204];
/* position of the objective row in J */
    lx0 = iw[298];
/* x0(nnL) = FP base point */
    objrow = iw[103];
/* Objective row of user-defined F */
    signobj = (doublereal) minimize;
    needf = 0;
    needg = 0;
    needh = 0;
    if (*modefg == 0 || *modefg == 2) {
	needf = 1;
    }
    if (*modefg == 1 || *modefg == 2) {
	needg = 1;
    }
    if (*modefg == 4) {
	needh = 1;
    }
/* ----------------------------------------------------------------- */
/* Save x in natural order. */
/* For safety, zero out the linear components of x, just in case */
/* they are used by  userfg. */
/* ----------------------------------------------------------------- */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	jn = kxn[j];
	if (j <= *nnh) {
	    xn[jn] = x[j];
	} else {
	    xn[jn] = 0.;
	}
    }
    if (needh > 0) {
/* -------------------------------------------------------------- */
/* The Hessian of the Lagrangian is requested. */
/* Expand the multipliers into Fmul. */
/* -------------------------------------------------------------- */
	dload_(nf, &c_b26, &fmul[1], &c__1);
	if (objrow > 0) {
	    fmul[objrow] = 1.;
	}
	i__1 = *nncon;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    j = *n + i__;
	    in = kxn[j];
	    fmul[in] = -signobj * ycon[i__];
	}
    }
/* ----------------------------------------------------------------- */
/* Compute the user-defined functions. */
/* ----------------------------------------------------------------- */
    if (lvltim >= 2) {
	s1time_(&c__4, &c__0, &iw[1], leniw, &rw[1], lenrw);
    }
    (*userfg)(status, n, &xn[1], &needf, nf, &f[1], &needg, leng, &g[1], cu + 
	    8, lencu, &iu[1], leniu, &ru[1], lenru, (ftnlen)8);
    if (lvltim >= 2) {
	s1time_(&c_n4, &c__0, &iw[1], leniw, &rw[1], lenrw);
    }
/* ----------------------------------------------------------------- */
/* Gather the nonlinear elements of F. */
/* ----------------------------------------------------------------- */
    if (needf > 0) {
	i__1 = *nncon;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    j = *n + i__;
	    in = kxn[j];
	    fcon[i__] = f[in];
	}
	if (*nnobj == 0) {
	    *fobj = 0.;
	} else if (objrow == 0) {
	    dcopy_(nnobj, &x[1], &c__1, &gobj[1], &c__1);
	    daxpy_(nnobj, &c_b34, &rw[lx0], &c__1, &gobj[1], &c__1);
	    *fobj = proxweight * .5 * ddot_(nnobj, &gobj[1], &c__1, &gobj[1], 
		    &c__1);
	    dscal_(nnobj, &proxweight, &gobj[1], &c__1);
	} else {
	    *fobj = f[objrow];
	}
    }
    if (needg > 0) {
/* Extract the derivatives from G. */
	i__1 = *neg;
	for (k = 1; k <= i__1; ++k) {
	    i__ = igfun[k];
	    j = jgvar[k];
	    if (i__ == iobj) {
		gobj[j] = g[k];
	    } else {
		nextg = locg[j];
		indg[nextg] = i__;
		gcon[nextg] = g[k];
		locg[j] = nextg + 1;
	    }
	}
    } else {
/* Keep locG in synch ready for any constant elements. */
	i__1 = *neg;
	for (k = 1; k <= i__1; ++k) {
	    i__ = igfun[k];
	    j = jgvar[k];
	    if (i__ != iobj) {
		++locg[j];
	    }
	}
    }
    if (needf > 0) {
/* Add to fCon the linear term associated with every constant */
/* element in the nonlinear part of Jcol. */
	i__1 = *nnjac;
	for (j = 1; j <= i__1; ++j) {
	    nextg = locg[j];
	    i__2 = nglin[j];
	    for (k = 1; k <= i__2; ++k) {
		i__ = indg[nextg];
		fcon[i__] += gconu[nextg] * x[j];
		++nextg;
	    }
	}
    }
/* locG(j) points to the first linear element in column j. */
/* Update it so that it points to the start of column j. */
    for (j = *nnjac; j >= 2; --j) {
	locg[j - 1] += nglin[j - 1];
	locg[j] = locg[j - 1];
    }
    locg[1] = 1;
    return 0;
} /* s0fga1_ */

