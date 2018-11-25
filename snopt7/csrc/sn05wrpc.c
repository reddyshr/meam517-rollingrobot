/* ./src/sn05wrpc.f -- translated by f2c (version 20100827).
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
static integer c__4 = 4;
static integer c__0 = 0;
static doublereal c_b16 = -1.;
static integer c_n4 = -4;
static integer c__3 = 3;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*     File  sn05wrpc.f  --- user-function interfaces for snOptC. */

/*     s0fgC */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Subroutine */ int s0fgc_(integer *iexit, integer *modefg, logical *
	nonlinearcon, logical *nonlinearobj, integer *n, integer *negcon, 
	integer *nncon0, integer *nncon, integer *nnjac, integer *nnh, 
	integer *nnobj0, integer *nnobj, S_fp userfun, U_fp dummyf, U_fp 
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
	    " of\002,i8,\002   constraint gradients.\002)";
    static char fmt_2010[] = "(\002 SnOptC  will define \002,i8,\002   gradi"
	    "ents for the \002,\002 FP objective.\002)";
    static char fmt_2000[] = "(\002 The user has defined\002,i8,\002   out"
	    " of\002,i8,\002   objective  gradients.\002)";
    static char fmt_2100[] = "(\002 XXX  Some objective  derivatives are mis"
	    "sing ---\002,\002 derivative level reduced to\002,i3)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static integer lxscaled, lvlscale, supplied, l, lg, lx0;
    static doublereal proxweight;
    static integer statususer;
    static char str[80];
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer mode;
    extern /* Subroutine */ int s8callstatus_(integer *, integer *, integer *)
	    , dload_(integer *, doublereal *, doublereal *, integer *), 
	    dscal_(integer *, doublereal *, doublereal *, integer *), ddscl_(
	    integer *, doublereal *, integer *, doublereal *, integer *), 
	    dddiv_(integer *, doublereal *, integer *, doublereal *, integer *
	    ), dcopy_(integer *, doublereal *, integer *, doublereal *, 
	    integer *), daxpy_(integer *, doublereal *, doublereal *, integer 
	    *, doublereal *, integer *), s1time_(integer *, integer *, 
	    integer *, integer *, doublereal *, integer *);
    static logical scaled;
    static integer lgsave, lgobju, nnglin, lgconu, minmax, nnobju;
    static doublereal gdummy;
    static logical fponly;
    static integer lvltim, status;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);
    static integer lscales;
    extern /* Subroutine */ int s8scaleg_(integer *, doublereal *, doublereal 
	    *, doublereal *, integer *), s8scalej_(integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, integer *, integer 
	    *, integer *, doublereal *, doublereal *, integer *);

    /* Fortran I/O blocks */
    static icilist io___22 = { 0, str, 0, fmt_1100, 80, 1 };
    static icilist io___24 = { 0, str, 0, fmt_2010, 80, 1 };
    static icilist io___25 = { 0, str, 0, fmt_2000, 80, 1 };
    static icilist io___26 = { 0, str, 0, fmt_2100, 80, 1 };


/*     ================================================================== */
/*     s0fgC   is a funwrapper interface that calls the user-written */
/*     routine  userfun  to evaluate the problem functions and possibly */
/*     their gradients. */

/*     Argument  userfun  is called using modefg to control */
/*     the gradients as follows: */

/*     modefg        Task */
/*     ------        ---- */
/*       2     Assign fCon, fObj and all known elements of gCon and gObj. */
/*       1     Assign all known elements of gCon and gObj. */
/*             (fObj and fCon are ignored). */
/*       0     Assign fObj, fCon.  (gCon and gObj are ignored). */

/*     If s0fgC is called with minmax = 0 (feasible point only) then */
/*     nnObj = max(nnJac,nnObj)  and the user objective is not used. */

/*     31 Oct 1998: First version based on snwrap in SNOPT 5.3-4. */
/*     30 Dec 2000: Housekeeping for Status = 1 included. */
/*     03 Aug 2003: snEXIT and snPRNT adopted. */
/*     04 Jan 2005: v, Hv added. */
/*     16 Jun 2008: Call-status implemented correctly. */
/*     15 Nov 2010: Call-status removed from the argument list. */
/*     11 Sep 2014: nnObjU and nnObj separated for FP mode. */
/*     04 Nov 2014: neH, indH, locH, Hcol added as arguments. */
/*     07 Nov 2014: yCon added as arguments. */
/*     13 Nov 2014: Scales applied for lvlScale > 0. */
/*     14 Feb 2015: proxWeight added. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/* = 0,1,2 or 3, deriv. level */
/* > 0 => some differences needed */
/* > 0 => some exact derivs */
/* > 0 => constant Jacob elements */
/* calls to fCon: all modes */
/* calls to fCon: modes = 1, 2 */
/* calls to fObj: all modes */
/* calls to fObj: modes = 1, 2 */
/*     ------------------------------------------------------------------ */
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
    nnobju = iw[22];
/* # of objective variables */
    lvlscale = iw[75];
/* scale option */
    minmax = iw[87];
/* 1, 0, -1  => MIN, FP, MAX */
    lvltim = iw[182];
/* Timing level */
    lscales = iw[296];
/* scales(nb) = row and column scales */
    lx0 = iw[298];
/* x0(nnH)    = Feasible starting point */
    lxscaled = iw[302];
/* xScaled(n) = copy of scaled  x */
    lgconu = iw[319];
/* record of unknown derivatives and constants */
    lgobju = iw[323];
/* record of unknown derivatives */
    lgsave = iw[339];
/* gSav(nnObjU) holds user-defined gObj */
    gdummy = rw[69];
/* definition of an 'unset' value */
    proxweight = rw[91];
/* Proximal-point weight */
    *iexit = 0;
    fponly = minmax == 0;
    scaled = lvlscale > 0;
    mode = *modefg;
/*     Determine the status of this call. */
    s8callstatus_(&status, &iw[1], leniw);
    if (status == 1) {
/* -------------------------------------------------------------- */
/* First evaluation of the problem functions in snOptC */
/* On entry, lvlScale = 0. */
/* -------------------------------------------------------------- */
	iw[183] = 0;
/* > 0 => some differences needed */
	iw[184] = 0;
/* > 0 => some exact derivatives provided */
	iw[185] = 0;
/* > 0 => constant Jacobian elements provided */
	snprnt_(&c__13, " ", &iw[1], leniw, (ftnlen)1);
	dload_(negcon, &gdummy, &gcon[1], &c__1);
	dload_(nnobj, &gdummy, &gobj[1], &c__1);
    }
/* ----------------------------------------------------------------- */
/* Unscale x (never required for Status = 1) */
/* ----------------------------------------------------------------- */
    if (scaled) {
	dcopy_(nnh, &x[1], &c__1, &rw[lxscaled], &c__1);
	ddscl_(nnh, &rw[lscales], &c__1, &x[1], &c__1);
/* If the Jacobian has some constant elements, they are wrecked */
/* by the scaling.  Restore them from gConu. */
	if (*nonlinearcon) {
	    if (*modefg > 0 && iw[185] > 0) {
		dcopy_(negcon, &rw[lgconu], &c__1, &gcon[1], &c__1);
	    }
	}
    }
/* ================================================================= */
/* Compute the user-defined functions and derivatives. */
/* ================================================================= */
    statususer = status;
/* In case userfun alters Status */
    if (lvltim >= 2) {
	s1time_(&c__4, &c__0, &iw[1], leniw, &rw[1], lenrw);
    }
    if (fponly) {
	(*userfun)(&mode, &nnobju, nncon, nnjac, nnh, negcon, &x[1], fobj, &
		rw[lgsave], &fcon[1], &gcon[1], &statususer, cu + 8, lencu, &
		iu[1], leniu, &ru[1], lenru, (ftnlen)8);
	dcopy_(nnobj, &x[1], &c__1, &gobj[1], &c__1);
	daxpy_(nnobj, &c_b16, &rw[lx0], &c__1, &gobj[1], &c__1);
	*fobj = proxweight * .5 * ddot_(nnobj, &gobj[1], &c__1, &gobj[1], &
		c__1);
	dscal_(nnobj, &proxweight, &gobj[1], &c__1);
    } else {
/* nnObj = nnObjU */
	(*userfun)(&mode, &nnobju, nncon, nnjac, nnh, negcon, &x[1], fobj, &
		gobj[1], &fcon[1], &gcon[1], &statususer, cu + 8, lencu, &iu[
		1], leniu, &ru[1], lenru, (ftnlen)8);
    }
    if (lvltim >= 2) {
	s1time_(&c_n4, &c__0, &iw[1], leniw, &rw[1], lenrw);
    }
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
	dcopy_(nnh, &rw[lxscaled], &c__1, &x[1], &c__1);
	if (*nonlinearcon && mode >= 0) {
	    dddiv_(nncon, &rw[lscales + *n], &c__1, &fcon[1], &c__1);
	    if (*modefg > 0 && iw[184] > 0) {
		s8scalej_(nncon, nnjac, negcon, n, &rw[lscales], nej, nlocj, &
			locj[1], &indj[1], &gcon[1], &rw[1], lenrw);
	    }
	}
	if (*nonlinearobj && mode >= 0) {
	    if (*modefg > 0 && iw[184] > 0) {
		s8scaleg_(nnobj, &rw[lscales], &gobj[1], &rw[1], lenrw);
	    }
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
/* Do some housekeeping on the first entry. */
/* ================================================================= */
    if (status == 1 && *iexit == 0) {
	if (*nonlinearcon) {
/* ----------------------------------------------------------- */
/* Count how many Jacobian elements are provided. */
/* ----------------------------------------------------------- */
	    nnglin = 0;
	    supplied = 0;
	    i__1 = *negcon;
	    for (l = 1; l <= i__1; ++l) {
		if (gcon[l] != gdummy) {
		    ++supplied;
		}
	    }
	    s_wsfi(&io___22);
	    do_fio(&c__1, (char *)&supplied, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*negcon), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)80);
	    if (supplied < *negcon) {
/* Some Jacobian elements are missing. */
		if (iw[70] >= 2) {
/* ----------------------------------------------------- */
/* All the Jacobian is known.  Any undefined elements */
/* are assumed constant, and are restored from gConu. */
/* ----------------------------------------------------- */
		    snprnt_(&c__3, " ==>  Some constraint derivatives are mi"
			    "ssing,  assumed constant.", &iw[1], leniw, (
			    ftnlen)65);
		    snprnt_(&c__3, " ", &iw[1], leniw, (ftnlen)1);
		    lg = lgconu;
		    i__1 = *negcon;
		    for (l = 1; l <= i__1; ++l) {
			if (gcon[l] == gdummy) {
			    gcon[l] = rw[lg];
			    ++nnglin;
			}
			++lg;
		    }
		} else {
/* ----------------------------------------------------- */
/* Save a permanent copy of gCon in gConu so that we */
/* know which derivatives must be estimated. */
/* ----------------------------------------------------- */
		    dcopy_(negcon, &gcon[1], &c__1, &rw[lgconu], &c__1);
		}
	    }
/* supplied < negCon */
	    if (supplied + nnglin < *negcon) {
		iw[183] = 1;
	    }
	    if (supplied > 0) {
		iw[184] = 1;
	    }
	    if (nnglin > 0) {
		iw[185] = 1;
	    }
	}
	if (*nonlinearobj) {
/* ----------------------------------------------------------- */
/* Count how many working gradient elements are known. */
/* (These may be the gradients of the FP objective.) */
/* ----------------------------------------------------------- */
	    if (fponly) {
		supplied = *nnobj;
		iw[184] = 1;
		s_wsfi(&io___24);
		do_fio(&c__1, (char *)&(*nnobj), (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)80);
	    } else {
		supplied = 0;
		i__1 = nnobju;
		for (l = 1; l <= i__1; ++l) {
		    if (gobj[l] != gdummy) {
			++supplied;
		    }
		}
		s_wsfi(&io___25);
		do_fio(&c__1, (char *)&supplied, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&nnobju, (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)80);
		if (supplied > 0) {
		    iw[184] = 1;
		}
		if (supplied < nnobju) {
/* Some objective gradients are missing. */
		    iw[183] = 1;
		    if (iw[70] == 1 || iw[70] == 3) {
/* -------------------------------------------------- */
/* The objective gradient was meant to be known. */
/* -------------------------------------------------- */
			--iw[70];
			s_wsfi(&io___26);
			do_fio(&c__1, (char *)&iw[70], (ftnlen)sizeof(integer)
				);
			e_wsfi();
			snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)80);
		    }
		}
	    }
/* -------------------------------------------------------- */
/* Copy gObj into gObju. */
/* -------------------------------------------------------- */
	    dcopy_(nnobj, &gobj[1], &c__1, &rw[lgobju], &c__1);
	}
    }
    return 0;
} /* s0fgc_ */

