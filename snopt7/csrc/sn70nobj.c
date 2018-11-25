/* ../snopt7/src/sn70nobj.f -- translated by f2c (version 20100827).
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

static doublereal c_b3 = .8;
static doublereal c_b4 = .2;
static integer c__0 = 0;
static integer c__1 = 1;
static integer c__11 = 11;
static doublereal c_b29 = 0.;
static doublereal c_b235 = -1.;
static doublereal c_b242 = 1.;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*     File  sn70nobj.f */

/*     s7checkA   s7checkG   s7checkHv   s7checkp   s7fixX   s7Jac */
/*     s7pert     s7step */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Subroutine */ int s7checka_(integer *iexit, integer *lvlchk, S_fp userfg, 
	integer *nf, integer *n, integer *igmax, doublereal *egmax, integer *
	elem, integer *set, integer *igfun, integer *jgvar, integer *leng, 
	integer *neg, doublereal *x, doublereal *f, doublereal *g, doublereal 
	*x1, doublereal *f1, doublereal *g1, doublereal *w, doublereal *w1, 
	doublereal *y, integer *iw, integer *leniw, doublereal *rw, integer *
	lenrw, char *cu, integer *lencu, integer *iu, integer *leniu, 
	doublereal *ru, integer *lenru, ftnlen cu_len)
{
    /* Initialized data */

    static char lwrong[4] = "bad?";
    static char lright[4] = "ok  ";

    /* Format strings */
    static char fmt_1601[] = "(\002 -->  The largest discrepancy was\002,1p,"
	    "e12.2,\002  in row\002,i6)";
    static char fmt_2001[] = "(\002 Column       x(j)        dx(j)\002,4x"
	    ",\002 Element\002,7x,\002Row\002,7x,\002 Derivative    Differenc"
	    "e approxn\002)";
    static char fmt_2101[] = "(i7,1p,e16.8,e10.2,2i10,2e18.8,2x,a4)";
    static char fmt_2201[] = "(33x,2i10,1pe18.8,e18.8,2x,a4)";
    static char fmt_2301[] = "(i7,2x,\002Nonzero not in sparse structure "
	    "\002,\002??\002,4x,i6,18x,1p,e18.8,2x,a4)";
    static char fmt_2302[] = "(9x,\002Nonzero not in sparse structure \002"
	    ",\002??\002,4x,i6,18x,1p,e18.8,2x,a4)";
    static char fmt_2501[] = "(\002 All\002,i7,\002 assigned derivatives see"
	    "m to be OK.\002)";
    static char fmt_2601[] = "(\002 XXX  There seem to be\002,i6,\002 incorr"
	    "ect derivatives.\002)";
    static char fmt_2701[] = "(\002 -->  The largest relative error was\002,"
	    "1p,e12.2,\002   in row\002,i6,\002,  column\002,i6,\002, elemen"
	    "t\002,i6)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), sqrt(doublereal);
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__, j, k, jj, ng;
    static doublereal dx, xj, yj, gfd, gij, eps;
    static char key[4];
    static doublereal err;
    static char str[120];
    static doublereal eps0, eps5, gmax;
    static integer nset;
    extern /* Subroutine */ int s8callstatus_(integer *, integer *, integer *)
	    ;
    extern doublereal s1eps_(void);
    static integer needf, needg;
    extern /* Subroutine */ int dload_(integer *, doublereal *, doublereal *, 
	    integer *), iload_(integer *, integer *, integer *, integer *);
    static integer jgmax, kgmax;
    static logical first;
    static doublereal fdint1;
    extern integer idamax_(integer *, doublereal *, integer *);
    static doublereal damper;
    static integer rightg;
    extern doublereal dnormi_(integer *, doublereal *, integer *);
    static doublereal gdummy;
    static integer wrongg, status;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen), s6userf_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, S_fp, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, char *, integer 
	    *, integer *, integer *, doublereal *, integer *, integer *, 
	    integer *, ftnlen);
    static logical badgijs;

    /* Fortran I/O blocks */
    static icilist io___21 = { 0, str, 0, fmt_1601, 120, 1 };
    static icilist io___22 = { 0, str, 0, fmt_2001, 120, 1 };
    static icilist io___35 = { 0, str, 0, fmt_2101, 120, 1 };
    static icilist io___36 = { 0, str, 0, fmt_2201, 120, 1 };
    static icilist io___37 = { 0, str, 0, fmt_2301, 120, 1 };
    static icilist io___38 = { 0, str, 0, fmt_2302, 120, 1 };
    static icilist io___39 = { 0, str, 0, fmt_2501, 120, 1 };
    static icilist io___40 = { 0, str, 0, fmt_2601, 120, 1 };
    static icilist io___41 = { 0, str, 0, fmt_2701, 120, 1 };


/*     ================================================================== */
/*     s7checkA   does the work for sncheckA, the stand-alone derivative */
/*     checker for snOptA. */

/*     01 Jun 2006: First version of s7checkA. */
/*     22 Apr 2007: Some features added and exit messages revised. */
/*     15 Jun 2008: Call-status implemented correctly. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --f1;
    --f;
    --elem;
    --y;
    --w1;
    --w;
    --x1;
    --x;
    --set;
    --g1;
    --g;
    --jgvar;
    --igfun;
    --iw;
    --rw;
    cu -= 8;
    --iu;
    --ru;

    /* Function Body */
/*     ------------------------------------------------------------------ */
    *iexit = 0;
    nset = 0;
/* number of derivatives to be checked */
    badgijs = FALSE_;
/* so far so good ... */
    if (*lvlchk < 0) {
	goto L900;
    }
/* Determine the status of this call. */
    s8callstatus_(&status, &iw[1], leniw);
    needf = 1;
/* Get both F and G */
    needg = 1;
    (*userfg)(&status, n, &x[1], &needf, nf, &f[1], &needg, leng, &g[1], cu + 
	    8, lencu, &iu[1], leniu, &ru[1], lenru, (ftnlen)8);
    if (status < 0) {
	goto L900;
    }
/*     ------------------------------------------------------------------ */
/*     Assign some constants. */
/*     ------------------------------------------------------------------ */
    gdummy = rw[69];
/* definition of an 'unset' value */
    eps = s1eps_();
/* machine precision */
    eps0 = pow_dd(&eps, &c_b3);
    eps5 = pow_dd(&eps, &c_b4);
/* eps**(1/5) */
    fdint1 = sqrt(eps);
    *igmax = 0;
    *egmax = 0.;
/*     ------------------------------------------------------------------ */
/*     Cheap test. */
/*     ------------------------------------------------------------------ */
    yj = 1. / *n;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	y[j] = yj;
	yj = -yj * .99999;
    }
/*     Thumb through G and count the assigned elements in each column. */
/*     y(j) is set to zero if the jth column has some unknown elements. */
    nset = 0;
    iload_(n, &c__0, &set[1], &c__1);
    i__1 = *neg;
    for (k = 1; k <= i__1; ++k) {
	j = jgvar[k];
	if (g[k] == gdummy) {
/* G(k) was not set */
	    y[j] = 0.;
/* x(j) is not perturbed */
	} else {
	    ++set[j];
	    ++nset;
	}
    }
    if (nset == 0) {
	snprnt_(&c__11, " No derivatives to be checked", &iw[1], leniw, (
		ftnlen)29);
    } else {
	if (*lvlchk == 0) {
	    snprnt_(&c__11, " Performing a cheap test of the derivatives...", 
		    &iw[1], leniw, (ftnlen)46);
	} else if (*lvlchk > 0) {
	    snprnt_(&c__11, " Performing a full test of the derivatives...", &
		    iw[1], leniw, (ftnlen)45);
	}
/*        --------------------------------------------------------------- */
/*        Compute F at  x1 = x + dx*y. */
/*        s6Userf reduces the step until F(x1) is well-defined. */
/*        --------------------------------------------------------------- */
	dx = fdint1 * (dnormi_(n, &x[1], &c__1) + 1.);
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    x1[j] = x[j] + dx * y[j];
	}
	needf = 1;
	needg = 0;
/* G isn't needed */
	s6userf_(&status, n, &x[1], &x1[1], &damper, (S_fp)userfg, &needf, nf,
		 &f1[1], &needg, leng, &g1[1], cu + 8, lencu, &iu[1], leniu, &
		ru[1], lenru, &iw[1], leniw, (ftnlen)8);
	if (status < 0) {
	    goto L900;
	}
/*        Set   w1 = (F1 - F)/dx - G*y.  This should be small. */
	dx = damper * dx;
	i__1 = *nf;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    w1[i__] = (f1[i__] - f[i__]) / dx;
	}
	i__1 = *neg;
	for (k = 1; k <= i__1; ++k) {
	    i__ = igfun[k];
	    j = jgvar[k];
	    w1[i__] -= g[k] * y[j];
	}
	*igmax = idamax_(nf, &w1[1], &c__1);
	gmax = (f1[*igmax] - f[*igmax]) / dx;
	*egmax = (d__1 = w1[*igmax], abs(d__1)) / (abs(gmax) + 1.);
	if (*egmax <= .1) {
	    snprnt_(&c__11, " All the assigned derivatives seem to be OK.", &
		    iw[1], leniw, (ftnlen)44);
	} else {
	    snprnt_(&c__11, " XXX  Some derivatives seem to be incorrect.", &
		    iw[1], leniw, (ftnlen)44);
	}
	s_wsfi(&io___21);
	do_fio(&c__1, (char *)&(*egmax), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&(*igmax), (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__11, str, &iw[1], leniw, (ftnlen)120);
	snprnt_(&c__1, " ", &iw[1], leniw, (ftnlen)1);
	badgijs = *egmax > .1;
    }
    if (*lvlchk == 0) {
	goto L900;
    }
/*     ------------------------------------------------------------------ */
/*     Proceed with column-wise verification. */
/*     ------------------------------------------------------------------ */
    s_wsfi(&io___22);
    e_wsfi();
    snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)120);
    ng = 0;
    rightg = 0;
    wrongg = 0;
    *egmax = -1.;
    *igmax = 0;
    jgmax = 0;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*        ------------------------------------------------------------ */
/*        Estimate the jth column of G. */
/*        Check the estimate against the user-supplied values. */
/*        Don't bother printing a line for an exact zero. */
/*        Look for nonzeros that are not in the sparse data structure. */
/*        ------------------------------------------------------------ */
	if (set[j] > 0) {
	    xj = x[j];
	    dx = fdint1 * (abs(xj) + 1.);
	    x[j] = xj + dx;
	    status = 0;
/* no need to check the call-status */
	    (*userfg)(&status, n, &x[1], &needf, nf, &f1[1], &needg, leng, &
		    g1[1], cu + 8, lencu, &iu[1], leniu, &ru[1], lenru, (
		    ftnlen)8);
	    if (status < 0) {
		goto L900;
	    }
/* Find the elements in this column */
	    dload_(nf, &c_b29, &w[1], &c__1);
	    i__2 = *neg;
	    for (k = 1; k <= i__2; ++k) {
		jj = jgvar[k];
		if (jj == j) {
		    i__ = igfun[k];
		    w[i__] = g[k];
		    elem[i__] = k;
		}
	    }
	    first = TRUE_;
	    i__2 = *nf;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		gij = w[i__];
		gfd = (f1[i__] - f[i__]) / dx;
		w1[i__] = gfd;
		k = elem[i__];
		if (gij != gdummy) {
		    ++ng;
/* Computing MAX */
		    d__3 = (d__1 = gfd - gij, abs(d__1)) / (abs(gij) + 1.), 
			    d__4 = (d__2 = gfd - gij, abs(d__2)) / (abs(gfd) 
			    + 1.);
		    err = max(d__3,d__4);
		    if (*egmax < err) {
			*egmax = err;
			*igmax = i__;
			jgmax = j;
			kgmax = k;
		    }
		    if (err <= eps5) {
			s_copy(key, lright, (ftnlen)4, (ftnlen)4);
			++rightg;
		    } else {
			s_copy(key, lwrong, (ftnlen)4, (ftnlen)4);
			++wrongg;
		    }
		    if (abs(gij) + err > eps0) {
			if (first) {
			    s_wsfi(&io___35);
			    do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer))
				    ;
			    do_fio(&c__1, (char *)&xj, (ftnlen)sizeof(
				    doublereal));
			    do_fio(&c__1, (char *)&dx, (ftnlen)sizeof(
				    doublereal));
			    do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer))
				    ;
			    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(
				    integer));
			    do_fio(&c__1, (char *)&gij, (ftnlen)sizeof(
				    doublereal));
			    do_fio(&c__1, (char *)&gfd, (ftnlen)sizeof(
				    doublereal));
			    do_fio(&c__1, key, (ftnlen)4);
			    e_wsfi();
			    snprnt_(&c__11, str, &iw[1], leniw, (ftnlen)120);
			    first = FALSE_;
			} else {
			    s_wsfi(&io___36);
			    do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer))
				    ;
			    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(
				    integer));
			    do_fio(&c__1, (char *)&gij, (ftnlen)sizeof(
				    doublereal));
			    do_fio(&c__1, (char *)&gfd, (ftnlen)sizeof(
				    doublereal));
			    do_fio(&c__1, key, (ftnlen)4);
			    e_wsfi();
			    snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)120);
			}
		    }
		}
		w1[i__] = gdummy;
/* Done with this row */
	    }
/*           Check that elements remaining in w are zero. */
	    i__2 = *nf;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		if (w1[i__] != gdummy) {
/* w1(i) not yet processed */
		    err = (d__1 = w1[i__], abs(d__1));
		    if (err > eps0) {
			++wrongg;
			if (first) {
			    s_wsfi(&io___37);
			    do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer))
				    ;
			    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(
				    integer));
			    do_fio(&c__1, (char *)&w1[i__], (ftnlen)sizeof(
				    doublereal));
			    do_fio(&c__1, lwrong, (ftnlen)4);
			    e_wsfi();
			    snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)120);
			    first = FALSE_;
			} else {
			    s_wsfi(&io___38);
			    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(
				    integer));
			    do_fio(&c__1, (char *)&w1[i__], (ftnlen)sizeof(
				    doublereal));
			    do_fio(&c__1, lwrong, (ftnlen)4);
			    e_wsfi();
			    snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)120);
			}
			if (*egmax < err) {
			    *egmax = err;
			    *igmax = i__;
			    jgmax = j;
			}
		    }
		}
	    }
	    x[j] = xj;
	}
/* set */
    }
/*     ------------------------------------------------------------------ */
/*     Final tally of the good, the bad and the ugly. */
/*     ------------------------------------------------------------------ */
    if (wrongg == 0) {
	s_wsfi(&io___39);
	do_fio(&c__1, (char *)&rightg, (ftnlen)sizeof(integer));
	e_wsfi();
    } else {
	s_wsfi(&io___40);
	do_fio(&c__1, (char *)&wrongg, (ftnlen)sizeof(integer));
	e_wsfi();
    }
    snprnt_(&c__11, str, &iw[1], leniw, (ftnlen)120);
    s_wsfi(&io___41);
    do_fio(&c__1, (char *)&(*egmax), (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&(*igmax), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&jgmax, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&kgmax, (ftnlen)sizeof(integer));
    e_wsfi();
    snprnt_(&c__11, str, &iw[1], leniw, (ftnlen)120);
/*      call snPRNT(  1, ' ', iw, leniw ) */
    badgijs = badgijs || wrongg > 0;
/*     ------------------------------------------------------------------ */
/*     Set iExit for exceptions. */
/*     ------------------------------------------------------------------ */
L900:
    if (*lvlchk < 0 || nset == 0) {
	*iexit = 106;
/* no derivatives were checked */
    } else if (badgijs) {
	*iexit = 55;
/* Some bad derivatives */
    } else if (status < 0) {
	if (status == -1) {
	    *iexit = 62;
/* undefined function at the first point */
	} else {
	    *iexit = 71;
/* terminated during function evaluation */
	}
    }
    return 0;
} /* s7checka_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s7checkA */
/* Subroutine */ int s7checkg_(integer *iexit, integer *n, integer *nncon0, 
	integer *nncon, integer *nnjac, integer *nnl, integer *nnobj0, 
	integer *nnobj, S_fp funwrapper, U_fp fgcon, U_fp fgobj, U_fp userhv, 
	doublereal *x, doublereal *x1, doublereal *bl, doublereal *bu, 
	doublereal *fobj, doublereal *gobj, doublereal *ycon, integer *nej, 
	integer *nlocj, integer *locj, integer *indj, integer *neh, integer *
	nloch, integer *loch, integer *indh, doublereal *hcol, integer *
	negcon, integer *nlocg, integer *locg, doublereal *fcon, doublereal *
	gcon, doublereal *gobj2, doublereal *fcon2, doublereal *gcon2, 
	doublereal *y, doublereal *y1, char *cu, integer *lencu, integer *iu, 
	integer *leniu, doublereal *ru, integer *lenru, char *cw, integer *
	lencw, integer *iw, integer *leniw, doublereal *rw, integer *lenrw, 
	ftnlen cu_len, ftnlen cw_len)
{
    /* Initialized data */

    static char lwrong[4] = "bad?";
    static char lright[4] = "ok  ";

    /* Format strings */
    static char fmt_1601[] = "(\002 -->  The largest discrepancy was\002,1p,"
	    "e12.2,\002  in constraint\002,i6)";
    static char fmt_1602[] = "(\002 Gradient projected in one direction\002,"
	    "1p,e20.11)";
    static char fmt_1603[] = "(\002 Difference approximation           \002,"
	    "1p,e20.11)";
    static char fmt_2001[] = "(\002 Column       x(j)        dx(j)\002,3x"
	    ",\002 Element no.    Row        Derivative    Difference approx"
	    "n\002)";
    static char fmt_2002[] = "(6x,\002j\002,7x,\002x(j)\002,8x,\002dx(j)\002"
	    ",11x,\002g(j)\002,9x,\002Difference approxn\002)";
    static char fmt_2101[] = "(i7,1p,e16.8,e10.2,2i10,2e18.8,2x,a4)";
    static char fmt_2201[] = "(33x,2i10,1pe18.8,e18.8,2x,a4)";
    static char fmt_2301[] = "(i7,2x,\002Nonzero not in sparse structure "
	    "\002,\002??\002,4x,i6,18x,1p,e18.8,2x,a4)";
    static char fmt_2302[] = "(9x,\002Nonzero not in sparse structure \002"
	    ",\002??\002,4x,i6,18x,1p,e18.8,2x,a4)";
    static char fmt_2102[] = "(i7,1p,e16.8,e10.2,10x,\002 Objective\002,2e18"
	    ".8,2x,a4)";
    static char fmt_2103[] = "(i7,1p,e16.8,e10.2,2e18.8,2x,a4)";
    static char fmt_2501[] = "(i7,\002  Jacobian elements in cols \002,i6"
	    ",\002  thru\002,i6,\002  seem to be OK.\002)";
    static char fmt_2601[] = "(\002 XXX  There seem to be\002,i6,\002  incor"
	    "rect Jacobian elements in cols\002,i6,\002  thru\002,i6)";
    static char fmt_2701[] = "(\002 -->  The largest relative error was\002,"
	    "1p,e12.2,\002   in row\002,i6,\002,  column\002,i6)";
    static char fmt_2502[] = "(i7,\002  objective gradients out of\002,i6"
	    ",\002  thru\002,i6,\002  seem to be OK.\002)";
    static char fmt_2602[] = "(\002 XXX  There seem to be\002,i6,\002  incor"
	    "rect objective gradients in cols\002,i6,\002  thru\002,i6)";
    static char fmt_2702[] = "(\002 -->  The largest relative error was\002,"
	    "1p,e12.2,\002   in column\002,i6)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__, j, k, l, j1, j2, j3, j4, k1, k2;
    static doublereal gj;
    static integer ng, jn, ir, nj;
    static doublereal dx, gp, xj, yj, gij;
    static char key[4];
    static integer irn;
    static doublereal err;
    static char str[120];
    static doublereal eps0, eps5, agij;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static logical done;
    static integer kmax;
    static doublereal fobj2;
    static logical nonlinearobj, nonlinearcon;
    static doublereal agdif, gdiff;
    static logical cheap;
    static integer nfeas;
    static doublereal emaxg, emaxj;
    extern doublereal dasum_(integer *, doublereal *, integer *);
    static integer imaxj;
    static logical found;
    static integer jlast, jmaxg, jmaxj;
    static doublereal gmaxg, gmaxj;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical first;
    static doublereal fdint1;
    extern integer s2coln_(integer *, integer *, integer *), s2rown_(integer *
	    , integer *, integer *);
    static integer modefg;
    extern integer idamax_(integer *, doublereal *, integer *);
    static integer rightg, rightj;
    extern doublereal dnormi_(integer *, doublereal *, integer *);
    static integer jfirst;
    static doublereal gdummy;
    static integer wrongg, wrongj, lvlver, nknown;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen), s7checkp_(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *);

    /* Fortran I/O blocks */
    static icilist io___77 = { 0, str, 0, fmt_1601, 120, 1 };
    static icilist io___81 = { 0, str, 0, fmt_1602, 120, 1 };
    static icilist io___82 = { 0, str, 0, fmt_1603, 120, 1 };
    static icilist io___83 = { 0, str, 0, fmt_2001, 120, 1 };
    static icilist io___84 = { 0, str, 0, fmt_2002, 120, 1 };
    static icilist io___104 = { 0, str, 0, fmt_2101, 120, 1 };
    static icilist io___105 = { 0, str, 0, fmt_2201, 120, 1 };
    static icilist io___106 = { 0, str, 0, fmt_2301, 120, 1 };
    static icilist io___107 = { 0, str, 0, fmt_2302, 120, 1 };
    static icilist io___109 = { 0, str, 0, fmt_2102, 120, 1 };
    static icilist io___110 = { 0, str, 0, fmt_2103, 120, 1 };
    static icilist io___111 = { 0, str, 0, fmt_2501, 120, 1 };
    static icilist io___112 = { 0, str, 0, fmt_2601, 120, 1 };
    static icilist io___113 = { 0, str, 0, fmt_2701, 120, 1 };
    static icilist io___114 = { 0, str, 0, fmt_2502, 120, 1 };
    static icilist io___115 = { 0, str, 0, fmt_2602, 120, 1 };
    static icilist io___116 = { 0, str, 0, fmt_2702, 120, 1 };


/*     ================================================================== */
/*     s7checkG  verifies the objective and constraint */
/*     gradients using finite differences. */

/*     First, a cheap heuristic test is performed, as in */
/*     subroutine chkgrd by the following authors: */
/*     Philip E. Gill, Walter Murray, Susan M. Picken and Hazel M. Barber */
/*     DNAC, National Physical Laboratory, England  (circa 1975). */

/*     Next, a more reliable test is performed on each component of the */
/*     gradient, for indices in the range  jverif(1)  thru  jverif(2). */

/*     lvlVer is the verify level, which has the following meaning: */

/*     -1         do not perform any check. */
/*      0         do the cheap test only. */
/*      1 or 3    do both cheap and full test on objective gradients. */
/*      2 or 3    do both cheap and full test on the Jacobian. */

/*     10 Oct 1998: First version based on combining s7chkJ and s7chkg. */
/*     27 Dec 2000: Permutations of the natural order included. */
/*     02 Aug 2003: snEXIT and snPRNT adopted. */
/*     20 Mar 2006: Relative error based on exact and estimated values. */
/*     16 Jun 2008: Call-status implemented correctly. */
/*     15 Nov 2010: Call-status removed from argument list. */
/*     19 Oct 2014: Added user-defined Hessian. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --x1;
    --x;
    --y1;
    --fcon2;
    --fcon;
    --ycon;
    --y;
    --bu;
    --bl;
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
/*     ------------------------------------------------------------------ */
    lvlver = iw[78];
/* Verify level */
    if (lvlver < 0) {
	return 0;
    }
    j1 = iw[110];
/* start col for obj.    derivative checking */
    j2 = iw[111];
/* stop  col for obj.    derivative checking */
    j3 = iw[112];
/* start col for constr. derivative checking */
    j4 = iw[113];
/* stop  col for constr. derivative checking */
    eps0 = rw[2];
/* eps**(4/5) */
    eps5 = rw[7];
/* eps**(1/5) */
    fdint1 = rw[76];
/* (1) forwrd diff. interval */
    gdummy = rw[69];

    *iexit = 0;
    j1 = max(j1,1);
    j2 = min(j2,*n);
    j3 = max(j3,1);
    j4 = min(j4,*n);
    jfirst = min(j1,j3);
    jlast = max(j2,j4);
/*     The problem functions are called to provide functions only. */
    modefg = 0;
/*     Cheap = do the cheap test only */
    cheap = lvlver == 0 || jfirst > jlast;
/*     ------------------------------------------------------------------ */
/*     Cheap test. */
/*     ------------------------------------------------------------------ */
/*     Generate a direction in which to perturb  x. */
    dcopy_(n, &x[1], &c__1, &x1[1], &c__1);
    yj = 1. / *nnl;
    i__1 = *nnl;
    for (j = 1; j <= i__1; ++j) {
	y[j] = yj;
	x1[j] = yj;
/* x1 = y */
	yj = -yj * .99999;
    }
/*     If needed, alter y to ensure that it will be a feasible direction. */
/*     If this gives zero, go back to original y and forget feasibility. */
    dx = fdint1 * (dnormi_(nnl, &x[1], &c__1) + 1.);
    s7checkp_(nnl, &bl[1], &bu[1], &x[1], &dx, &y[1], &nfeas);
    if (nfeas == 0) {
	dcopy_(nnl, &x1[1], &c__1, &y[1], &c__1);
    }
/*     ------------------------------------------------------------------ */
/*     Do not perturb x(j) if the jth column contains unknown elements. */
/*     ------------------------------------------------------------------ */
    nknown = 0;
    l = 0;
    i__1 = *nnl;
    for (j = 1; j <= i__1; ++j) {
/*        Do not perturb x(j) if g(j) is unknown. */
	if (j <= *nnobj) {
	    if (gobj[j] == gdummy) {
		y[j] = 0.;
	    }
	}
	if (j <= *nnjac) {
	    k1 = locj[j];
	    k2 = locj[j + 1] - 1;
	    i__2 = k2;
	    for (k = k1; k <= i__2; ++k) {
		ir = indj[k];
		if (ir > *nncon) {
		    goto L130;
		}
		++l;
		if (gcon[l] == gdummy) {
		    y[j] = 0.;
		}
	    }
	}
L130:
	if (y[j] != 0.) {
	    ++nknown;
	}
    }
    if (nknown > 0) {
	if (cheap) {
	    snprnt_(&c__11, " Cheap test of user-supplied problem derivative"
		    "s...", &iw[1], leniw, (ftnlen)51);
	} else {
	    snprnt_(&c__11, " Verification of user-supplied problem derivati"
		    "ves.", &iw[1], leniw, (ftnlen)51);
	}
/*        --------------------------------------------------------------- */
/*        Compute functions at a short step along  y. */
/*        --------------------------------------------------------------- */
	dx = fdint1 * (dasum_(nnl, &x[1], &c__1) + 1.);
	i__1 = *nnl;
	for (j = 1; j <= i__1; ++j) {
	    x1[j] = x[j] + dx * y[j];
	}
	nonlinearcon = *nnjac > 0;
	nonlinearobj = *nnobj > 0;
	(*funwrapper)(iexit, &modefg, &nonlinearcon, &nonlinearobj, n, negcon,
		 nncon0, nncon, nnjac, nnl, nnobj0, nnobj, (U_fp)fgcon, (U_fp)
		fgobj, (U_fp)userhv, &x1[1], &ycon[1], nej, nlocj, &locj[1], &
		indj[1], neh, nloch, &loch[1], &indh[1], &hcol[1], &fcon2[1], 
		&fobj2, &gcon2[1], &gobj2[1], cu + 8, lencu, &iu[1], leniu, &
		ru[1], lenru, cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (
		ftnlen)8, (ftnlen)8);
	if (*iexit != 0) {
	    goto L900;
	}
/*        --------------------------------------------------------------- */
/*        Cheap test for the constraint Jacobian. */
/*        --------------------------------------------------------------- */
/*        Set   y1 = (fCon2 - fCon)/dx - gCon*y.  This should be small. */
	if (nonlinearcon) {
	    i__1 = *nncon;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		y1[i__] = (fcon2[i__] - fcon[i__]) / dx;
	    }
	    l = 0;
	    i__1 = *nnjac;
	    for (j = 1; j <= i__1; ++j) {
		yj = y[j];
		k = locj[j];
		kmax = locj[j + 1] - 1;
		done = FALSE_;
/* +             --------------------------------------------------------- */
/* +             while (k .le. kmax  .and.  .not. Done) do */
L180:
		if (k <= kmax && ! done) {
		    ir = indj[k];
		    if (ir > *nncon) {
			done = TRUE_;
		    } else {
			++l;
			y1[ir] -= gcon[l] * yj;
		    }
		    ++k;
		    goto L180;
/* +             end while */
/* +             --------------------------------------------------------- */
		}
	    }
	    imaxj = idamax_(nncon, &y1[1], &c__1);
	    gmaxj = (fcon2[imaxj] - fcon[imaxj]) / dx;
	    emaxj = (d__1 = y1[imaxj], abs(d__1)) / (abs(gmaxj) + 1.);
	    if (emaxj <= .1) {
		snprnt_(&c__11, " The constraint gradients seem to be OK.", &
			iw[1], leniw, (ftnlen)40);
	    } else {
		snprnt_(&c__11, " XXX  The constraint gradients seem to be i"
			"ncorrect.", &iw[1], leniw, (ftnlen)52);
	    }
	    s_wsfi(&io___77);
	    do_fio(&c__1, (char *)&emaxj, (ftnlen)sizeof(doublereal));
	    i__1 = *n + s2rown_(&imaxj, leniw, &iw[1]);
	    do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__11, str, &iw[1], leniw, (ftnlen)120);
	    snprnt_(&c__1, " ", &iw[1], leniw, (ftnlen)1);
	}
/*        --------------------------------------------------------------- */
/*        Cheap test for the objective gradient. */
/*        --------------------------------------------------------------- */
	if (nonlinearobj) {
	    gp = ddot_(nnobj, &gobj[1], &c__1, &y[1], &c__1);
	    gmaxg = (fobj2 - *fobj) / dx;
/* Computing MAX */
	    d__3 = (d__1 = gmaxg - gp, abs(d__1)) / (abs(gmaxg) + 1.), d__4 = 
		    (d__2 = gmaxg - gp, abs(d__2)) / (abs(gp) + 1.);
	    emaxg = max(d__3,d__4);
/*           Set an error indicator if emaxG is too large. */
	    if (emaxg <= .1) {
		snprnt_(&c__11, " The objective  gradients seem to be OK.", &
			iw[1], leniw, (ftnlen)40);
	    } else {
		snprnt_(&c__11, " XXX  The objective  gradients seem to be i"
			"ncorrect.", &iw[1], leniw, (ftnlen)52);
	    }
	    s_wsfi(&io___81);
	    do_fio(&c__1, (char *)&gp, (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    snprnt_(&c__11, str, &iw[1], leniw, (ftnlen)120);
	    s_wsfi(&io___82);
	    do_fio(&c__1, (char *)&gmaxg, (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)120);
	}
    }
    if (cheap) {
	goto L900;
    }
/*     ------------------------------------------------------------------ */
/*     Proceed with the verification of column elements. */
/*     Evaluate columns  jFirst thru jLast  of the problem derivatives. */
/*     ------------------------------------------------------------------ */
    snprnt_(&c__11, " ", &iw[1], leniw, (ftnlen)1);
    if (j3 <= j4) {
	s_wsfi(&io___83);
	e_wsfi();
	snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)120);
    } else {
	s_wsfi(&io___84);
	e_wsfi();
	snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)120);
	snprnt_(&c__1, " ", &iw[1], leniw, (ftnlen)1);
    }
    wrongj = 0;
    rightj = 0;
    wrongg = 0;
    rightg = 0;
    jmaxj = 0;
    jmaxg = 0;
    nj = 0;
    ng = 0;
    emaxj = -1.;
    emaxg = -1.;
    i__1 = *nnl;
    for (j = 1; j <= i__1; ++j) {
	jn = s2coln_(&j, leniw, &iw[1]);
	nonlinearobj = j <= *nnobj && jn >= j1 && jn <= j2;
	nonlinearcon = j <= *nnjac && jn >= j3 && jn <= j4;
	if (nonlinearcon || nonlinearobj) {
/*           See if there are any known gradients in this column. */
	    found = FALSE_;
	    if (nonlinearobj) {
		if (gobj[j] != gdummy) {
		    found = TRUE_;
		}
	    }
	    if (nonlinearcon) {
/*              See if there are any known gradients in this column. */
		l = locg[j];
		k = locj[j];
		kmax = locj[j + 1] - 1;
		done = FALSE_;
/* +             --------------------------------------------------------- */
/* +             while (k .le. kmax  .and. .not.(Found  .or.  Done)) do */
L200:
		if (k <= kmax && ! (found || done)) {
		    ir = indj[k];
		    if (ir > *nncon) {
			done = TRUE_;
		    } else {
			if (gcon[l] != gdummy) {
			    found = TRUE_;
			}
			++l;
		    }
		    ++k;
		    goto L200;
		}
/* +             end while */
/* +             --------------------------------------------------------- */
	    }
	    if (found) {
/*              --------------------------------------------------------- */
/*              gObj(j) or an element of the jth column of J is known. */
/*              --------------------------------------------------------- */
		xj = x[j];
		dx = fdint1 * (abs(xj) + 1.);
		if (bl[j] < bu[j] && xj >= bu[j]) {
		    dx = -dx;
		}
		x[j] = xj + dx;
		(*funwrapper)(iexit, &modefg, &nonlinearcon, &nonlinearobj, n,
			 negcon, nncon0, nncon, nnjac, nnl, nnobj0, nnobj, (
			U_fp)fgcon, (U_fp)fgobj, (U_fp)userhv, &x[1], &ycon[1]
			, nej, nlocj, &locj[1], &indj[1], neh, nloch, &loch[1]
			, &indh[1], &hcol[1], &fcon2[1], &fobj2, &gcon2[1], &
			gobj2[1], cu + 8, lencu, &iu[1], leniu, &ru[1], lenru,
			 cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)
			8, (ftnlen)8);
		if (*iexit != 0) {
		    goto L900;
		}
/*              --------------------------------------------------------- */
/*              Estimate the jth column of the Jacobian. */
/*              Check the estimate against the user-supplied values. */
/*              Don't bother printing a line for an exact zero. */
/*              Look for nonzeros not in the sparse data structure. */
/*              --------------------------------------------------------- */
		if (nonlinearcon) {
		    i__2 = *nncon;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			y1[i__] = (fcon2[i__] - fcon[i__]) / dx;
		    }
		    l = locg[j];
		    k = locj[j];
		    first = TRUE_;
		    done = FALSE_;
/* +                ------------------------------------------------------ */
/* +                while (k .le. kmax  .and.  .not. Done) do */
L260:
		    if (k <= kmax && ! done) {
			ir = indj[k];
			irn = s2rown_(&ir, leniw, &iw[1]);
			if (ir > *nncon) {
			    done = TRUE_;
			} else {
			    gij = gcon[l];
			    agij = abs(gij);
			    if (gij != gdummy) {
				++nj;
				gdiff = y1[ir];
				agdif = abs(gdiff);
/* Computing MAX */
				d__3 = (d__1 = gdiff - gij, abs(d__1)) / (
					agij + 1.), d__4 = (d__2 = gdiff - 
					gij, abs(d__2)) / (agdif + 1.);
				err = max(d__3,d__4);
				if (emaxj < err) {
				    emaxj = err;
				    imaxj = irn;
				    jmaxj = jn;
				}
				if (err <= eps5) {
				    s_copy(key, lright, (ftnlen)4, (ftnlen)4);
				    ++rightj;
				} else {
				    s_copy(key, lwrong, (ftnlen)4, (ftnlen)4);
				    ++wrongj;
				}
				if (agij + err > eps0) {
				    if (first) {
					s_wsfi(&io___104);
					do_fio(&c__1, (char *)&jn, (ftnlen)
						sizeof(integer));
					do_fio(&c__1, (char *)&xj, (ftnlen)
						sizeof(doublereal));
					do_fio(&c__1, (char *)&dx, (ftnlen)
						sizeof(doublereal));
					do_fio(&c__1, (char *)&l, (ftnlen)
						sizeof(integer));
					do_fio(&c__1, (char *)&irn, (ftnlen)
						sizeof(integer));
					do_fio(&c__1, (char *)&gij, (ftnlen)
						sizeof(doublereal));
					do_fio(&c__1, (char *)&gdiff, (ftnlen)
						sizeof(doublereal));
					do_fio(&c__1, key, (ftnlen)4);
					e_wsfi();
					snprnt_(&c__11, str, &iw[1], leniw, (
						ftnlen)120);
					first = FALSE_;
				    } else {
					s_wsfi(&io___105);
					do_fio(&c__1, (char *)&l, (ftnlen)
						sizeof(integer));
					do_fio(&c__1, (char *)&irn, (ftnlen)
						sizeof(integer));
					do_fio(&c__1, (char *)&gij, (ftnlen)
						sizeof(doublereal));
					do_fio(&c__1, (char *)&gdiff, (ftnlen)
						sizeof(doublereal));
					do_fio(&c__1, key, (ftnlen)4);
					e_wsfi();
					snprnt_(&c__1, str, &iw[1], leniw, (
						ftnlen)120);
				    }
				}
			    }
			    ++l;
/*                       Mark this row as being in the sparse structure. */
			    y1[ir] = gdummy;
			}
			++k;
			goto L260;
		    }
/* +                end while */
/* +                ------------------------------------------------------ */
/*                 Check that column elements not in gCon  are zero. */
/*                 These are the unmarked elements of  y1. */
		    i__2 = *nncon;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			if (y1[i__] != gdummy) {
			    irn = s2rown_(&i__, leniw, &iw[1]);
			    err = (d__1 = y1[i__], abs(d__1));
			    if (err > eps0) {
				++wrongj;
				if (first) {
				    s_wsfi(&io___106);
				    do_fio(&c__1, (char *)&jn, (ftnlen)sizeof(
					    integer));
				    do_fio(&c__1, (char *)&irn, (ftnlen)
					    sizeof(integer));
				    do_fio(&c__1, (char *)&y1[i__], (ftnlen)
					    sizeof(doublereal));
				    do_fio(&c__1, lwrong, (ftnlen)4);
				    e_wsfi();
				    snprnt_(&c__1, str, &iw[1], leniw, (
					    ftnlen)120);
				    first = FALSE_;
				} else {
				    s_wsfi(&io___107);
				    do_fio(&c__1, (char *)&irn, (ftnlen)
					    sizeof(integer));
				    do_fio(&c__1, (char *)&y1[i__], (ftnlen)
					    sizeof(doublereal));
				    do_fio(&c__1, lwrong, (ftnlen)4);
				    e_wsfi();
				    snprnt_(&c__1, str, &iw[1], leniw, (
					    ftnlen)120);
				}
				if (emaxj < err) {
				    emaxj = err;
				    imaxj = irn;
				    jmaxj = jn;
				}
			    }
			}
		    }
		}
/*              --------------------------------------------------------- */
/*              Estimate gObj(j) */
/*              --------------------------------------------------------- */
/* NonlinearCon */
		if (nonlinearobj) {
		    gj = gobj[j];
		    if (gj != gdummy) {
			++ng;
			gdiff = (fobj2 - *fobj) / dx;
/* Computing MAX */
			d__3 = (d__1 = gdiff - gj, abs(d__1)) / (abs(gj) + 1.)
				, d__4 = (d__2 = gdiff - gj, abs(d__2)) / (
				abs(gdiff) + 1.);
			err = max(d__3,d__4);
			if (err > emaxg) {
			    emaxg = err;
			    jmaxg = jn;
			}
			if (err <= eps5) {
			    s_copy(key, lright, (ftnlen)4, (ftnlen)4);
			    ++rightg;
			} else {
			    s_copy(key, lwrong, (ftnlen)4, (ftnlen)4);
			    ++wrongg;
			}
			if (abs(gj) + err > eps0) {
			    if (j3 <= j4) {
				s_wsfi(&io___109);
				do_fio(&c__1, (char *)&jn, (ftnlen)sizeof(
					integer));
				do_fio(&c__1, (char *)&xj, (ftnlen)sizeof(
					doublereal));
				do_fio(&c__1, (char *)&dx, (ftnlen)sizeof(
					doublereal));
				do_fio(&c__1, (char *)&gj, (ftnlen)sizeof(
					doublereal));
				do_fio(&c__1, (char *)&gdiff, (ftnlen)sizeof(
					doublereal));
				do_fio(&c__1, key, (ftnlen)4);
				e_wsfi();
			    } else {
				s_wsfi(&io___110);
				do_fio(&c__1, (char *)&jn, (ftnlen)sizeof(
					integer));
				do_fio(&c__1, (char *)&xj, (ftnlen)sizeof(
					doublereal));
				do_fio(&c__1, (char *)&dx, (ftnlen)sizeof(
					doublereal));
				do_fio(&c__1, (char *)&gj, (ftnlen)sizeof(
					doublereal));
				do_fio(&c__1, (char *)&gdiff, (ftnlen)sizeof(
					doublereal));
				do_fio(&c__1, key, (ftnlen)4);
				e_wsfi();
			    }
			    snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)120);
			}
		    }
		}
/* NonlinearObj */
		x[j] = xj;
	    }
/* Found */
	}
    }
/*     ------------------------------------------------------------------ */
/*     Final tally of the good, the bad and the ugly. */
/*     ------------------------------------------------------------------ */
    if (j3 <= j4 && nj > 0) {
	if (wrongj == 0) {
	    s_wsfi(&io___111);
	    do_fio(&c__1, (char *)&rightj, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&j3, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&j4, (ftnlen)sizeof(integer));
	    e_wsfi();
	} else {
	    s_wsfi(&io___112);
	    do_fio(&c__1, (char *)&wrongj, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&j3, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&j4, (ftnlen)sizeof(integer));
	    e_wsfi();
	}
	snprnt_(&c__11, str, &iw[1], leniw, (ftnlen)120);
	s_wsfi(&io___113);
	do_fio(&c__1, (char *)&emaxj, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&imaxj, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&jmaxj, (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__11, str, &iw[1], leniw, (ftnlen)120);
	snprnt_(&c__1, " ", &iw[1], leniw, (ftnlen)1);
	if (emaxj >= 1.) {
	    *iexit = 52;
/* Bad constraint gradients. */
	}
    }
    if (j1 <= j2 && ng > 0) {
	if (wrongg == 0) {
	    s_wsfi(&io___114);
	    do_fio(&c__1, (char *)&rightg, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&j1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&j2, (ftnlen)sizeof(integer));
	    e_wsfi();
	} else {
	    s_wsfi(&io___115);
	    do_fio(&c__1, (char *)&wrongg, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&j1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&j2, (ftnlen)sizeof(integer));
	    e_wsfi();
	}
	snprnt_(&c__11, str, &iw[1], leniw, (ftnlen)120);
	s_wsfi(&io___116);
	do_fio(&c__1, (char *)&emaxg, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&jmaxg, (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__11, str, &iw[1], leniw, (ftnlen)120);
	snprnt_(&c__1, " ", &iw[1], leniw, (ftnlen)1);
	if (emaxg >= 1.) {
	    *iexit = 51;
/* Bad objective gradients. */
	}
    }
/*     ------------------------------------------------------------------ */
/*     Print a message if we had to abandon the check. */
/*     ------------------------------------------------------------------ */
L900:
    if (*iexit < 0) {
	snprnt_(&c__11, " XXX  Unable to complete derivative check.", &iw[1], 
		leniw, (ftnlen)42);
	*iexit = 0;
    }
    return 0;
} /* s7checkg_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s7checkG */
/* Subroutine */ int s7checkhv_(integer *iexit, integer *n, integer *nncon0, 
	integer *nncon, integer *nnjac, integer *nnh, integer *nnobj0, 
	integer *nnobj, S_fp funwrapper, U_fp funcon, U_fp funobj, U_fp 
	userhv, doublereal *x, doublereal *x1, doublereal *bl, doublereal *bu,
	 integer *nej, integer *nlocj, integer *locj, integer *indj, integer *
	negcon, integer *nlocg, integer *locg, integer *neh, integer *nloch, 
	integer *loch, integer *indh, doublereal *hcol, doublereal *hv, 
	doublereal *gcon, doublereal *gobj, doublereal *fcon2, doublereal *
	gcon2, doublereal *gobj2, doublereal *ycon, doublereal *w, doublereal 
	*w1, char *cu, integer *lencu, integer *iu, integer *leniu, 
	doublereal *ru, integer *lenru, char *cw, integer *lencw, integer *iw,
	 integer *leniw, doublereal *rw, integer *lenrw, ftnlen cu_len, 
	ftnlen cw_len)
{
    /* Initialized data */

    static char lwrong[4] = "bad?";
    static char lright[4] = "ok  ";

    /* Format strings */
    static char fmt_1600[] = "(\002 -->  The largest discrepancy was\002,1p,"
	    "e12.2,\002  in row and column\002,i6)";
    static char fmt_2001[] = "(\002 Column       x(j)        dx(j)\002,3x"
	    ",\002     Row        Derivative    Difference approxn\002)";
    static char fmt_2101[] = "(i7,1p,e16.8,e10.2,i10,2e18.8,2x,a4)";
    static char fmt_2201[] = "(33x,i10,1pe18.8,e18.8,2x,a4)";
    static char fmt_2501[] = "(\002 All\002,i7,\002  nonzero Hessian element"
	    "s in cols \002,i6,\002  thru\002,i6,\002  seem to be OK.\002)";
    static char fmt_2601[] = "(\002 XXX  There seem to be\002,i6,\002  incor"
	    "rect Hessian elements in cols\002,i6,\002  thru\002,i6)";
    static char fmt_2701[] = "(\002 -->  The largest relative error was\002,"
	    "1p,e12.2,\002   in row\002,i6,\002,  column\002,i6)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__, j, j5, j6, jn;
    static doublereal dx, wj, xj, hij;
    static char key[4];
    static integer irn;
    static doublereal err;
    static char str[120];
    static doublereal eps0, eps5, maxh, fobj2;
    static logical nonlinearobj, nonlinearcon;
    static doublereal abhij;
    static logical cheap;
    extern /* Subroutine */ int dload_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal hijfd;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer nfeas;
    static doublereal emaxh;
    extern doublereal dasum_(integer *, doublereal *, integer *);
    static integer imaxh, jmaxh;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    static logical first;
    static integer nzero;
    static doublereal fdint1;
    extern integer s2coln_(integer *, integer *, integer *);
    extern /* Subroutine */ int s8hxnp_(U_fp, integer *, integer *, integer *,
	     integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, char *, integer *, integer *, integer *, doublereal *, 
	    integer *, char *, integer *, integer *, integer *, doublereal *, 
	    integer *, ftnlen, ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
    static integer righth;
    extern doublereal dnormi_(integer *, doublereal *, integer *);
    static integer wrongh, lvlver, status;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen), s8gprod_(integer *, doublereal *, integer *, integer *
	    , integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, integer *);
    static integer funmode;
    extern /* Subroutine */ int s7checkp_(integer *, doublereal *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, integer *);

    /* Fortran I/O blocks */
    static icilist io___140 = { 0, str, 0, fmt_1600, 120, 1 };
    static icilist io___141 = { 0, str, 0, fmt_2001, 120, 1 };
    static icilist io___155 = { 0, str, 0, fmt_2101, 120, 1 };
    static icilist io___156 = { 0, str, 0, fmt_2201, 120, 1 };
    static icilist io___157 = { 0, str, 0, fmt_2501, 120, 1 };
    static icilist io___158 = { 0, str, 0, fmt_2601, 120, 1 };
    static icilist io___159 = { 0, str, 0, fmt_2701, 120, 1 };


/* ================================================================= */
/* s7checkHv  verifies the Hessian of the Lagrangian using */
/* finite differences. */

/* 10 May 2003: First version of s7checkHv. */
/* 29 Dec 2008: Estimate the number of Hessian nonzeros. */
/* 10 Oct 2010: Fixed bug in bad/ok printing */
/* ================================================================= */
/* ----------------------------------------------------------------- */
    /* Parameter adjustments */
    --x1;
    --x;
    --w1;
    --ycon;
    --fcon2;
    --w;
    --hv;
    --bu;
    --bl;
    --gobj2;
    --gobj;
    --indj;
    --locj;
    --gcon2;
    --gcon;
    --locg;
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
/* ----------------------------------------------------------------- */
    lvlver = iw[78];
/* Verify level */
    if (lvlver < 0) {
	return 0;
    }
    j5 = iw[114];
/* start col for Hessian checking */
    j6 = iw[115];
/* stop  col for Hessian checking */
    eps0 = rw[2];
/* eps**(4/5) */
    eps5 = rw[7];
/* eps**(1/5) */
    fdint1 = rw[76];
/* The problem functions are called to provide gradients only. */
/* (1) forwrd diff. interval */
    *iexit = 0;
/* Cheap  => do the Cheap test only */
    cheap = lvlver == 0;
    if (cheap) {
	snprnt_(&c__11, " Cheap test of user-supplied  Lagrangian Hessian...",
		 &iw[1], leniw, (ftnlen)51);
    } else {
	snprnt_(&c__11, " Verification of user-supplied second derivatives.", 
		&iw[1], leniw, (ftnlen)50);
    }
/* ----------------------------------------------------------------- */
/* Compute functions at a short step along  w. */
/* ----------------------------------------------------------------- */
    dcopy_(n, &x[1], &c__1, &x1[1], &c__1);
    wj = 1. / *nnh;
    i__1 = *nnh;
    for (j = 1; j <= i__1; ++j) {
	w[j] = wj;
	x1[j] = wj;
	wj = -wj * .99999;
    }
/* If needed, fix w to give a feasible direction. */
/* If this gives 0, go back to original w and forget feasibility. */
    dx = fdint1 * (dnormi_(nnh, &x[1], &c__1) + 1.);
    s7checkp_(nnh, &bl[1], &bu[1], &x[1], &dx, &w[1], &nfeas);
    if (nfeas == 0) {
	dcopy_(nnh, &x1[1], &c__1, &w[1], &c__1);
    }
/* ----------------------------------------------------------------- */
/* Compute functions at a short step along  w. */
/* ----------------------------------------------------------------- */
    dx = fdint1 * (dasum_(nnh, &x[1], &c__1) + 1.);
    dcopy_(nnh, &x[1], &c__1, &x1[1], &c__1);
    daxpy_(nnh, &dx, &w[1], &c__1, &x1[1], &c__1);
    nonlinearcon = *nnjac > 0;
    nonlinearobj = *nnobj > 0;
    funmode = 2;
    (*funwrapper)(iexit, &funmode, &nonlinearcon, &nonlinearobj, n, negcon, 
	    nncon0, nncon, nnjac, nnh, nnobj0, nnobj, (U_fp)funcon, (U_fp)
	    funobj, (U_fp)userhv, &x1[1], &ycon[1], nej, nlocj, &locj[1], &
	    indj[1], neh, nloch, &loch[1], &indh[1], &hcol[1], &fcon2[1], &
	    fobj2, &gcon2[1], &gobj2[1], &ycon[1], cu + 8, lencu, &iu[1], 
	    leniu, &ru[1], lenru, cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw,
	     (ftnlen)8, (ftnlen)8);
    if (*iexit != 0) {
	goto L900;
    }
    if (*nnobj > 0) {
/* Compute gL2 - gL */
	dcopy_(nnobj, &gobj2[1], &c__1, &w1[1], &c__1);
	daxpy_(nnobj, &c_b235, &gobj[1], &c__1, &w1[1], &c__1);
    }
    nzero = *nnh - *nnobj;
    if (nzero > 0) {
	dload_(&nzero, &c_b29, &w1[*nnobj + 1], &c__1);
    }
    if (*nncon > 0) {
	s8gprod_(&c__1, &eps0, nej, nlocj, &locj[1], &indj[1], negcon, nlocg, 
		&locg[1], &gcon2[1], &c_b235, &ycon[1], nncon, &c_b242, &w1[1]
		, nnjac);
	s8gprod_(&c__1, &eps0, nej, nlocj, &locj[1], &indj[1], negcon, nlocg, 
		&locg[1], &gcon[1], &c_b242, &ycon[1], nncon, &c_b242, &w1[1],
		 nnjac);
    }
/* ----------------------------------------------------------------- */
/* Cheap test for the Lagrangian Hessian. */
/* ----------------------------------------------------------------- */
/* Set   Hv =  H*w */
    status = 0;
    s8hxnp_((U_fp)userhv, nnh, neh, nloch, &loch[1], &indh[1], &hcol[1], &w[1]
	    , &hv[1], &status, cu + 8, lencu, &iu[1], leniu, &ru[1], lenru, 
	    cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)8, (ftnlen)8)
	    ;
/* Set Hv =  (gL2 - gL)'w/dx - H*w.  This should be small. */
    i__1 = *nnh;
    for (j = 1; j <= i__1; ++j) {
	hv[j] = w1[j] / dx - hv[j];
    }
    jmaxh = idamax_(nnh, &hv[1], &c__1);
    maxh = w1[jmaxh] / dx;
    emaxh = (d__1 = hv[jmaxh], abs(d__1)) / (abs(maxh) + 1.);
    if (emaxh <= .1) {
	snprnt_(&c__11, " The Lagrangian Hessian seems to be OK.", &iw[1], 
		leniw, (ftnlen)39);
    } else {
	snprnt_(&c__11, " XXX  The Lagrangian Hessian seems to be incorrect.",
		 &iw[1], leniw, (ftnlen)51);
    }
    s_wsfi(&io___140);
    do_fio(&c__1, (char *)&emaxh, (ftnlen)sizeof(doublereal));
    i__1 = s2coln_(&jmaxh, leniw, &iw[1]);
    do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
    e_wsfi();
    snprnt_(&c__11, str, &iw[1], leniw, (ftnlen)120);
    snprnt_(&c__1, " ", &iw[1], leniw, (ftnlen)1);
    if (cheap) {
	goto L900;
    }
/* ----------------------------------------------------------------- */
/* Proceed with the verification of individual columns. */
/* ----------------------------------------------------------------- */
    snprnt_(&c__11, " ", &iw[1], leniw, (ftnlen)1);
    if (j5 <= j6) {
	s_wsfi(&io___141);
	e_wsfi();
	snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)120);
    }
    wrongh = 0;
    righth = 0;
    emaxh = -1.;
    i__1 = *nnh;
    for (j = 1; j <= i__1; ++j) {
	jn = s2coln_(&j, leniw, &iw[1]);
	if (jn >= j5 && jn <= j6) {
	    dload_(nnh, &c_b29, &w[1], &c__1);
	    xj = x[j];
	    dx = fdint1 * (abs(xj) + 1.);
	    w[j] = 1.;
	    if (bl[j] < bu[j] && xj >= bu[j]) {
		dx = -dx;
	    }
	    x[j] = xj + dx;
	    funmode = 2;
	    (*funwrapper)(iexit, &funmode, &nonlinearcon, &nonlinearobj, n, 
		    negcon, nncon0, nncon, nnjac, nnh, nnobj0, nnobj, (U_fp)
		    funcon, (U_fp)funobj, (U_fp)userhv, &x[1], &ycon[1], nej, 
		    nlocj, &locj[1], &indj[1], neh, nloch, &loch[1], &indh[1],
		     &hcol[1], &fcon2[1], &fobj2, &gcon2[1], &gobj2[1], &ycon[
		    1], cu + 8, lencu, &iu[1], leniu, &ru[1], lenru, cw + 8, 
		    lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)8, (ftnlen)8)
		    ;
	    if (*iexit != 0) {
		goto L900;
	    }
	    if (*nnobj > 0) {
/* Compute gL2 - gL */
		dcopy_(nnobj, &gobj2[1], &c__1, &w1[1], &c__1);
		daxpy_(nnobj, &c_b235, &gobj[1], &c__1, &w1[1], &c__1);
	    }
	    nzero = *nnh - *nnobj;
	    if (nzero > 0) {
		dload_(&nzero, &c_b29, &w1[*nnobj + 1], &c__1);
	    }
	    if (*nncon > 0) {
		s8gprod_(&c__1, &eps0, nej, nlocj, &locj[1], &indj[1], negcon,
			 nlocg, &locg[1], &gcon2[1], &c_b235, &ycon[1], nncon,
			 &c_b242, &w1[1], nnjac);
		s8gprod_(&c__1, &eps0, nej, nlocj, &locj[1], &indj[1], negcon,
			 nlocg, &locg[1], &gcon[1], &c_b242, &ycon[1], nncon, 
			&c_b242, &w1[1], nnjac);
	    }
/* Set w1 =  (gL2 - gL)/dx (which is approximately H*w). */
	    d__1 = 1. / dx;
	    dscal_(nnh, &d__1, &w1[1], &c__1);
	    x[j] = xj;
/* Set   Hv =  H*w */
	    status = 0;
	    s8hxnp_((U_fp)userhv, nnh, neh, nloch, &loch[1], &indh[1], &hcol[
		    1], &w[1], &hv[1], &status, cu + 8, lencu, &iu[1], leniu, 
		    &ru[1], lenru, cw + 8, lencw, &iw[1], leniw, &rw[1], 
		    lenrw, (ftnlen)8, (ftnlen)8);
	    first = TRUE_;
	    i__2 = *nnh;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		irn = s2coln_(&i__, leniw, &iw[1]);
		hij = hv[i__];
		hijfd = w1[i__];
		abhij = abs(hij);
		err = (d__1 = hijfd - hij, abs(d__1)) / (abhij + 1.);
		if (emaxh < err) {
		    emaxh = err;
		    imaxh = irn;
		    jmaxh = jn;
		}
		if (err <= eps5) {
		    s_copy(key, lright, (ftnlen)4, (ftnlen)4);
		    ++righth;
		} else {
		    s_copy(key, lwrong, (ftnlen)4, (ftnlen)4);
		    ++wrongh;
		}
		if (abhij + err > eps0) {
		    if (first) {
			s_wsfi(&io___155);
			do_fio(&c__1, (char *)&jn, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&xj, (ftnlen)sizeof(doublereal))
				;
			do_fio(&c__1, (char *)&dx, (ftnlen)sizeof(doublereal))
				;
			do_fio(&c__1, (char *)&irn, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&hij, (ftnlen)sizeof(doublereal)
				);
			do_fio(&c__1, (char *)&hijfd, (ftnlen)sizeof(
				doublereal));
			do_fio(&c__1, key, (ftnlen)4);
			e_wsfi();
			snprnt_(&c__11, str, &iw[1], leniw, (ftnlen)120);
			first = FALSE_;
		    } else {
			s_wsfi(&io___156);
			do_fio(&c__1, (char *)&irn, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&hij, (ftnlen)sizeof(doublereal)
				);
			do_fio(&c__1, (char *)&hijfd, (ftnlen)sizeof(
				doublereal));
			do_fio(&c__1, key, (ftnlen)4);
			e_wsfi();
			snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)120);
		    }
		}
	    }
	}
/* jN .ge. j5  .and.  jN .le. j6 */
    }
/* ----------------------------------------------------------------- */
/* Print a summary. */
/* ----------------------------------------------------------------- */
    if (j5 <= j6) {
	if (wrongh == 0) {
	    s_wsfi(&io___157);
	    do_fio(&c__1, (char *)&righth, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&j5, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&j6, (ftnlen)sizeof(integer));
	    e_wsfi();
	} else {
	    s_wsfi(&io___158);
	    do_fio(&c__1, (char *)&wrongh, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&j5, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&j6, (ftnlen)sizeof(integer));
	    e_wsfi();
	}
	snprnt_(&c__11, str, &iw[1], leniw, (ftnlen)120);
	s_wsfi(&io___159);
	do_fio(&c__1, (char *)&emaxh, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&imaxh, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&jmaxh, (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__11, str, &iw[1], leniw, (ftnlen)120);
	snprnt_(&c__1, " ", &iw[1], leniw, (ftnlen)1);
	if (emaxh >= 1.) {
	    *iexit = 54;
/* Bad second derivatives. */
	}
    }
/* ----------------------------------------------------------------- */
/* Print a message if the problem functions are undefined at x. */
/* ----------------------------------------------------------------- */
L900:
    if (*iexit < 0) {
	snprnt_(&c__11, " XXX  Unable to complete Hessian check.", &iw[1], 
		leniw, (ftnlen)39);
	*iexit = 0;
    }
    return 0;
} /* s7checkhv_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* of s7checkHv */
/* Subroutine */ int s7checkp_(integer *n, doublereal *bl, doublereal *bu, 
	doublereal *x, doublereal *dx, doublereal *p, integer *nfeas)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j;
    static doublereal b1, b2, pj, xj, xnew;

/*     ================================================================== */
/*     s7checkp  checks that x + dx*p is feasible. */
/*     It is used by s7checkG for the cheap gradient checks. */

/*     Original:    Looked at the sign of p for variables on a bound. */
/*     13 Mar 1992: dx added as a parameter to make certain that */
/*                  x + dx*p does not lie outside the bounds. */
/*                  p may be altered to achieve this. */
/*     09 Aug 1992: First version based on Minos routine m7chkg. */
/*     27 Feb 2000: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --p;
    --x;
    --bu;
    --bl;

    /* Function Body */
    *nfeas = 0;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	xj = x[j];
	b1 = bl[j];
	b2 = bu[j];
	if (b1 == b2) {
	    p[j] = 0.;
	}
	if (p[j] != 0.) {
/*           x(j) is not fixed, so there is room to move. */
/*           If xj + dx*pj is beyond one bound, reverse pj */
/*           and make sure it is not beyond the other. */
/*           Give up and use set pj = zero if both bounds are too close. */
	    pj = p[j];
	    xnew = xj + *dx * pj;
	    if (pj > 0.) {
		if (xnew > b2) {
		    pj = -pj;
		    xnew = xj + *dx * pj;
		    if (xnew < b1) {
			pj = 0.;
		    }
		}
	    } else {
		if (xnew < b1) {
		    pj = -pj;
		    xnew = xj + *dx * pj;
		    if (xnew > b2) {
			pj = 0.;
		    }
		}
	    }
	    p[j] = pj;
	    if (pj != 0.) {
		++(*nfeas);
	    }
	}
    }
    return 0;
} /* s7checkp_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s7checkp */
/* Subroutine */ int s7fixx_(integer *n, integer *nfixed, logical *perturball,
	 doublereal *boundpert, doublereal *bl, doublereal *bu, doublereal *x)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static integer j;
    static doublereal absxj;

/*     ================================================================== */
/*     s7fixX  ensures that variables satisfy their simple bounds. */

/*     If PerturbAll = true, the bounds for the fixed variables are */
/*     relaxed boundPert. */

/*     On exit, nFixed = the number of fixed variables. */

/*     27 Sep 2003: First version of s7fixX. */
/*     27 Sep 2003: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --x;
    --bu;
    --bl;

    /* Function Body */
    *nfixed = 0;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	d__1 = x[j], d__2 = bl[j];
	x[j] = max(d__1,d__2);
/* Computing MIN */
	d__1 = x[j], d__2 = bu[j];
	x[j] = min(d__1,d__2);
	absxj = (d__1 = x[j], abs(d__1)) + 1.;
	if (bl[j] == bu[j]) {
	    ++(*nfixed);
	    if (*perturball) {
		bl[j] -= *boundpert * absxj;
		bu[j] += *boundpert * absxj;
	    }
	}
    }
    return 0;
} /* s7fixx_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s7fixX */
/* Subroutine */ int s7jac_(integer *iexit, S_fp userfg, logical *perturball, 
	integer *nf, integer *n, doublereal *infbnd, integer *imaxj, 
	doublereal *emaxj, integer *iafun, integer *javar, integer *lena, 
	integer *nea, doublereal *a, integer *igfun, integer *jgvar, integer *
	leng, integer *neg, integer *rowtype, integer *coltype, doublereal *x,
	 doublereal *xlow, doublereal *xupp, doublereal *f, doublereal *g, 
	doublereal *w, doublereal *y, doublereal *z__, doublereal *fw, 
	doublereal *fy, doublereal *fz, doublereal *gcolw, doublereal *gcoly, 
	doublereal *gcolz, integer *iw, integer *leniw, char *cu, integer *
	lencu, integer *iu, integer *leniu, doublereal *ru, integer *lenru, 
	ftnlen cu_len)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k, k0;
    static doublereal wj, yj, zj, dwj, dyj, dzj, eps, tol;
    static integer nea0;
    static doublereal dyj0, eps0, jwij, jyij, jzij;
    extern /* Subroutine */ int s8callstatus_(integer *, integer *, integer *)
	    ;
    extern doublereal s1eps_(void);
    static integer needf, needg;
    extern /* Subroutine */ int dload_(integer *, doublereal *, doublereal *, 
	    integer *), iload_(integer *, integer *, integer *, integer *);
    static integer seeds[3];
    static doublereal gmaxj, xnorm, ynorm, fdint1;
    extern /* Subroutine */ int s7pert_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), s7step_(
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *), s7fixx_(integer *, integer *, 
	    logical *, doublereal *, doublereal *, doublereal *, doublereal *)
	    ;
    extern integer idamax_(integer *, doublereal *, integer *);
    static doublereal damper;
    static integer nfixed;
    extern doublereal dnormi_(integer *, doublereal *, integer *);
    static integer status;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen), s6userf_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, S_fp, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, char *, integer 
	    *, integer *, integer *, doublereal *, integer *, integer *, 
	    integer *, ftnlen);

/*     ================================================================== */
/*     s7Jac   does the work for snJac */

/*     26 Oct 2002: First version. */
/*     27 Sep 2003: More thorough check for feasibility. */
/*     22 Apr 2007: Arguments of s6Userf changed. */
/*     15 Jun 2008: Call-status implemented correctly. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --gcolz;
    --gcoly;
    --gcolw;
    --fz;
    --fy;
    --fw;
    --f;
    --rowtype;
    --z__;
    --y;
    --w;
    --xupp;
    --xlow;
    --x;
    --coltype;
    --a;
    --javar;
    --iafun;
    --g;
    --jgvar;
    --igfun;
    --iw;
    cu -= 8;
    --iu;
    --ru;

    /* Function Body */
    eps = s1eps_();
/* machine precision */
    eps0 = pow_dd(&eps, &c_b3);
    fdint1 = sqrt(eps);
/* Computing MAX */
    d__1 = 1., d__2 = dnormi_(n, &x[1], &c__1);
    xnorm = max(d__1,d__2);
    tol = eps0 * xnorm;
    *neg = 0;
    *nea = 0;
    *imaxj = 0;
    *emaxj = 0.;
/*     Move x inside its bounds. */
/*     If requested, perturb bounds for fixed variables. */
    s7fixx_(n, &nfixed, perturball, &c_b242, &xlow[1], &xupp[1], &x[1]);
/* Determine the status of this call. */
    s8callstatus_(&status, &iw[1], leniw);
    needf = 1;
    needg = 0;
    (*userfg)(&status, n, &x[1], &needf, nf, &f[1], &needg, leng, &g[1], cu + 
	    8, lencu, &iu[1], leniu, &ru[1], lenru, (ftnlen)8);
    if (status < 0) {
	goto L999;
    }
/*     Define random points  w, y and z near the user-supplied x. */
    seeds[0] = 5872;
    seeds[1] = 8584;
    seeds[2] = 4879;
/*     Create a random feasible perturbation of the variables. */
    s7pert_(n, seeds, &x[1], &xlow[1], &xupp[1], &xnorm, &w[1]);
    s7pert_(n, seeds, &x[1], &xlow[1], &xupp[1], &xnorm, &y[1]);
    s7pert_(n, seeds, &x[1], &xlow[1], &xupp[1], &xnorm, &z__[1]);
    iload_(n, &c__0, &coltype[1], &c__1);
    iload_(nf, &c__0, &rowtype[1], &c__1);
/*     Evaluate F at  w, y and z. */
    s6userf_(&status, n, &x[1], &w[1], &damper, (S_fp)userfg, &needf, nf, &fw[
	    1], &needg, leng, &g[1], cu + 8, lencu, &iu[1], leniu, &ru[1], 
	    lenru, &iw[1], leniw, (ftnlen)8);
    if (status < 0) {
	goto L999;
    }
    s6userf_(&status, n, &x[1], &y[1], &damper, (S_fp)userfg, &needf, nf, &fy[
	    1], &needg, leng, &g[1], cu + 8, lencu, &iu[1], leniu, &ru[1], 
	    lenru, &iw[1], leniw, (ftnlen)8);
    if (status < 0) {
	goto L999;
    }
    s6userf_(&status, n, &x[1], &z__[1], &damper, (S_fp)userfg, &needf, nf, &
	    fz[1], &needg, leng, &g[1], cu + 8, lencu, &iu[1], leniu, &ru[1], 
	    lenru, &iw[1], leniw, (ftnlen)8);
    if (status < 0) {
	goto L999;
    }
    status = 0;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*        --------------------------------------------------------------- */
/*        Compute functions at  w + dwj*e_j, y + dyj*e_j  and  z + dzj*e_j. */
/*        --------------------------------------------------------------- */
	wj = w[j];
	yj = y[j];
	zj = z__[j];
/*        Evaluate F at w + dwj*e_j. */
	d__1 = abs(wj);
	s7step_(infbnd, &xlow[j], &xupp[j], &wj, &d__1, &dwj);
	if (dwj == 0.) {
	    dload_(nf, &c_b29, &gcolw[1], &c__1);
	} else {
	    w[j] = wj + dwj;
	    (*userfg)(&status, n, &w[1], &needf, nf, &f[1], &needg, leng, &g[
		    1], cu + 8, lencu, &iu[1], leniu, &ru[1], lenru, (ftnlen)
		    8);
	    if (status < 0) {
		goto L999;
	    }
	    i__2 = *nf;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		gcolw[i__] = (f[i__] - fw[i__]) / dwj;
	    }
	}
/*        Evaluate F at y + dyj*e_j. */
	d__1 = abs(yj) * 2.;
	s7step_(infbnd, &xlow[j], &xupp[j], &yj, &d__1, &dyj);
	if (dyj == 0.) {
	    dload_(nf, &c_b29, &gcoly[1], &c__1);
	} else {
	    y[j] = yj + dyj;
	    (*userfg)(&status, n, &y[1], &needf, nf, &f[1], &needg, leng, &g[
		    1], cu + 8, lencu, &iu[1], leniu, &ru[1], lenru, (ftnlen)
		    8);
	    if (status < 0) {
		goto L999;
	    }
	    i__2 = *nf;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		gcoly[i__] = (f[i__] - fy[i__]) / dyj;
	    }
	}
/*        Evaluate F at z + dzj*e_j. */
	d__1 = abs(zj) * 3.;
	s7step_(infbnd, &xlow[j], &xupp[j], &zj, &d__1, &dzj);
	if (dzj == 0.) {
	    dload_(nf, &c_b29, &gcolz[1], &c__1);
	} else {
	    z__[j] = zj + dzj;
	    (*userfg)(&status, n, &z__[1], &needf, nf, &f[1], &needg, leng, &
		    g[1], cu + 8, lencu, &iu[1], leniu, &ru[1], lenru, (
		    ftnlen)8);
	    if (status < 0) {
		goto L999;
	    }
	    i__2 = *nf;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		gcolz[i__] = (f[i__] - fz[i__]) / dzj;
	    }
	}
/*        --------------------------------------------------------------- */
/*        Figure out which elements derivatives are zero or constant */
/*        --------------------------------------------------------------- */
	i__2 = *nf;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    jwij = gcolw[i__];
	    jyij = gcoly[i__];
	    jzij = gcolz[i__];
	    if (abs(jwij) + abs(jyij) + abs(jzij) > tol) {
/*              Nonzero Jacobian element, to be stored in A or G. */
		if ((d__1 = jwij - jyij, abs(d__1)) <= tol && (d__2 = jwij - 
			jzij, abs(d__2)) <= tol) {
/*                 Constant  Jij */
		    ++(*nea);
		    if (*nea > *lena) {
			snprnt_(&c__1, " Increase lenA.", &iw[1], leniw, (
				ftnlen)15);
			*iexit = 91;
			goto L999;
		    }
		    iafun[*nea] = i__;
		    javar[*nea] = j;
		    a[*nea] = (jwij + jzij) / 2.;
		} else {
/*                 Nonlinear Jij */
		    ++(*neg);
		    if (*neg > *leng) {
			snprnt_(&c__1, " Increase lenG.", &iw[1], leniw, (
				ftnlen)15);
			*iexit = 91;
			goto L999;
		    }
		    igfun[*neg] = i__;
		    jgvar[*neg] = j;
		    rowtype[i__] = 1;
		    coltype[j] = 1;
		}
	    }
	}
	w[j] = wj;
	y[j] = yj;
	z__[j] = zj;
    }
/*     ----------------------------------------------------------------- */
/*     Constant elements in nonlinear rows and columns must be treated */
/*     as nonlinear. */
/*     ----------------------------------------------------------------- */
    k = 1;
    nea0 = *nea;
    i__1 = nea0;
    for (k0 = 1; k0 <= i__1; ++k0) {
	i__ = iafun[k0];
	j = javar[k0];
	if (rowtype[i__] == 1 && coltype[j] == 1) {
	    --(*nea);
	    ++(*neg);
	    if (*neg > *leng) {
		snprnt_(&c__1, " Increase lenG.", &iw[1], leniw, (ftnlen)15);
		*iexit = 91;
		goto L999;
	    }
	    igfun[*neg] = i__;
	    jgvar[*neg] = j;
	} else {
	    if (k < k0) {
		iafun[k] = iafun[k0];
		javar[k] = javar[k0];
		a[k] = a[k0];
	    }
	    ++k;
	}
    }
    if (*neg < *n * *nf) {
/*        =============================================================== */
/*        J has some constant elements.  Better check we got it right. */

/*        Compare J*p with (G + A)*p, everything computed at y */
/*        =============================================================== */
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
/*           ------------------------------------------------------------ */
/*           Compute functions at  y + dyj*e_j. */
/*           ------------------------------------------------------------ */
	    yj = y[j];
/*           Evaluate F at y + dyj*e_j. */
	    dyj0 = abs(yj) * fdint1;
	    s7step_(infbnd, &xlow[j], &xupp[j], &yj, &dyj0, &dyj);
	    if (dyj == 0.) {
		dload_(nf, &c_b29, &gcoly[1], &c__1);
	    } else {
		y[j] = yj + dyj;
		(*userfg)(&status, n, &y[1], &needf, nf, &f[1], &needg, leng, 
			&g[1], cu + 8, lencu, &iu[1], leniu, &ru[1], lenru, (
			ftnlen)8);
		if (status < 0) {
		    goto L999;
		}
		i__2 = *nf;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    gcoly[i__] = (f[i__] - fy[i__]) / dyj;
		}
	    }
	    i__2 = *neg;
	    for (k = 1; k <= i__2; ++k) {
		if (jgvar[k] == j) {
		    i__ = igfun[k];
		    g[k] = gcoly[i__];
		}
	    }
	    y[j] = yj;
	}
/*        --------------------------------------------------------------- */
/*        Compute a new random feasible perturbation z and compute */
/*        the functions at  w = y + z. */
/*        --------------------------------------------------------------- */
	ynorm = dnormi_(n, &y[1], &c__1) + 1.;
	d__1 = fdint1 * ynorm;
	s7pert_(n, seeds, &y[1], &xlow[1], &xupp[1], &d__1, &w[1]);
/*        Evaluate F at w + z. */
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    z__[j] = w[j] - y[j];
	}
	(*userfg)(&status, n, &w[1], &needf, nf, &f[1], &needg, leng, &g[1], 
		cu + 8, lencu, &iu[1], leniu, &ru[1], lenru, (ftnlen)8);
	if (status < 0) {
	    goto L999;
	}
/*        --------------------------------------------------------------- */
/*        Cheap test for the  Jacobian. */
/*        --------------------------------------------------------------- */
/*        Compute  F - (Fy + (G + A)*z).  This should be small. */
	i__1 = *nf;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    fw[i__] = f[i__] - fy[i__];
	}
	i__1 = *neg;
	for (k = 1; k <= i__1; ++k) {
	    i__ = igfun[k];
	    j = jgvar[k];
	    fw[i__] -= g[k] * z__[j];
	}
	i__1 = *nea;
	for (k = 1; k <= i__1; ++k) {
	    i__ = iafun[k];
	    j = javar[k];
	    fw[i__] -= a[k] * z__[j];
	}
	*imaxj = idamax_(nf, &fw[1], &c__1);
	gmaxj = f[*imaxj] - fy[*imaxj];
	*emaxj = (d__1 = fw[*imaxj], abs(d__1)) / (abs(gmaxj) + 1.);
    }
    return 0;
L999:
    if (status == -1) {
	*iexit = 63;
/* unable to proceed into undefined regio */
    } else if (status <= -2) {
	*iexit = 71;
/* terminated during function evaluation */
    }
    return 0;
} /* s7jac_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s7Jac */
/* Subroutine */ int s7pert_(integer *n, integer *seeds, doublereal *x, 
	doublereal *bl, doublereal *bu, doublereal *order, doublereal *xpert)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static integer j, direction;
    static doublereal theta;
    extern /* Subroutine */ int ddrand_(integer *, doublereal *, integer *, 
	    integer *);

/*     ================================================================== */
/*     s7pert finds a feasible perturbation xPert of a feasible point x. */
/*     The order of the perturbation will be at most "order" and the */
/*     elements of xPert  will lie between bl and bu for bl < bu. */

/*     seeds(1:3) are seed values for the random number generator. */

/*     26 Oct 2002: First version based on snadiopt routine */
/*     27 Sep 2003: Current version of s7pert. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
/*     Get a vector of pseudo-random numbers between 0 and 1 */
    /* Parameter adjustments */
    --xpert;
    --bu;
    --bl;
    --x;
    --seeds;

    /* Function Body */
    ddrand_(n, &xpert[1], &c__1, &seeds[1]);
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*        If direction .eq. free, then x(j) is not at it's bound. */
	direction = 2;
	if (x[j] == bl[j]) {
/*           x(j) is on its lower bound, we must perturb up. */
	    direction = 1;
	}
	if (x[j] == bu[j]) {
/*           x(j) is on its upper bound... */
	    if (direction == 2) {
/*              ... but not on its lower bound.  Perturb down. */
		direction = -1;
	    } else {
/*              ... and is also on its lower bound, and so it must be */
/*              a fixed variable.  It is not perturbed. */
		direction = 0;
	    }
	}
	if (direction == 2) {
/*           If the variable is not on either bound, choose randomly */
/*           whether to go up or down. */
	    if (xpert[j] > .5) {
		xpert[j] = (xpert[j] - .5) * 2.;
		direction = 1;
	    } else {
		xpert[j] *= 2.;
		direction = -1;
	    }
	}
	theta = xpert[j] * .5;
	if (direction == -1) {
/*           Perturb between 1/8 and 1/4 of the way to the lower bound, */
/*           but no more than "order". */
/* Computing MAX */
	    d__1 = (x[j] * 3. + bl[j]) / 4., d__2 = x[j] - *order;
	    xpert[j] = theta * x[j] + (1. - theta) * max(d__1,d__2);
	} else if (direction == 1) {
/*           Perturb between 1/8 and 1/4 of the way to the upper bound, */
/*           but no more than "order". */
/* Computing MIN */
	    d__1 = (x[j] * 3. + bu[j]) / 4., d__2 = x[j] + *order;
	    xpert[j] = theta * x[j] + (1. - theta) * min(d__1,d__2);
	} else {
/*           direction .eq. fixed */
/*           xPert(j) = theta*x(j) + (one - theta)*order */
	    xpert[j] = x[j];
	}
    }
    return 0;
} /* s7pert_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s7pert */
/* Subroutine */ int s7step_(doublereal *infbnd, doublereal *blj, doublereal *
	buj, doublereal *xj, doublereal *dxj0, doublereal *dxj)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Local variables */
    static integer move;

/*     ================================================================== */
/*     s7step finds a feasible direction of perturbation dxj for xj. */

/*     dxj0  is a suggested value for the step. */

/*     26 Oct 2002: First version. */
/*     27 Sep 2003: Current version of s7step. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
/*     If move .eq. free, then xj is not at one of it's bounds. */
    move = 2;
    if (*xj == *blj) {
/*        xj is on its lower bound, we must perturb up. */
	move = 1;
    }
    if (*xj == *buj) {
/*        xj is on its upper bound... */
	if (move == 2) {
/*           ... but not on its lower bound.  Perturb down. */
	    move = -1;
	} else {
/*           ... and is also on its lower bound, and so it must be */
/*           a fixed variable.  It is not perturbed. */
	    move = 0;
	}
    }
    if (move == 2) {
/*        If the variable is not on either bound, choose randomly */
/*        whether to go up or down. */
	if (*xj > .5) {
	    move = 1;
	} else {
	    move = -1;
	}
    }
    *dxj = *dxj0;
    if (move == 1) {
	if (*buj < *infbnd) {
/* Computing MIN */
	    d__1 = *dxj, d__2 = *buj - *xj;
	    *dxj = min(d__1,d__2);
	}
    } else if (move == -1) {
	*dxj = -(*dxj);
	if (*blj > -(*infbnd)) {
/* Computing MAX */
	    d__1 = *dxj, d__2 = *blj - *xj;
	    *dxj = max(d__1,d__2);
	}
    } else if (move == 0) {
	*dxj = 0.;
    }
    return 0;
} /* s7step_ */

