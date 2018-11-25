/****************************************************************
Copyright (C) 1997, 1998, 2001 Lucent Technologies
All Rights Reserved

Permission to use, copy, modify, and distribute this software and
its documentation for any purpose and without fee is hereby
granted, provided that the above copyright notice appear in all
copies and that both that the copyright notice and this
permission notice and warranty disclaimer appear in supporting
documentation, and that the name of Lucent or any of its entities
not be used in advertising or publicity pertaining to
distribution of the software without specific, written prior
permission.

LUCENT DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS.
IN NO EVENT SHALL LUCENT OR ANY OF ITS ENTITIES BE LIABLE FOR ANY
SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER
IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION,
ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
THIS SOFTWARE.
****************************************************************/

#include "signal.h"
#include "ctype.h"
#include "getstub.h"
#include "obj_adj.h"

/* For (defective) Linux on Alpha chips where infinity arithmetic fails: */
/* compile, e.g., with -DFinite_Plinfy=1e154 */
#ifdef Finite_Plinfy
#undef Infinity
#define Infinity Finite_Plinfy
#endif

#define N_FNAMES 100
static char *fnames[N_FNAMES];

#ifdef __cplusplus
extern "C" {
#endif

typedef void fgCon ANSI((fint *MODE, fint *M, fint *N, fint *NJAC,
		real *X, real *F, real *G, fint *NSTATE,
		char *CU, fint *LENCU, fint *IU, fint *LENIU,
		real *RU, fint *LENRU, ftnlen len_cu));

typedef void fgObj ANSI((fint *MODE, fint *N,
		real *X, real *F, real *G, fint *NSTATE,
		char *CU, fint *LENCU, fint *IU, fint *LENIU,
		real *RU, fint *LENRU, ftnlen len_cu));

typedef void fQp ANSI((fint*, real*, real*, fint*, char*, fint*,
		fint*, fint*, real*, fint*, ftnlen));

 static fgCon funcon;
 static fgObj funobj;
 extern void s1file_ ANSI((fint *, fint*, fint*));
 extern int snmemb_ ANSI((fint *info, fint *m, fint *n, fint *ne, fint *negcon,
	fint *nncon, fint *nnjac, fint *nnobj, fint *mincw, fint *miniw,
	fint *minrw, char *cw, fint *lencw, fint *iw, fint *leniw,
	real *rw, fint *lenrw, ftnlen cw_len));
 extern void sninit_  ANSI((fint *iprint, fint *isumm, char *cw,
		fint *lencw, fint *iw, fint *leniw, real *rw,
		fint *lenrw, ftnlen cw_len));
 extern void snopt_ ANSI((char *start, fint *m, fint *n, fint *ne,
		fint *nname, fint *nncon, fint *nnobj, fint *nnjac,
		fint *iobj, real *objadd, char *prob, fgCon *fgcon,
		fgObj *fgobj, real *jcol, fint *indj, fint *locj, real *bl,
		real *bu, char *names, fint *hs, real *x, real *pi,
		real *rc, fint *info, fint *mincw, fint *miniw,
		fint *minrw, fint *ns, fint *ninf, real *sinf, real *obj,
		char *cu, fint *lencu, fint *iu, fint *leniu, real *ru,
		fint *lenru, char *cw, fint *lencw, fint *iw, fint *leniw,
		real *rw, fint *lenrw, ftnlen start_len, ftnlen prob_len,
		ftnlen names_len, ftnlen cu_len, ftnlen cw_len));
 extern void snset_(char *buffer, fint *iprint, fint *isumm, fint *errors,
		char *cw, fint *lencw, fint *iw, fint *leniw, real *rw,
		fint *lenrw, ftnlen buffer_len, ftnlen cw_len);
 extern void sntitle_ ANSI((char*, ftnlen));
 extern void sqopt_ ANSI((char *start, fQp*, fint *m, fint *n,
	fint *ne, fint *nname, fint *lenc, fint *ncolh,
	fint *iobj, real *objadd, char *prob,
	real *acol, fint *inda, fint *loca, real *bl, real *bu,
		real *c, char *names,
	fint *helast, fint *hs, real *x, real *pi, real *rc,
	fint *inform, fint *mincw, fint *miniw, fint *minrw,
	fint *ns, fint *ninf, real *sinf, real *obj,
	char *cu, fint *lencu, fint *iu, fint *leniu, real *ru, fint *lenru,
	char *cw, fint *lencw, fint *iw, fint *leniw, real *rw, fint *lenrw,
	ftnlen start_len, ftnlen prob_len, ftnlen names_len, ftnlen cu_len,
	ftnlen cw_len));

#define Aijtol	rw[94]
#define Itn	iw[420]
#define Objsen	iw[86]
#define Objtru	rw[420]
#define Plinfy	rw[69]
#define iPrint	iw[11]
#define iSumm	iw[12]
#define itnlim	iw[88]
#define lDenJ	iw[104]
#define lenalu	iw[212]
#define lprSol	iw[83]
#define lvlScl	iw[74]
#define lvlTim	iw[181]
#define maxcw	iw[6]
#define maxiw	iw[4]
#define maxrw	iw[2]
#define maxS	iw[52]
#define minlen	iw[162]
#define nMajor	iw[421]
#define nQNmod	iw[380]
#define neGcon	iw[19]
#define nnCon	iw[22]
#define nnJac	iw[20]
#define nnObj	iw[21]
#define objno	iw[107]
#define stkyOp	iw[115]

static fint Objno, iObj;
extern char **xargv;

 typedef struct
Winfo {
	char *cw;
	fint *iw;
	real *rw;
	fint lencw, leniw, lenrw;
	ftnlen cw_len;
	} Winfo;

 static Jmp_buf Jb;
 static Winfo winfo;
 static char Emsgbuf[96], *Errmsg = Emsgbuf;
 static int Errnum = 599; /* should be assigned if needed */
 static int functimes, wantfuncs;
 static int nfuncon, nfungrd, nfunjac, nfunobj;
 static int objrep = 3, outlev = 1, qpcheck = 1;
 static fint i_exit, last_state, memadj, meminc, timing;
 static real Times[6];
 static char Title[32];
#define asl cur_ASL

 static char *usage_msg[] = {
 "  where  stub  is from  `ampl -obstub`  or  `ampl -ogstub`.  Assignments",
 "  have the form  u=filename  or  spec=value  in which the integers u are",
 "  file unit numbers that appear in spec assignments (don't use 5 or 6",
 "  for u) and  spec  is a SNOPT SPECS file keyword phrase, with keywords",
 "  in the phrase separated by _ (underscore).  Use",
 "	outlev=0 for no options echoed on stdout,",
 "	outlev=1 for neither log nor summary file on stdout (default),",
 "	outlev=2 to see summary output on stdout,",
 "	outlev=3 to see detailed (log file) output on stdout, and",
 "	outlev=4 to get log file plus solution on stdout.",
 "  For outlev <= 2, use 7=logfilename to get the log file.  Assignments",
 "  may also be given in $snopt_options .",
		0 };

 typedef struct Phrase Phrase;
 struct
Phrase {
	Phrase *next;
	ftnlen len;
	char val[1]; /*really [len+1]*/
	};

 static Phrase *firstphrase, **lastnext;
 static char *phbuf;
 static size_t phblen;

 static fint
mkey(char *phrase, ftnlen len)
{
	Phrase **ln, *p;
	char *s;
	size_t L, L1;

	ln = lastnext;
	L1 = (len + sizeof(Phrase) + sizeof(char*)) & ~(sizeof(char*)-1);
	if (L1 > phblen) {
		L = L1;
		if (L < 4096)
			L = 4096;
		phbuf = (char*)M1alloc(phblen = L);
		if (!ln)
			ln = &firstphrase;
		}
	*ln = p = (Phrase*)phbuf;
	phbuf += L1;
	lastnext = &p->next;
	p->next = 0;
	memcpy(s = p->val, phrase, (size_t)(p->len = len));
	s[len] = 0;
	return 0;
	}

 static fint
nkey(fint *np, char *fname, ftnlen L)
{
	fint n = *np;
	char *s;

	if (n < 0 || n >= N_FNAMES) {
		printf("Unit number > %d\n", N_FNAMES-1);
		return 1;
		}
	fnames[n] = s = (char*)M1alloc(L+1);
	memcpy(s, fname, L);
	s[L] = 0;
	return 0;
	}

 static char *
set_outlev(Option_Info *oi, keyword *kw, char *s)
{
	char *rv = I_val(oi,kw,s);
	if (!*(int *)kw->info)
		oi->option_echo = 0;
	return rv;
	}

 static char qpcheck_msg[] =
	"whether to pass LPs and QPs to SQOPT rather than SNOPT:\n\
		0 = no\n\
		1 = yes (default)";

 static char objrep_desc[] = "Whether to replace\n\t\t\t\tminimize obj: v;\n\t\t\twith\n\t\t\t\t"
				"minimize obj: f(x)\n\t\t\twhen variable v appears linearly\n\t"
				"\t\tin exactly one constraint of the form\n\t\t\t\t"
				"s.t. c: v >= f(x);\n\t\t\tor\n\t\t\t\ts.t. c: v == f(x);\n"
				"\t\t\tPossible objrep values:\n\t\t\t0 = no\n"
				"\t\t\t1 = yes for v >= f(x)\n"
				"\t\t\t2 = yes for v == f(x)\n"
				"\t\t\t3 = yes in both cases (default)";

 static keyword keywds[] = {
	KW("ftimes",  I_val,  &functimes, "report function eval. times"),
	KW("maxfwd",  IA_val, voffset_of(ASL,p.maxfwd_), "max vars in fwd AD of common exprs (default 5)"),
	KW("meminc",  L_val,  &meminc, "increment to minimum memory allocation"),
	KW("objno",   L_val,  &Objno, "objective number: 0 = none, 1 = first (default)"),
	KW("objrep",  I_val, &objrep, objrep_desc),
	KW("outlev",  set_outlev, &outlev, "output level; 1 = default"),
	KW("qpcheck", I_val, &qpcheck, qpcheck_msg),
	KW("timing",  I_val, &timing,  "report I/O and solution times: 1 = stdout, 2 = stderr, 3 = both"),
	KW("version", Ver_val, 0, "report version"),
	KW("wantsol", WS_val, 0, WS_desc_ASL+5)
	};

 static keyword options[] = {
	KW("f", IK1_val, &wantfuncs, "list available user-defined functions"),
	KW("t", IK1_val, &functimes, "time function evaluations")
	};

 static Option_Info Oinfo = {
	"snopt", Title, "snopt_options", keywds, nkeywds, 1, 0,
	usage_msg, mkey, nkey, options, sizeof(options)/sizeof(keyword),
	20151025 };

 static void
get_Title(VOID)
{
	/* Get title from sntitl and omit excessive blanks. */

	char *s, *se, *t, *t1;
	size_t L, L1;
	int keepblank = 0;

	sntitle_(Title, sizeof(Title));
	s = t = Title;
	se = s + sizeof(Title);
	t1 = 0;
	while(s < se) {
		while(*s > ' ') {
			if (*s == '('/*)*/ && !t1)
				t1 = t;
			*t++ = *s++;
			if (s >= se)
				goto done;
			}
		if (!*s || s+1 >= se)
			break;
		if (keepblank || s[1] == ' ') {
			*t++ = *s++;
			keepblank = 0;
			}
		while(*s <= ' ')
			if (!*s || ++s >= se)
				goto done;
		if (*s == '(')/*)*/
			keepblank = 1;
		}
 done:
	if (!t1)
		t1 = t;
	*t = 0;
	L = t - Title;
	L1 = t1 - Title;
	Oinfo.bsname = M1alloc(L + L1 + 2);
	Oinfo.version = se = Oinfo.bsname + L1 + 1;
	strncpy(Oinfo.bsname, Title, L1);
	Oinfo.bsname[L1] = 0;
	memcpy(se, Title, L);
	se[L] = 0;
	}

 static SufDecl
suftab[] = {
	{ "sstatus", 0, ASL_Sufkind_var, 1 },
	{ "sstatus", 0, ASL_Sufkind_con, 1 }
	};

 static int
envopt(void)
{
	Phrase *p;
	fint inform, *iw;
	int rv;
	static fint lprint, lsumm = 6;

	iw = winfo.iw;
	if (fnames[7])
		iPrint = 7;
	if (fnames[0] && !freopen(fnames[0], "w", Stderr))
		printf("Can't redirect Stderr to %s\n", fnames[0]);
	rv = 0;
	for(p = firstphrase; p; p = p->next) {
		inform = 0;
		snset_(p->val, &lprint, &lsumm, &inform,
		winfo.cw, &winfo.lencw,
		winfo.iw, &winfo.leniw,
		winfo.rw, &winfo.lenrw,
		p->len, (ftnlen)winfo.lencw);
		if (inform)
			++rv;
		}
	fflush(stdout);
	return rv;
	}

#undef asl

 static fint
objmunge(fint M, fint *mp, fint N, fint NZ, fint *nzp, fint *hs,
		real *lb, real *ub, real *A, fint *ha, fint *ja,
		real *objadj, real *x)
{
	ASL *asl;
	Objrep *od, **pod;
	cgrad *cg, *cg1, **cgp, **cgx, *ncg;
	char *h, *he;
	fint *ha1, *ha2, *hae, *iw, *ja1, *kadj, *kadj1;
	fint i, j, k, na0, ne, nlclim, nlvlim, nz, rv, si, vi;
	int *cm, *vmi;
	ograd *og, **ogp;
	real *a1, *a2, *lbe, *rw;
	real aijtol, ninf, pinf, t, t1;

	asl = cur_ASL;
	iw = winfo.iw;
	rw = winfo.rw;
	rv = ne = na0 = 0;
	i = N;
	cgp = 0;
	ogp = 0;
	cm = vmi = 0;
	od = 0;
	if (Objno >= 0 && (pod = asl->i.Or) && (od = pod[Objno])) {
		if (!asl->i.cmap)
			i += asl->i.n_con0;
		if (!asl->i.vmap)
			i += asl->i.n_var0;
		}
	kadj = (fint *)Malloc(i *= sizeof(fint));
	memset((char *)kadj, 0, N*sizeof(fint));
	if (od) {
		vmi = (int*)kadj + N;
		if (!(cm = asl->i.cmap)) {
			cm = vmi;
			vmi += j = asl->i.n_con0;
			while(j > 0) {
				--j;
				cm[j] = j;
				}
			}
		if (asl->i.vmap)
			vmi = get_vminv_ASL(asl);
		else {
			j = asl->i.n_var0;
			while(j > 0) {
				--j;
				vmi[j] = j;
				}
			}
		}
	aijtol = Aijtol;
	nz = NZ;
	nlvlim = nlvc;
	nlclim = nlc;

	/* omit tiny components of A */

	if (aijtol > 0) {
		a1 = a2 = A;
		ha1 = ha2 = ha;
		ja1 = ja;
		for(i = 0; i < N; i++) {
			j = *++ja1;
			hae = ha + j - 1;
			while(ha1 < hae) {
				t = *a1++;
				if ((si = *ha1++) > nlclim
				 && i >= nlvlim
				 && (t < 0 ? -t : t) < aijtol) {
					na0++;
					continue;
					}
				*a2++ = t;
				*ha2++ = si;
				}
			*ja1 = j - na0;
			}
		nz -= na0;
		}

	/* find objective, count gradient components */

	*objadj = 0;
	if (Objno >= 0) {
		if (Objsen == 2)
			Objsen = objtype[Objno] ? -1 : 1;
		if (nl_obj(Objno))
			rv = od ? nlvc : nlvo;
		else
			*objadj = objconst(Objno);
		if (od) {
			if (!(cgp = asl->i.Cgrad0))
				cgp = asl->i.Cgrad_;
			cgp += od->ico;
			for(cg = *cgp; cg; cg = cg->next) {
				if ((i = vmi[cg->varno]) < rv)
					goto ckeep;
				if ((t = cg->coef) < 0)
					t = -t;
				if (t >= aijtol) {
				 ckeep:
					kadj[i] = 1;
					ne++;
					}
				else
					na0++;
				}
			}
		else {
			ogp = &Ograd[Objno];
			for(og = *ogp; og; og = og->next) {
				if (og->varno < rv)
					goto keep;
				else {
					if ((t = og->coef) < 0)
						t = -t;
					if (t >= aijtol) {
					 keep:
						kadj[og->varno] = 1;
						ne++;
						}
					else
						na0++;
					}
				}
			}
		}
	else
		Objsen = 0;
	h = havex0;
	he = h + N;

	pinf = Plinfy;
	ninf = -pinf;

	while(h < he) {
		si = 2;
		if ((t1 = *lb) <= ninf)
			t1 = *lb = ninf;
		else
			si = 0;
		if ((t = *ub) >= pinf)
			t = *ub = pinf;
		else
			si = t <= 0;
		lb++;
		ub++;
		if (*h++ && si != 2) {
			if (*x <= t1) {
				*x = t1;
				si = 4;
				}
			else if (*x >= t) {
				*x = t;
				si = 5;
				}
			}
		*hs++ = si;
		x++;
		}

	/* adjust bounds */

	lbe = lb + M;
	while(lb < lbe) {
		si = 2;
		if ((t = *lb) <= ninf)
			t = ninf;
		else
			si = 0;
		if ((t1 = *ub) >= pinf)
			t1 = pinf;
		else
			si = 0;
		*lb++ = t;
		*ub++ = t1;
		*hs++ = si;
		}

	/* insert objective */

	si = M;
	if (ne) {
		iObj = *mp = ++si;
		*lb = ninf;
		*ub = pinf;
		*hs = 0;
		a1 = A + nz;
		ha1 = ha + nz;
		k = nz;
		nz += ne;
		a2 = A + nz;
		ha2 = ha + nz;
		ja1 = ja + N;
		kadj1 = kadj + N;
		for(;;) {
			*ja1-- += ne;
			if (*--kadj1) {
				*kadj1 = a2 - A;
				*--a2 = 0;
				*--ha2 = si;
				if (!--ne)
					break;
				}
			for(j = *ja1; k >= j; k--) {
				*--a2 = *--a1;
				*--ha2 = *--ha1;
				}
			}
		if (cgp) {
			od->cg0 = cg = *cgp;
			for(k = 0; cg; cg = cg->next) {
				if ((i = kadj[vi = vmi[cg->varno]]) && vi >= rv)
					A[i - 1] = -cg->coef;
				else
					++k;
				}
			if (k) {
				cg1 = (cgrad*)M1alloc(k*sizeof(cgrad));
				for(cg = *cgp; cg; cg = cg->next) {
					if (!kadj[vi = vmi[cg->varno]] || vi < rv) {
						*cg1 = *cg;
						*cgp = cg1;
						cgp = &cg1->next;
						++cg1;
						}
					}
				}
			*cgp = 0;
			}
		else {
			while((og = *ogp)) {
				if ((i = kadj[og->varno]) && og->varno >= rv) {
					A[i - 1] = og->coef;
					*ogp = og->next;
					}
				else
					ogp = &og->next;
				}
			}
		}
	else
		*mp = M;
	*nzp = nz;

	/* adjust for computing Jacobian */

	neGcon = 0;
	if (nlvlim) {
		hae = ha + (ja[nlvlim] - 1);
		for(ha1 = ha, k = 0; ha1 < hae;)
			if (*ha1++ <= nlclim)
				k++;
		neGcon = k;
		if ((cgx = asl->i.Cgrad0)) {
			memset(kadj, 0, N*sizeof(fint));
			j = 0;
			if (cm)
			    for(; j < nlclim; ++j) {
				for(cg = cgx[cm[j]]; cg; cg = cg->next)
					++kadj[vmi[cg->varno]];
				}
			else
			    for(; j < nlclim; ++j) {
				for(cg = cgx[j]; cg; cg = cg->next)
					++kadj[cg->varno];
				}
			for(j = k = 0; j < N; ++j) {
				i = kadj[j];
				kadj[j] += k;
				k += i;
				}
			j = nlclim;
			if (cm)
			    while(j > 0) {
				for(cgp = &cgx[cm[--j]]; (cg = *cgp); cgp = &cg->next) {
					i = vmi[cg->varno];
					if (i >= nlvlim) {
						*cgp = 0;
						break;
						}
					cg->goff = --kadj[i];
					}
				}
			else
			    while(j > 0) {
				for(cgp = &cgx[--j]; (cg = *cgp); cgp = &cg->next) {
					i = cg->varno;
					if (i >= nlvlim) {
						*cgp = 0;
						break;
						}
					cg->goff = --kadj[i];
					}
				}
			for(a1 = A, ha1 = ha; ha1 < hae; a1++) {
				if (*ha1++ <= nlclim)
					*a1 = 0.;
				}
			}
		else {
			free(kadj);
			kadj = 0;
			Cgrad = (cgrad **)M1alloc(nlclim*sizeof(cgrad *));
			memset((char *)Cgrad, 0, nlclim*sizeof(cgrad *));
			for(ha1 = ha, k = 0; ha1 < hae;)
				if (*ha1++ <= nlclim)
					k++;
			ncg = (cgrad *)M1alloc(k*sizeof(cgrad));
			i = k = 0;
			ja1 = ja + 1;
			for(a1 = A, ha1 = ha; ha1 < hae; a1++) {
				if ((j = *ha1++) <= nlclim) {
					cgx = Cgrad + j - 1;
					cg = ncg++;
					cg->next = *cgx;
					*cgx = cg;
					cg->goff = k++;
					j = ha1 - ha;
					while(j >= *ja1) {
						ja1++;
						i++;
						}
					cg->varno = i;
					cg->coef = *a1;
					*a1 = 0.;
					}
				}
			}
		}
	if (kadj)
		free(kadj);
	if (na0) {
		printf("%d coefficients less than aij = %g treated as 0.\n",
			na0, aijtol);
		need_nl = 0;
		}
	c_vars = nlvlim;
	return o_vars = rv;
	}

 static void
time_out(FILE *f)
{
	fprintf(f," SNOPT times:\n read: %10.2f\n solve: %9.2f",
		Times[1] - Times[0], Times[2] - Times[1]);
	fprintf(f, "\n write: %9.2f\n total: %9.2f\n",
		Times[3] - Times[2], Times[3] - Times[0]);
	if (functimes) {
		fprintf(f, "\n");
		if (nfuncon)
			fprintf(f,
	" constraints: %9.2f for %d evaluations, %d Jacobians\n",
				Times[4], nfuncon, nfunjac);
		if (nfunobj)
			fprintf(f,
	" objective: %11.2f for %d evaluations, %d gradients\n",
				Times[5], nfunobj, nfungrd);
		}
	}

 static void
negate(fint M, register real *x)
{
	register real *xe;
	for(xe = x + M; x < xe; x++)
		*x = -*x;
	}

 static fint
vtrans[] = { 0, 3, 2, 0, 1, 0, 0 },
ctrans[] = { 0, 3, 2, 1, 0, 0, 0 };

 static void
hs1_adjust(fint N, int *s, fint *ss, fint *trans)
{
	fint i;

	for(i = 0; i < N; i++)
		ss[i] = trans[s[i]];
	}

 static char *
hs_adjust(ASL *asl, fint N, fint M, fint m1, fint *vss, int *vs, SufDesc *vsd, SufDesc *csd, real *A, fint *ha, fint *ja)
{
	fint *jae, nerror;
	fint *hae, i1, i2, is;
	int nlin;
	real t, *x, *xe;

	if (!(vsd->kind & ASL_Sufkind_input)
	 || !(csd->kind & ASL_Sufkind_input))
		return "Cold";
	hs1_adjust(N, vs, vss, vtrans);
	hs1_adjust(M, vs+N, vss+N, ctrans);
	if (m1 > M)
		vss[M+N] = 3;
	/* Why can't minos compute the initial slacks? */
	/* Then we could eliminate the following mess... */
	if (nlc) {
		nerror = 0;
		conval(X0, x = X0+N, &nerror);
		if (nerror)
			memset(X0+N, 0, nlc*sizeof(real));
		else
			for(xe = x + nlc; x < xe; x++)
				*x = -*x;
		}
	if ((nlin = n_con - nlc) > 0) {
		x = X0 + N + nlc;
		i1 = nlc + 1;
		i2 = n_con;
		memset(x, 0, nlin*sizeof(real));
		hae = ha;
		xe = X0;
		for(jae = ja + N; ja < jae; ja++) {
			hae += ja[1] - ja[0];
			for(t = *xe++; ha < hae; A++)
				if ((is = *ha++) >= i1 && is <= i2)
					x[is - i1] -= t**A;
			}
		}
	return "Warm";
	}

 static void
send_status(fint N, fint *ss, int *s, fint *trans, real *L, real *U)
{
	fint i;

	for(i = 0; i < N; i++)
		if ((((s[i] = trans[ss[i]]) + 1) & ~1) == 4 && L[i] == U[i])
			s[i] = 5;
	}

#undef scream
#undef asl

 static void
scream(char *fmt, ...)
{
	va_list ap;
	va_start(ap, fmt);
	fprintf(Stderr, "snopt: ");
	vfprintf(Stderr, fmt, ap);
	va_end(ap);
	longjmp(Jb.jb, 1);
	}

 static fint *colq, ncolh, *rowq;
 static real *delsq;

 static void
qpHx(fint *ncolH, real *x, real *Hx, fint *nState,
	char *cu, fint *lencu, fint *iu, fint *leniu, real *ru, fint *lenru,
	ftnlen culen)
{
	fint *cq, i, n, *rq, *rqe;
	real *h, t, *xe;

	n = *ncolH;
	if (n != ncolh) {
		scream("qpHx expected ncolH = %d but got %d\n",
			(int)n, (int)ncolh);
		/* suppress warning about unused variables: */
		Not_Used(nState);
		Not_Used(cu);
		Not_Used(lencu);
		Not_Used(iu);
		Not_Used(leniu);
		Not_Used(ru);
		Not_Used(lenru);
		Not_Used(culen);
		}
	for(i = 0; i < n; i++)
		Hx[i] = 0;
	xe = x + n;
	cq = colq;
	rq = rowq;
	h = delsq;
	while(x < xe) {
		t = *x++;
		rqe = rowq + *++cq;
		while(rq < rqe)
			Hx[*rq++] += t**h++;
		}
	}

typedef struct { fint *itnlimp; fint zapped; } Itnzapper;

 static Itnzapper Itnzap;

 static void
itnlimzap(int sig, void *v)
{
	Itnzapper *I = (Itnzapper*)v;
	if (*I->itnlimp > 0) {
		*I->itnlimp = 0;
		I->zapped = 1;
		}
	}

 void
MAIN__(void)
{
	ASL *asl;
	FILE *nl;
	SufDesc *csd, *vsd;
	char **av, *cw, *cw0, *msg1, *stub, *s, *start;
	char Prob[8], buf[64], msg[400];
	fint *ha, *helast, *hs,*iw, *iw0, *ja;
	fint INFORM, M, MXROW, MXCOL, N, N1, NB, NO, NS, NZ;
	fint i, m1, mincw, mincw0, miniw, miniw0, minrw, minrw0, nInf, nint, nz, nzh;
	int flags, qpsol, *varstat;
	real *A, *lb, *pi, *rc, *ub, *y, *z;
	real c0, nlobj, objadj, sInf, *rw, *rw0;
	size_t L0, L2;
	void **vw;
	typedef struct { char *msg; int code, wantobj; } Sol_info;
	Sol_info *SI;
	static int	objrep_flags[4] = { 0, ASL_obj_replace_ineq, ASL_obj_replace_eq,
						ASL_obj_replace_ineq | ASL_obj_replace_eq };
	static Sol_info
		info_optimal	= {/*1*/ "Optimal solution found", 0, 1},
		info_feaspt	= {/*2*/ "Feasible point found", 1, 1},
		info_inacc	= {/*3*/ "Requested accuracy could not be achieved", 102, 1},
		info_infeas_lin = {/*11*/ "Infeasible linear constraints", 202, 1},
		info_infeas_leq = {/*12*/ "Infeasible linear equalities", 203, 1},
		info_infeas_nlm = {/*13*/ "Nonlinear infeasibilities minimized", 204, 1},
		info_infeas_min = {/*14*/ "Infeasibilities minimized", 205, 1},
		info_unbdd_obj  = {/*21*/ "Unbounded objective", 301, 1},
		info_unbdd_con  = {/*22*/ "Constraint violation limit reached", 302, 1},
		info_limit_iter = {/*31*/ "Too many iterations", 400, 1},
		info_limit_maji = {/*32*/ "Major iteration limit reached", 401, 1},
		info_limit_supb = {/*33*/ "The superbasics limit is too small", 402, 1},
		info_cannot_imp = {/*41*/ "The current point cannot be improved", 101, 1},
		info_sing_basis = {/*42*/ "Singular basis", 508, 1},
		info_infeasgenc = {/*43*/ "Cannot satisfy nonlinear constraints", 206, 1},
		info_ill_condnb = {/*44*/ "Ill-conditioned null-space basis", 522, 1},
		info_bad_objgrd = {/*51*/ "Incorrect gradients from funobj", 502, 0},
		info_bad_congrd = {/*52*/ "Incorrect gradients from funcon", 503, 0},
		info_indef_hes  = {/*53*/ "Indefinite Hessian", 530, 0},
		info_bad_obj0	= {	  "Error evaluating initial objective", 515, 0},
		info_bad_init	= {/*61*/ "Error evaluating initial constraint values", 514, 0},
		info_bad_feas	= {/*62*/ "Error evaluating objective at 1st feas. pt.", 513, 1},
		info_bad_prog	= {/*63*/ "Cannot proceed into undefined region", 523, 1},
		info_short_work	= {/*81*/ "Work arrays too small", 550, 0},
		info_short_char = {/*82*/ "Not enough character storage", 551, 0},
		info_short_int  = {/*83*/ "Not enough integer storage", 552, 0},
		info_short_real = {/*84*/ "Not enough real storage", 553, 0},
		info_bad_input  = {/*91*/ "Invalid input argument", 554, 0},
		info_bad_basis  = {/*92*/ "Wrong basis file dimentions", 524, 0},
		info_bad_numbv	= {/*141*/ "Wrong number of basic variables", 555, 0},
		info_basis_err  = {/*142*/ "Error in basis package", 556, 0},
		info_surprise	= { 0, 560, 0},
		info_surprisr	= { "Bug: \"scream\" called", 561, 0 };
	static char names[8] = "        ";
	static fint I0 = 0, I1 = 1, openF = 1;
	static fint	sctrans[] = { 4, 3, 2, 1 },
			svtrans[] = { 3, 4, 2, 1 };

	Times[0] = xectim_();

	asl = ASL_alloc(ASL_read_fg);
	get_Title();
	Times[4] = Times[5] = 0;
	winfo.lencw = winfo.leniw = winfo.lenrw = 500;
	L0 = 500*(sizeof(fint) + sizeof(real) + 8);
	winfo.rw = rw0 = rw = (real*)M1alloc(L0);
	winfo.iw = iw0 = iw = (fint*)(rw + 500);
	winfo.cw = cw0 = cw = (char*)(iw + 500);
	winfo.cw_len = 8;
	memset(rw, 0, L0);
	sninit_(&I0, &I0, cw, &winfo.lencw, iw, &winfo.leniw,
		rw, &winfo.lenrw, winfo.cw_len);
	Plinfy = Infinity;
	lDenJ = 2;
	lvlTim = 1;
	asl->i.congrd_mode = 2;	/* sparse Jacobians */
	av = xargv;
	stub = getstub(&av, &Oinfo);
	if (!stub)
		usage_ASL(&Oinfo, 1);

	nl = jacdim(stub, &M, &N, &NO, &NZ, &MXROW, &MXCOL, (fint)strlen(stub));
	Objno = n_obj > 0;
	meminc = 20*(M + N);
	if (getopts(av,&Oinfo))
		exit(1);
	if (wantfuncs) {
		show_funcs();
		exit(0);
		}

	if (N <= 0) {
		sprintf(Errmsg = Emsgbuf, "%s has no variables\n", filename);
		Errnum = 516;
 bailout:
		solve_result_num = Errnum;
		write_sol(Errmsg, 0, 0, &Oinfo);
		exit(0);
		}

	suf_declare(suftab, sizeof(suftab)/sizeof(SufDecl));
	itnlim = 99999999;
	nnCon = n_conjac[1] = nlc;
	nnJac = nlvc;
	N1 = N + 1;
	m1 = M + 1;
	NB = N + m1;
	nz = NZ + N;
	Fortran = 1;

	LUv = lb = (real *)M1alloc((4*NB + m1 + nz)*sizeof(real)
					+ (N1 + NB + nz)*sizeof(fint) + N);
	LUrhs = lb + N;
	Uvx = ub = LUrhs + m1;
	Urhsx = ub + N;
	X0 = Urhsx + m1;
	A = A_vals = X0 + NB;
	memset((char *)(X0 + N), 0, m1*sizeof(real));
	rc = A + nz;
	pi = pi0 = rc + NB;
	ja = (fint *)(pi + m1);
	ha = ja + N1;
	hs = ha + nz;
	varstat = (int*)M1alloc(NB*sizeof(int));
	vsd = suf_iput("sstatus", ASL_Sufkind_var, varstat);
	csd = suf_iput("sstatus", ASL_Sufkind_con, varstat + N);
	havex0 = (char *)(hs + NB);
	if (sizeof(int) == sizeof(fint))
		A_rownos = (int *)ha;
	else
		A_rownos = (int *)M1alloc(NZ*sizeof(int));
	if (sizeof(int) == sizeof(fint))
		A_colstarts = (int *)ja;
	else
		A_colstarts = (int *)M1alloc(N1*sizeof(int));
	asl->i.nlvog = nlvo;
	Objsen = 2;
	if (meminc <= 0)
		meminc = 20*(M + N);
	nQNmod = 20;
	lvlTim = -999;
	if (setjmp(Jb.jb))
		goto bailout;
	if (envopt()) {
		Errmsg = "Bad options";
		Errnum = 517;
		goto bailout;
		}
	if (Objno > n_obj || Objno < 0) {
		sprintf(Errmsg = Emsgbuf, "Objno = %ld must be in [0, %d]",
			(long)Objno, n_obj);
		Errnum = 518;
		goto bailout;
		}
	--Objno;
	objno = obj_no = Objno;
	flags = objrep_flags[objrep & 3];
	if ((qpsol = qpcheck && nlc == 0) && nlo > 0) {
		qp_read(nl, flags);
		nzh = nqpcheck((int)Objno, &rowq, &colq, &delsq);
		ncolh = 0;	/* not necessary -- it's file static */
		if (nzh > 0) {
			/* undo Fortran numbering and compute ncolh */
			for(i = 1; i <= N; i++)
				if (colq[i-1]-- < colq[i])
					ncolh = i;
			i = --colq[N];
			while(--i >= 0)
				rowq[i]--;
			}
		else if (nzh < 0)
			switch(nzh + 2) {
			 case 0: /* unlikely */
				Errmsg = "Objective involves division by 0";
				Errnum = 519;
				goto bailout;

			 case 1: /* nonlinear */
				qpsol = 0;
				qp_opify();
				break;

			 default:
				sprintf(Errmsg = Emsgbuf,
					"Unexpected return %ld from nqpcheck\n",
					(long)nzh);
				Errnum = 521;
				goto bailout;
			 }
		}
	else
		fg_read(nl,flags);

	if (n_con < M) {
		nnCon = n_conjac[1] = nlc;
		M = n_con;
		N = n_var;
		N1 = N + 1;
		m1 = M + 1;
		NB = N + m1;
		NZ = nzc;
		nz = NZ + N;
		y = lb + N;
		z = LUrhs;
		for(i = 0; i < M; ++i)
			y[i] = z[i];
		LUrhs = y;
		y = ub + N;
		z = Urhsx;
		for(i = 0; i < M; ++i)
			y[i] = z[i];
		Urhsx = y;
		}

	if (sizeof(int) != sizeof(fint)) {
		fint *ja1 = ja;
		fint *ja1e = ja + N1;
		int *ja2 = A_colstarts;
		while(ja1 < ja1e)
			*ja1++ = *ja2++;
		}

	if ((nint = nlogv + niv + nlvbi + nlvci + nlvoi)) {
		printf("ignoring integrality of %ld variables\n", (long)nint);
		need_nl = 0;
		}

	nnObj = objmunge(M,&m1,N,NZ,&nz,hs,lb,ub,A,ha,ja,&objadj,X0);
	helast = 0;
	if (qpsol) {
		if (nnObj < ncolh)
			nnObj = ncolh;
		i = N + m1;
		helast = (fint*)M1zapalloc(i*sizeof(fint));
		while(--i >= N)
			helast[i] = 3;
		}
	start = hs_adjust(asl, N, M, m1, hs, varstat, vsd, csd, A, ha, ja);
	if ((i = nnObj) < nnJac)
		i = nnJac;
	if (lvlScl < 0)
		lvlScl = i > 0 ? 1 : 2;

	if (outlev > 1) {
		if (lvlTim == -999)
			lvlTim = 3;
		if (outlev == 2)
			iSumm = 6;
		else {
			iPrint = 6;
			lprSol = 0;
			if (outlev > 3)
				lprSol = 2;
			}
		}
	else if (lvlTim == -999)
		lvlTim = iPrint == 0 && iSumm == 0 ? 0 : 3;

	if ((iSumm == 6 || iPrint == 6) && need_nl) {
		printf("\n");
		need_nl = 0;
		}
	s1file_(&openF, iw, &winfo.leniw);
	mincw0 = miniw0 = minrw0 = 501;
	i = 0;
	snmemb_(&i, &m1, &N, &nz, &neGcon, &nnCon, &nnJac, &nnObj,
		&mincw0, &miniw0, &minrw0, cw, &winfo.lencw,
		iw, &winfo.leniw, rw, &winfo.lenrw, winfo.cw_len);
	if (i != 104)
		scream("Unexpected return %ld from snMemb.", (long)i);
	stkyOp = 1; /* retain lenalu and minlen */
	vw = 0;
 tryagain:
	maxiw = winfo.leniw = miniw0 + 2*meminc;
	maxrw = winfo.lenrw = minrw0 + meminc;
	maxcw = winfo.lencw = mincw0;
	mincw = miniw = minrw = 501;
	L2 = winfo.leniw*sizeof(fint) + winfo.lenrw*sizeof(real) + winfo.lencw*8;
	if (!vw)
		vw = M1record(Malloc(L2));
	else
		*vw = Realloc(*vw, L2);
	winfo.rw = (real*)*vw;
	winfo.iw = (fint*)(winfo.rw + winfo.lenrw);
	winfo.cw = (char*)(winfo.iw + winfo.leniw);
	memcpy(winfo.rw, rw0, 500*sizeof(real));
	memcpy(winfo.iw, iw0, 500*sizeof(fint));
	memcpy(winfo.cw, cw0, 500*8);
	rw = winfo.rw;
	iw = winfo.iw;
	cw = winfo.cw;
	*stub_end = 0;
	i = strlen(s = basename(stub));
	if (i > 8)
		i = 8;
	memcpy(Prob, s, i);
	while(i < 8)
		Prob[i++] = ' ';
	INFORM = NS = 0;
	Itnzap.itnlimp = &itnlim;
	Itnzap.zapped = 0;
	int_catch(itnlimzap, &Itnzap);
	err_jmp1 = &Jb;
	SI = 0;
	if (setjmp(Jb.jb))
 longjumped:
		SI = &info_surprisr;
	else {
		if (setjmp(fpe_jmpbuf)) {
			report_where(asl);
			printf("\nFloating point error.\n");
			fflush(stdout);
			need_nl = 0;
			goto longjumped;
			}
#ifdef SYMANTEC
		if (_8087)
#endif
		signal(SIGFPE, fpecatch);
		if (Objsen < 0)
			negate(M,pi);
		if (m1 > M)
			pi[M] = 0;

		Times[1] = xectim_();
		if (qpsol)
		  sqopt_(start, qpHx, &m1,
			&N, &nz, &I1, &I0, &ncolh,
			&iObj, &objadj, Prob,
			A, ha, ja, lb, ub, &c0, names,
			helast, hs, X0, pi, rc,
			&INFORM, &mincw, &miniw, &minrw,
			&NS, &nInf, &sInf, &nlobj,
			0, &I0, 0, &I0, 0, &I0,
			cw, &winfo.lencw, iw, &winfo.leniw, rw, &winfo.lenrw,
			(ftnlen)4, (ftnlen)8, (ftnlen)8, (ftnlen)0,
			winfo.cw_len);
		else
		  snopt_(start, &m1, &N, &nz, &I1,
			&nnCon, &nnObj, &nnJac,
			&iObj, &objadj, Prob,
			funcon, funobj,
			A, ha, ja, lb, ub, names,
			hs, X0, pi, rc,
			&INFORM, &mincw, &miniw, &minrw,
			&NS, &nInf, &sInf, &nlobj,
			0, &I0, 0, &I0, 0, &I0,
			cw, &winfo.lencw, iw, &winfo.leniw, rw, &winfo.lenrw,
			(ftnlen)4, (ftnlen)8, (ftnlen)8, (ftnlen)0,
			winfo.cw_len);
		i = INFORM;
		if (i == 84 && !memadj && minlen > lenalu) {
			memadj = meminc += minlen - lenalu + 10;
			goto tryagain;
			}
		}
	Times[2] = xectim_();
	if (functimes && !timing)
		timing = 1;
	if (need_nl && timing & 1) {
		printf("\n");
		need_nl = 0;
		}
	msg1 = msg + Sprintf(msg, "%s: ", Oinfo.bsname);
	if (!SI)
	    switch(i) {
		case 1:		SI = &info_optimal;	break;
		case 2:		SI = &info_feaspt;	break;
		case 3:		SI = &info_inacc;	break;
		case 11:	SI = &info_infeas_lin;	break;
		case 12:	SI = &info_infeas_leq;	break;
		case 13:	SI = &info_infeas_nlm;	break;
		case 14:	SI = &info_infeas_min;	break;
		case 21:	SI = &info_unbdd_obj;	break;
		case 22:	SI = &info_unbdd_con;	break;
		case 31:	SI = &info_limit_iter;	break;
		case 32:	SI = &info_limit_maji;	break;
		case 33:	SI = &info_limit_supb;	break;
		case 41:	SI = &info_cannot_imp;	break;
		case 42:	SI = &info_sing_basis;	break;
		case 43:	SI = &info_infeasgenc;	break;
		case 44:	SI = &info_ill_condnb;	break;
		case 51:	SI = &info_bad_objgrd;	break;
		case 53:	SI = &info_indef_hes;	break;
		case 52:	SI = &info_bad_congrd;	break;
		case 61:	SI = &info_bad_init;
				if (last_state && i_exit == 20)
					SI = &info_bad_obj0;
				break;
		case 62:	SI = &info_bad_feas;	break;
		case 63:	SI = &info_bad_prog;	break;
		case 81:	SI = &info_short_work;	break;
		case 82:	SI = &info_short_char;	break;
		case 83:	SI = &info_short_int;	break;
		case 84:	SI = &info_short_real;	break;
		case 91:	SI = &info_bad_input;	break;
		case 92:	SI = &info_bad_basis;	break;
		case 141:	SI = &info_bad_numbv;	break;
		case 142:	SI = &info_basis_err;	break;
		default:	SI = &info_surprise;
				Sprintf(SI->msg = buf,
					"Surprise INFO = %ld from snOptB", i);
		}
	solve_result_num = SI->code;
	msg1 += Sprintf(msg1, SI->msg, (long)INFORM);
	msg1 += Sprintf(msg1, ".\n%ld iterations", Itn);
	if (SI->wantobj) {
		g_fmtop(buf, Objtru);
		msg1 += Sprintf(msg1, ", objective %s", buf);
		}
	else if (i == 1) {
		g_fmtop(buf, sInf);
		msg1 += Sprintf(msg1, ", infeasibility sum %s", buf);
		}
	if (nfuncon + nfunobj) {
		msg1 += Sprintf(msg1, "\nNonlin evals: ");
		if (nfunobj)
			msg1 += Sprintf(msg1, "obj = %d, grad = %d%s",
				nfunobj, nfungrd, nfuncon ? ", " : "");
		if (nfuncon)
			msg1 += Sprintf(msg1, "constrs = %d, Jac = %d",
				nfuncon, nfunjac);
		msg1 += Sprintf(msg1, ".");
		}
	if (memadj)
		Sprintf(msg1, "\nAssumed meminc = %ld.", (long)memadj);	
	send_status(N, hs, varstat, svtrans, LUv, Uvx);
	send_status(M, hs+N, varstat+N, sctrans, LUrhs, Urhsx);
	write_sol(msg, X0, pi, &Oinfo);
	if (timing) {
		Times[3] = xectim_();
		if (timing & 1)
			time_out(stdout);
		if (timing & 2) {
			fflush(stdout);
			time_out(Stderr);
			}
		}
	ASL_free(&asl);	/* for Purify */
	}

 void
getfnm_(fint *i0, char *fname, fint *last, ftnlen fname_len)
{
	int i = *i0, L;
	char *s;

	if (i >= N_FNAMES  || i < 0) {
		sprintf(Errmsg = Emsgbuf, "File assignment nn=filename with nn = %d not in [0, %d]", i, N_FNAMES-1);
		Errnum = 517;
		longjmp(Jb.jb, 1);
		}
	if ((s = fnames[i]))
		L = Sprintf(fname, "%s", s);
	else
		L = Sprintf(fname, "fort.%d", i);
	*last = L;
	while (L < fname_len)
		fname[L++] = ' ';
	}

 static void
funcon(fint *MODE, fint *M, fint *N, fint *NJAC,
		real *X, real *F, real *G, fint *NSTATE,
		char *CU, fint *LENCU, fint *IU, fint *LENIU,
		real *RU, fint *LENRU, ftnlen len_cu)
{
	real t0;
	int wantder;
	ASL *asl;
	fint nerror = 0;

	if (Itnzap.zapped) {
		*MODE = -2;
		return;
		}
	t0 = 0.;
	if (functimes)
		t0 = xectim_();

	asl = cur_ASL;
	if ((last_state = *NSTATE)) {
		if (*NSTATE != 1)
			return;
		if (*N != c_vars) {
			scream("funcon expected N = %d but got %d\n",
				c_vars, (int)*N);
			/* suppress warning about unused vars: */
			Not_Used(M);
			Not_Used(NJAC);
			Not_Used(CU);
			Not_Used(LENCU);
			Not_Used(IU);
			Not_Used(LENIU);
			Not_Used(RU);
			Not_Used(LENRU);
			Not_Used(len_cu);
			}
		}
	nfuncon++;
	wantder = 0;
	if ((int)*MODE & 2) {
		wantder = 1;
		nfunjac++;
		}
	want_deriv = wantder;
	xknowne(X, &nerror);
	if (nerror <= 0)
		conval(X, F, &nerror);
	if (nerror > 0) {
		i_exit = 21;
		*MODE = -1;
		return;
		}
	if (wantder)
		jacval(X, G, 0);
	xunknown();
	if (functimes)
		Times[4] += xectim_() - t0;
	}

 static void
funobj(fint *MODE, fint *N,
		real *X, real *F, real *G, fint *NSTATE,
		char *CU, fint *LENCU, fint *IU, fint *LENIU,
		real *RU, fint *LENRU, ftnlen len_cu)
{
	int i, wantder;
	real t0;
	fint nerror = 0;
	ASL *asl = cur_ASL;

	if (Itnzap.zapped) {
		*MODE = -2;
		return;
		}
	t0 = 0.;
	if (functimes)
		t0 = xectim_();
	if ((last_state = *NSTATE)) {
		if (*NSTATE != 1)
			return;
		if (*N != o_vars) {
			scream("funobj expected N = %d but got %d\n",
				o_vars, (int)*N);
			/* suppress warning about unused vars: */
			Not_Used(CU);
			Not_Used(LENCU);
			Not_Used(IU);
			Not_Used(LENIU);
			Not_Used(RU);
			Not_Used(LENRU);
			Not_Used(len_cu);
			}
		}
	nfunobj++;
	wantder = 0;
	if ((int)*MODE & 2) {
		nfungrd++;
		wantder = 1;
		}
	want_deriv = wantder;
	i = Objno;
	if (i < 0 || i >= n_obj) {
		*F = 0.;
		if (wantder)
			memset(G, 0, *N*sizeof(real));
		}
	else {
		xknowne(X, &nerror);
		if (nerror >= 0)
			*F = objval(i,X,&nerror);
		if (nerror > 0) {
			i_exit = 20;
			*MODE = -1;
			return;
			}
		if (wantder)
			objgrd(i,X,G,0);
		xunknown();
		}
	if (functimes)
		Times[5] += xectim_() - t0;
	}

#ifdef __cplusplus
}
#endif
