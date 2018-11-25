*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     File  dad.f
*
*     Module to support dense Adifor-style automatic diffentiation of
*     an optimization problem.
*
*     nwddmx   defnfg   inddgn   decl2x   deivtb   spcopy   decsp
*     declsp   readsp   comprC
*
*     Written by Mike Gertz - 18-Jun-99
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine nwddmx
     &   ( ia, lenia, isize, morei )

      implicit
     &     none
      integer
     &     lenia, ia(lenia), isize, morei

*     ==================================================================
*     nwddmx
*
*     07 Mar 2000: First version of nwddmx.
*     18 Jun 2000: Current version.
*     ==================================================================
      integer            kind,     lext
      parameter         (kind = 1, lext = 0)
*     ------------------------------------------------------------------

      call n2admx
     &   ( kind, lext, ia, lenia, isize, morei )

      end ! subroutine nwddmx

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine defnfg
     &   ( nState,
     &     v, g_v, nnvar, nv, isizes,
     &     needF, lenF, F, g_cf,
     &     needG, irowG, neG, jcolG, G,
     &     itcon, mncon, jtvar,
     &     g_funcf, g_bndf,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     nState, needF, nnvar, nv, isizes(*), needG, lenF,
     &     neG, irowG(neG), jcolG(neG), mncon, itcon(lenF),
     &     lencu, leniu, lenru, lencw, leniw, lenrw,
     &     jtvar(nv), iu(leniu), iw(leniw)
      double precision
     &     v(nv), g_v(nnvar,nv), F(lenF), g_cf(nnvar,lenF), G(neG),
     &     ru(lenru), rw(lenrw)
      character*8
     &     cu(lencu), cw(lencw)
      external
     &     g_bndf, g_funcf

*     ==================================================================
*     defnfg calls an Adifor-generated routine to find the function and
*     gradient values.
*
*     This subroutine should be called from within the bndfg routine.
*     Most of the parameters to this routine will have the same meaning
*     as in bndfg. The additional parameters are needed specifically for
*     Adifor-generated code. The new parameters are:
*
*     itcon(lenF) - constraint k is linear if and only if itcon(k) is
*                   zero.
*
*     mncon       - The number of nonlinear constraints.
*
*     jtvar       - The variable v(k) occurs nonlinearly if and only if
*                   jtvar(k) .ne. 0. If the variable v(k) occurs
*                   nonlinearly, then g_v(jtvar)k),:) is the
*                   differential corresponding to v(k).
*                   (See below for a definition of g_v).
*
*     nnvar       - The number of variables that occur nonlinearly.
*
*     g_v(nnvar,nv)
*                 - holds the differentials of the nonlinear independent
*                   variables. Because this is used by Adifor, it is
*                   transposed.  Note that the shape of g_v is
*                   g_v(nnvar,nv), where nv = 2*(nx + 1) + np.
*                   Because nnvar may be zero, g_v may be empty!
*
*     g_cf(nnvar,lenF)
*                 - holds the differentials of the dependent variables.
*                   The matrix g_cf will be filled by a call to Adifor-
*                   generated code.
*
*     g_bndf      - The Adifor-generated routine to be called.
*
*     07 Mar 2000: First version of defnfg.
*     18 Jun 2000: Current version.
*     ==================================================================
*     By the time this function is called, the structure of the
*     derivatives have been initialized, so the user's
*     function has already been called. SNOPT may not know this, and
*     may set nState .eq. 1. We trap this, and pass nState .eq. 0 down
*     to the user

      integer
     &     nnStat
*     ------------------------------------------------------------------
      if ( nState .eq. 1 ) then
         nnStat = 0
      else
         nnStat = nState
      end if

      call zrlin
     &   ( jtvar, nv, v )

      call g_funcf
     &   ( nnStat,
     &     v, g_v, nnvar, nv, isizes,
     &     needF, lenF, F, g_cf,
     &     g_bndf,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      call spcopy
     &   ( g_cf, nnvar, lenF, nnvar,
     &     irowG, neG, jcolG, G,
     &     jtvar, nv )

      call zrlin
     &   ( itcon, lenF, F )

      end ! subroutine defnfg

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine inddgn
     &   ( v0, nv, isizes, blv, buv,
     &     bl, lenF, bu,
     &     funfg, g_funf, adtlc, ladtlc,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     ia, lenia, ra, lenra, lusdi, lusdr,
     &     morei, morer )

      implicit
     &     none
      integer
     &     nv, lenF, ladtlc, lusdi, lusdr, morei, morer, lencu, lencw,
     &     leniu, leniw, lenia, lenru, lenrw, lenra, isizes(*),
     &     iu(leniu), iw(leniw), ia(lenia)
      double precision
     &     v0(nv), blv(nv), buv(nv), bl(lenF), bu(lenF),
     &     ru(lenru), rw(lenrw), ra(lenra)
      character*8
     &     cu(lencu), cw(lencw)
      external
     &     funfg, g_funf, adtlc, ladtlc

*     ==================================================================
*     inddgn
*
*     07 Mar 2000: First version of inddgn.
*     18 Jun 2000: Current version.
*     ==================================================================
      integer
     &     iirowC, neC, ijcolC, ic, iirowG, neG, ijcolG, iitcon, mncon,
     &     ijtvar, nnvar, ig_cf, ig_cf2, ig_v, ifc, iv
*     ------------------------------------------------------------------
      double precision   zero,          one
      parameter         (zero = 0.0d+0, one = 1.0d+0)
*     ------------------------------------------------------------------
      morei = 0
      morer = 0

*     Now, lusdi is correct (stadmx knows the size of the integer
*     scalars, this routine does not.) Allocate the type arrays.

      iitcon = lusdi + 1
      ijtvar = iitcon + lenF
      lusdi  = ijtvar + nv - 1
      if ( lusdi .gt. lenia ) goto 600

      ig_cf  = lusdr  + 1
      ig_cf2 = ig_cf  + nv*lenF
      ig_v   = ig_cf2 + nv*lenF
      ifc    = ig_v   + nv*nv
      iv     = ifc    + 2 *lenF
      lusdr  = iv     + nv - 1


      if ( lusdr .gt. lenra ) goto 600

*     Allocation was successful, pass the new arrays to an auxiliary
*     routine that calls the user's function twice.
*     First, zero out the differentials.

      call dload
     &   ( (ifc-ig_cf), zero, ra(ig_cf), 1 )

      call decl2x
     &   ( v0, nv, blv, buv, isizes, ra(iv), ra(ig_v),
     &     lenF, ra(ifc), ra(ig_cf),
     &     funfg, g_funf,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      call decsp
     &   ( nv, lenF, nv,
     &     iirowC, ijcolC, neC, ra(ig_cf),
     &     iirowG, ijcolG, neG,
     &     ra(ifc), ra(iv),
     &     ia(iitcon), ia(ijtvar), mncon, nnvar,
     &     adtlc, ladtlc,
     &     ia, lenia, lusdi, morei )

      if ( morei .ne. 0 ) goto 600

*     Adjust the bounds to reflect the rhs of the affine constraints.

      call daxpy
     &   ( lenF, (-one), ra(ifc), 1, bl, 1 )
      call daxpy
     &   ( lenF, (-one), ra(ifc), 1, bu, 1 )

*     Put C in the same location as g_cf, which is more than big enough.

      ic    = ig_cf

*     Deallocate the real scratch arrays, and make C permanent.

      lusdr = iC + neC - 1

*     Now allocate the arrays of differentals (reallocating g_cf).

      ig_cf  = lusdr + 1
      ig_v   = ig_cf + nnvar*lenF
      lusdr  = ig_v  + nnvar*nv   - 1
      if ( lusdr .gt. lenra ) goto 600

*     Initialize the inverse table and differentials.

      if ( nnvar. ne. 0 ) then
         call deivtb
     &      ( ia(ijtvar), nv, ra(ig_v), nnvar, nnvar )
      end if

*     Now set all the globals.

      call setmtx
     &   ( iirowC, neC, ijcolC, ic,
     &     iirowG, neG, ijcolG, lenF,
     &     ia, lenia )

      call stadmx
     &   ( iitcon, mncon, ijtvar, nnvar,
     &     ig_v, ig_cf,
     &     ia, lenia )

      return

*     Handle memory errors.

  600 morei = max( morei, lusdi - lenia )
      morer = max( morer, lusdr - lenra )

      lusdi = 0
      lusdr = 0

      end ! subroutine inddgn

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine decl2x
     &   ( v0, nv, blv, buv, isizes, v, g_v,
     &     lenF, F, g_f, funfg, g_funf,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     nv, lenF, lencu, lencw, leniu, leniw, lenru, lenrw,
     &     isizes(*), iu(leniu), iw(leniw)
      double precision
     &     v0(nv), blv(nv), buv(nv), v(nv), g_v(nv,nv),
     &     F(2*lenF), g_f(nv, 2*lenF), ru(lenru), rw(lenrw)
      external
     &     funfg, g_funf
      character*8
     &     cu(lencu), cw(lencw)

*     ==================================================================
*     decl2x  calls an Adifor-generated "dense" subroutine twice at
*     randomly generated points. decl2x is designed to be an auxiliary
*     routine for inddgn and should not be called by any other
*     routine.
*
*     07 Mar 2000: First version of decl2x.
*     18 Jun 2000: Current version.
*     ==================================================================
      integer
     &     k, iseed(3)
*     ------------------------------------------------------------------
      double precision   zero,          one
      parameter         (zero = 0.0d+0, one = 1.0d+0)
      integer            nState,     needF
      parameter         (needF = 1)
*     ------------------------------------------------------------------
*     Get a new random seed.

      call gtrdsd
     &   ( iseed(1), iseed(2), iseed(3) )

*     Form the identity in the differentials of the dependent variables.
*     They were already set to zero by the calling routine.
*     "zero out" F, so that the returned value is well-defined even if
*     the user doesn't set it.

      call dload ( nv, one, g_v(1,1), nv+1 )
      call dload ( 2*lenF, zero, F, 1 )

      do k = 1, 0, -1

*        Create a random perturbation of the variables.

         call rndprb( nv, v0, blv, buv, one, v, iseed )

*        nState should be set to one for the first call and
*        zero for subsequent calls.
         nState = k
         call funfg
     &      ( nState,
     &        v, g_v, nv, nv, isizes,
     &        needF, lenF, F(1+k*lenF), g_f(1,1+k*lenF),
     &        g_funf,
     &        cu, lencu, iu, leniu, ru, lenru,
     &        cw, lencw, iw, leniw, rw, lenrw )

*     end do k = 0,1
      end do

      end ! subroutine decl2x

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine deivtb
     &   ( jtvar, nv, g_v, ldg_v, nnvar )

      implicit
     &     none
      integer
     &     nv, ldg_v, jtvar(nv), nnvar
      double precision
     &     g_v(ldg_v, nv)

*     ==================================================================
*     deivtb initialize the "inverse table" and the differentials of the
*     independent variables. This routine is an auxiliarly routine for
*     inddgn, and should not be called by any other subroutine.
*
*     Because some of the independent variables occur only in linear
*     terms, it is not necessary to compute these (necessarily constant)
*     terms at each iteration.  We omit this computation by not passing
*     the the differential corresponding to those terms to the
*     Adifor generated code. Thus, g_v is not the identity matrix, but
*     rather a matrix which only includes the differentials of those
*     variables with occur nonlinearly. (I would say it is the identity
*     matrix, with the columns corresponding to the linear variables
*     omitted.  Recall, however, that g_v is transposed, so the rows,
*     are omitted not the columms.)
*
*     jtvar(nv) - (input/output) On entry, v(k) occurs nonlinearly
*                 if and only if jtvar(k) .ne. 0.
*                 On exit, if v(k) occurs nonlinearly, then
*                 g_v(jtvar(k), :) is the differential corresponding
*                 to v(k).  If v(k) occurs only linearly, then there is
*                 no corresponding differential and we set jtvar(k)=0.
*
*     Notice that on entry and exit, v(k) occurs nonlinearly if and
*     only if jtvar(k) .ne. 0.  The value of jtvar(k), however, may
*     change.
*
*     07 Mar 2000: First version of deivtb.
*     18 Jun 2000: Current version.
*     ==================================================================
      integer
     &     i, j
*     ------------------------------------------------------------------
      integer            linear
      parameter         (linear = 0)
      double precision   zero,         one
      parameter         (zero = 0.0d0, one = 1.0d0)
*     ------------------------------------------------------------------

      call dload
     &   ( nnvar*nv, zero, g_v, 1 )

      j = 1
      do i = 1, nv
         if ( jtvar(i) .ne. linear ) then
            g_v(j,i) = one
            jtvar(i) = j
            j = j + 1
         end if
      end do

      end ! subroutine deivtb

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine spcopy
     &   ( gf, ldgf, m, n, irowG, neG, jcolG, G, jtvar, nv )

      implicit
     &     none
      integer
     &     n, m, ldgf, neG, nv, jtvar(nv), irowG(neG), jcolG(neG)
      double precision
     &     gf(ldgf, m), G(neG)

*     ==================================================================
*     spcopy copies the contents of the Adifor-style dense matrix  gf
*     into the sparse structure (irowG, jcolG, G).  As with all
*     Adifor-style matrices,  gf  is transposed.
*
*     gf(ldfg, m)    - is a dense matrix
*     irowG(lenG), etc.
*                    - is a sparse matrix. The kth element of G is given
*                      by the triple (irowG(k), jcolG(k), G(k)).
*
*     jtvar(nv)      - is an array that maps the columns numbers of G to
*                      the column numbers of gf.  In other words:
*                        G(k) = gf(jtvar(jcolG(k)), irowG(k))
*
*     07 Mar 2000: First version of spcopy.
*     18 Jun 2000: Current version.
*     ==================================================================
      integer
     &     i, j, k
*     ------------------------------------------------------------------

      do k = 1, neG
         i = irowG(k)
         j = jcolG(k)

         G(k) = gf(jtvar(j), i)
      end do

      end ! subroutine spcopy

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine  decsp
     &   ( ldg, m, n,
     &     iirowC, ijcolC, neC, C,
     &     iirowG, ijcolG, neG,
     &     F, v,
     &     itcon, jtvar, mncon, nnvar,
     &     adtlc, ladtl, ia, lenia, lusdia, mormem )

      implicit
     &     none
      integer
     &     ldg, m, n, mormem, ladtl, nnvar, mncon, jtvar(n), itcon(m),
     &     lusdia, lenia,  ia(lenia), iirowC, ijcolC, neC,
     &     iirowG, ijcolG, neG
      double precision
     &     C(2*ldg*m), F(m), v(n)
      external
     &     adtlc, ladtl

*     ==================================================================
*     decsp  finds the sparsity pattern and constant elements of the
*     Jacobian.
*
*     Matrices G1 and G2 represent the (transpose of the) m by n
*     Jacobian matrix evaluated at two distinct points. It is expected
*     that those elements that remain constant in both G1 and G2 are
*     truly constant, and that those elements that are zero in both G1
*     and G2 are identically zero for all x. It is recommended that G1
*     and G2 be evaluated at a pseudo-random perturbation from an
*     initial guess.
*
*     G1 and G2 are the transpose of the desired Jacobian because the
*     transpose is computed by Adifor-generated code.
*     However, we still think of G1(j,i) as the Jacobian element at row
*     i, column j. To avoid confusion, we will refer to G1(j,i) as
*     the Jacobian element for constraint i, variable j.
*
*     ia(lenia)    - an array whose elements describe the sparsity
*                    pattern and pattern of constant elements in the
*                    Jacobian
*
*     (ia(iirowC+k-1), ia(ijcolC+k-1))
*                  - is the (constraint, variable) index of the kth
*                    nonzero constant element of the Jacobian.
*                    There are neC such constant elements, stored in
*                    row-major order.
*
*     (ia(iirowG+k-1), ia(ijcolG+k-1))
*                  - is the (constraint, variable) index of the kth
*                    nonzero non-constant element of the Jacobian.
*                    There are neG such elements stored in row-major
*                    order.
*
*     itcon(m)     - Each element, itcon(k), tells whether the given
*                    constraint is linear (0) or nonlinear (non-zero).
*                    A constraint is linear if every Jacobian element
*                    in that constraint is constant.
*
*     itvar(n)     - Each element, itvar(k), tells whether the given
*                    variable is linear (0) or nonlinear(non-zero).
*                    A variable is linear if every Jacobian element in
*                    the variable is constant.
*
*     lusdia       - On entry, lusdia + 1 is the first free element
*                    of ia.
*                    On exit, lusdia is the total number of elements
*                    of ia that are allocated (i.e. lusdia+1 is the new
*                    first free element of ia, if lusdia+1 .le. lenia.)
*
*     An element of the Jacobian is considered to be constant if it is
*     the same in G1 and G2 AND either it occurs in a linear constraint
*     or it occurs in a linear variable. Formally, the rule for
*     accepting an element as constant is as follows. Scan to find which
*     elements remain the same in G1 and G2 and determine which
*     variables are linear/nonlinear using this information. Then choose
*     those elements which remain the same in G1 and G2 AND either occur
*     in a linear constraint or a linear variable to be constant.
*
*     mormem       - if mormem > 0, then the subroutine failed because
*                    ia was to small. The length of ia must be increase
*                    by at least mormem for the subroutine to find the
*                    sparsity structure.
*
*     Local variables:
*
*     lenC         - is the actual allocated lengths of the arrays at
*                    ia(iirowC), iaia(ijcolC).  On completion of decsp,
*                    lenC = neC.
*     lenG         - is similar to lenC.
*     nz           - is the number of nonzero elements.

*     07 Mar 2000: First version of decsp.
*     18 Jun 2000: Current version.
*     ==================================================================
      integer
     &     lenC, lenG, nz, i, j, k, iold, iG2, necadl
*     ------------------------------------------------------------------
      double precision   zero
      parameter         (zero = 0.0d0)
*     ------------------------------------------------------------------
      mormem = 0
      iG2    = 1 + ldg * m

*     find out which variable
      call declsp
     &   ( C(1), ldg, m, n, C(iG2),
     &     itcon, jtvar, mncon, nnvar, nz )

*     Set lenG and lenC to their maximium possible sizes

      lenG = min( nnvar * mncon, nz )
      lenC = nz - max(nnvar, mncon)

      iirowG = lusdia + 1
      ijcolG = iirowG + lenG

      iirowC = ijcolG + lenG
      ijcolC = iirowC + lenC

*     Make sure we have the memory

      lusdia = ijcolC + lenC - 1
      if ( ijcolC + lenC - 1 .gt. lenia ) goto 900

*     Fill in the sparsity pattern

      call readsp
     &   ( C(1), ldg, m, n, C(iG2), ldg,
     &     jtvar, itcon,
     &     ia(iirowG), lenG, ia(ijcolG), neG,
     &     ia(iirowC), lenC, ia(ijcolC), neC )

*     Compress the results

      lenG   = neG
      iold   = ijcolG
      ijcolG = iirowG + lenG

      call icopy
     &   ( lenG, ia(iold), 1, ia(ijcolG), 1 )

      lenC   = neC
      iold   = iirowC
      iirowC = ijcolG + lenG
      call icopy
     &   ( lenC, ia(iold), 1, ia(iiRowC), 1 )

      iold   = ijcolC
      ijcolC = iirowC + lenC
      call icopy
     &   ( lenC, ia(iold), 1, ia(ijcolC), 1 )

      lusdia = ijcolC + lenC - 1
      call cmprC
     &   ( C(iG2), ldg, m, n, ia(iirowC), neC, ia(ijcolC), C )

*     Store the constant of affine constraints in F.

      do k = 0, neC-1
         i = ia( iirowC + k )
         if ( itcon(i) .eq. 0 ) then

*           Affine constraint.

            j    = ia( ijcolC + k )
            F(i) = F(i)  - C(k+1) * v(j)
         end if
      end do

      do i = 1, m
         if ( itcon(i) .ne. 0 ) then

*           Nonlinear constraint.  No meaningful constant.

            F(i) = zero
         end if
      end do

*     Get any additional constants that the user may have defined.  The
*     variable neC is at most n*m, where length(C) is 2*m*n, so
*     C(neC + 1) exists.

      necadl = ladtl()
      if ( necadl .gt. 0 ) then
         lusdia = iirowC + 2 * (neC + necadl) - 1
         if ( lusdia .gt. lenia ) goto 900
*        Move ijcolC
         do i = neC-1, 0, -1
            ia(iirowC+neC+necadl+i) = ia(ijcolC+i)
         end do
         ijcolC = iirowC + neC + necadl
         call adtlc
     &      ( ia(iirowC+neC), necadl, ia(ijcolC+neC), C(neC+1))
         neC = neC + necadl
         call lexsrt
     &      ( ia(iirowC), neC, ia(ijcolC), C )
      end if

      return

 900  continue

*     Handle memory errors. Unallocate all arrays and set mormem

      mormem = lusdia - lenia
      iirowC = 0
      ijcolC = 0
      iirowG = 0
      ijcolG = 0
      lusdia = 0

      end ! subroutine decsp

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine declsp
     &   ( G1, ldg, m, n, G2, itcon, jtvar, mncon, nnvar, nz)

      implicit
     &     none
      integer
     &     ldg, m, n, i, j, nz, itcon(m), jtvar(n), mncon, nnvar
      double precision
     &     G1(ldg,m), G2(ldg,m)

*     ==================================================================
*     declsp
*
*     07 Mar 2000: First version of declsp.
*     18 Jun 2000: Current version.
*     ==================================================================
      integer            linear,     nonlin
      parameter         (linear = 0, nonlin = 1)
      double precision   zero
      parameter         (zero   = 0.0d0)
*     ------------------------------------------------------------------
      do i = 1, m
         itcon(i) = linear
      end do

*     Find which variables/constraints contain nonlinear elements.

      nz = 0
      do j = 1, n
         jtvar(j) = linear
         do i = 1, m
            if ( G1(j, i) .ne. zero  .or.  G2(j, i) .ne. zero ) then
               nz = nz + 1
            end if

            if ( G1(j, i) .ne. G2(j,i) ) then
*              nonlinear element
               jtvar(j) = nonlin
               itcon(i) = nonlin
            end if
         end do
      end do

*     Count the number of variables that occur as nonlinear elts.
      nnvar = 0
      do j = 1, n
         if ( jtvar(j) .eq. nonlin ) then
            nnvar = nnvar + 1
         end if
      end do

*     Count the number of nonlinear constraints.

      mncon = 0
      do i = 1, m
         if ( itcon(i) .eq. nonlin ) then
            mncon = mncon + 1
         end if
      end do

      end ! subroutine declsp

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine readsp
     &   ( G1, ldG1, m, n, G2, ldG2, jtvar, itcon,
     &     irowG, lenG, jcolG, neG, irowC, lenC, jcolC, neC )

      implicit
     &     none
      integer
     &     ldG1, m, n, ldG2, jtvar(n), itcon(m), lenG, neG, lenC, neC,
     &     irowG(lenG), jcolG(lenG), irowC(lenC), jcolC(lenC)
      double precision
     &     G1(ldG1, m), G2(ldG2, m)

*     ==================================================================
*     readsp finds the sparsity pattern and constant elements of a Jacobian
*     matrix. This is an auxilary routine to decsp, so most of the
*     parameters have the same meaning.
*
*     The parameters irowG, jcolG, irowC and jcolC are modified by this
*     routine to contain the actual sparsity patterns. The arrays must
*     be long enough to hold these patterns.
*
*     On entry, itcon and jtvar must correctly identify rows/cols as
*     linear/nonlinear. This routine adds elements to irowG/jcolG if
*     both the row and column of the element is nonlinear. Otherwise,
*     a non-zero element is added to irowC/jcolC.  irowC/jcolC as
*     appropriate.  If itcon(i) .eq. 0, then row i is linear.
*     Otherwise, it is nonlinear. If jtvar(j) .eq. 0, then col j is
*     linear. Otherwise, it is nonlinear.
*
*     An element is considered to be non-zero if it is non-zero in
*     EITHER G1 or G2.
*
*     Remember that ADIFOR matrices are tranposed!
*
*     07 Mar 2000: First version of readsp.
*     18 Jun 2000: Current version.
*     ==================================================================
      integer
     &     i, j
*     ------------------------------------------------------------------
      double precision   zero
      parameter         (zero = 0.0d0)
      integer            linear
      parameter         (linear = 0)
*     ------------------------------------------------------------------
      neG = 0
      neC = 0
      do i = 1, m
         do j = 1, n
            if ( (G1(j,i) .ne. zero) .or. (G2(j,i) .ne. zero) ) then
               if((jtvar(j) .ne. linear) .and. (itcon(i) .ne. linear))
     &         then
*                 Nonlinear element.
                  neG        = neG + 1
                  irowG(neG) = i
                  jcolG(neG) = j
               else
*                 Linear element.
                  neC        = neC + 1
                  irowC(neC) = i
                  jcolC(neC) = j
               end if
            end if
         end do
      end do

      end ! subroutine readsp

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine cmprC
     &   ( G, ldG, lenF, n, irowC, lenC, jcolC, C )

      implicit
     &     none
      integer
     &     n, lenF, ldG, lenC, irowC(lenC), jcolC(lenC)
      double precision
     &     G( ldG, lenF ), C(lenC)

*     ==================================================================
*     cmprC copies the data from the dense matrix G into the sparse
*     matrix C.
*
*     On entry, irowC and jcolC contain the (row, col) locations
*     of the non-zero elements of the transpose of  G. On exit, the
*     triple (irowC, jcolC, C) will be a sparse representation of
*     the transpose of G.
*
*     ADIFOR, unfortunately, forces us to work with G transpose, rather
*     than G.
*
*     If the elements of C are given in row-major order, then it is safe
*     for C and G to map to the same location in memory. Thus, this
*     routine "compresses" G into C.
*
*     07 Mar 2000: First version of cmprC.
*     18 Jun 2000: Current version.
*     ==================================================================
      integer
     &     k
*     ------------------------------------------------------------------
      do k = 1, lenC
         C(k) = G( jcolC(k), irowC(k) ) !  G is transposed
      end do

      end ! subroutine cmprC
