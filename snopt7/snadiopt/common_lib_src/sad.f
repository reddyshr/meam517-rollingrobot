*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     File  sad.f
*
*     Module to suport sparse Adifor-style automatic diffentiation of an
*     optimization problem.
*
*     inadfr   nwsdmx   spfnfg   spunwrp   speye   insdgn   spcl2x
*     spivtb   spovrl   spcsp    spreed    spccnt  spcmrg
*
*     Written by Mike Gertz - 03-Mar-00
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine inadfr
      integer       iadifr
      save          iadifr
      data          iadifr / 0 /

      if (iadifr .eq. 0) then
         iadifr = 1
         call xspini
      end if

      end ! subroutine inadfr

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine nwsdmx( ia, lenia, isize, morei )
      implicit      none

      integer       lenia, ia(lenia), isize, morei
      integer       kind, lext
      parameter     ( kind = 2, lext = 0 )

      call inadfr
      call n2admx( kind, lext, ia, lenia, isize, morei )

      end ! subroutine nwsdmx

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine spfnfg( nState,
     &     v, j_v, nv, isizes,
     &     needF, mCon, f, j_cf,
     &     needG, irowg, neg, jcolg, g,
     &     itcon, mncon, jtvar, nnvar,
     &     iwork, rwork,
     &     funcf, g_funcf,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit            none

      integer             nState, nv, j_v(nv), isizes(*)
      double precision    v(nv)
      integer             needF, mCon, j_cf(mCon)
      double precision    f(mCon)
      integer             needG, neg, irowg(neg), jcolg(neg)
      double precision    g(neg)
      integer             itcon(mCon), mncon, jtvar(nv), nnvar
      integer             iwork(nv)
      double precision    rwork(nv)
      external            funcf, g_funcf
      integer             lencu, lencw
      integer             leniu, leniw
      integer             lenru, lenrw
      character*8         cu(lencu), cw(lencw)
      integer             iu(leniu), iw(leniw)
      double precision    ru(lenru), rw(lenrw)

      integer             info

      integer             nnStat

*     By the time this function is called, the structure of the
*     derivatives have definitely been initialized, so the user's
*     function has already been called. SNOPT may not know this, and
*     may set nState .eq. 1. We trap this, and pass nState .eq. 0 down
*     to the user

      if ( nState .eq. 1 ) then
         nnStat = 0
      else
         nnStat = nState
      end if

      call zrlin( jtvar, nv, v )

      call g_funcf( nnStat,
     &     v, j_v, nv, isizes,
     &     needF, mCon, f, j_cf,
     &     funcf,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      call zrlin( itcon, mCon, f )

      call spovrl( jtvar, nv, nnvar, j_cf, mCon,
     &     irowg, neg, jcolg, g,
     &     iwork, rwork, info )

      end ! subroutine spfnfg

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine spunwrp( ig, n, irowc, lenc, nec, ishift, jcolc, c )
*      implicit            none
      integer              n, ig(n), lenc, irowc(lenc), jcolc(lenc)
      double precision     c(lenc)

      integer              ic

      ic = 1
      do i = 1, n
         call dspxsq( jcolc(ic), c(ic), lenc - ic + 1,
     &        ig(i), isize, info )

         do k = ic, ic + isize - 1
            irowc(k) = i + ishift
         end do
         ic = ic + isize
      end do

      nec = ic - 1

      end ! subroutine spunwrp

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine speye ( ig, n, ishift )
      implicit none

      integer             n, ig(n), ishift

      double precision    one
      parameter           ( one = 1.0d0 )
      integer             k

      integer             iwork(1)
      double precision    dwork(1)

      dwork(1) = one

      do k = 1, n
         iwork(1) = k + ishift
         call dspsd( ig(k), iwork, dwork, 1 )
      end do

      end ! subroutine speye

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine insdgn( v, nv, isizes, blv, buv,
     &     mCon, bl, bu,
     &     g_funfg, g_funf, adtlc, ladtl,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     ia, lenia, ra, lenra, lusdi, lusdr,
     &     morei, morer )

      implicit             none

      integer              nv, isizes(*)
      double precision     v(nv), blv(nv), buv(nv)
      integer              mCon
      double precision     bl(mCon), bu(mCon)
      integer              ladtl
      external             g_funfg, g_funf, adtlc, ladtl
      integer              lencu, lencw
      integer              leniu, leniw, lenia
      integer              lenru, lenrw, lenra
      character*8          cu(lencu), cw(lencw)
      integer              iu(leniu), iw(leniw), ia(lenia)
      double precision     ru(lenru), rw(lenrw), ra(lenra)
      integer              lusdi, lusdr, morei, morer

      integer              iirowc, nec, ijcolc, ic
      integer              iirowg, neg, ijcolg
      integer              iitcon, mncon, ijtvar, nnvar

      integer              ig_v, ig_cf, ig_cf2, ifc, ifc2, iv

      integer              indvec(1)
      double precision     valvec(1)

      integer              i, j

*     Now, lusdi is correct (stadmx knows the size of the integer
*     scalars, this routine does not.) Allocate the type arrays.
      iitcon = lusdi + 1
      ijtvar = iitcon + mCon
      lusdi = ijtvar + nv - 1
      if ( lusdi .gt. lenia ) goto 600

      ig_v      = lusdi  + 1
      ig_cf     = ig_v   + nv
      ig_cf2    = ig_cf  + mCon
      lusdi     = ig_cf2 + mCon - 1

      ifc       = lusdr  + 1
      ifc2      = ifc    + mCon
      iv        = ifc2   + mCon
      lusdr     = iv     + nv    - 1
      if ( lusdr .gt. lenra .or. lusdi .gt. lenia) goto 600
*     Allocation was successful, pass the new arrays to an auxiliary
*     subroutine that calls the user's function twice.
*     But first, zero out all the differentials
      do j = ig_v, lusdi
         ia(j) = 0
      end do

      call spcl2x(  v, ia(ig_v), nv, isizes, blv, buv,
     &     mCon, ra(ifc), ia(ig_cf),
     &     g_funfg, g_funf,
     &     ra(iv),
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

         call spcsp( ia(ig_cf), mCon, nv,
     &     ia(ig_cf2),
     &     iirowc, ijcolc, nec, ic,
     &     iirowg, ijcolg, neg,
     &     ra(ifc), ra(iv),
     &     ia(iitcon), ia(ijtvar), mncon, nnvar,
     &     adtlc, ladtl,
     &     ia, lenia, lusdi,
     &     ra, lenra, lusdr, morei, morer )


      if ( morei .ne. 0 .or. morer .ne. 0 ) goto 600

      do i = 1, mCon
         bl(i) = bl(i) - ra(ifc + i - 1)
         bu(i) = bu(i) - ra(ifc + i - 1)
      end do
*     deallocate ifc and iv: move ic
      do i = 0, nec - 1
         ra(ifc + i) = ra(ic + i)
      end do
      ic = ifc
      lusdr = ic + nec - 1

*     Deallocate ig_cf2
      indvec(1) = 1
      valvec(1) = 1d0
      do i = 0, mCon - 1
         call dspsd(ia(ig_cf2 + i), indvec(1), valvec(1), 0)
      end do

*     Deallocate ig_cf2 in ia
      call mvi( ia, ig_cf2 + mCon,
     &     ig_cf2, lusdi - (ig_cf2 + mCon) + 1)
      ig_cf2 = 0
      lusdi  = lusdi - mCon

*     Correct the variables effected by the shift
      iirowc = iirowc - mCon
      ijcolc = ijcolc - mCon
      iirowg = iirowg - mCon
      ijcolg = ijcolg - mCon

*     Initailize the inverse table and differentials in an auxiliary routine
      call spivtb( ia(ijtvar), nv, ia(ig_v) )

*     Now set all the globals
      call setmtx( iirowc, nec, ijcolc, ic,
     &     iirowg, neg, ijcolg, mCon,
     &     ia, lenia )

      call stadmx( iitcon, mncon, ijtvar, nnvar,
     &     ig_v, ig_cf,
     &     ia, lenia )

      return
 600  continue
*     Handle memory errors.
      morei = max( morei, lusdi - lenia )
      morer = max( morer, lusdr - lenra )

      lusdi = 0
      lusdr = 0

      end ! subroutine insdgn

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine spcl2x(  v0, j_v, nv, isizes, blv, buv,
     &     mCon, f, j_f,
     &     g_funfg, g_funf,
     &     v,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit            none

      integer             nv, j_v(nv), isizes(*)
      double precision    v0(nv), blv(nv), buv(nv)
      integer             mCon, j_f(2*mCon)
      double precision    f(2*mCon)
      external            g_funfg, g_funf
      double precision    v(nv)
      integer             lencu, lencw
      integer             leniu, leniw
      integer             lenru, lenrw
      character*8         cu(lencu), cw(lencw)
      integer             iu(leniu), iw(leniw)
      double precision    ru(lenru), rw(lenrw)

      integer             iseed(3)

      integer             nState,     needF
      parameter         ( needF = 1 )
      double precision    one
      parameter         ( one = 1.0d0 )

      integer             k

*     Get a new random seed.
      call gtrdsd( iseed(1), iseed(2), iseed(3) )

*     Form the identity in the differentials of the dependent variables.
*     They were already set to zero by the calling routine.
      call speye( j_v, nv, 0 )

*     Set f to zero before sending to the user routines.
      do k = 1, 2 * mCon
         f(k) = 0.0d0
      end do

      do k = 1, 0, -1
*        Create a random perturbation of the input variables.
         call rndprb( nv, v0, blv, buv, one, v, iseed )

*        nState should be set to one for the first call, and zero
*        for subsequent calls.
         nState = k
         call g_funfg( nState,
     &        v, j_v, nv, isizes,
     &        needF, mCon, f(1 + k*mCon), j_f(1 + k*mCon),
     &        g_funf,
     &        cu, lencu, iu, leniu, ru, lenru,
     &        cw, lencw, iw, leniw, rw, lenrw )

*     end do k = 0,1
      end do

      end  ! subroutine spcl2x

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine spivtb( jtvar, nv, g_v )

*     Initialize the "inverse table" and the differentials of the
*     independent variables. This routine is meant to be an auxiliarly
*     routine for iniadgen, and should not be called by any other
*     subroutine.
*
*     Because some of the independent varibles occur only in linear
*     terms, it is not necessary to compute these (necessarily constant)
*     terms at each iteration. We omit this computation by not passing
*     the the differential corresponding to those terms to the
*     Adifor generated code. Thus, g_v is not the identity matrix, but
*     rather a matrix which only includes the differentials of those
*     variables with occur non-linearly. (I would say it is the identity
*     matrix, with the columns corresponding to the linear variables
*     ommited. Recall, however, that g_v is transposed, so we omit rows,
*     not columms.)
*
*     jtvar(nv) - (input/output) On entry, v(k) occurs non-linearly
*          if and only if jtvar(k) .ne. 0. On exit, if v(k) occurs
*          non-linearly, then g_v(jtvar(k), :) is the differential
*          corresponding to v(k).  If v(k) occurs only linearly, then
*          there is no corresponding differential and we let
*          jtvar(k) .eq. 0.
*
*     Notice that on entry and on exit, v(k) occurs non-linearly if and
*     only if jtvar(k) .ne. 0. The value of jtvar(k), however, may
*     change.


      implicit            none

      integer             nv, jtvar(nv)
      integer             g_v(nv), indvec(1)

      double precision    one(1)

      integer             linear
      parameter           ( linear = 0 )

      integer             i, j

      one(1) = 1.0d0

      do i = 1, nv
         g_v(i) = 0
      end do

      j = 1
      do i = 1, nv
         if( jtvar(i) .ne. linear ) then
*           g_v(j,i) = one
            indvec(1) = j
            call dspsd( g_v(i), indvec, one, 1 )
            jtvar(i) = j
            j = j + 1
         else
            call dspsd( g_v(i), indvec, one, 0 )
         end if
      end do

      end ! subroutine spivtb

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine spovrl( jtvar, nv, nnvar, idg, m,
     &     irowg, neg, jcolg, g,
     &     iwork, dwork, info)

      implicit            none
*      implicit double precision (a-h, o-z)

      integer             nv, jtvar(nv), nnvar, m, idg(m)
      integer             neg, irowg(neg), jcolg(neg)
      double precision    g(neg)

      integer             iwork(nnvar)
      double precision    dwork(nnvar)

      double precision    zero
      parameter           ( zero = 0.0d0 )

      integer             info
      integer             i, jg, jtg, jtidg, kg, kidg
      integer             ksize, icurrent, info2

      info  = 1
      info2 = 0

      kg = 1
      icurrent = 0

      do while (kg .le. neg)
*        processing a  constraint
         i = irowg(kg)
         icurrent = i

*        get the current constraint from idg
         call dspxsq(iwork, dwork, nnvar,
     &        idg(i), ksize, info2)
         if (info2 .ne. 0) then
            info = 1
            return
         end if

*        jg is the variable index in the context of g, jtg is the varible index
*        in the context of idg.
         jg     = jcolg(kg)
         jtg    = jtvar(jg)
         kidg   = 1

         do while ( i .eq. icurrent )
            if ( kidg .gt. ksize ) then
*              There are no more elements left in idg, but there are
*              still elements in g
               g(kg) = zero
            else
*              idg is not empty
               jtidg = iwork(kidg)
               if ( jtg .eq. jtidg ) then
*                 The element in g and the element in idg match up exactly.
                  g(kg)  = dwork(kidg)
                  kidg   = kidg + 1
               else if ( jtg .lt. jtidg ) then
*                 This element of g is not in idg.
                  g(kg) = zero
               else
*                 jtg .gt. jtidg.  This element of idg is not in g. The
*                 sparsity pattern of idg is not a subset of the
*                 sparsity pattern of g.
                  info = 1
                  return
               end if
            end if
            kg  = kg + 1
            if (kg .le. neg) then
               i   = irowg(kg)
               jg  = jcolg(kg)
               jtg = jtvar(jg)
            else
               i   = icurrent + 1
            end if
         end do

         if ( kidg .lt. ksize ) then
*           There are still elements left in idg for this
*           constraint. The sparsity pattern of idg is not a subset of
*           the sparsity pattern of g.
            info = 1
            return
         end if
*     next constraint
      end do

      info = 0

      end ! subroutine spovrl

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine  spcsp( ig1, m, n, ig2,
     &             iirowc, ijcolc, nec, ic,
     &             iirowg, ijcolg, neg, f, v,
     &             itcon, jtvar, mncon, nnvar,
     &             adtlc, ladtl,
     &             ia, lenia, lusdia,
     &             ra, lenra, lusdra,
     &             morei, morer )

      implicit            none

      integer             m, n, ig1(m), ig2(m)
      double precision    f(m), v(n)
      integer             jtvar(n), itcon(m), nnvar, mncon
      integer             lenia,  ia( lenia ), lusdia
      integer             lenra, lusdra
      double precision    ra( lenra )

      integer             iirowc, ijcolc, nec, ic
      integer             iirowg, ijcolg, neg
      integer             morei, morer

*     Compute the sparsity pattern from the sparse matrices in ig1 and ig2

      integer             linear
      parameter           ( linear = 0 )

      integer             sentnl
      parameter           ( sentnl = -1 )

      integer i, j, k
      integer ijcol1, nrg1, irg1
      integer ijcol2, nrg2, irg2
      integer isppat, igpat
      integer igpatk, ispptk, info, nptdel, isppto, igpato
      integer lpat, nepat, necadl

      external            adtlc, ladtl
      integer             ladtl

      morei = 0
      morer = 0
*     Everything linear by default, will be changed later as needed
      do i = 1, m
         itcon(i) = linear
      end do
      do j = 1, n
         jtvar(j) = linear
      end do

      ijcol1 = lusdia + 1
      ijcol2 = ijcol1 + n
      irg1   = lusdra + 1
      irg2   = irg1   + n
      isppat  = ijcol2 + n
*     Check memory allocation
      if ( isppat .gt. lenia .or. irg2 + n .gt. lenra ) then
         morei = max( ijcol2 + n - lenia, 0 )
         morer = max( irg2   + n - lenra, 0 )
         goto 900
      end if

*     The elements ia(isppat:end) will represent the sparsity pattern of
*     the matrix. Each row in the matrix will be represented as an
*     acending sequence of column numbers, terminated by a sentinal
*     value (currently -1) in other words, the sequence 1, 3 6, -1, 2,
*     3, -1 represents the pattern (1,1), (1,3), (1,6), (2, 2), (2,3).
*     ispptk is the first unused element in ia(isppat:end)
      ispptk = isppat

*     ra(igpat:end) are the values corresponding to the indices in
*     ia(isppat:end)
      igpat     = irg2   + n
      igpatk    = igpat
*     igpat and isppat must be the same size. We call that size lpat
      lpat = min( lenia - isppat + 1, lenra - igpat + 1)
*     nepat is the number of elements actually used
      nepat = 0

      do i = 1, m
         call dspxsq( ia(ijcol1), ra(irg1), n, ig1(i), nrg1, info )
         call dspxsq( ia(ijcol2), ra(irg2), n, ig2(i), nrg2, info )
*        Merge the sparsity patterns in jcol1 and jcol2 into jcolc

         if ( lpat - nepat .ge. max( nrg1, nrg2 ) + 1 ) then
*           We have enough memory to attempt the merge
            call spcmrg( ia(ijcol1), nrg1, ra(irg1),
     &           ia(ijcol2), nrg2, ra(irg2),
     &           ia(ispptk), lpat - nepat, nptdel,
     &           ra(igpatk), itcon(i), jtvar, n,
     &           info )

            if ( info .eq. 0 ) then
*              We had enough elements in the pattern
               nepat = nepat + nptdel + 1
            else
*              We need at least two more elements in the pattern, one to
*              hold a column index, and one to hold the sentinal
               nepat = lpat + 2
            end if
         else
*           We don't have enough memory to complete the merge
            nepat = nepat + max( nrg1, nrg2 ) + 1
         end if

         if ( nepat .le. lpat ) then
*           spcmrg was successful and we have enough memory to complete
*           processing this row. Make the changes  permanent. First
*           add a sentinal value (in this case -1) to mark the end of the
*           row
            ia(ispptk + nptdel) = sentnl
*           Move ispptk and igpatk to the next unused element.
            ispptk = ispptk + nptdel + 1
            igpatk = igpatk + nptdel + 1
         else
*           We don't have enough memory.
            morei = max( 0, isppat + nepat - 1 - lenia )
            morer = max( 0, igpat  + nepat - 1 - lenra )
            goto 900
         end if
      end do
*     Count the true number of constant and non-linear elements of the
*     matrix.
      call spccnt ( itcon, m, jtvar, n, ia(isppat), nepat,
     &     nec, neg )

*     Save the old location of isppat
      isppto = isppat
*     Allocate irowg, jcolg, and irowc (overwriting ijcol1 and ijcol2)
      iirowg = lusdia + 1
      ijcolg = iirowg + neg
      iirowc = ijcolg + neg
*     Move isppat to a new location
      isppat = iirowc + nec
*     but first check for memory
      if ( isppat + nepat - 1 .gt. lenia ) then
         morei = isppat + nepat - 1 - lenia
         goto 900
      end if
      call mvi    ( ia, isppto, isppat, nepat )

*     Move igpat to the first unused element in ra (overwriting irg1 and irg2)
*     This is necessarily a shift left (i.e. igpato > igpat)
      igpato = igpat
      igpat  = lusdra + 1
      do i = 0, nepat - 1
         ra(igpat + i) = ra(igpato + i)
      end do

*     At this point, all of the sparsity structure is still encoded in
*     isppat and igpat. Call spreed to fill in appropriate values in irowg,
*     jcolg and irowc.
*     Set ic and ijcolc to point at the sparsity pattern of the matrix
      ijcolc = isppat
      ic     = igpat
*     on exit from spreed, jcolc and c will contain the correct values
*     (i.e. (irowc, jcolc, c) will correctly represent the sparse
*     constants.)
      call spreed ( itcon, m, jtvar, n,
     &     ia(iirowc), nec, ia(ijcolc), nepat, ra(ic),
     &     ia(iirowg), neg, ia(ijcolg) )

*     Count the number of non-linear variables and non-linear
*     constraints
      mncon = 0
      do i = 1, m
         if ( itcon(i) .ne. linear ) then
            mncon = mncon + 1
         end if
      end do

      nnvar = 0
      do i = 1, n
         if ( jtvar(i) .ne. linear ) then
            nnvar = nnvar + 1
         end if
      end do

      lusdia = ijcolc + nec - 1
      lusdra = ic     + nec - 1


*     Put the constants of the affine constraints into f
      do i = 1, m
         if ( itcon(i) .ne. 0 ) then
*           This is not an affine constraint
            f(i) = 0.0d0
         end if
      end do

      do k = 0, nec - 1
         i = ia( iirowc + k )
         if ( itcon(i) .eq. 0 ) then
*           This is an affine constraint
            j = ia( ijcolc + k )
            f(i) = f(i) - v(j) * ra(ic + k)
         end if
      end do


*     Get any additional constants that the user may have defined.  The
*     variable nec is at most n * m, where length(c) is 2 * m * n, so
*     c(nec + 1) exists.
      necadl = ladtl()
      if ( necadl .gt. 0 ) then
         lusdia = iirowc + 2 * (nec + necadl) - 1
         if ( lusdia .gt. lenia ) goto 900
*        Move ijcolc
         do i = nec - 1, 0, -1
            ia( iirowc + nec + necadl + i ) = ia( ijcolc + i)
         end do
         ijcolc = iirowc + nec + necadl
         call adtlc( ia( iirowc + nec ), necadl, ia( ijcolc + nec),
     &        ra(ic + nec) )
         nec = nec + necadl
         call lexsrt( ia(iirowc), nec, ia(ijcolc), ra(ic) )
      end if

      return
 900  continue

      end ! subroutine spcsp

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine spreed ( itcon, m,   jtvar, n,
     &                    irowc, nec, jcolc, lenc, c,
     &                    irowg, neg, jcolg )

*     Split a sparse matrix into linear and non-linear
*     parts.
*
*     On exit, the triple (irowc, jcolc, c) will contain the linear
*     (constant) part of the matrix and (irowg, jcolg) will contain the
*     sparsity pattern of the linear part.
*
*     On entry, jcolc and c contain a quite different structure. The
*     elements of jcolc contain the column indicies of each row of c,
*     sorted in ascending order. Each row of c is separated by a
*     sentinal value, which is currently -1. For instancet he sequence
*     1, 3 6, -1, 2, 3, -1 represents the pattern (1,1), (1,3), (1,6),
*     (2, 2), (2,3). The element c(k) is the value corresponding to the
*     index represented by jcolc(k). T
*
*     On entry, itcon(i) is zero if the i^th function contains only
*     linear elements, and jtvar(j) is zero if the j^th variable only
*     occurs linearly.
*
*     nec and neg are the number of linear and non-linear elements,
*     respectively. These may be counted using spccnt.

      implicit            none

      integer             m, itcon(m), n, jtvar(n)
      integer             nec, irowc(nec), lenc, jcolc(lenc)
      double precision    c(lenc)
      integer             neg, irowg(neg), jcolg(neg)

      integer             linear
      parameter           ( linear = 0 )

      integer             i, j, k, kg, kc

*     kg - the first unused element of ( irowg, jcolg )
      kg = 1
*     To save memory, we overwrite the input sparsity pattern
*     ( jcolc, c ) with the column and value information of the triple
*     ( irowc, jcolc, c ). This is safe because we are always able to
*     write the current output element in ( jcolc, c ) to a lower index (kc)
*     than the index in ( jcolc, c ) from which we read the current
*     input element (k). This is a little confusing, but works quite
*     well. Just remember ( or convince yourself ), than kc <= k.
      kc = 1
      i  = 1

      do k = 1, lenc
         j = jcolc(k)
         if ( j .lt. 0 ) then
*           We must start a new row
            i = i + 1
         else
            if ( jtvar(j) .eq. linear .or. itcon(i) .eq. linear ) then
*              This is a linear (constant) element.
               irowc(kc) = i
               jcolc(kc) = j
               c(kc)     = c(k)
               kc        = kc + 1
*               print *, 'Const', irowc(kc-1), jcolc(kc-1), c(kc-1)
            else
*              This is a non-linear element.
               irowg(kg) = i
               jcolg(kg) = j
               kg        = kg + 1
*               print *, 'Spars', irowg(kg-1), jcolg(kg-1)
            end if
         end if
      end do

      end ! subroutine spreed

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      subroutine spccnt( itcon, m, jtvar, n, jcolc, lenc,
     &     nec, neg )

      implicit            none

      integer              m, itcon(m), n, jtvar(n), lenc, jcolc(lenc)
      integer              nec, neg

*     Count the number of linear and non-linear elements

      integer              linear
      parameter            ( linear = 0 )
      integer              i, j, k

      i = 1
      nec = 0
      neg = 0
      do k = 1, lenc
         j = jcolc(k)
         if ( j .lt. 0 ) then
            i = i + 1
         else
            if ( jtvar(j) .eq. linear .or. itcon(i) .eq. linear ) then
               nec = nec + 1
            else
               neg = neg + 1
            end if
         end if
      end do

      end ! subroutine spccnt

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine spcmrg( jcol1, nec1, g1,
     &     jcol2, nec2, g2,
     &     jcolc, lenc, nec, c,
     &     itconi, jtvar, n,
     &     info )

*     Form the union of the sparsity pattern in (jcol1, g1) and (jcol2,
*     g2). The results of the merge are stored in (jcolc, c).
*
*     The number of elements in (jcol1, g1) is nec1 and is known on
*     entry. Similarly nec2 is known on entry, but nec, the number of
*     elements in (jcolc, c), is computed within this routine.
*
*     jtvar(n) - if element k of the merged vector is found to be
*     non-linear, then jtvar(jcolc(k)) will be set to 1. Otherwise
*     jtvar(jcolc(k) will be left untouched.
*
*     itconi - if any element of the merged vector is found to be
*     non-linear then itconi will be set to 1. Otherwise it will be
*     untouched.
*
*     If info .eq. 0 then the merge was successful. Otherwise,
*     completing the merge would have required more than lenc elements
*     in jcolc and c.

      implicit            none

      integer             nec1, jcol1(nec1)
      double precision    g1(nec1)

      integer             nec2, jcol2(nec2)
      double precision    g2(nec2)

      integer             lenc, nec, jcolc(lenc)
      double precision    c(lenc)

      integer             itconi, n, jtvar(n)

      integer             info

      integer             nonlin
      parameter           ( nonlin = 1 )
      double precision    zero
      parameter           ( zero = 0.0d0 )

      integer             k, k1, k2

*     Set info so that by default the routine fails.
      info = 1
      nec  = 0
      k1 = 1
      k2 = 1
      do k = 1, lenc
*           Loop over all available space in (jcolc, c). If we exhaust,
*           (jcol1, c1) and (jcol2, g2) before running out of space in
*           (jcolc, c) then we quit the loop early.
         if ( k1 .gt. nec1 .and. k2 .gt. nec2 ) then
*           Sucess! We've exhausted all the elements in (jcol1, g1) and
*           (jcol2, g2).
            info = 0
            goto 600
         else if ( k1 .gt. nec1 ) then
*           We've exhausted all the elements in jcol1. Pull off a
*           (necessarily non-linear) element from jcol2
            jcolc(k)         = jcol2(k2)

            jtvar(jcol2(k2)) = nonlin
            itconi           = nonlin
            c(k)             = zero

            k2               = k2 + 1
         else if ( k2 .gt. nec2 ) then
*           We've exhasted all the elements in jcol2. Pull off a
*           (necessarily non-linear) element from jcol1
            jcolc(k)         = jcol1(k1)

            jtvar(jcol1(k1)) = nonlin
            itconi           = nonlin
            c(k)             = zero

            k1               = k1 + 1
         else if ( jcol1(k1) .eq. jcol2(k2) ) then
*           jcol1 and jcol2 both contain an element with the same index
            jcolc(k) = jcol1(k1)
            if( g1(k1) .eq. g2(k2) ) then
*              The element is linear
               c(k)  = g1(k1)
            else
*              The element is non-linear
               jtvar(jcol1(k1)) = nonlin
               itconi           = nonlin
               c(k)             = zero
            end if
            k1 = k1 + 1
            k2 = k2 + 1
         else if ( jcol1(k1) .lt. jcol2(k2) ) then
*           The next element in the merge is in jcol1. It is necessarily
*           non-linear because there is no corresponding element in
*           jcol2.
            jcolc(k)         = jcol1(k1)

            jtvar(jcol1(k1)) = nonlin
            itconi           = nonlin
            c(k)             = zero

            k1               = k1 + 1
         else
*           jcol2(k2) .lt. jcol1(k1)
*           The next element in the merge is in jcol2. It is necessarily
*           non-linear because there is no corresponding element in
*           jcol1.
            jcolc(k)         = jcol2(k2)

            jtvar(jcol2(k2)) = nonlin
            itconi           = nonlin
            c(k)             = zero

            k2               = k2 + 1
         end if

         nec = nec + 1
      end do
      if ( k1 .gt. nec1 .and. k2 .gt. nec2 ) then
         info = 0
      end if
 600  continue

      end ! subroutine spcmrg












