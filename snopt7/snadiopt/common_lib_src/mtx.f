*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*     mtx.f
*     Mike Gertz 27-Jul-99
*
*     Contains the routines needed to access the "mtx" data type, as
*     well as routines that are closely related to that type.
*
*     subroutines:
*
*     newmtx _______ n2wmtx
*     getmtx
*     kndmtx
*     setmtx
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine newmtx( ia, lenia, isize, moreia )
      integer            lenia, ia(lenia), isize, moreia

*     ==================================================================
*     Allocate a new instance of the "mtx" data type in the array
*     ia(lenia).  This type represents the constant elements and
*     sparsity pattern of a Jacobian matrix.
*
*     type. If you are calling this routine directly, you will typically
*     set lext to zero.
*
*     isize  - The total number of elements allocated.
*     moreia - Zero if the allocation was successful. Otherwise at least
*              moreia more elements of ia are needed.
*     ==================================================================
      integer            kind,     lext
      parameter         (kind = 0, lext = 0)

      call n2wmtx( kind, lext, ia, lenia, isize, moreia )

      end ! subroutine newmtx

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine n2wmtx( kind, lext, ia, lenia, isize, moreia )
      implicit           none

*     ==================================================================
*     This routine is similar to the newmtx routine. The only difference
*     is that it allows the user to set the value of kind and lext.
*     Unless you are writing an extension to the mtx datatype, you
*     should not need to call this routine.
*
*     lext - reserve lext elements of ia for a extension to the
*     datatype.  For example, if an automatic differentiation scheme is
*     used, one might reserve lext elements to hold the data needed by
*     the automatic differentiation routines.
*
*     kind - Specify what "kind" of mtx this is. Every extension to the
*     mtx datatype will have a different value for kind.
*     ==================================================================
      integer            kind, lext, lenia, ia(lenia), isize, moreia
      integer            lmtx
      parameter         (lmtx = 11)

      isize = lmtx + lext
      if ( isize .gt. lenia ) then
         moreia = isize - lenia
      else
*        ia(1) is iiext
         ia(1) = lmtx + 1
         ia(2) = lext
         ia(3) = kind
         moreia = 0
      end if

      end ! subroutine n2wmtx

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine getmtx( iirowc, nec,   ijcolc, ic,
     &     iirowg, neg, ijcolg, m, ia, lenia )

      implicit
     &     none

      integer
     &     iirowc, nec, ijcolc, ic,
     &     iirowg, neg, ijcolg, m,
     &     lenia, ia(lenia)

*     ==================================================================
*     Gets all the components of an instance of the "mtx" data type.
*     This type represents the constant elements and sparsity pattern
*     of a Jacobian matrix.
*
*     iirowc - The index in ia of the integer array irowc(nec).
*     ijcolc - The index in ia of the integer array jcolc(nec).
*     ic - The index in a double precision array of the array c(nec).
*
*     The triple (irowc(k), jcolc(k), c(k)) represents the k^th constant
*     element of the Jacobian.
*
*     iirowg - The index in ia of the integer array irowg(neg).
*     ijcolg - The index in ia of the integer array jcolg(neg).
*
*     The pair (irowg(k), jcolg(k)) represents the location of the k^th
*     non-zero, non-constant element of the Jacobian matrix.
*
*     m - the highest index row of the Jacobian matrix.
*
*     ia(lenia) - An array containing an instance of this data type
*     starting at ia(1).
*
*
*     Note that iiext/lext are not set by setmtx, but by newmtx (i.e.
*     the space is reserved during the allocation of the "mtx".)
*     The value of kind is also set when the mtx is created.
*     ==================================================================


*     iiext   is ia(1)
*     lext    is ia(2)
*     kind    is ia(3)
      iirowc  =  ia(4)
      nec     =  ia(5)
      ijcolc  =  ia(6)
      ic      =  ia(7)
      iirowg  =  ia(8)
      neg     =  ia(9)
      ijcolg  =  ia(10)
      m       =  ia(11)

      end ! subroutine getmtx

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      function kndmtx( ia, lenia )

      implicit     none
      integer      kndmtx, lenia, ia(lenia)

*     ==================================================================
*     kndmtx - What "kind" of mtx is this.
*     ==================================================================
      kndmtx = ia(3)

      end ! function kndmtx

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine setmtx( iirowc, nec, ijcolc, ic,
     &                   iirowg, neg, ijcolg, m,
     &                   ia, lenia )

      implicit           none

      integer            iirowc, nec, ijcolc, ic
      integer            iirowg, neg, ijcolg, m
      integer            lenia, ia(lenia)

*     ==================================================================
*     Sets all the components of an instance of the "mtx" data type.
*     This type represents the constant elements and sparsity pattern
*     of a Jacobian matrix.
*
*     See getmtx for a description of this data type.
*     ==================================================================
*     ia(1)   is iiext
*     ia(2)   is lext
*     ia(3)   is kind
      ia(4)   =  iirowc
      ia(5)   =  nec
      ia(6)   =  ijcolc
      ia(7)   =  ic
      ia(8)   =  iirowg
      ia(9)   =  neg
      ia(10)  =  ijcolg
      ia(11)  =  m

      end ! subroutine setmtx

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c$$$      subroutine cntfg( nState, mode,
c$$$     &     x, nx, p, np, t, ndot,
c$$$     &     f, mCon,
c$$$     &     irowg, neg, jcolg, g,
c$$$     &     v,
c$$$     &     cu, lencu, iu, leniu, ru, lenru,
c$$$     &     cw, lencw, iw, leniw, rw, lenrw,
c$$$     &     ia, lenia, ra, lenra,
c$$$     &     icnt, licnt, rcnt, lrcnt,
c$$$     &     iiext, lext, kind )
c$$$
c$$$*     Compute the function values and derivatives of the interior ("cnt"
c$$$*     for continuous interior) constraints.
c$$$*
c$$$*     nState, mode - flags set by snopt which are passed to user code.
c$$$*     See snopt documentation for details.
c$$$*
c$$$*     x(nx) - State and control variables. These variables have values
c$$$*     vary with each time step.
c$$$*
c$$$*     p(np) - "Parameters". These are variables used in the optimization
c$$$*     which are constant in time.
c$$$*
c$$$*     t     - The time value
c$$$*
c$$$*     ndot - The number of elements of x which are derivatives. The
c$$$*     variables x(1:ndot) are the derivatives of the variables
c$$$*     x(ndot+1:2*ndot).
c$$$*
c$$$*     f(0:mCon) - (output) On exit, f will contain the function values.
c$$$*
c$$$*     irowg(neg), jcolg(neg) - (input) The triple (irowg(k), jcolg(k),
c$$$*     g(k)) is the k^th element of the Jacobian matrix g. On entry, the
c$$$*     vectors irowg and jcolg contain the correct sparsity patern, which
c$$$*     must not be changed by cntfg.
c$$$*
c$$$*     g(neg) - (output) The Jacobian elements. These elemenst are to be
c$$$*     computed by cntfg.
c$$$*
c$$$*     v - A scratch array that is at least as long as the
c$$$*     concatenation of x, p and t.
c$$$*
c$$$*     [cir][uw](len[cir][uw]) - Work arrays that must be passed to the
c$$$*     user's computational routine.
c$$$*
c$$$*     [ir]cnt(l[ir]cnt) - Global variables and work space needed by this
c$$$*     routine.
c$$$
c$$$      implicit none
c$$$
c$$$      integer             nState, mode
c$$$      integer             nx, np, ndot
c$$$      double precision    x(nx), p(0:np), t
c$$$      integer             mCon, neg, irowg(neg), jcolg(neg)
c$$$      double precision    f(0:mCon),  g(neg)
c$$$      double precision    v( 2*(nx + 1) + np )
c$$$
c$$$      integer             lencu, lencw
c$$$      integer             leniu, leniw, lenia, licnt
c$$$      integer             lenru, lenrw, lenra, lrcnt
c$$$      character*8         cu(lencu), cw(lencw)
c$$$      integer             iu(leniu), iw(leniw), ia(lenia), icnt(licnt)
c$$$      double precision    ru(lenru), rw(lenrw), ra(lenra), rcnt(lrcnt)
c$$$
c$$$      integer             iirwcc, necc, ijclcc, icc
c$$$      integer             iirwgc, negc, ijclgc, mcnt
c$$$      integer             kind, iiext, lext
c$$$
c$$$      if ( kind .eq. 1 ) then
c$$$         call dectfg( nState, mode,
c$$$     &        x, nx, p, np, t, ndot,
c$$$     &        f, mCon,
c$$$     &        irowg, neg, jcolg, g,
c$$$     &        v,
c$$$     &        cu, lencu, iu, leniu, ru, lenru,
c$$$     &        cw, lencw, iw, leniw, rw, lenrw,
c$$$     &        iirwcc, necc, ijclcc, icc,
c$$$     &        iirwgc, negc, ijclgc, mcnt,
c$$$     &        ia, lenia, ra, lenra,
c$$$     &        icnt, licnt, rcnt, lrcnt, iiext, lext )
c$$$      else
c$$$         print *, 'Internal error in cntfg: unknown kind of problem'
c$$$         stop 1
c$$$      end if
c$$$
c$$$*     dispatch to the appropriate
c$$$
c$$$      end
c$$$
c$$$****************************************************************
c$$$      subroutine bndfg ( nState, mode,
c$$$     &     x, nx, xf, p, np, t0, dT, ndot,
c$$$     &     f,  m, irowg, neg, jcolg, g,
c$$$     &     v,
c$$$     &     cu, lencu, iu, leniu, ru, lenru,
c$$$     &     cw, lencw, iw, leniw, rw, lenrw,
c$$$     &     ia, lenia, ra, lenra,
c$$$     &     ibnd, libnd, rbnd, lrbnd )
c$$$*     Compute the function values and derivatives of the boundary
c$$$*     ("bnd") constraints.
c$$$*
c$$$*     nState, mode - flags set by snopt which are passed to user code.
c$$$*     See snopt documentation for details.
c$$$*
c$$$*     x(nx) - State and control variables at the initial time (left-hand
c$$$*     endpoint of the time interval.)
c$$$*
c$$$*     xf(nx) - State and control variables at the final time (right-hand
c$$$*     endpoint of the time interval.)
c$$$*
c$$$*     p(np) - "Parameters". These are variables used in the optimization
c$$$*     which are constant in time.
c$$$*
c$$$*     t0    - The initial time (left-hand endpoint of the time
c$$$*     interval.)
c$$$*
c$$$*     dT    - The width of the time interval.
c$$$*
c$$$*     ndot - The number of elements of x which are derivatives. The
c$$$*     variables x(1:ndot) are the derivatives of the variables
c$$$*     x(ndot+1:2*ndot).
c$$$*
c$$$*     f(0:m) - (output) On exit, f will contain the function values.
c$$$*
c$$$*     irowg(neg), jcolg(neg) - (input) The triple (irowg(k), jcolg(k),
c$$$*     g(k)) is the k^th element of the Jacobian matrix g. On entry, the
c$$$*     vectors irowg and jcolg contain the correct sparsity patern, which
c$$$*     must not be changed by cntfg.
c$$$*
c$$$*     g(neg) - (output) The Jacobian elements. These elemenst are to be
c$$$*     computed by cntfg.
c$$$*
c$$$*     v - A scratch array that is at least as long as the
c$$$*     concatenation of x, p and t.
c$$$*
c$$$*     [cir][uw](len[cir][uw]) - Work arrays that must be passed to the
c$$$*     user's computational routine.
c$$$*
c$$$*     [ir]bnd(l[ir]bnd) - Global variables and work space needed by this
c$$$*     routine.
c$$$
c$$$      implicit            none
c$$$
c$$$      integer             nState, mode
c$$$      integer             nx, np, ndot
c$$$      double precision    x(nx), xf(nx), p(np), t0, dT
c$$$      double precision    v( 2*(nx + 1) + np )
c$$$      integer             m, neg, irowg(neg), jcolg(neg)
c$$$      double precision    f(0:m), g(neg)
c$$$      integer             lencu, leniu, lenru
c$$$      integer             lencw, leniw, lenrw
c$$$      integer             lenia, lenra
c$$$      integer             libnd, lrbnd
c$$$      character*8         cu(lencu), cw(lencw)
c$$$      integer             iu(leniu), iw(leniw), ia(lenia), ibnd(libnd)
c$$$      double precision    ru(lenru), rw(lenrw), ra(lenra), rbnd(lrbnd)
c$$$
c$$$      integer             i1, i2, i3, i4, i5, i6, i7, i8, isize
c$$$      integer             iitcon, mncon, ijtvar, nnvar
c$$$      integer             ig_v, ig_cf
c$$$
c$$$      end
