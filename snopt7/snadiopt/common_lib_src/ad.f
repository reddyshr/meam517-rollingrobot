*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     File  ad.f
*     Routines needed by both sparse and dense versions of  Adifor.
*
*     n2admx   stadmx   gtadmx   zrlin    rndprb   gtrdsd   lexsrt
*     mvi      mvr
*
*     Written by Mike Gertz - 28-Jul-99
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine n2admx
     &   ( kind, lext, ia, lenia, isize, morei )

      implicit
     &     none
      integer
     &     kind, lext, lenia, ia(lenia), isize, morei

*     ==================================================================
*     Create a new admx, which contains info about a Jacobian matrix,
*     along with data structures needed for the automatic
*     differentiation program ADIFOR.
*
*     As ADIFOR can do both sparse and dense derivatives, this
*     routine must be supplied a value for the ikind field. It is
*     recommended, however, that a user *not* call this routine
*     directly, but rather call nwddmx to use dense ADIFOR or nwsdmx to
*     use sparse ADIFOR. These routines know the appropriate value of
*     the ikind parameter, and can do any other needed initialization.
*
*     The admx datatype is created as an extension to the mtx datatype.
*     Therefore, this routine will create an "mtx" in ia with the
*     additional fields as its extension. See gtadmx for a description
*     of how to get its fields.
*     ==================================================================
      integer
     &     iiself
*     ------------------------------------------------------------------
      integer       mysize
      parameter    (mysize = 8)
*     ------------------------------------------------------------------

      call n2wmtx
     &   ( kind, lext+mysize, ia, lenia, isize, morei )

      if ( morei .eq. 0 ) then
         iiself = isize - (lext + mysize)
*        ia(iiself + 1) is  iiext
         ia(iiself + 1) =   mysize + 1
         ia(iiself + 2) =   lext
      end if

      end ! subroutine n2admx

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine stadmx
     &   ( iitcon, mncon, ijtvar, nnvar, ig_v, ig_cf, ia, lenia )

      implicit
     &     none
      integer
     &     iitcon, mncon, ijtvar, nnvar,
     &     ig_v, ig_cf, lenia, ia(lenia)

*     ==================================================================
*     Set the fields of an instance of the admx datatype. See gtadmx for
*     a description of the arguments of this subroutine.
*
*     The admx datatype is an extension of the mtx datatype. ia should
*     be an array which was previously passed to n3admx, nwddmx, nwsdmx
*     or some similar allocation subroutine.
*     ==================================================================
      integer
     &     iiadmx
*     ------------------------------------------------------------------

      iiadmx =   ia(1)
*     ia(iiadmx + 0)  is  iiext - which is only set in nwadmx
*     ia(iiadmx + 1)  is  lext  - which is only set in nwadmx
      ia(iiadmx + 2)  =  ig_v
      ia(iiadmx + 3)  =  ig_cf
      ia(iiadmx + 4)  =  iitcon
      ia(iiadmx + 5)  =  mncon
      ia(iiadmx + 6)  =  ijtvar
      ia(iiadmx + 7)  =  nnvar

      end ! subroutine stadmx

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine gtadmx
     &   ( iitcon, mncon, ijtvar, nnvar, ig_x, ig_cf, ia, lenia )

      implicit
     &     none
      integer
     &     iitcon, mncon, ijtvar, nnvar, ig_x, ig_cf, lenia, ia(lenia)

*     ==================================================================
*     Get all elements of the admx datatype.
*
*     iitcon    - The index in ia of the integer array itcon(0:m).
*                 Constraint k is linear if and only if itcon(k) is zero.
*
*     mncon     - The number of non-linear constraints.
*
*     ijtvar    - The index in ia of the integer array jtvar(nv), where nv
*                 is the number of variables. The variable v(k) occurs only
*                 linearly if and only if jtvar(k) .eq. 0. If the variable v(k)
*                 occurs non-linearly, then g_v(jtvar(k),:) is the differential
*                 corresponding to v(k). (See below for a definition of g_v).
*
*     nnvar     - The number of variables that occur non-linearly.
*
*     ig_v      - The index (pointer into a double precision work array)
*                 of the differentials of the non-linear independent
*                 variables. Because this matrix is to be used bu Adifor, it is
*                 transposed. The shape of g_v is g_v(nnvar, nv) where nv is
*                 the number of variables.
*
*     ig_cf     - The index (pointer into a double precision work array) of
*                 the differentials of the dependent variables.  The matrix
*                 g_cf will be filled by a call to Adifor-generated code.  The
*                 shape of g_cf is g_cf(nnvar, 0:m).
*
*     ia(lenia) - An array containing the fields of an admx datatype.
*
*     The admx datatype is an extension of the mtx datatype. ia should
*     be an array which was previously passed to n3admx, nwddmx, nwsdmx
*     or some similar allocation subroutine.
*     ==================================================================
      integer
     &     iiadmx
*     ------------------------------------------------------------------

      iiadmx  = ia(1)

*     iiext is  ia(iiadmx + 0)
*     lext  is  ia(iiadmx + 1)
      ig_x    = ia(iiadmx + 2)
      ig_cf   = ia(iiadmx + 3)
      iitcon  = ia(iiadmx + 4)
      mncon   = ia(iiadmx + 5)
      ijtvar  = ia(iiadmx + 6)
      nnvar   = ia(iiadmx + 7)

      end ! subroutine gtadmx

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine zrlin
     &   ( jtvar, nv, v )

      implicit
     &     none
      integer
     &     nv, jtvar(nv)
      double precision
     &     v(nv)
*     ==================================================================
*     Makes v(j) = 0  if  jtvar(j) = 0.  In the context in which this
*     routine is called this "zeros" out the "linear" variables of v.
*     ==================================================================
      integer
     &     j
*     ------------------------------------------------------------------
      integer            linear
      parameter         (linear = 0)
      double precision   zero
      parameter         (zero   = 0.0d0)
*     ------------------------------------------------------------------

      do j = 1, nv
         if ( jtvar(j) .eq. linear ) then
            v(j) = zero
         end if
      end do

      end ! subroutine zrlin

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine rndprb
     &   ( n, x, bl, bu, order, x2, iseed )

      implicit
     &     none
      integer
     &     n, iseed(3)
      double precision
     &     x(n), x2(n), bl(n), bu(n), order

*     ==================================================================
*     Let x2 be a random perturbation of x.  The order of the
*     perturbation will be at most "order" and the elements of x2 will
*     lie between bl and bu.
*
*     iseed1, iseed2, iseed2 are seed values for the random number
*     generator.
*     ==================================================================
      integer
     &     j, direction
      double precision
     &     theta
*     ------------------------------------------------------------------
      integer            down,      fixed,     up,     free
      parameter         (down = -1, fixed = 0, up = 1, free = 2)
      double precision   half        , two
      parameter         (half = 0.5d0, two = 2.0d0)
*     ------------------------------------------------------------------

*     Get a vector of pseudo-random numbers between 0 and 1

      call ddrand
     &   ( n, x2, 1, iseed )

      do j = 1, n

*        If direction .eq. free, then x(j) is not at one of it's bounds.

         direction = free
         if ( x(j) .eq. bl(j) ) then

*           x(j) is on its lower bound, we must perturb up.

            direction = up
         end if

         if ( x(j) .eq. bu(j) ) then

*           x(j) is on its upper bound.

            if ( direction .eq. free ) then

*              but is not on its lower bound, so we perturb down

               direction = down
            else

*              ... and is also on its lower bound, and so it must be
*              a fixed variable.  It is not perturbed.

               direction = fixed
            end if
         end if

         if ( direction .eq. free ) then

*           If the variable is not on either bound, choose randomly
*           whether to go up or down.

            if ( x2(j) .gt. half ) then
               x2(j) = two * ( x2(j) - half )
               direction = up
            else
               x2(j) = two * x2(j)
               direction = down
            end if
         end if

         theta = half * x2(j)
         if ( direction .eq. down ) then

*           Perturb between 1/8 and 1/4 of the way to the lower bound,
*           but no more than "order".

            x2(j) = theta * x(j) +
     &           (1 - theta) * max( ( 3*x(j) + bl(j) )/4, x(j) - order )
         else if (direction .eq. up ) then

*           Perturb between 1/8 and 1/4 of the way to the upper bound,
*           but no more than "order".

            x2(j) = theta * x(j) +
     &           (1 - theta) * min( ( 3*x(j) + bu(j) )/4, x(j) + order )
         else

*           direction .eq. fixed

            x2(j) = theta * x(j) + (1 - theta) * order
         end if

      end do

      end ! subroutine rndprb

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine gtrdsd
     &   ( iseed1, iseed2, iseed3 )

      implicit
     &     none
      integer
     &     iseed1, iseed2, iseed3

*     ==================================================================
*     Get seed values for a pseudo-random number generator.
*
*     REPLACE ME - This version just returns the same values
*     every time.
*
*     ==================================================================

      iseed1 = 1547
      iseed2 = 2671
      iseed3 = 3770

      end ! subroutine gtrdsd

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lexsrt
     &   ( irowC, neC, jcolC, C )

      implicit
     &     none
      integer
     &     neC, irowC(neC), jcolC(neC)
      double precision
     &     C(neC)

*     ==================================================================
*     Do a lexicographic sort on (irowC x jcolC).
*
*     Implemented using Shell's sort with Sedgewick's intervals.
*
*     (Other reasonable (or unreasonable) choices of intervals are 2^k-1
*     (Hibbard's intervals) and neC/2, neC/4, etc. (Shell's intervals)
*     ==================================================================
      double precision
     &     Ctemp
      integer
     &     inc, kinc, k, j, itemp, jtemp
*     ------------------------------------------------------------------
*     Set incs, the majic numbers.
      integer            nincs
      parameter         (nincs = 12)
      integer            incs(nincs)
      data               incs /1, 5, 19, 41, 109, 209, 505,
     &                         929, 2161, 3905, 8929, 16001/
      save               incs
*     ------------------------------------------------------------------
      inc = 1

      do k = 1, nincs
         kinc = k
         if ( incs(kinc) .gt. neC/2 ) then
            kinc = kinc - 1
            goto 10
         end if
      end do
 10   continue

*     incs(kinc) is the greatest value in the sequence that is also less than
*     or equal to neC/2

      do while (kinc .gt. 0)
*        Loop over all increments
         inc = incs(kinc)
         do k = inc + 1, neC

*           Loop over all subarrays defined by the current increment

            itemp = irowC(k)
            jtemp = jcolC(k)
            Ctemp = C(k)
            j     = k

*           Insert element k into the sorted subarray.

            do while (j .gt. inc)

*              Loop over the elements in the current subarray

               if (itemp .lt. irowC(j-inc) .or.
     &              (itemp .eq. irowC(j-inc) .and.
     &               jtemp .lt. jcolC(j-inc) ) ) then

*                 Swap element j and j - inc (implicitly use the fact that
*                 tmp holds element j to avoid having to assign to element
*                 j-inc)

                  irowC(j) = irowC(j-inc)
                  jcolC(j) = jcolC(j-inc)
                  C(j)     =     C(j-inc)
               else
*                 There are no more elements in this sorted subarray
*                 which are less than element j
                  goto 100
               end if
               j = j - inc
*           end loop over the elements in the current subarray
            end do
 100        continue

            if ( j .ne. k ) then

*              Move element j out of temporary storage

               irowC(j) = itemp
               jcolC(j) = jtemp
               C(j)     = Ctemp
            end if
*        end loop over all subarrays defined by the current increment
         end do
         kinc = kinc - 1
*     end loop over all increments
      end do

      end ! subroutine lexsrt

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine mvi
     &     ( ia, iold, inew, n )

      implicit
     &     none
      integer
     &     ia(*), iold, inew, n

*     ==================================================================
*     Move n contiguous integers from ia(iold) to ia(inew).
*     ==================================================================
      integer
     &     i
*     ------------------------------------------------------------------

      if ( iold .lt. inew ) then
         do i = n - 1, 0, -1
            ia( inew + i ) = ia( iold + i )
         end do
      else if ( inew .lt. iold ) then
         do i = 0, n - 1
            ia( inew + i ) = ia( iold + i )
         end do
      end if

      end ! subroutine mvi

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine mvr
     &   ( ra, iold, inew, n )

      implicit
     &     none
      integer
     &     iold, inew, n
      double precision
     &     ra(*)

*     ==================================================================
*     Move n contiguous doubles from ra(iold) to ra(inew).
*     ==================================================================
      integer
     &     i
*     ------------------------------------------------------------------
      if ( iold .lt. inew ) then
         do i = n-1, 0, -1
            ra( inew + i ) = ra( iold + i )
         end do
      else if ( inew .lt. iold ) then
         do i = 0, n - 1
            ra( inew + i ) = ra( iold + i )
         end do
      end if

      end ! subroutine mvr
