!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     File  sn05wrpb.f  --- user function wrapper for snOptB (aka snOpt)
!
!     s0fgB
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s0fgB
     &   ( iExit,
     &     modefg, NonlinearCon, NonlinearObj,
     &     n, negCon, nnCon0, nnCon,
     &     nnJac, nnH, nnObj0, nnObj,
     &     fgcon, fgobj, dummyH,
     &     x, yCon,
     &     neJ, nlocJ, locJ, indJ,
     &     neH, nlocH, locH, indH, Hcol,
     &     fCon, fObj, gCon, gObj,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     fgcon, fgobj, dummyH
      logical
     &     NonlinearCon, NonlinearObj
      integer
     &     iExit, lencu, leniu, lenru, lencw, leniw, lenrw, modefg, n,
     &     negCon, neH, neJ, nlocH, nlocJ, nnCon0, nnCon, nnJac, nnH,
     &     nnObj0, nnObj, indH(neH), indJ(neJ), locH(nlocH),
     &     locJ(nlocJ),  iu(leniu), iw(leniw)
      double precision
     &     fObj, fCon(nnCon0), gObj(nnObj0), gCon(negCon),
     &     Hcol(neH), x(n), yCon(nnCon0), ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

!     ==================================================================
!     s0fgB  is an instance of funwrapper that calls the user-written
!     routines  fgcon  and  fgobj  to evaluate the problem functions
!     and possibly their gradients.
!
!     Arguments  fgcon  and  fgobj  are called using modefg to control
!     the gradients as follows:
!
!     modefg        Task
!     ------        ----
!       2     Assign fCon, fObj and all known elements of gCon and gObj.
!       1     Assign all known elements of gCon and gObj.
!             (fObj and fCon are ignored).
!       0     Assign fObj, fCon.  (gCon and gObj are ignored).
!
!     If s0fgB is called with minmax = 0 (feasible point only) then
!     nnObj = max(nnJac,nnObj)  and the user objective is not used.
!
!     09-Jan 1992: First version based on Minos 5.4 routine m6fobj.
!     28-Sep 1993: Test for mode = -1 exit from funobj.
!     22-Dec 1994: Simplified to conform to npfobj in npopt.
!     16-Jan 1997: fgobj added as argument.
!     19-Mar 1997: First instantiated version for SNOPT 5.1.
!     17-Jul 1997: First thread-safe version.
!     02-Oct 1997: Character workspace added.
!     11 Oct 1998: Logicals NonlinearCon and NonlinearObj added.
!     30 Dec 2000: Housekeeping for Status = 1 included.
!     03 Aug 2003: snEXIT and snPRNT adopted.
!     16 Jun 2008: Call-status implemented correctly.
!     15 Nov 2010: Call-status removed from the argument list.
!     11 Sep 2014: nnObjU and nnObj separated for FP mode.
!     07 Nov 2014: yCon added as argument.
!     13 Nov 2014: Scales applied for lvlScale > 0.
!     14 Feb 2015: proxWeight added.
!     ==================================================================
      character
     &     str*80
      external
     &     ddot
      logical
     &     FPonly, Scaled
      integer
     &     gotFD, gotG, conJij, l, lG, lgConU, lgObjU, lvlScale, lvlTim,
     &     lvlDer, lscales, lx0, lxScaled, minmax, modeC, modeF,
     &     nfCon1, nfCon2, nfObj1, nfObj2, nnGlin, nnObjU,
     &     supplied, status, statusUser
      double precision
     &     ddot, gDummy, proxWeight
!     ------------------------------------------------------------------
      parameter         (lvlDer =  70) ! = 0,1,2 or 3, deriv. level
      parameter         (gotFD  = 183) ! > 0 => some differences needed
      parameter         (gotG   = 184) ! > 0 => some exact derivs
      parameter         (conJij = 185) ! > 0 => constant Jacob elements
      parameter         (nfCon1 = 189) ! calls to fCon: mode = 0
      parameter         (nfCon2 = 190) ! calls to fCon  mode > 0
      parameter         (nfObj1 = 194) ! calls to fObj: mode = 0
      parameter         (nfObj2 = 195) ! calls to fObj: mode > 0

      double precision   half,            one
      parameter         (half   = 0.5d+0, one   = 1.0d+0)
!     ------------------------------------------------------------------
      nnObjU     = iw( 22) ! # of user-defined objective variables
      lvlScale   = iw( 75) ! scale option
      minmax     = iw( 87) ! 1, 0, -1  => MIN, FP, MAX
      lvlTim     = iw(182) ! Timing level

      lscales    = iw(296) ! scales(nb) = row and column scales
      lx0        = iw(298) ! x0(nnL)    = FP base point
      lxScaled   = iw(302) ! xScaled(n) = copy of scaled  x
      lgConU     = iw(319) ! record of unknown derivatives and constants
      lgObjU     = iw(323) ! record of unknown derivatives

      gDummy     = rw( 69) ! definition of an 'unset' value
      proxWeight = rw( 91) ! Proximal-point weight

      iExit      = 0

      FPonly     = minmax   .eq. 0
      Scaled     = lvlScale .gt. 0

      modeC      = modefg
      modeF      = modefg

      ! Determine the user-function call-status.

      call s8callStatus( status, iw, leniw )

      if (status .eq. 1) then
         !--------------------------------------------------------------
         ! First evaluation of the problem functions in snOptB
         ! On entry, lvlScale = 0.
         !--------------------------------------------------------------
         iw(gotFD)  =  0 ! > 0 => some differences needed
         iw(gotG)   =  0 ! > 0 => some exact derivatives provided
         iw(conJij) =  0 ! > 0 => constant Jacobian elements provided
         call snPRNT( 13, ' ', iw, leniw )
         call dload ( negCon, gDummy, gCon, 1 )
         call dload ( nnObj , gDummy, gObj, 1 )
      end if

      !-----------------------------------------------------------------
      ! Unscale x (never required for status = 1)
      !-----------------------------------------------------------------
      if (Scaled) then
         call dcopy ( nnH,           x, 1, rw(lxScaled), 1 )
         call ddscl ( nnH, rw(lscales), 1,            x, 1 )

         ! If the Jacobian has some constant elements, they are wrecked
         ! by the scaling.  Restore them from gConU.

         if (NonlinearCon) then
            if (modefg .gt. 0  .and.  iw(conJij) .gt. 0) then
               call dcopy ( negCon, rw(lgConU), 1, gCon, 1 )
            end if
         end if
      end if

      !=================================================================
      ! Compute the constraint functions.
      !=================================================================
      if (NonlinearCon) then
         statusUser = status    ! In case fgCon alters status
         if (lvlTim .ge. 2) call s1time( 4, 0, iw, leniw, rw, lenrw )
         call fgcon
     &      ( modeC, nnCon, nnJac, negCon,
     &        x, fCon, gCon, statusUser,
     &        cu, lencu, iu, leniu, ru, lenru )
         if (lvlTim .ge. 2) call s1time(-4, 0, iw, leniw, rw, lenrw )
         iw(nfCon1) = iw(nfCon1) + 1
         if (modefg .gt. 0)
     &        iw(nfCon2) = iw(nfCon2) + 1
      end if

      !=================================================================
      ! Compute the objective function.
      !=================================================================
      if (NonlinearObj  .and.  modeC .ge. 0) then
         if (FPonly) then
            call dcopy ( nnObj, x, 1, gObj, 1 )
            call daxpy ( nnObj, (-one), rw(lx0), 1, gObj, 1 )
            fObj = half*proxWeight*ddot ( nnObj, gObj, 1, gObj, 1 )
            call dscal ( nnObj, proxWeight, gObj, 1 )
         else ! nnObj = nnObjU
            statusUser = status ! In case fgobj alters status
            if (lvlTim .ge. 2) call s1time( 5, 0, iw, leniw, rw, lenrw )
            call fgobj
     &         ( modeF, nnObjU, x, fObj, gObj, statusUser,
     &           cu, lencu, iu, leniu, ru, lenru )
            if (lvlTim .ge. 2) call s1time(-5, 0, iw, leniw, rw, lenrw )
            iw(nfObj1) = iw(nfObj1) + 1
            if (modefg .gt. 0)
     &           iw(nfObj2) = iw(nfObj2) + 1
         end if
      end if

      !-----------------------------------------------------------------
      ! Scale  x and the derivatives.
      !-----------------------------------------------------------------
      if (Scaled) then
         call dcopy ( nnH, rw(lxScaled), 1, x, 1 )

         if (NonlinearCon) then
            call dddiv ( nnCon, rw(lscales+n), 1, fCon, 1 )
            if (modefg .gt. 0  .and.  iw(gotG) .gt. 0) then
               call s8scaleJ
     &            ( nnCon, nnJac, negCon, n, rw(lscales),
     &              neJ, nlocJ, locJ, indJ, gCon, rw, lenrw )
            end if
         end if

         if (NonlinearObj  .and.  modeC .ge. 0) then
            if (modefg .gt. 0  .and.  iw(gotG) .gt. 0) then
               call s8scaleg
     &            ( nnObj, rw(lscales), gObj, rw, lenrw )
            end if
         end if
      end if

      if (modeC .lt. 0  .or.  modeF .lt. 0) then
         !--------------------------------------------------------------
         ! The user may be saying the function is undefined (mode = -1)
         ! or may just want to stop                         (mode < -1).
         !--------------------------------------------------------------
         if (modeC .eq. -1  .or.  modeF .eq. -1) then
            iExit = -1
         else
            if (modeC .lt. 0) then
               iExit = 72
            else
               iExit = 73
            end if
         end if
      end if

      !=================================================================
      ! Do some housekeeping on the first snOptB entry.
      !=================================================================
      if (status .eq. 1  .and.  iExit .eq. 0) then
         if (NonlinearCon) then
            !-----------------------------------------------------------
            ! Count how many Jacobian elements are provided.
            !-----------------------------------------------------------
            nnGlin   = 0
            supplied = 0
            do l = 1, negCon
               if (gCon(l) .ne. gDummy) supplied = supplied + 1
            end do

            write(str, 1100) supplied, negCon
            call snPRNT( 3, str, iw, leniw )

            if (supplied .lt. negCon) then

               ! Some Jacobian elements are missing.

               if (iw(lvlDer) .ge. 2) then
                  !-----------------------------------------------------
                  ! All the Jacobian is known.  Any undefined elements
                  ! are assumed constant, and are restored from gConU.
                  !-----------------------------------------------------
                  call snPRNT( 3,
     &              ' ==>  Some constraint derivatives are missing, '
     &            //' assumed constant.', iw, leniw )
                  call snPRNT( 3, ' ', iw, leniw )

                  lG  = lgConU
                  do l  = 1, negCon
                     if (gCon(l) .eq. gDummy) then
                        gCon(l) = rw(lG)
                        nnGlin  = nnGlin + 1
                     end if
                     lG = lG + 1
                  end do
               else
                  !-----------------------------------------------------
                  ! Save a permanent copy of gCon in gConU so that we know
                  ! which derivatives must be estimated.
                  !-----------------------------------------------------
                  call dcopy ( negCon, gCon, 1, rw(lgConU), 1 )
               end if
            end if ! supplied < negCon
            if (supplied + nnGlin .lt. negCon) iw(gotFD ) = 1
            if (supplied          .gt.      0) iw(gotG  ) = 1
            if (nnGlin            .gt.      0) iw(conJij) = 1
         end if

         if (NonlinearObj) then
            !-----------------------------------------------------------
            ! Count how many working gradient elements are known.
            ! (These may be the gradients of the FP objective.)
            !-----------------------------------------------------------
            if (FPonly) then
               supplied = nnObj
               iw(gotG) = 1
               write(str, 2010) nnObj
               call snPRNT( 3, str, iw, leniw )

            else
               supplied = 0
               do l = 1, nnObjU
                  if (gObj(l) .ne. gDummy) then
                     supplied = supplied + 1
                  end if
               end do

               write(str, 2000) supplied, nnObjU
               call snPRNT( 3, str, iw, leniw )

               if (supplied .gt. 0) then
                  iw(gotG) = 1
               end if

               if (supplied .lt. nnObjU) then

                  ! Some objective gradients are missing.

                  iw(gotFD) = 1

                  if (iw(lvlDer) .eq. 1  .or.  iw(lvlDer) .eq. 3) then
                     !--------------------------------------------------
                     ! The objective gradient was meant to be known.
                     !--------------------------------------------------
                     iw(lvlDer) = iw(lvlDer) - 1
                     write(str, 2100) iw(lvlDer)
                     call snPRNT( 3, str, iw, leniw )
                  end if
               end if
            end if

            !--------------------------------------------------------
            ! Copy gObj into gObjU.
            !--------------------------------------------------------
            call dcopy ( nnObj, gObj, 1, rw(lgObjU), 1 )

         end if
      end if

      return

 1100 format(' The user has defined', i8, '   out of', i8,
     &       '   constraint gradients.')
 2000 format(' The user has defined', i8, '   out of', i8,
     &       '   objective  gradients.')
 2010 format(' SnOptB  will define ', i8, '   gradients for the ',
     &       ' FP objective.')
 2100 format(' XXX  Some objective  derivatives are missing ---',
     &       ' derivative level reduced to', i3)

      end ! subroutine s0fgB

