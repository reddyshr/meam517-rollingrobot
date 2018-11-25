*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     File  snmain.f
*
*     usrfun   usrini
*
*     User-defined routines for running SnadiOpt.
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine usrfun
     &   ( Status, mode,
     &     neF, n, x, F,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     Status, mode, neF, n, lencu, leniu, lenru, lencw, leniw,
     &     lenrw, iu(leniu), iw(leniw)
      double precision
     &     F(neF), x(n), ru(lenru), rw(lenrw)
      character*8
     &     cu(lencu), cw(lencw)

*     ==================================================================
*     Usrfun computes the objective and constraint functions for the
*     Hexagon problem featured in the Snopt users guide.
*
*     minimize  -x_2 x_6 + x_1 x_7 - x_3 x_7 - x_5 x_8 + x_4 x_9
*                                            + x_3 x_8
*     subject to
*          0 <=   -x_1 + x_2                            <= 1
*          0 <=         -x_2 + x_3                      <= 1
*          0 <=                x_3 - x_4                <= 1
*          0 <=                      x_4 - x_5          <= 1
*                  x_1^2 + x_6^2                        <= 1
*                 (x_2   - x_1)^2  +  (x_7 - x_6)^2     <= 1
*                 (x_3   - x_1)^2  +   x_6^2            <= 1
*                 (x_1   - x_4)^2  +  (x_6 - x_8)^2     <= 1
*                 (x_1   - x_5)^2  +  (x_6 - x_9)^2     <= 1
*                  x_2^2 + x_7^2                        <= 1
*                 (x_3   - x_2)^2  +   x_7^2            <= 1
*                 (x_4   - x_2)^2  +  (x_8 - x_7)^2     <= 1
*                 (x_2   - x_5)^2  +  (x_7 - x_9)^2     <= 1
*                 (x_4   - x_3)^2  +   x_8^2            <= 1
*                 (x_5   - x_3)^2  +   x_9^2            <= 1
*                  x_4^2 +  x_8^2                       <= 1
*                 (x_4   - x_5)^2 + (x_9 - x_8)^2       <= 1
*                  x_5^2 + x_9^2                        <= 1
*
*     ==================================================================
      integer
     &     Obj
*     ------------------------------------------------------------------
      Obj    = 19                  ! The objective row

*     Constraint functions

      F( 1)  =   -x(1)    + x(2)
      F( 2)  =             -x(2) + x(3)
      F( 3)  =                     x(3) - x(4)
      F( 4)  =                            x(4) - x(5)
      F( 5)  =    x(1)**2 + x(6)**2
      F( 6)  =   (x(2)   - x(1))**2  +  (x(7) - x(6))**2
      F( 7)  =   (x(3)   - x(1))**2  +   x(6)**2
      F( 8)  =   (x(1)   - x(4))**2  +  (x(6) - x(8))**2
      F( 9)  =   (x(1)   - x(5))**2  +  (x(6) - x(9))**2
      F(10)  =    x(2)**2 + x(7)**2
      F(11)  =   (x(3)   - x(2))**2  +   x(7)**2
      F(12)  =   (x(4)   - x(2))**2  +  (x(8) - x(7))**2
      F(13)  =   (x(2)   - x(5))**2  +  (x(7) - x(9))**2
      F(14)  =   (x(4)   - x(3))**2  +   x(8)**2
      F(15)  =   (x(5)   - x(3))**2  +   x(9)**2
      F(16)  =    x(4)**2 +  x(8)**2
      F(17)  =   (x(4)   - x(5))**2 + (x(9) - x(8))**2
      F(18)  =    x(5)**2 + x(9)**2

*     Objective (minimized).

      F(Obj) =  -x(2)*x(6) + x(1)*x(7) - x(3)*x(7) - x(5)*x(8)
     &                                 + x(4)*x(9) + x(3)*x(8)

      end ! subroutine usrfun

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine usrini
     &   ( ObjAdd, ObjRow, Prob,
     &     x, xlow, xupp, xstate, Names,
     &     Fmul, Flow, Fupp, Fstate, FNames,
     &     iSpecs, iPrint, iSumm, iErr,
     &     cu,  iu,  ru,
     &     cw,  iw,  rw )

      implicit
     &     none
      integer
     &     n, neF, nName, nFnames, ObjRow,
     &     lencw, leniw, lenrw, lencu, leniu, lenru
      parameter        ( lencw = 501, leniw = 10000, lenrw = 20000 )
      parameter        ( lencu =   1, leniu =     1, lenru =     1 )
      parameter        ( n     =   9, neF   =    19 )
      parameter        ( nName =   1, nFnames =   1 )
      integer
     &     iSpecs, iPrint, iSumm, iErr, xstate(n), Fstate(neF),
     &     iu(leniu), iw(leniw)
      double precision
     &     ObjAdd, x(n), xlow(n), xupp(n), Flow(neF), Fupp(neF),
     &     Fmul(neF), ru(lenru), rw(lenrw)
      character*8
     &     Prob, Names(nName), FNames(nFnames), cu(lencu), cw(lencw)

*     ==================================================================
*     usrini defines input data for the hexagon problem discussed in the
*     snopt  Users Guide.
*
*     ==================================================================
      integer
     &     i, inform, iPrt, iSum, j
      logical
     &     byname
      character*20
     &     lfile
*     ------------------------------------------------------------------
      double precision   InfBnd
      parameter        ( InfBnd = 1.0d+20 )
*     ------------------------------------------------------------------
*     Give the problem a name.

      Prob   = 'Hexagon '

      ObjAdd = 0.0d+0
      ObjRow = 19

*     iPrt   = 0
*     iSum   = 0
*     call snset ( 'Maximize   ', iPrt, iSum, inform,
*    &     cw, lencw, iw, leniw, rw, lenrw )

*     ----------------
*     Initial x.
*     ----------------
      x(1)   =  .1d+0
      x(2)   =  .125d+0
      x(3)   =  .666666d+0
      x(4)   =  .142857d+0
      x(5)   =  .111111d+0
      x(6)   =  .2d+0
      x(7)   =  .25d+0
      x(8)   = -.2d+0
      x(9)   = -.25d+0

*     ------------------------------------------------------------------
*     Ranges for the constraint functions.
*     ------------------------------------------------------------------
*     Linear constraints.

      do i = 1, 4
         Flow(i) =  0.0d+0
         Fupp(i) =  InfBnd
      end do

*     Nonlinear constraints.

      do i = 5, 18
         Flow(i) = -InfBnd
         Fupp(i) =  1.0d+0
      end do

*     ------------------------------------------------------------------
*     Ranges on the variables
*     ------------------------------------------------------------------
      do j = 1, n
         xlow(j)   = -InfBnd
         xupp(j)   =  InfBnd
         xstate(j) = 0
      end do

      xlow(1) =  0.0d+0
      xlow(3) = -1.0d+0
      xlow(5) =  0.0d+0
      xlow(6) =  0.0d+0
      xlow(7) =  0.0d+0

      xupp(3) =  1.0d+0
      xupp(8) =  0.0d+0
      xupp(9) =  0.0d+0

      do i = 1, neF
         Fmul(i) = 0.0d+0
      end do

      iSpecs = 4
      iPrint = 15

      lfile = 'snmain.spc'
      open( iSpecs, file=lfile, status='OLD',     err=800 )

      lfile = 'snmain.out'
      open( iPrint, file=lfile, status='UNKNOWN', err=800 )

      return

 800  write(iErr, 4000) 'Error while opening file', lfile

 4000 format(/  a, 2x, a  )

      end ! subroutine usrini

