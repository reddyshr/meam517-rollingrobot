*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     File  banana.f  -- an unconstrained test problem.
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
*     Usrfun computes the objective function for the unconstrained
*     problem banana (aka Rosenbrock's function).
*
*     minimize  100*(x_2 - x_1**2)**2  + (1 - x_1)**2
*
*     subject to    -10 <= x_1, x_2 <= 10
*
*     ==================================================================

      F(1) = 100.0d+0*(x(2) - x(1)**2)**2  + (1.0d+0 - x(1))**2

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
      parameter        ( n     =   2, neF   =     1 )
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
*     usrini defines input data for the unconstrained problem banana.
*
*     ==================================================================
      integer
     &     i, inform, iPrt, iSum, j
      logical
     &     byname
      character*20
     &     lfile
*     ------------------------------------------------------------------
*     Give the problem a name.

      Prob   = 'Banana '

      ObjRow = 1

*     ----------------
*     Initial x.
*     ----------------
      x(1) = - 1.2d+0
      x(2) =   1.0d+0

*     ------------------------------------------------------------------
*     Ranges on the variables
*     ------------------------------------------------------------------
      do j = 1, n
         xlow(j) = -10.0d+0
         xupp(j) =  10.0d+0
      end do

      iSpecs = 4
      iPrint = 15

      lfile = 'banana.spc'
      open( iSpecs, file=lfile, status='OLD',     err=800 )

      lfile = 'banana.out'
      open( iPrint, file=lfile, status='UNKNOWN', err=800 )

      return

 800  write(iErr, 4000) 'Error while opening file', lfile

 4000 format(/  a, 2x, a  )

      end ! subroutine usrini

