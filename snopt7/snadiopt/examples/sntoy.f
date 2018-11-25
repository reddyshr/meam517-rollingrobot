*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     File  sntoy.f
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
*     Usrfun computes the objective and constraint functions for the toy
*     problem featured in the SnadiOpt users guide.
*
*     Minimize      3*x(1) + (x(1) + x(2) + x(3))**2 + 5*x(4)
*
*     subject to             4*x(2)    + 2*x(3)               >= 0
*                     x(1) +   x(2)**2 +   x(3)**2             = 2
*                              x(2)**4 +   x(3)**4   +   x(4)  = 4
*
*                     x(1) >= 0,                         x(4) >= 0.
*
*
*     ==================================================================
      integer
     &     Obj
*     ------------------------------------------------------------------
      Obj  = 1                  ! The objective row

      F(Obj) = 3.0d+0*x(1) + (x(1) + x(2) + x(3))**2 + 5.0d+0*x(4)

*     Nonlinear components of the constraints.

      F(2)   =              4.0d+0*x(2)    + 2.0d+0*x(3)
      F(3)   =       x(1) +        x(2)**2 +        x(3)**2
      F(4)   =                     x(2)**4 +        x(3)**4  + x(4)

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
      parameter        ( n     =   4, neF   =     4 )
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
*     datToy defines input data for the toy problem discussed in the
*     snopt Users Guide.
*
*     Minimize      3*x(1) + (x(1) + x(2) + x(3))**2 + 5*x(4)
*
*     subject to             4*x(2)    + 2*x(3)               >= 0
*                     x(1) +   x(2)**2 +   x(3)**2             = 2
*                              x(2)**4 +   x(3)**4   +   x(4)  = 4
*
*                     x(1) >= 0,                         x(4) >= 0.
*
*     ==================================================================
      integer
     &     i
      logical
     &     byname
      character*20
     &     lfile
*     ------------------------------------------------------------------
      double precision   plInfy
      parameter         ( plInfy = 1.0d+20 )
*     ------------------------------------------------------------------
*     Give the problem a name.

      Prob   = 'Toy prob'

      ObjAdd = 0.0d+0
      ObjRow = 1

*     ----------------
*     Initial x.
*     ----------------
      x(1)   =  1.0d+0
      x(2)   =  1.0d+0
      x(3)   =  1.0d+0
      x(4)   =  1.0d+0

      do i = 1, n
         xlow(i)   = -plInfy
         xupp(i)   =  plInfy
         xstate(i) =  0
      end do

      xlow(1) = 0.0d+0
      xlow(4) = 0.0d+0

*     Impose bounds on the constraint rows.

      Flow(2)   =  0.0d+0
      Fupp(2)   =  plInfy

      Flow(3)   =  2.0d+0       ! Equality constraint
      Fupp(3)   =  2.0d+0

      Flow(4)   =  4.0d+0       ! Equality constraint
      Fupp(4)   =  4.0d+0

      do i = 1, neF
         Fmul(i) = 0.0d+0
      end do

      iSpecs = 4
      iPrint = 15

      lfile = 'sntoy.spc'
      open( iSpecs, file=lfile, status='OLD',     err=800 )

      lfile = 'sntoy.out'
      open( iPrint, file=lfile, status='UNKNOWN', err=800 )

      return

 800  write(iErr, 4000) 'Error while opening file', lfile

 4000 format(/  a, 2x, a  )

      end ! subroutine usrini

