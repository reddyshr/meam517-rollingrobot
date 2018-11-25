*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     File  hs47.f
*
*     usrfun   usrini
*
*     User-defined routines for running snadifor.
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
*     HS problem 47.
*     Minimize    (x(1)-x(2))**2 + (x(2)-x(3))**3 + (x(3)-x(4))**4
*                                                 + (x(4)-x(5))**4
*
*     subject to                        x(3)    + x(4) + x(5) >= 3
*                 x(1)      + x(2)**2 + x(3)**3                = 3
*                 x(1)      + x(2)                            >= 1
*                             x(2)    - x(3)**2 + x(4)         = 1
*                 x(1)*x(5)                                    = 1
*     ==================================================================
      integer
     &     Obj
*     ------------------------------------------------------------------

      Obj  = 1

      F(Obj) = (x(1)-x(2))**2 + (x(2)-x(3))**3 +
     &         (x(3)-x(4))**4 + (x(4)-x(5))**4

      F(2)   =  x(3) + x(4) + x(5)
      F(3)   =  x(1) + x(2)*x(2) + x(3)*x(3)*x(3)
      F(4)   =  x(1) + x(2)
      F(5)   =  x(2) - x(3)*x(3) + x(4)
      F(6)   =  x(1)*x(5)

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
     &     n, neF, nName, nFnames, ObjRow, lencw, leniw, lenrw, lencu,
     &     leniu, lenru
      parameter        ( lencw = 501, leniw = 10000, lenrw = 20000 )
      parameter        ( lencu =   1, leniu =     1, lenru =     1 )
      parameter        ( n     =   5, neF   =     6 )
      parameter        ( nName =   1, nFnames =   1 )
      integer
     &     iSpecs, iPrint, iSumm, iErr, xstate(n), Fstate(neF),
     &     iu(leniu), iw(leniw)
      double precision
     &     ObjAdd, x(n), xlow(n), xupp(n), Flow(neF), Fupp(neF),
     &     Fmul(neF), ru(lenru), rw(lenrw)
      character*8
     &     Prob, Names(nName), FNames(nFnames), cu(lencu), cw(lencw)

*     ------------------------------------------------------------------
*
*     Minimize    (x(1)-x(2))**2 + (x(2)-x(3))**3 + (x(3)-x(4))**4
*                                                 + (x(4)-x(5))**4
*     subject to  x(1)      + x(2)**2 + x(3)**3            = 3
*                             x(2)    - x(3)**2 + x(4)     = 1
*                 x(1)*x(5)                                = 1
*
*     ------------------------------------------------------------------
      double precision   zero, pt5, one, two, three, plInfy
      parameter          ( zero  = 0.0d+0, pt5    = 0.5d+0  )
      parameter          ( one   = 1.0d+0, two    = 2.0d+0  )
      parameter          ( three = 3.0d+0, plInfy = 1.0d+20 )

      integer            i
      logical            byname
      character*20       lfile
*     ------------------------------------------------------------------
*     Give the problem a name.

      Prob   = 'HS47    '

      ObjAdd = zero

*     ----------------
*     Initial x.
*     ----------------
      x(1)   =  two
      x(2)   =  sqrt(two) - one
      x(3)   =  x(2)
      x(4)   =  two
      x(5)   =  pt5

      do i = 1, n
         xlow(i) = -plInfy
         xupp(i) =  plInfy
      end do

      Flow(2)  = three
      Fupp(2)  = plInfy
      Flow(3)  = three
      Fupp(3)  = three
      Flow(4)  = one
      Fupp(4)  = plInfy
      Flow(5)  = one
      Fupp(5)  = one
      Flow(6)  = one
      Fupp(6)  = one

      do i = 1, neF
         Fmul(i) = zero
      end do

      iSpecs = 4
      iPrint = 15

      lfile = 'hs47.spc'
      open( iSpecs, file=lfile, status='OLD',     err=800 )

      lfile = 'hs47.out'
      open( iPrint, file=lfile, status='UNKNOWN', err=800 )

      return

 800  write(iErr, 4000) 'Error while opening file', lfile

 4000 format(/  a, 2x, a  )

      end ! subroutine usrini

