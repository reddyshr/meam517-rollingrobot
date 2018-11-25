*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     File  train.f
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

      implicit           none
      integer            Status, mode, neF, n
      integer            lencu, leniu, lenru, lencw, leniw, lenrw
      double precision   F(neF), x(n)
      character*8        cu(lencu), cw(lencw)
      integer            iu(leniu), iw(leniw)
      double precision   ru(lenru), rw(lenrw)

      integer            eNN
      parameter          ( eNN = 2002 )
      call trainf( x(1), eNN,
     &     x(    eNN + 1),    x(2 * eNN + 1),
     &     x(3 * eNN + 1),    x(4 * eNN + 1),
     &     F(1),  F(2), F( eNN + 2 ) )

      end ! subroutine usrfun

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine trainf
     &   ( x, N, v, a, ua, ub, Obj, F, C )

      implicit           none
      integer            N
      double precision   x(N), v(N), a(N), ua(N), ub(N)
      double precision   Obj, F(N), C( 2 * (N - 1) )

      double precision   sum
      integer            i, j

      double precision   h, time
      parameter          ( time = 4.8d0 )
      integer            ns
      parameter          ( ns = 3 )
      double precision   s(ns), z(ns - 1)

      double precision   aa, bb, cc, eps, pi
      parameter          ( aa  = 0.3d0,  bb = 0.14d0, cc = 0.16d0 )
      parameter          ( eps = 0.05d0, pi = 3.14159d0 )

      data            s / 2.0, 0.0, -2.0 /
      data            z / 2.0, 4.0 /

      h = time / ( 1.0d0 * ( N - 1 ) )

      Obj = 0.0d0
      do i = 1, N
         Obj = Obj + ua(i) * v(i) * h
      end do

      do i = 1, N
         sum = 0.0d0
         do j = 1, ns - 1
            sum = sum + (s(j+1) - s(j)) * atan( (x(i) - z(j)) /eps )/ pi
         end do

         F(i) = a(i) + sum + aa + bb * v(i) + cc * v(i) * v(i) -
     &        ua(i) + ub(i)
      end do

      do i = 1, N - 1
         C(2 * i - 1) = x(i + 1) - x(i) -
     &        h * ( v(i + 1) + v(i) ) / 2
         C(2 * i)     = v(i + 1) - v(i) -
     &        h * ( a(i + 1) + a(i) ) / 2
      end do

      end  ! subroutine trainf

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine usrini
     &   ( ObjAdd, ObjRow, Prob,
     &     x, xlow, xupp, xstate, Names,
     &     Fmul, Flow, Fupp, Fstate, FNames,
     &     iSpecs, iPrint, iSumm, iErr,
     &     cu,  iu,  ru,
     &     cw,  iw,  rw )

      implicit           none
      double precision   ObjAdd
      integer            ObjRow

      integer            n, neF, nName, nFnames
      integer            lencw, leniw, lenrw, lencu, leniu, lenru

      parameter        ( lencw = 501, leniw =20000000, lenrw =50000000 )
      parameter        ( lencu =   1, leniu =       1, lenru =       1 )
      parameter        ( n     = 10010,   neF  = 6005 )
      parameter        ( nName = 1, nFnames = 1 )

      character*8        Prob, Names(nName), FNames(nFnames)
      double precision   x(n), xlow(n), xupp(n)
      double precision   Flow(neF), Fupp(neF), Fmul(neF)

      integer            iSpecs, iPrint, iSumm, iErr
      integer            xstate(n),  Fstate(neF)

      integer            iu(leniu), iw(leniw)
      character*8        cu(lencu), cw(lencw)
      double precision   ru(lenru), rw(lenrw)

      double precision   zero, pt5, one, two, three, plInfy
      parameter          ( zero  = 0.0d+0, pt5    = 0.5d+0  )
      parameter          ( one   = 1.0d+0, two    = 2.0d+0  )
      parameter          ( three = 3.0d+0, plInfy = 1.0d+20 )

      integer            i
      logical            byname
      character*20       lfile

      integer            eNN
      parameter          ( eNN = 2002 )
*     n = eNN * 5
*     neF = 1 + eNN + 2 * (eNN - 1)
*     ------------------------------------------------------------------
*     Give the problem a name.

      Prob   = 'train   '

      ObjAdd = zero

      iSpecs = 4
      iPrint = 15

      lfile = 'train.spc'
      open( iSpecs, file=lfile, status='OLD',     err=800 )

      lfile = 'train.out'
      open( iPrint, file=lfile, status='UNKNOWN', err=800 )

      call trainini( x(1), eNN,
     &     x(    eNN + 1),    x(2 * eNN + 1),
     &     x(3 * eNN + 1),    x(4 * eNN + 1),
     &     xlow(1),
     &     xlow(    eNN + 1), xlow(2 * eNN + 1),
     &     xlow(3 * eNN + 1), xlow( 4 * eNN + 1 ),
     &     xupp(1),
     &     xupp(    eNN + 1), xupp(2 * eNN + 1),
     &     xupp(3 * eNN + 1), xupp(4 * eNN + 1),
     &     Flow, neF, Fupp )

      return

 800  write(iErr, 4000) 'Error while opening file', lfile

 4000 format(/  a, 2x, a  )

      end ! subroutine usrini

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine trainini
     &   ( x, N, v, a, ua, ub,
     &     xlow, vlow, alow, ualow, ublow,
     &     xupp, vupp, aupp, uaupp, ubupp,
     &     Flow, neF, Fupp )

      implicit none

      integer               N
      double precision      x(N), v(N), a(N), ua(N), ub(N)
      double precision      xlow(N), vlow(N), alow(N)
      double precision      ualow(N), ublow(N)
      double precision      xupp(N), vupp(N), aupp(N)
      double precision      uaupp(N), ubupp(N)

      integer               neF
      double precision      Flow(neF), Fupp(neF)

      integer               i

      do i = 1, N
         xlow(i)  = 0.0d0
         vlow(i)  = 0.0d0
         ualow(i) = 0.0d0
         ublow(i) = 0.0d0

         uaupp(i) = 10.0d0
         ubupp(i) = 2.0d0
      end do

*     boundary values
      x(1)    = 0.0d0
      xlow(1) = 0.0d0
      xupp(1) = 0.0d0

      x(N)    = 6.0d0
      xlow(N) = 6.0d0
      xupp(N) = 6.0d0


      v(1)    = 0.0d0
      vlow(1) = 0.0d0
      vupp(1) = 0.0d0

      v(N)    = 0.0d0
      vlow(N) = 0.0d0
      vupp(N) = 0.0d0

*     there are only equality constraints
      do i = 1, neF
         Flow(i) = 0.0d0
         Fupp(i) = 0.0d0
      end do

      end ! subroutine trainini

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine dumpout
     &   ( x, N, v, a, ua, ub )

      implicit
     &     none
      integer
     &     n
      double precision
     &     x(N), v(N), a(N), ua(N), ub(N)
      integer
     &     i

      do i = 1, N
         print *, x(i), v(i), a(i), ua(i), ub(i)
      end do

      end ! subroutine dumpout
