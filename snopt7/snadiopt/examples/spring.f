*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     File  spring.f   -- an optimal control problem
*
*     usrfun   usrini
*
*     User-defined routines for use with  SnadiOpt.
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

*     ==================================================================
*     The problem size depends on a parameter T.  There are
*     2T constraints and 3T + 2 variables, as well as bounds
*     on the variables.  The first T constraints are quadratic in
*     T + 1 variables, and the objective function is quadratic in
*     T + 1 other variables.
*
*     The control problem models a spring, mass and damper system.
*     It is of the form
*
*   --------------------------------------------------------------------
*   | minimize    1/2 sum x(t)**2   (t = 0 to T)                       |
*   |                                                                  |
*   | subject to                                                       |
*   |     x(t+1)  =  x(t)  +  0.2  y(t),                               |
*   |                                                                  |
*   |     y(t+1)  =  y(t)  -  0.01 y(t)**2  -  0.004 x(t)  +  0.2 u(t) |
*   |                                                                  |
*   |     y(t)   >=  -1,     -0.2  <=  u(t)  <=  0.2,                  |
*   |                                                                  |
*   |                (all for t = 0 to T-1)                            |
*   | and                                                              |
*   |     y(0)    =   0,      y(T)  =  0,       x(0) = 10.             |
*   --------------------------------------------------------------------
*
*     For large enough T (e.g. T >= 90), the optimal objective value
*     is about 1186.382.
*
*     This model with T = 100 was used as test problem 5.11 in
*     B. A. Murtagh and M. A. Saunders (1982), A projected Lagrangian
*     algorithm and its implementation for sparse nonlinear constraints,
*     Mathematical Programming Study 16, 84--117.
*
*     14 Nov 1994: First version of spring.f, derived from manne.f.
*     02 Jun 2001: Updated for SNOPT 6.1
*     ==================================================================
      double precision
     &     FObj, ut, xt, xtp1, yt, ytp1
      integer
     &     T, jt, jx0, jx, jy0, jy, ju0, ju, lin, nln, ObjRow
*     ------------------------------------------------------------------
      T      = (n-2)/3

      ObjRow = 1
      lin    = ObjRow + 1       ! points to linear    constraints in F.
      nln    = lin    + T       ! points to nonlinear constraints in F.

      jy0  = 1                  ! points to state   y(0)  in x.
      jx0  = jy0 + T + 1        ! points to state   x(0)  in x.
      ju0  = jx0 + T + 1        ! points to control u(0)  in x

      FObj = 0.0d0

      do jt = 0, T-1
         jx = jx0 + jt
         jy = jy0 + jt
         ju = ju0 + jt

         xt   = x(jx)
         xtp1 = x(jx+1)

         yt   = x(jy)
         ytp1 = x(jy+1)

         ut   = x(ju)

         F(nln) = 1.0d-2*yt**2 - yt + ytp1 + 4.0d-3*xt       - 0.2d0*ut
         F(lin) = -0.2d0*yt                -        xt + xtp1
         FObj   = FObj +  xt**2

         nln    = nln + 1
         lin    = lin + 1
      end do

*     Set the objective row.

      F(ObjRow) = (FObj + xtp1**2)/2.0d+0

      end ! of usrfun

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
      integer            n, neF, nName, nFnames, ObjRow
      integer            lencw, leniw, lenrw, lencu, leniu, lenru

      parameter        ( lencw = 501, leniw = 1000000, lenrw = 200000 )
      parameter        ( lencu =   1, leniu =       1, lenru =      1 )
*                        n     =3*T+2   neF     =2*T+1
      parameter        ( n     =  302,  neF     =  201 )
      parameter        ( nName =    1,  nFnames =    1 )

      character*8        Prob, Names(nName)
      double precision   x(n), xlow(n), xupp(n)
      double precision   Flow(neF), Fupp(neF), Fmul(neF)
      character*8        FNames(nFnames)

      integer            iSpecs, iPrint, iSumm, iErr
      integer            xstate(n),  Fstate(neF)

      integer            iu(leniu), iw(leniw)
      character*8        cu(lencu), cw(lencw)
      double precision   ru(lenru), rw(lenrw)

*   --------------------------------------------------------------------
*   | minimize    1/2 sum x(t)**2   (t = 0 to T)                       |
*   |                                                                  |
*   | subject to                                                       |
*   |     x(t+1)  =  x(t)  +  0.2  y(t),                               |
*   |                                                                  |
*   |     y(t+1)  =  y(t)  -  0.01 y(t)**2  -  0.004 x(t)  +  0.2 u(t) |
*   |                                                                  |
*   |     y(t)   >=  -1,     -0.2  <=  u(t)  <=  0.2,                  |
*   |                                                                  |
*   |                (all for t = 0 to T-1)                            |
*   | and                                                              |
*   |     y(0)    =   0,      y(T)  =  0,       x(0) = 10.             |
*   --------------------------------------------------------------------
      integer
     &     T, i, jt, jx0, jx, jy0, jy, ju0, ju, k, m, nCon
      logical
     &     byname
      character*20
     &     lfile
*   --------------------------------------------------------------------
      double precision   zero,          one,          plInfy
      parameter         (zero = 0.0d+0, one = 1.0d+0, plInfy = 1.0d+20)
*   --------------------------------------------------------------------

*     Given n, compute the interval length T.

      T    = (n - 2)/3

*     Write T into the problem name.

      write(prob, '(i8)') T
      if      (T .lt.  100) then
         prob(1:6) = 'Spring'
      else if (T .lt. 1000) then
         prob(1:5) = 'Sprin'
      else
         prob(1:3) = 'Spr'
      end if

      write(iSumm, *) 'Problem SPRING.    T =', T

      iSpecs =  4
      iPrint = 15
      iSumm  =  6

      lfile = 'spring.spc'
      open( iSpecs, file=lfile, status='OLD',     err=800 )

      lfile = 'spring.out'
      open( iPrint, file=lfile, status='UNKNOWN', err=800 )

      ObjAdd = zero
      nCon   = 2*T

*     ObjRow = 2*T+1
      ObjRow = 1

*     The variables are ordered as follows:
*     x Variables    1: T+1  states
*     y            T+2:2T+2  states
*     u           2T+2:3T+2  controls
*     jx, jy, ju are the base indices for the corresponding variables

      jy0  = 1                  ! points to state   y(0)  in x.
      jx0  = jy0 + T + 1        ! points to state   x(0)  in x.
      ju0  = jx0 + T + 1        ! points to control u(0)  in x

      do jt = 0, T
         jx = jx0 + jt
         jy = jy0 + jt

*        Initialize the bounds and values of the states x.

         xlow(jx)   = -plInfy
         xupp(jx)   =  plInfy
         x(jx)      =  zero
         xstate(jx) =  3

*        Initialize the bounds and values of the states y.

         xlow(jy)   = -one
         xupp(jy)   =  plInfy
         x(jy)      = -one
         xstate(jy) =  0
      end do

      do jt = 0, T-1
         ju = ju0 + jt

*        Initialize the bounds and values of the constrols u.

         xlow(ju)   = -0.2d0
         xupp(ju)   =  0.2d0
         x(ju)      =  zero
         xstate(ju) =  3
      end do

*     Fix the boundary conditions.

      xlow(jx0)   = 10.0d0
      xupp(jx0)   = 10.0d0
      x(jx0)      = 10.0d0

      xlow(jy0)   = zero
      xupp(jy0)   = zero

      xlow(jy0+T) = zero
      xupp(jy0+T) = zero
      x(jy0+T)    = zero

*     Bounds on F
*     Set the objective first

      Fmul(ObjRow) =  zero
      Flow(ObjRow) = -plInfy
      Fupp(ObjRow) =  plInfy

      k = ObjRow + 1
      do i = 1, nCon
         Flow(k) = zero
         Fupp(k) = zero
         Fmul(k) = zero
         k = k + 1
      end do

      return

 800  write(iErr, 4000) 'Error while opening file', lfile

 4000 format(/  a, 2x, a  )

      end ! subroutine usrini
