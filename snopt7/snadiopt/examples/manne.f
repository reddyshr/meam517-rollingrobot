*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     File  manne.f
*
*     usrfun   usrini
*
*     User-defined routines for running snadifor.
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine usrfun
     &   ( State, mode,
     &     neF, n, x, F,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     State, mode, neF, n,
     &     lencu, leniu, lenru, lencw, leniw, lenrw
      double precision
     &     F(neF), x(n)
      character*8
     &     cu(lencu), cw(lencw)
      integer
     &     iu(leniu), iw(leniw)
      double precision
     &     ru(lenru), rw(lenrw)

*     ==================================================================
*     This is usrfun for problem manne.
*
*     The data bt(*) is computed by usrfun  on its first entry.
*
*     ==================================================================
      integer
     &     j, jc, ji, m, ObjRow, T
      double precision
     &     FObj, a, grow, beta, gfac, xc0, xCon, xi0, xk0
*     ------------------------------------------------------------------
      double precision   growth
      parameter         (growth = .03d+0)
      double precision   zero,          one
      parameter         (zero = 0.0d+0, one  = 1.0d+0)

      double precision   b, at(365), bt(365)
      common    /manne / b, at     , bt
*     ------------------------------------------------------------------

      T      = n/3
      m      = 2*T + 1
      ObjRow = m

*     ---------------------------------------
*     First entry.  Define b, at(*) and bt(*)
*     for this and all subsequent entries.
*     ---------------------------------------
      if (State .eq. 1) then
         grow   = 0.03d+0
         beta   = 0.95d+0
         xk0    = 3.0d+0
         xc0    = 0.95d+0
         xi0    = 0.05d+0
         b      = 0.25d+0

         a      = (xc0 + xi0) / xk0**b
         gfac   = (one + grow)**(one - b)
         at(1)  = a*gfac
         bt(1)  = beta

         do j  = 2, T
            at(j) = at(j-1)*gfac
            bt(j) = bt(j-1)*beta
         end do

         bt(T) = bt(T) / (one - beta)
      end if

*     -------------
*     Normal entry.
*     -------------
*     jC and jI are base indices for the Ct and It variables.

      jC    =   T
      jI    = 2*T

*     The first T components are linear, and the next T are nonlinear.

      do j = 1, T-1
         F(j) = - x(j) + x(j+1)       - x(jI+j)
      end do

*     The next linear component is special.

      F(T) = growth*x(T)   - x(jI+T)

      do j = 1, T
         F(T+j)   = at(j)*x(j)**b - x(jC+j) - x(jI+j)
      end do

*     Set the objective element.

      FObj = zero

      do j = 1, T
         xCon = x(T+j)
         FObj = FObj  +  bt(j) * log(xCon)
*Min     FObj = FObj  -  bt(j) * log(xCon)
      end do

      F(ObjRow) = FObj

      end ! subroutine usrfun

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine usrini
     &   ( ObjAdd, ObjRow, Prob,
     &     x, xlow, xupp, xstate, Names,
     &     Fmul, Flow, Fupp, Fstate, FNames,
     &     iSpecs, iPrint, iSumm, iErr,
     &     cu,  iu,  ru,
     &     cw,  iw,  rw )

*     ==================================================================
*     This is problem Manne.
*
*     ==================================================================
      implicit           none
      double precision   ObjAdd
      integer            n, neF, nName, nFnames, ObjRow
      integer            lencw, leniw, lenrw, lencu, leniu, lenru

      parameter        ( lencw = 501, leniw = 10000, lenrw = 20000 )
      parameter        ( lencu =   1, leniu =     1, lenru =     1 )
*                        n     =3*T,  neF   = 2*T+1  (max T = 365)
      parameter        ( n     = 90,  neF   =    61 )
      parameter        ( nName = 1,   nFnames =   1 )

      character*8        Prob, Names(nName)
      double precision   x(n), xlow(n), xupp(n)
      double precision   Flow(neF), Fupp(neF), Fmul(neF)
      character*8        FNames(nFnames)

      integer            iSpecs, iPrint, iSumm, iErr
      integer            xstate(n),  Fstate(neF)

      integer            iu(leniu), iw(leniw)
      character*8        cu(lencu), cw(lencw)
      double precision   ru(lenru), rw(lenrw)
*     ------------------------------------------------------------------
      integer            i, inform, iP, iS, jc, ji, jm, jy, k, m, T
      double precision   scale
      logical            byname
      character*20       lfile

      double precision   bplus,             bminus
      parameter         (bplus  = 1.0d+20,  bminus = - bplus)
      double precision   zero,              one
      parameter         (zero  = 0.0d+0,    one    =   1.0d+0)
*     ------------------------------------------------------------------

      T      = n/3
      m      = 2*T + 1          ! = neF

*     Write T into the problem name.

      write(prob, '(i8)') T
      if      (T .lt.  1000) then
         prob(1:5) = 'Manne'
      else if (T .lt. 10000) then
         prob(1:4) = 'Mann'
      else
         prob(1:3) = 'Man'
      end if

      ObjAdd = zero

      iP  = 0
      iS  = 0
      call snset ( 'Maximize', iP, iS, inform,
     &     cw, lencw, iw, leniw, rw, lenrw )

*     The variables are ordered as follows:
*     Kt Variables    1: T    represent Kapital     (nonlinear)
*     Ct            T+1:2T    represent Consumption
*     It           2T+1:3T    represent Investment
*     jC and jI are base indices for the Ct and It variables.

      jC     =   T
      jI     = 2*T

*     Set lower and upper bounds for Kt, Ct, It.
*     Also initial values and initial states for all variables.
*     The nonlinear variables are the most important.
*     We make them all superbasic.
*     The rest are ok nonbasic.
*     For test purposes, we want the initial x to be infeasible
*     with respect to the linear constraints.
*     Try setting the last Kapital too high.

      do k = 1, T
         xlow(   k) = 3.05d+0
         xupp(   k) = bplus
         xlow(jC+k) = 0.95d+0
         xupp(jC+k) = bplus
         xlow(jI+k) = 0.05d+0
         xupp(jI+k) = bplus

         x(   k)  = 3.0d+0 + (k - 1)/10.0d+0
         x(jC+k)  = xlow(jC+k)
         x(jI+k)  = xlow(jI+k)

*-->     xstate(   k) = 2
         xstate(   k) = 0
         xstate(jC+k) = 0
         xstate(jI+k) = 0

         if (k .eq. T) then
            x(k)      = 1.0d+3
            xstate(k) = 2
         end if
      end do

*     The first Capital is fixed.
*     The last three Investments are bounded.
*     Fudge them to be the normal ones for T = 10.

      scale        = T / 10.0d+0
      xupp(1)      = xlow(1)
      x(1)         = xlow(1)
      xstate(1)    = 0
      xupp(jI+T-2) = 0.112d+0 * scale
      xupp(jI+T-1) = 0.114d+0 * scale
      xupp(jI+T  ) = 0.116d+0 * scale

*     Set the bounds on F.
*     The T    linear (Capacity) components are <=.
*     The T nonlinear (Money)    components are >=.
*     We no longer need to set initial values and states for slacks
*     (assuming SNOPT does a cold start).

      jM     = T
      jY     = 0

      do k = 1, T
         Flow(jY+k) = bminus
         Fupp(jY+k) = zero
         Flow(jM+k) = zero
         Fupp(jM+k) = bplus

*-       x (jM+k) = zero
*-       x (jY+k) = zero
*-       Fstate(jM+k) = 0
*-       Fstate(jY+k) = 0
      end do

*     The last Money and Capacity components have a Range.

      Fupp(jM+T) =   10.0d+0
      Flow(jY+T) = - 20.0d+0

*     Set the objective and its bounds.

      ObjRow  = m
      Fmul(ObjRow) = zero
      Flow(ObjRow) = bminus
      Fupp(ObjRow) = bplus

*     Initialize pi.

      do i = 1, T
         Fmul(jM+i) = - one
         Fmul(jY+i) = + one
      end do

      iSpecs =  4
      iPrint = 15
      iSumm  =  6

      lfile = 'manne.spc'
      open( iSpecs, file=lfile, status='OLD',     err=800 )

      lfile = 'manne.out'
      open( iPrint, file=lfile, status='UNKNOWN', err=800 )

      if (iSumm .gt. 0)
     &   write(iSumm , *) 'Problem MANNE.    T =', T
      if (iPrint .gt. 0)
     &   write(iPrint, *) 'Problem MANNE.    T =', T

      return

  800 write(iErr, 4000) 'Error while opening file', lfile
 4000 format(/  a, 2x, a  )

      end ! subroutine usrini

