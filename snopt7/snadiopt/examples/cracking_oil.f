*     The cracking oil problem
*     Mike Gertz 20-Jul-2000
*
*     \documentclass{article}
*     \bibliographystyle{plain}
*
*     \newcommand\nT{N}
*     \newcommand\minim{\mathop{\hbox{\rm minimize}}}
*     \newcommand\subject{\hbox{\rm subject to}}
*     \newcommand\T{^T\!}
*
*     \begin{document}
*
*     Let $0 = t_0 < t_1 < \cdots < t_N$ be $N + 1$ given points
*     dividing the interval $[0, t_N]$ into $N$ subintervals. Let $z_1$
*     and $z_2$ be given, constant vectors with $N$ elements. The
*     cracking oil problem is to find the vector $\theta$ that minimizes
*     \begin{equation}
*       \label{objective}
*          \sum_{k = 1}^{20}
*             \left\{ ( y_1(t_k) - z_1^{(k)} ) ^ 2  +
*                     ( y_2(t_k) - z_2^{(k)} ) ^ 2 \right\}
*     \end{equation}
*
*     subject to the constraint that the function $y(t)$ satisfies the
*     ODE
*     \begin{equation}
*       \label{ode}
*       \begin{array}{lcl}
*         \dot{y}_1 + (\theta_1 + \theta_3) * y_1^2     & =    & 0 \\
*         \dot{y}_2 - \theta_1 * y_1^2 + \theta_2 * y_2 & =    & 0 \\
*                                                y(t_0) & =    & (1, 0)\T
*       \end{array}
*     \end{equation}
*     and that the variables satisfy the bounds
*     \begin{equation}
*       \label{bounds}
*       \begin{array}{rcccl}
*         0 & \leq & y(t)   & \leq & 1 \\
*         0 & \leq & \theta & \leq & 20
*       \end{array}
*     \end{equation}
*
*     We follow \cite{Averick1992, Tjoa1991} and attack this problem
*     using a collocation scheme that approximates the functions
*     $y_1(t)$ and $y_2(t)$ by continuous, piecewise fourth-order
*     polynomials. Four collocation points are chosen non-uniformly in
*     each interval in a manner that this collocation scheme equivalent
*     to Gaussian-Legendre quadrature (note that none of the
*     collocation points are the endpoints of the intervals). We
*     require the bounds to be satisfied at the nodes, $t_k, \ k = 1,
*     \ldots, N$ {\em and} the collocation points.
*
*     \bibliography{cracking_oil}
*
*     \end{document}
*
*     @TechReport{Averick1992,
*       author =		 {B.M. Averick and R.G. Carter and J.J. Mor\'{e}
*                       and G.L. Xue},
*       title =		 {The MINPACK-2 test problem collection},
*       institution =	 {Mathematics and Computer Science Division,
*                       Argonne National Laboratory},
*       year =		 1992,
*       type =         Preprint,
*       number =		 {MCS-P153-0694},
*       address =		 {Argonne, Illinois}
*     }
*
*     @Article{Tjoa1991,
*       author =		 {I.B. Tjoa and L.T. Biegler},
*       title =		 {Simultaneous solution and optimization
*                       strategies for parameter estimation of
*                       differential-algebraic equations systems },
*       journal =		 {Ind. Eng. Chem. Res.},
*       year =		 1991,
*       volume =		 30,
*       pages =		 {376--385}
*     }


****************************************************************
      subroutine crack1
     &   ( v1, w1, y1, v2, w2, y2, theta, h, f )

      implicit
     &     none

      integer              np
      parameter          ( np = 4 )

      double precision     v1, w1(np), y1(np + 1)
      double precision     v2, w2(np), y2(np + 1)
      double precision     theta(3)
      double precision     h, f(2 * np)

      double precision     y1dot(np), y2dot(np)
      integer              k

      double precision     rho(np + 1)
      data rho /
     &     .0694318413734436035d0, .330009490251541138d0,
     &     .669990539550781250d0,  .930568158626556396d0,
     &     1 /


*     Compute the values of an order $\np$ polynomial at $\np + 1$
*     points $h * \rho_k, \ k = 1, ..., np + 1$. The first np values in
*     $\rho$ are the points needed for Gaussian quadrature.

      call yanddot( v1, w1, np, rho, np + 1, h, y1, y1dot, np )
      call yanddot( v2, w2, np, rho, np + 1, h, y2, y2dot, np )

*     The first np points are collocation points. We must satisfy
*     the ODE for the cracking-oil problem:
*     \begin{eqnarray}
*         \dot{y}_1 + (theta_1 + theta_3) * y_1^2     & = & 0 \\
*         \dot{y}_2 - theta_1 * y_1^2 + theta_2 * y_2 & = & 0
*     \end{eqnarray}
*     at those points.
*
*     The last point will be used in continuity constraints.
      do k = 1, np
         f(2 * k - 1) = y1dot(k) + ( theta(1) + theta(3) ) * y1(k)**2
         f(2 * k )    = y2dot(k) - theta(1) * y1(k)**2 +
     &        theta(2) * y2(k)
      end do

      end

****************************************************************

      subroutine yanddot( v, w, np, rho, my, h, y, ydot, mydot )
      implicit    none

*     Find approximate values and derivative values at $\my$ points using
*     interpolation of a polynomial of degree $\np$.  $v$ and $w$ are the
*     coeficients of this polynomial, which is represented by the
*     formula:
*
*     \begin{equation}
*     y(t) = v + \sum_{j=1}^{np} \frac{ t^j } { j! h^{j-1} } w_j
*     \end{equation}
*
*     where $t = \rho(k) * h$. The variables $y$ and $\ydot$ hold the
*     approximate function and derivative values. Since we may want
*     fewer derivative values than function values, we only
*     compute $\mydot \leq \my$ derivatives.

      integer              np
      double precision     v, w(np)
      integer              my
      double precision     rho(my), h
      integer              mydot
      double precision     y(my),   ydot(mydot)

      integer              j, k
      integer              lastdt
      double precision     acc

      lastdt = min( my, mydot )
*     acc - holds rho^j/j!
      do k = 1, lastdt
         y(k)    = v
         ydot(k) = 0.0d0
         acc     = 1.0d0
         do j = 1, np
            ydot(k)     = ydot(k) + acc * w(j)
            acc         = acc * rho(k)/j
            y(k)        = y(k)    + h * acc * w(j)
         end do
      end do
      do k = lastdt, my
         y(k)    = v
         acc     = 1.0d0
         do j = 1, np
            acc         = acc *  rho(k) / j
            y(k)        = y(k)     + h * acc * w(j)
         end do
      end do

      end

****************************************************************

      subroutine usrfun( Status, mode,
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

*     nT - the number of time intervals
*     np - the degree of polynomial interpolation
      integer        nT, np
      parameter    ( nT = 20, np = 4 )

      double precision     y1(np+1), y2(np+1)
      integer              mcnt
      parameter          ( mcnt = 4 * np + 2 )

*     itheta - the index of theta in x. "theta" comes after all the
*     state variables There np + 1 parameters are needed to represent
*     each state variable in each time block. There are two state
*     variables and nT time blocks. Thus, theta starts
*     at 2 * nT * (np + 1) + 1

      integer              itheta
      parameter          ( itheta = 2 * nT * (np + 1) + 1 )
      double precision     z1dat(nT+1), z2dat(nT+1), taudat(nT+1)

*     z1dat, z2dat - tabular data. We are trying to find a least squares
*     fit of y1 to z1data and y2 to z2dat

      data z1dat / 1.0000, 0.8105, 0.6208, 0.5258, 0.4345, 0.3903,
     &     0.3342, 0.3034, 0.2735, 0.2405, 0.2283, 0.2071, 0.1669,
     &     0.1530, 0.1339, 0.1265, 0.1200, 0.0990, 0.0870, 0.0770,
     &     0.0690 /

      data z2dat / 0, 0.2000, 0.2886, 0.3010, 0.3215, 0.3123, 0.2716,
     &     0.2551, 0.2258, 0.1959, 0.1789, 0.1457, 0.1198, 0.0909,
     &     0.0719, 0.0561, 0.0460, 0.0280, 0.0190, 0.0140, 0.0100 /

*     taudat - the time points at which the tabular data lies.
      data taudat / 0.0,
     &     0.025, 0.05, 0.075, 0.10, 0.125, 0.150, 0.175,
     &     0.20, 0.225, 0.250, 0.30, 0.35, 0.40, 0.45,
     &     0.50, 0.55, 0.65, 0.75, 0.85, 0.95 /

      save z1dat, z2dat, taudat

      call usrblk( taudat, nT, z1dat, z2dat,
     &     x, np, x(itheta), F(1), F(2), y1, y2 )

      end


****************************************************************
      subroutine usrblk( taudat, nT, z1dat, z2dat,
     &     x, np, theta, phi, F, y1, y2 )
      implicit           none
*
*     Evaluate the objective and constraints over all the time blocks
*
*     nT - the number of time interval
*
*     np - the order of polynomal interpolation in each interval
*
*     z1dat, z2dat - tabular data. We are trying to find a least squares
*         fit of y1 to z1data and y2 to z2dat.
*
*     taudat - the time points at which the tabular data lies.
*
*     x(np + 1, 2, nT) - the variables. There are np + 1 parameters
*         used to represent each state variable, 2 state variables and
*         nt time itervals
*
*     theta - the parameters
*
*     phi - the objective value
*
*     F( 4 * np + 2, nT ) - the constraints
*     There are
*         2 * np collocation constraints
*         2 * np bounds on the variables at the collocation points
*         2 continuity constraints.
*     which yeilds 4 * np + 2 constraints over nT time intervals
*
*     y1(np + 1), y2(np + 1) - Work arrays. For each time iterval, these
*          arrays will hold the computed value of the state variables
*          at the collocation points and at the next grid point

      integer            nT
      double precision   taudat(nT + 1), z1dat(nT + 1), z2dat(nT + 1)
      integer            np
      double precision   x( np + 1, 2, nT )
      double precision   theta(3)
      double precision   phi

      double precision   F( 4 * np + 2, nT )
      double precision   y1(np + 1), y2(np + 1)

      integer            j, k

      phi = 0.0d0
      do k = 1, nT
*        Set the first 2 * np elts of the current f
         call crack1(
     &        x(1, 1, k),  x(2, 1, k), y1,
     &        x(1, 2, k),  x(2, 2, k), y2,
     &        theta,
     &        taudat(k+1) - taudat(k),
     &        f(1, k) )


*        The rest of the y's must be between zero and one
         do j = 1, np
            f( 2 * np + 2 * j - 1, k ) = y1( j )
            f( 2 * np + 2 * j,     k ) = y2( j )
         end do

         if ( k .lt. nT ) then
*           Set the continuity constraints
            f( 4 * np + 1, k ) = y1(np + 1) - x(1, 1, k + 1)
            f( 4 * np + 2, k ) = y2(np + 1) - x(1, 2, k + 1)
         else
*           There is no continuity constraint for the last grid
*           point, but the computed y is constrained to lie between
*           zero and one
            f( 4 * np + 1, k ) = y1(np + 1)
            f( 4 * np + 2, k ) = y2(np + 1)
         end if
*        objective - the square distance between the computed y's at the
*        grid points and the tabular data.
         phi = phi + (x(1, 1, k) - z1dat(k))**2
     &        + (x(1, 2, k) - z2dat(k))**2
      end do

*     Add the square distance of the last computed y to the objective
*     function.
      phi = phi + (y1(np + 1) - z1dat(nT+1))**2
     &     + (y2(np + 1) - z2dat(nT + 1))**2

      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine usrini( ObjAdd, ObjRow, Prob,
     &     x, xlow, xupp, xstate,
     &     Names, Fmul, Flow, Fupp, Fstate, FNames,
     &     iSpecs, iPrint, iSumm, iErr,
     &     cu,  iu,  ru,
     &     cw,  iw,  rw )

      implicit           none

      integer            ObjRow
      double precision   ObjAdd
      integer            n, neF, nName, nFnames
      integer            lencw, leniw, lenrw, lencu, leniu, lenru

      parameter        ( lencw = 501, leniw = 100000, lenrw = 200000 )
      parameter        ( lencu =   1, leniu =      1, lenru =      1 )
      parameter        ( n     = 203, neF   =    361 )
      parameter        ( nName =   1, nFnames =    1 )

      character*8        Prob, Names(nName)
      double precision   x(n), xlow(n), xupp(n)
      double precision   Flow(neF), Fupp(neF), Fmul(neF)
      integer            xstate(n), Fstate(n)
      character*8        FNames(nFnames)

      integer            iSpecs, iPrint, iSumm, iErr

      integer            iu(leniu), iw(leniw)
      character*8        cu(lencu), cw(lencw)
      double precision   ru(lenru), rw(lenrw)

      character*32        lfile

*     nT - the number of time intervals
*     np - the degree of polynomial interpolation
      integer        nT, np
      parameter    ( nT = 20, np = 4 )

*     itheta - the index of theta in x. "theta" comes after all the
*     state variables There np + 1 parameters are needed to represent
*     each state variable in each time block. There are two state
*     variables and nT time blocks. Thus, theta starts
*     at 2 * nT * (np + 1) + 1
      integer      itheta
      parameter  ( itheta = 2 * nT * (np + 1) + 1 )

      double precision z1dat(nT + 1), z2dat(nT + 1)

      data z1dat / 1.0000, 0.8105, 0.6208, 0.5258, 0.4345, 0.3903,
     &     0.3342, 0.3034, 0.2735, 0.2405, 0.2283, 0.2071, 0.1669,
     &     0.1530, 0.1339, 0.1265, 0.1200, 0.0990, 0.0870, 0.0770,
     &     0.0690 /

      data z2dat / 0, 0.2000, 0.2886, 0.3010, 0.3215, 0.3123, 0.2716,
     &     0.2551, 0.2258, 0.1959, 0.1789, 0.1457, 0.1198, 0.0909,
     &     0.0719, 0.0561, 0.0460, 0.0280, 0.0190, 0.0140, 0.0100 /

*     ------------------------------------------------------------------
*     Give the problem a name.

      Prob   = 'crackoil'
      ObjRow = 1

*     Initialize the variables and constraints for the various time
*     blocks

      call iniblk( z1dat, nT, z2dat,
     &     x, np, xlow, xupp,
     &     Flow(2), Fupp(2),
     &     x(itheta), xlow(itheta), xupp(itheta) )

      iSpecs = 4
      iPrint = 15

      lfile = 'cracking_oil.spc'
      open( iSpecs, file=lfile, status='OLD',     err=800 )

      lfile = 'cracking_oil.out'
      open( iPrint, file=lfile, status='UNKNOWN', err=800 )

      return

 800  write(iErr, 4000) 'Error while opening file', lfile

 4000 format(/  a, 2x, a  )

      end
****************************************************************
      subroutine iniblk( z1dat, nT, z2dat,
     &     x, np, xlow, xupp,
     &     Flow, Fupp,
     &     theta, thlow, thupp )
*
*    Initialize the variables and constriants by time block
*
*     nT - the number of time interval
*     np - the order of polynomal interpolation in each interval
*
*     z1dat, z2dat - tabular data. We are trying to find a least squares
*         fit of the state varibles to z1dat and z2dat.
*
*     x(np + 1, 2, nT), xlow(,,), xupp(,,) - the variables and bounds
*         on the variables. There are np + 1 parameters
*         used to represent each state variable, 2 state variables and
*         nT  time itervals
*
*     Flow( 4*np +  2, nT), Fupp(,) - bounds on the constraint functions
*     There are
*         2 * np collocation constraints
*         2 * np bounds on the variables at the collocation points
*         2 continuity constraints.
*     which yeilds 4 * np + 2 constraints over nT time intervals
*
*     theta(3), thlow(3), thupp(3) - the "parameters" and bounds on the
*         parameters. The parameters are values that we are allowed to
*         change to influence the ODE
      implicit none


      integer              nT
      double precision     z1dat(nT), z2dat(nT)

      integer              np
*     We are doing polynomial interpolation of the state variables. Our
*     polynomials are of order np, which means that there are nblock =
*     np + 1 variables for each block of state variables. Since there
*     are two state variables, there are nv = 2 * nblock variables for
*     each interval

*     Note, 2*(np + 1) .eq. nv
      double precision     x   ( np + 1, 2, nT )
      double precision     xlow( np + 1, 2, nT )
      double precision     xupp( np + 1, 2, nT )

*     mcnt - the number of constraints in a normal time interval (in
*     other words, any interval but the last.)
*     There are
*         2 * np collocation constraints
*         2 * np bounds on the variables at the collocation points
*         2 continuity constraints.

*     Note, 4 * np + 2 .eq. mcnt
      double precision     Flow( 4 * np + 2, nT )
      double precision     Fupp( 4 * np + 2, nT )
      double precision     theta(3), thlow(3), thupp(3)

      integer i, j, k

      double precision   zero, pt5, one, two, three, plInfy
      parameter          ( zero  = 0.0d+0, pt5    = 0.5d+0  )
      parameter          ( one   = 1.0d+0, two    = 2.0d+0  )
      parameter          ( three = 3.0d+0, plInfy = 1.0d+20 )

      do k = 1, nT
*        x(1, 1, *) and x(1, 2, *) are the values of the states at the
*        gridpoints (which are not the collocation points). Bound them
*        between zero and one
         x   ( 1, 2, k ) = z1dat(k)
         xlow( 1, 2, k ) = zero
         xupp( 1, 2, k ) = one

         x   ( 1, 2, k ) = z2dat(k)
         xlow( 1, 2, k ) = zero
         xupp( 1, 2, k ) = one

*        Set the rest of the variables to 0.5
         do i = 1, np
            x( 1 + i, 1, k ) = pt5
            x( 1 + i, 2, k ) = pt5
         end do


      end do

*     Set the initial conditions on the ODE
      xlow( 1, 1, 1 ) = one
      xupp( 1, 1, 1 ) = one

      xlow( 1, 2, 1 ) = zero
      xupp( 1, 2, 1 ) = zero

      do k = 1, nT
*        Make the ODEs into equality constraints
         do i = 1, 2 * np
            Flow( i, k ) = zero
            Fupp( i, k ) = zero
         end do

*        Set bounds on the computed y's
         do i = 2 * np + 1, 4 * np
            Flow( i, k ) = zero
            Fupp( i, k ) = one
         end do
      end do

      do k = 1, nT - 1
*        Set bounds on the continuity constraints
         do i = 4 * np + 1, 4 * np + 2
            Flow( i, k ) = zero
            Fupp( i, k ) = zero
         end do
      end do
*     There is no continuity constraint on the last grid point, but
*     there is the constraint that the computed value lies between
*     zero and one
      do i = 4 * np + 1, 4 * np + 2
         Flow( i, nT ) = zero
         Fupp( i, nT ) = one
      end do

*     the elements of theta must lie between 0 and 20
      do k = 1, 3
         theta(k)       = zero
         thlow(k)       = zero
         thupp(k)       = 20d0
      end do

      end
