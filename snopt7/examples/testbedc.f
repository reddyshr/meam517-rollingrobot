!     ------------------------------------------------------------------
!     File another_bug.f
!
!     This solves an LP in fully nonlinear mode and partly linear mode
!     with two different answers.
!     ------------------------------------------------------------------
      program
     &     tester

      implicit
     &     none
      integer
     &     maxm, maxn, maxne, nName
      parameter
     &     ( maxm   = 1000,
     &       maxn   = 1000,
     &       maxne  = 3000,
     &       nName  = 1 )

      character
     &     Prob*8, Names(nName)*8
      integer
     &     indA(maxne) , hs(maxn+maxm)
      integer
     &     locA(maxn+1)
      double precision
     &     Acol(maxne) , bl(maxn+maxm), bu(maxn+maxm),
     &     x(maxn+maxm), pi(maxm)     , rc(maxn+maxm),
     &     xLin(maxn+maxm), xNonlin(maxn+maxm)
      integer
     &     lenrw, leniw, lencw
      !-----------------------------------------------------------------
      ! SNOPT workspace
      parameter          (  lenrw = 20000)
      double precision   rw(lenrw)
      parameter          (  leniw = 10000)
      integer            iw(leniw)
      parameter          (  lencw =   500)
      character*8        cw(lencw)
      !-----------------------------------------------------------------
      logical
     &     byname
      character
     &     lfile*20
      external
     &     nonlinearUsrfun,
     &     linearUsrfun
      integer
     &     Errors, i1, i2, INFO, iObj, iPrint, iSpecs, iSumm, itnlim, j,
     &     m, mincw, miniw, minrw, n, ne, nInf, nnCon, nnJac, nnObj,
     &     nOut, nS
      double precision
     &     Obj, ObjAdd, sInf
      !-----------------------------------------------------------------
      ! Specify some of the SNOPT files.
      ! iSpecs  is the Specs file   (0 if none).
      ! iPrint  is the Print file   (0 if none).
      ! iSumm   is the Summary file (0 if none).
      ! nOut    is an output file used here by snmain.

      iSpecs = 4
      iPrint = 9
      iSumm  = 6
      nOut   = 6

      byname = .true.

      if ( byname ) then
         ! Unix and DOS systems.  Open the Specs and print files.
         lfile = 'testbedc.spc'
         open( iSpecs, file=lfile, status='OLD',     err=800 )

         lfile = 'testbedc.out'
         open( iPrint, file=lfile, status='UNKNOWN', err=800 )
      end if

      !-----------------------------------------------------------------
      ! First,  snInit MUST be called to initialize optional parameters
      ! to their default values.
      !-----------------------------------------------------------------
      call snInit
     &   ( iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw )

      !-----------------------------------------------------------------
      ! Read a Specs file (Optional).
      !-----------------------------------------------------------------
      call snSpec
     &   ( iSpecs, INFO, cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .ne. 101  .and.  INFO .ne. 107) then
         go to 990
      end if

      ! Set up the data structure for the sparse Jacobian.
      ! Assign dummy values for the nonlinear elements.

       call NonlinearDat
     &    ( Prob, maxm, maxn, maxne, INFO,
     &     m, n, ne, nnCon, nnObj, nnJac, iObj, ObjAdd,
     &     Acol, indA, locA, bl, bu, hs, x, pi )
      !-----------------------------------------------------------------
      ! Specify any options not set in the Specs file.
      ! i1 and i2 may refer to the Print and Summary file respectively.
      ! Setting them to 0 suppresses printing.
      !-----------------------------------------------------------------
      Errors = 0

      itnlim = 250
      i1     =   0
      i2     =   0
      call snSeti
     &   ( 'Iterations        ', itnlim, i1, i2, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      !-----------------------------------------------------------------
      ! Solve
      !-----------------------------------------------------------------
      call snOptC
     &   ( 'Cold', m, n, ne, nName,
     &     nnCon, nnObj, nnJac,
     &     iObj, ObjAdd, Prob,
     &     nonlinearUsrfun,
     &     Acol, indA, locA, bl, bu, Names,
     &     hs, x, pi, rc,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf, Obj,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
         go to 910
      end if

      do j = 1, n
         xNonlin(j) = x(j)
      end do

      write(nOut, *) ' '
      write(nOut, *) 'snOptC finished. (nonlinear)'
      write(nOut, *) 'Input  errors =', Errors
      write(nOut, *) 'snOptC INFO   =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      if (iObj .gt. 0) then
         write(nOut, *) 'Obj           =', ObjAdd + x(n+iObj) + Obj
      else
         write(nOut, *) 'Obj           =', ObjAdd + Obj
      end if
      write(nOut, *) 'xopt          =', x(1), x(2)
      if (INFO .ge. 30) go to 910

      !-----------------------------------------------------------------
      ! Solve
      !-----------------------------------------------------------------
      call LinearDat
     &   ( Prob, maxm, maxn, maxne, INFO,
     &     m, n, ne, nnCon, nnObj, nnJac, iObj, ObjAdd,
     &     Acol, indA, locA, bl, bu, hs, x, pi )

c$$$      do j = 1, n
c$$$         x(j) = xNonlin(j)
c$$$      end do
c$$$      x(1) = -4.0d+0
c$$$      x(2) = -5.5917625347658149d+0

      call snOptC
     &   ( 'Cold', m, n, ne, nName,
     &     nnCon, nnObj, nnJac,
     &     iObj, ObjAdd, Prob,
     &     linearUsrfun,
     &     Acol, indA, locA, bl, bu, Names,
     &     hs, x, pi, rc,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf, Obj,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
         go to 910
      end if

      write(nOut, *) ' '
      write(nOut, *) 'snOptC finished. (linear)'
      write(nOut, *) 'Input  errors =', Errors
      write(nOut, *) 'snOptC INFO   =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      if (iObj .gt. 0) then
         write(nOut, *) 'Obj           =', ObjAdd + x(n+iObj) + Obj
      else
         write(nOut, *) 'Obj           =', ObjAdd + Obj
      end if
      write(nOut, *) 'xopt          =', x(1), x(2)
      if (INFO .ge. 30) go to 910
      stop

      !-----------------------------------------------------------------
      ! Error exit.
      !-----------------------------------------------------------------
  800 write(nOut, 4000) 'Error while opening file', lfile
      stop

  910 write(nOut, *) ' '
      write(nOut, *) 'STOPPING because of error condition'

  990 stop

 4000 format(/  a, 2x, a  )

      end ! program snoptc

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine nonlinearUsrfun
     &   ( mode, nnObj, nnCon, nnJac, nnL, neJac,
     &     x, fObj, gObj, fCon, gCon,
     &     nState, cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none
      integer
     &     mode, nnObj, nnCon, nnJac, nnL, neJac, nState,
     &     lencu, leniu, lenru, iu(leniu)
      double precision
     &     x(nnL)
      double precision
     &     fObj, gObj(nnObj)
      double precision
     &     fCon(nnCon), gCon(neJac), ru(lenru)
      character
     &     cu(lencu)*8

      fObj = 7*x(1) - 4*x(2)

      gObj(1) =  7.0d+0
      gObj(2) = -4.0d+0

      !-----------------------------------------------------------------
      ! Constraints.
      !-----------------------------------------------------------------
      fCon(1) = -8*x(1) + 6*x(2)

      gCon(1) = -8.0d+0
      gCon(2) =  6.0d+0

      end ! subroutine nonlinearUsrfun

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine NonlinearDat
     &   ( Prob, maxm, maxn, maxne, INFO,
     &     m, n, ne, nnCon, nnObj, nnJac, iObj, ObjAdd,
     &     Acol, indA, locA, bl, bu, hs, x, pi )

      implicit
     &     none
      integer
     &     maxm, maxn, maxne, INFO, m, n, ne, nnCon, nnObj, nnJac,
     &     iObj, indA(maxne) , hs(maxn+maxm), locA(maxn+1)
      double precision
     &     ObjAdd, Acol(maxne) , bl(maxn+maxm), bu(maxn+maxm),
     &     x(maxn+maxm), pi(maxm)
      character*8
     &     Prob
      integer
     &     i, j
      !-----------------------------------------------------------------
      ! Define a name for the Problem.

      Prob   = 'Nonlin  '

      ne     = 2
      n      = 2
      m      = 1

      nnCon  = 1
      nnJac  = 2
      nnObj  = 2

      iObj    = 0
      ObjAdd  = 0.0d+0

      ! Check if there is enough storage.

      INFO   = 0
      if (m     .gt. maxm ) INFO = 1
      if (n     .gt. maxn ) INFO = 1
      if (ne    .gt. maxne) INFO = 1
      if (INFO  .gt.   0  ) return

      !--------------------------------------------
      ! Set up the data structure for the Jacobian.
      ! The nonzeros are stored columnwise.
      ! Nonlinear elements must be set (0 is okay).
      !--------------------------------------------

      !===========================================
      ! Column  1
      !===========================================
      locA(1) = 1

      indA(1) = 1
      Acol(1) = 0.0d+0  ! 0 because nonlinear

      !===========================================
      ! Column  2
      !===========================================
      locA(2) = 2

      indA(2) = 1
      Acol(2) = 0.0d+0  ! 0 because nonlinear

      !===========================================
      ! Don't forget to finish off  locA.
      ! This is crucial.
      !===========================================
      locA(3) =  ne + 1 ! = 3

      !-----------------------------------------------------------------
      ! Constraint ranges
      !-----------------------------------------------------------------
      bl(n+1) = -4.0d+0
      bu(n+1) = 1.0d+0

      !-----------------------------------------------------------------
      ! Variable ranges
      !-----------------------------------------------------------------
      bl(1) = -4.0d+0
      bu(1) = -1.0d+0

      bl(2) = -7.0d+0
      bu(2) =  9.0d+0

      !-----------------------------------------------------------------
      ! Initialize x, hs and pi.
      !-----------------------------------------------------------------
      x(1)   =  0d+0
      x(2)   =  0d+0

      do j = 1, n
         hs(j)  = 0
      end do

      do i = 1, m
         pi(i)  = 0.0d+0
      end do

      end ! subroutine nonlinearDat


*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine linearUsrfun
     &   ( mode, nnObj, nnCon, nnJac, nnL, neJac,
     &     x, fObj, gObj, fCon, gCon,
     &     nState, cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none
      integer
     &     mode, nnObj, nnCon, nnJac, nnL, neJac, nState,
     &     lencu, leniu, lenru, iu(leniu)
      double precision
     &     x(nnL)
      double precision
     &     fObj, gObj(nnObj)
      double precision
     &     fCon(nnCon), gCon(neJac), ru(lenru)
      character
     &     cu(lencu)*8

      ! The objective is totally linear but we put nonlinear jacobian
      ! variables first, so set nonlinear objective and it's gradient
      ! to 0 as per the manual.

      fObj    = 0.0d+0
      gObj(1) = 0.0d+0  ! objective totally linear

      !-----------------------------------------------------------------
      ! Constraints.
      !-----------------------------------------------------------------
      fCon(1) = -8.0d+0*x(1)
      gCon(1) = -8.0d+0

      end ! subroutine linearUsrfun

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine LinearDat
     &   ( Prob, maxm, maxn, maxne, INFO,
     &     m, n, ne, nnCon, nnObj, nnJac, iObj, ObjAdd,
     &     Acol, indA, locA, bl, bu, hs, x, pi )

      implicit
     &     none
      integer
     &     maxm, maxn, maxne, INFO, m, n, ne, nnCon, nnObj, nnJac,
     &     iObj, indA(maxne) , hs(maxn+maxm), locA(maxn+1)
      double precision
     &     ObjAdd, Acol(maxne) , bl(maxn+maxm), bu(maxn+maxm),
     &     x(maxn+maxm), pi(maxm)
      character*8
     &     Prob
      integer
     &     i, j
      !-----------------------------------------------------------------
      ! Define a name for the Problem.

      Prob   = 'linear  '

      ne     = 4
      n      = 2
      m      = 2

      nnCon  = 1
      nnJac  = 1
      nnObj  = 1

      iObj    = 2
      ObjAdd  = 0.0d+0

      ! Check if there is enough storage.

      INFO   = 0
      if (m     .gt. maxm ) INFO = 1
      if (n     .gt. maxn ) INFO = 1
      if (ne    .gt. maxne) INFO = 1
      if (INFO  .gt.   0  ) return

      !--------------------------------------------
      ! Set up the data structure for the Jacobian.
      ! The nonzeros are stored columnwise.
      ! Nonlinear elements must be set (0 is okay).
      !--------------------------------------------

      !===========================================
      ! Column  1
      !===========================================
      locA(1) = 1

      indA(1) = 1
      Acol(1) = 0.0d+0  ! 0 because nonlinear

      indA(2) = 2
      Acol(2) = 7.0d+0  ! objective row

      !===========================================
      ! Column  2
      !===========================================
      locA(2) = 3

      indA(3) = 1
      Acol(3) = 6.0d+0

      indA(4) = 2
      Acol(4) = -4.0d+0  ! objective row

      !===========================================
      ! Don't forget to finish off  locA.
      ! This is crucial.
      !===========================================
      locA(3) =  ne + 1 ! = 5

      !-----------------------------------------------------------------
      ! Constraint ranges
      !-----------------------------------------------------------------
      bl(n+1) = -4.0d+0
      bu(n+1) =  1.0d+0

      bl(n+2) = -1.1d+20
      bu(n+2) =  1.1d+20

      !-----------------------------------------------------------------
      ! Variable ranges
      !-----------------------------------------------------------------
      bl(1) = -4.0d+0
      bu(1) = -1.0d+0

      bl(2) = -7.0d+0
      bu(2) =  9.0d+0

      !-----------------------------------------------------------------
      ! Initialize x, hs and pi.
      !-----------------------------------------------------------------
      x(1)   =  0d+0
      x(2)   =  0d+0

      do j = 1, n
         hs(j)  = 0
      end do

      do i = 1, m
         pi(i)  = 0.0d+0
      end do

      end ! subroutine linearDat
