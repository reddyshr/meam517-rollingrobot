!     ------------------------------------------------------------------
!     File testerA.f
!     This is a main program to test various weird problem formats for
!     the free format interface snOptA,  part of the SNOPT 7 package.
!
!     29 Dec 2002: First version for SNOPT 6
!     10 Jul 2004: Updated for SNOPT 7
!     09 Apr 2008: Current version.
!     ------------------------------------------------------------------
      program
     &     testerA
      implicit
     &     none
      integer
     &     maxF, maxn, nxname, nFname, lenA, lenG
      parameter
     &   ( maxF   = 30,
     &     maxn   = 10,
     &     lenA   = 50, lenG   = 100,
     &     nxname =  1, nFname =   1 )
      integer
     &     iAfun(lenA), jAvar(lenA), iGfun(lenG), jGvar(lenG),
     &     xstate(maxn), Fstate(maxF)
      character
     &     lfile*20, Prob*8, xnames(nxname)*8, Fnames(nFname)*8
      double precision
     &     ObjAdd, sInf, A(lenA), Flow(maxF), Fupp(maxF), F(maxF),
     &     xlow(maxn), xupp(maxn), x(maxn), Fmul(maxF), xmul(maxn)
      integer
     &     DerOpt, Errors, neA, neG, ObjRow, INFO,
     &     iPrt, iPrint, iSpecs, iSum, iSumm,
     &     Major, mincw, miniw, minrw, nF, n, nInf, nS, OutS, OutP
      external
     &     usrfun1, usrfun2, usrfun3, usrfun4
!     ------------------------------------------------------------------
!     SNOPT workspace

      integer               lenrw
      parameter          (  lenrw = 20000)
      double precision   rw(lenrw)
      integer               leniw
      parameter          (  leniw = 10000)
      integer            iw(leniw)
      integer               lencw
      parameter          (  lencw =   500)
      character          cw(lencw)*8

      integer             Cold,       Basis,      Warm
      parameter          (Cold   = 0, Basis  = 1, Warm  = 2)
!     ------------------------------------------------------------------
!     Specify some of the SNOPT files.
!     iSpecs  is the Specs file   (0 if none).
!     iPrint  is the Print file   (0 if none).
!     iSumm   is the Summary file (0 if none).
!
!     Out     is an output file used by the calling program.

      iSpecs =  4
      iSumm  =  6
      iPrint =  9
      OutS   =  6
      OutP   =  9

      lfile = 'testera.spc'
      open( iSpecs, file=lfile, status='OLD',     err=800 )
      lfile = 'testera.out'
      open( iPrint, file=lfile, status='UNKNOWN', err=800 )

!     ------------------------------------------------------------------
!     First,  snInit MUST be called to initialize optional parameters
!     to their default values.
!     ------------------------------------------------------------------
      call snInit
     &   ( iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw )

!     ------------------------------------------------------------------
!     Read a Specs file (Optional).
!     ------------------------------------------------------------------
      call snSpec
     &   ( iSpecs, INFO, cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .ne. 101  .and.  INFO .ne. 107) then
         go to 990
      end if

!     ------------------------------------------------------------------
!     Test 1
!     No constant Jacobian elements
!     ------------------------------------------------------------------
      write(OutS, *) ' '
      write(OutS, *) ' --------------------------------------'
      write(OutS, *) ' Test 1. No constant Jacobian elements.'
      write(OutS, *) ' --------------------------------------'
      write(OutP, *) ' '
      write(OutP, *) ' --------------------------------------'
      write(OutP, *) ' Test 1. No constant Jacobian elements.'
      write(OutP, *) ' --------------------------------------'

      Errors = 0

      call test1
     &   ( Errors, maxF, maxn,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     Prob, nF, n,
     &     ObjAdd, ObjRow, xlow, xupp, Flow, Fupp,
     &     x, xstate, xmul, F, Fstate, Fmul,
     &     cw, lencw, iw, leniw, rw, lenrw )
      if (Errors .gt. 0) go to 910

!     ------------------------------------------------------------------
!     Specify any options not set in the Specs file.
!     i1 and i2 may refer to the Print and Summary file respectively.
!     Setting them to 0 suppresses printing.
!     ------------------------------------------------------------------
      Major    = 250
      iPrt     =   0
      iSum     =   0

      call snseti
     &   ( 'Major Iteration limit', Major, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      call snOptA
     &   ( Cold, nF, n, nxname, nFname,
     &     ObjAdd, ObjRow, Prob, usrfun1,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     xlow, xupp, xnames, Flow, Fupp, Fnames,
     &     x, xstate, xmul, F, Fstate, Fmul,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
         go to 920
      end if

      write(OutS, *) ' '
      write(OutS, *) 'testera (1) finished.'
      write(OutS, *) 'Input errors  =', Errors
      write(OutS, *) 'snOptA INFO   =', INFO
      write(OutS, *) 'nInf          =', nInf
      write(OutS, *) 'sInf          =', sInf
      if (ObjRow .gt. 0)
     &write(OutS, *) 'Obj           =', ObjAdd + F(ObjRow)
      if (INFO .ge. 30) go to 920

!     ------------------------------------------------------------------
!     Test 2(a)
!     Normal run
!     ------------------------------------------------------------------
      write(OutS, *) ' '
      write(OutS, *) ' ----------------------'
      write(OutS, *) ' Test 2(a). Normal run.'
      write(OutS, *) ' ----------------------'
      write(OutP, *) ' '
      write(OutP, *) ' ----------------------'
      write(OutP, *) ' Test 2(a). Normal run.'
      write(OutP, *) ' ----------------------'

      Errors = 0

      call test2
     &   ( Errors, maxF, maxn,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     Prob, nF, n,
     &     ObjAdd, ObjRow, xlow, xupp, Flow, Fupp,
     &     x, xstate, xmul, F, Fstate, Fmul,
     &     cw, lencw, iw, leniw, rw, lenrw )
      if (Errors .gt. 0) go to 910

      call snOptA
     &   ( Cold, nF, n, nxname, nFname,
     &     ObjAdd, ObjRow, Prob, usrfun2,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     xlow, xupp, xnames, Flow, Fupp, Fnames,
     &     x, xstate, xmul, F, Fstate, Fmul,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
         go to 920
      end if

      write(OutS, *) ' '
      write(OutS, *) 'testera (2a) finished.'
      write(OutS, *) 'Input errors  =', Errors
      write(OutS, *) 'snOptA INFO   =', INFO
      write(OutS, *) 'nInf          =', nInf
      write(OutS, *) 'sInf          =', sInf
      if (ObjRow .gt. 0)
     &write(OutS, *) 'Obj           =', ObjAdd + F(ObjRow)
      if (INFO .ge. 30) go to 920

!     ------------------------------------------------------------------
!     Test 2(b)
!     Normal run derivative option not set.
!     ------------------------------------------------------------------
      write(OutS, *) ' '
      write(OutS, *) ' ---------------------------------------------'
      write(OutS, *) ' Test 2(b). Derivative level intead of option.'
      write(OutS, *) ' ---------------------------------------------'
      write(OutP, *) ' '
      write(OutP, *) ' ---------------------------------------------'
      write(OutP, *) ' Test 2(b). Derivative level intead of option.'
      write(OutP, *) ' ---------------------------------------------'

      Errors = 0

      call test2
     &   ( Errors, maxF, maxn,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     Prob, nF, n,
     &     ObjAdd, ObjRow, xlow, xupp, Flow, Fupp,
     &     x, xstate, xmul, F, Fstate, Fmul,
     &     cw, lencw, iw, leniw, rw, lenrw )
      if (Errors .gt. 0) go to 910

      DerOpt = 3
      call snSeti
     &   ( 'Derivative level', DerOpt, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      call snOptA
     &   ( Cold, nF, n, nxname, nFname,
     &     ObjAdd, ObjRow, Prob, usrfun2,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     xlow, xupp, xnames, Flow, Fupp, Fnames,
     &     x, xstate, xmul, F, Fstate, Fmul,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
         go to 920
      end if

      write(OutS, *) ' '
      write(OutS, *) 'testera (2b) finished.'
      write(OutS, *) 'Input errors  =', Errors
      write(OutS, *) 'snOptA INFO   =', INFO
      write(OutS, *) 'nInf          =', nInf
      write(OutS, *) 'sInf          =', sInf
      if (ObjRow .gt. 0)
     &write(OutS, *) 'Obj           =', ObjAdd + F(ObjRow)
      if (INFO .ge. 30) go to 920

!     ------------------------------------------------------------------
!     Test 3(a)
!     Zero rows and columns in the Jacobian
!     ------------------------------------------------------------------
      write(OutS, *) ' '
      write(OutS, *) ' ---------------------------------'
      write(OutS, *) ' Test 3(a). Zero rows and columns.'
      write(OutS, *) ' ---------------------------------'
      write(OutP, *) ' '
      write(OutP, *) ' ---------------------------------'
      write(OutP, *) ' Test 3(a). Zero rows and columns.'
      write(OutP, *) ' ---------------------------------'
      Errors = 0

      call test3
     &   ( Errors, maxF, maxn,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     Prob, nF, n,
     &     ObjAdd, ObjRow, xlow, xupp, Flow, Fupp,
     &     x, xstate, xmul, F, Fstate, Fmul,
     &     cw, lencw, iw, leniw, rw, lenrw )
      if (Errors .gt. 0) go to 910

      call snOptA
     &   ( Cold, nF, n, nxname, nFname,
     &     ObjAdd, ObjRow, Prob, usrfun3,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     xlow, xupp, xnames, Flow, Fupp, Fnames,
     &     x, xstate, xmul, F, Fstate, Fmul,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
         go to 920
      end if

      write(OutS, *) ' '
      write(OutS, *) 'testera (3a) finished.'
      write(OutS, *) 'Input errors  =', Errors
      write(OutS, *) 'snOptA INFO   =', INFO
      write(OutS, *) 'nInf          =', nInf
      write(OutS, *) 'sInf          =', sInf
      if (ObjRow .gt. 0)
     &write(OutS, *) 'Obj           =', ObjAdd + F(ObjRow)
      if (INFO .ge. 30) go to 920

!     ------------------------------------------------------------------
!     Test 3(b)
!     ObjRow = 0
!     ------------------------------------------------------------------
      write(OutS, *) ' '
      write(OutS, *) ' ----------------------'
      write(OutS, *) ' Test 3(b). ObjRow = 0.'
      write(OutS, *) ' ----------------------'
      write(OutP, *) ' '
      write(OutP, *) ' ----------------------'
      write(OutP, *) ' Test 3(b). ObjRow = 0.'
      write(OutP, *) ' ----------------------'

      Errors = 0

      call test2
     &   ( Errors, maxF, maxn,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     Prob, nF, n,
     &     ObjAdd, ObjRow, xlow, xupp, Flow, Fupp,
     &     x, xstate, xmul, F, Fstate, Fmul,
     &     cw, lencw, iw, leniw, rw, lenrw )
      if (Errors .gt. 0) go to 910

      ObjRow = 0                ! Find a feasible point

      call snOptA
     &   ( Cold, nF, n, nxname, nFname,
     &     ObjAdd, ObjRow, Prob, usrfun2,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     xlow, xupp, xnames, Flow, Fupp, Fnames,
     &     x, xstate, xmul, F, Fstate, Fmul,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
         go to 920
      end if

      write(OutS, *) ' '
      write(OutS, *) 'testera (3b) finished.'
      write(OutS, *) 'Input errors  =', Errors
      write(OutS, *) 'snOptA INFO   =', INFO
      write(OutS, *) 'nInf          =', nInf
      write(OutS, *) 'sInf          =', sInf
      if (ObjRow .gt. 0)
     &write(OutS, *) 'Obj           =', ObjAdd + F(ObjRow)
      if (INFO .ge. 30) go to 920

!     ------------------------------------------------------------------
!     Test 4
!     Linear program
!     ------------------------------------------------------------------
      write(OutS, *) ' '
      write(OutS, *) ' -----------------------'
      write(OutS, *) ' Test 4. Linear program.'
      write(OutS, *) ' -----------------------'
      write(OutP, *) ' '
      write(OutP, *) ' -----------------------'
      write(OutP, *) ' Test 4. Linear program.'
      write(OutP, *) ' -----------------------'

      Errors = 0

      call test4
     &   ( Errors, Prob, maxF, maxn, nF, n, neG,
     &     iAfun, jAvar, lenA, neA, A,
     &     ObjAdd, ObjRow, xlow, xupp, Flow, Fupp, x, xstate, Fmul,
     &     cw, lencw, iw, leniw, rw, lenrw )
      if (Errors .gt. 0) go to 910

      call snOptA
     &   ( Cold, nF, n, nxname, nFname,
     &     ObjAdd, ObjRow, Prob, usrfun4,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     xlow, xupp, xnames, Flow, Fupp, Fnames,
     &     x, xstate, xmul, F, Fstate, Fmul,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
         go to 920
      end if

      write(OutS, *) ' '
      write(OutS, *) 'testera (4) finished.'
      write(OutS, *) 'Input errors  =', Errors
      write(OutS, *) 'snOptA INFO   =', INFO
      write(OutS, *) 'nInf          =', nInf
      write(OutS, *) 'sInf          =', sInf
      if (ObjRow .gt. 0)
     &write(OutS, *) 'Obj           =', ObjAdd + F(ObjRow)
      if (INFO .ge. 30) go to 920

!     ------------------------------------------------------------------
!     Test 5
!     Linear program with no objective specified
!     ------------------------------------------------------------------
      write(OutS, *) ' '
      write(OutS, *) ' ----------------------------------------'
      write(OutS, *) ' Test 5. LP  with no objective specified.'
      write(OutS, *) ' ----------------------------------------'
      write(OutP, *) ' '
      write(OutP, *) ' ----------------------------------------'
      write(OutP, *) ' Test 5. LP  with no objective specified.'
      write(OutP, *) ' ----------------------------------------'

      Errors = 0

      call test4
     &   ( Errors, Prob, maxF, maxn, nF, n, neG,
     &     iAfun, jAvar, lenA, neA, A,
     &     ObjAdd, ObjRow, xlow, xupp, Flow, Fupp, x, xstate, Fmul,
     &     cw, lencw, iw, leniw, rw, lenrw )
      if (Errors .gt. 0) go to 910

      ObjRow = 0
      Prob   = 'testA 5 '

      call snOptA
     &   ( Cold, nF, n, nxname, nFname,
     &     ObjAdd, ObjRow, Prob, usrfun4,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     xlow, xupp, xnames, Flow, Fupp, Fnames,
     &     x, xstate, xmul, F, Fstate, Fmul,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
         go to 920
      end if

      write(OutS, *) ' '
      write(OutS, *) 'testera (5) finished.'
      write(OutS, *) 'Input errors  =', Errors
      write(OutS, *) 'snOptA INFO   =', INFO
      write(OutS, *) 'nInf          =', nInf
      write(OutS, *) 'sInf          =', sInf
      if (ObjRow .gt. 0)
     &write(OutS, *) 'Obj           =', ObjAdd + F(ObjRow)
      if (INFO .ge. 30) go to 920

      stop

!     ------------------------------------------------------------------
!     Error exit.
!     ------------------------------------------------------------------
  800 write(OutS, 4000) 'Error while opening file', lfile
      stop

  910 write(OutS, *) ' '
      write(OutS, *) 'Insufficient space to hold the problem'
      stop

  920 write(OutS, *) ' '
      write(OutS, *) 'STOPPING because of error condition'

  990 stop

 4000 format(/  a, 2x, a  )
 4010 format(/  a, 2x, i6 )

      end ! main program

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine test1
     &   ( Errors, maxF, maxn,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     Prob, nF, n,
     &     ObjAdd, ObjRow, xlow, xupp, Flow, Fupp,
     &     x, xstate, xmul, F, Fstate, Fmul,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     Errors, maxF, maxn, neA, lenA, neG, lenG, nF, n,
     &     ObjRow, lencw, leniw, lenrw, xstate(maxn), Fstate(maxF),
     &     iAfun(lenA), jAvar(lenA), iGfun(lenG), jGvar(lenG),
     &     iw(leniw)
      double precision
     &     ObjAdd,
     &     A(lenA), xlow(maxn), xupp(maxn), Flow(maxF), Fupp(maxF),
     &     x(maxn), F(maxF), xmul(maxn), Fmul(maxF), rw(lenrw)
      character
     &     Prob*8, cw(lencw)*8

!     ==================================================================
!     Test1 defines input data for the toy problem discussed in the
!     snoptA Users Guide.
!
!     Minimize      3*x(1) + 5*x(2) + (x(1) + x(3) + x(4))**2
!
!     subject to      x(1)         +   x(3)**2 +   x(4)**2     = 2
!                                    2*x(3)    + 4*x(4)       >= 0
!                             x(2)             +   x(4)**4     = 4
!                     x(1) >= 0,                         x(4) >= 0.
!
!     On exit:
!        nF  is the number of objective and constraint functions
!               (including linear and nonlinear)
!        n    is the number of variables.
!
!        (iGfun(k),jGvar(k)), k = 1,2,...,neG, define the coordinates
!             of the nonzero problem derivatives.
!             If (iGfun(k),jGvar(k)) = (i,j), G(k) is the ijth element
!             of the problem vector F(i), i = 0,1,2,...,nF,  with
!             objective function in position 0 and constraint functions
!             in positions  1  through  m.
!
!        (iAfun(k),jAvar(k),a(k)), k = 1,2,...,neA, are the coordinates
!             of the nonzero constant problem derivatives.
!
!     ==================================================================
      integer
     &     i, Obj
!     ------------------------------------------------------------------
      double precision     zero,         one ,         two
      parameter           (zero =0.0d+0, one  =1.0d+0, two    =2.0d+0 )
      double precision     four,         five,         plInfy
      parameter           (four =4.0d+0, five =5.0d+0, plInfy =1.0d+20)
!     ------------------------------------------------------------------
!     Give the problem a name.

      Prob   = 'testerA '

!     Assign the dimensions of the constraint Jacobian.

      nF     = 4
      Obj    = 1                ! Toy problem objective row
      ObjRow = 1                ! Could be 0

      n      = 4

!     Check that there is enough storage.

      Errors = 0
      if (nF     .gt. maxF ) Errors = 1
      if (n      .gt. maxn ) Errors = 1
      if (Errors .gt.     0) return

      neA        =  0
      neG        =  0

      neG        =  neG + 1
!     G(neG)     =  two*sum + three
      iGfun(neG) =  Obj
      jGvar(neG) =  1

      neG        =  neG + 1
!     G(neG)     =  five
      iGfun(neG) =  Obj
      jGvar(neG) =  2

      neG        =  neG + 1
!     G(neG)     =  two*sum
      iGfun(neG) =  Obj
      jGvar(neG) =  3

      neG        =  neG + 1
!     G(neG)     =  two*sum
      iGfun(neG) =  Obj
      jGvar(neG) =  4

!     Nonlinear constraints (derivatives by row)

      neG        =  neG + 1
!     G(neG)     =  one
      iGfun(neG) =  2
      jGvar(neG) =  1

      neG        =  neG + 1
!     G(neG)     =  two*x(3)
      iGfun(neG) =  2
      jGvar(neG) =  3

      neG        =  neG + 1
!     G(neG)     =  two*x(4)
      iGfun(neG) =  2
      jGvar(neG) =  4

      neG        =  neG + 1
!     G(neG)     =  two
      iGfun(neG) =  3
      jGvar(neG) =  3

      neG        =  neG + 1
!     G(neG)     =  four
      iGfun(neG) =  3
      jGvar(neG) =  4

      neG        =  neG + 1
!     G(neG)     =  one
      iGfun(neG) =  4
      jGvar(neG) =  2

      neG        =  neG + 1
!     G(neG)     =  four*x(4)**3
      iGfun(neG) =  4
      jGvar(neG) =  4

!     ----------------
!     Initial x.
!     ----------------
      ObjAdd = zero

      x(1)   =  one
      x(2)   =  one
      x(3)   =  one
      x(4)   =  one

      do i = 1, n
         xlow(i)   = -plInfy
         xupp(i)   =  plInfy
         xstate(i) =  0
      end do

      xlow(1) = zero
      xlow(2) = zero

!     The objective row is a free row.

      Flow(Obj) = -plInfy
      Fupp(Obj) =  plInfy

      Flow(2)   = two           ! Equality constraint
      Fupp(2)   = two

      Flow(3)   = zero
      Fupp(3)   = plInfy

      Flow(4)   = four          ! Equality constraint
      Fupp(4)   = four

      do i = 1, nF
         Fstate(i) = 0
         Fmul(i)   = zero
         F(i)      = zero
      end do

      end ! subroutine test1

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine test2
     &   ( Errors, maxF, maxn,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     Prob, nF, n,
     &     ObjAdd, ObjRow, xlow, xupp, Flow, Fupp,
     &     x, xstate, xmul, F, Fstate, Fmul,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     Errors, maxF, maxn, neA, lenA, neG, lenG, nF, n,
     &     ObjRow, lencw, leniw, lenrw, xstate(maxn), Fstate(maxF),
     &     iAfun(lenA), jAvar(lenA), iGfun(lenG), jGvar(lenG),
     &     iw(leniw)
      double precision
     &     ObjAdd,
     &     A(lenA), xlow(maxn), xupp(maxn), Flow(maxF), Fupp(maxF),
     &     x(maxn), F(maxF), xmul(maxn), Fmul(maxF), rw(lenrw)
      character
     &     Prob*8, cw(lencw)*8

!     ==================================================================
!     Test2 defines input data for the toy problem discussed in the
!     snoptA Users Guide.
!
!     Minimize      3*x(1) + 5*x(2) + (x(1) + x(3) + x(4))**2
!
!     subject to      x(1)         +   x(3)**2 +   x(4)**2     = 2
!                                    2*x(3)    + 4*x(4)       >= 0
!                             x(2)             +   x(4)**4     = 4
!                     x(1) >= 0,                         x(4) >= 0.
!
!     On exit:
!        nF  is the number of objective and constraint functions
!               (including linear and nonlinear)
!        n    is the number of variables.
!
!        (iGfun(k),jGvar(k)), k = 1,2,...,neG, define the coordinates
!             of the nonzero problem derivatives.
!             If (iGfun(k),jGvar(k)) = (i,j), G(k) is the ijth element
!             of the problem vector F(i), i = 0,1,2,...,nF,  with
!             objective function in position 0 and constraint functions
!             in positions  1  through  m.
!
!        (iAfun(k),jAvar(k),a(k)), k = 1,2,...,neA, are the coordinates
!             of the nonzero constant problem derivatives.
!
!     ==================================================================
      integer
     &     i, Obj
!     ------------------------------------------------------------------
      double precision     zero,         one ,         two
      parameter           (zero =0.0d+0, one  =1.0d+0, two    =2.0d+0 )
      double precision     four,         five,         plInfy
      parameter           (four =4.0d+0, five =5.0d+0, plInfy =1.0d+20)
!     ------------------------------------------------------------------
!     Give the problem a name.

      Prob   = 'testerA '

!     Assign the dimensions of the constraint Jacobian.

      nF     = 4
      Obj    = 1                ! Toy problem objective row
      ObjRow = 1                ! Could be 0
      n      = 4

!     Check that there is enough storage.

      Errors = 0
      if (nF     .gt. maxF ) Errors = 1
      if (n      .gt. maxn ) Errors = 1
      if (Errors .gt.     0) return

      neG        =  0

      neG        =  neG + 1
!     G(neG)     =  two*sum + three
      iGfun(neG) =  Obj
      jGvar(neG) =  1

      neG        =  neG + 1
!     G(neG)     =  two*sum
      iGfun(neG) =  Obj
      jGvar(neG) =  3

      neG        =  neG + 1
!     G(neG)     =  two*sum
      iGfun(neG) =  Obj
      jGvar(neG) =  4

!     Nonlinear constraints (derivatives by row)

      neG        =  neG + 1
!     G(neG)     =  two*x(3)
      iGfun(neG) =  2
      jGvar(neG) =  3

      neG        =  neG + 1
!     G(neG)     =  two*x(4)
      iGfun(neG) =  2
      jGvar(neG) =  4

      neG        =  neG + 1
!     G(neG)     =  four*x(4)**3
      iGfun(neG) =  4
      jGvar(neG) =  4

!     neG        = 7 derivatives in all

!     -------------------------------------------------------
!     Next we assign the list of constant derivative entries.
!     -------------------------------------------------------
      neA        =  0

      neA        =  neA + 1
      iAfun(neA) =  Obj
      jAvar(neA) =  2
      A(neA)     =  five

      neA        =  neA + 1
      iAfun(neA) =  2
      jAvar(neA) =  1
      A(neA)     =  one

      neA        =  neA + 1
      iAfun(neA) =  3
      jAvar(neA) =  3
      A(neA)     =  two

      neA        =  neA + 1
      iAfun(neA) =  3
      jAvar(neA) =  4
      A(neA)     =  four

      neA        =  neA + 1
      iAfun(neA) =  4
      jAvar(neA) =  2
      A(neA)     =  one

!     neA        =  5  derivatives in all

!     ----------------
!     Initial x.
!     ----------------
      ObjAdd = zero

      x(1)   =  one
      x(2)   =  one
      x(3)   =  one
      x(4)   =  one

      do i = 1, n
         xlow(i)   = -plInfy
         xupp(i)   =  plInfy
         xstate(i) =  0
      end do

      xlow(1) = zero
      xlow(2) = zero

!     The objective row is a free row.

      Flow(Obj) = -plInfy
      Fupp(Obj) =  plInfy

      Flow(2)   = two           ! Equality constraint
      Fupp(2)   = two

      Flow(3)   = zero
      Fupp(3)   = plInfy

      Flow(4)   = four          ! Equality constraint
      Fupp(4)   = four

      do i = 1, nF
         Fstate(i) = 0
         F(i)      = zero
         Fmul(i)   = zero
      end do

      end ! subroutine test2

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine test3
     &   ( Errors, maxF, maxn,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     Prob, nF, n,
     &     ObjAdd, ObjRow, xlow, xupp, Flow, Fupp,
     &     x, xstate, xmul, F, Fstate, Fmul,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     Errors, maxF, maxn, neA, lenA, neG, lenG, nF, n,
     &     ObjRow, lencw, leniw, lenrw, xstate(maxn), Fstate(maxF),
     &     iAfun(lenA), jAvar(lenA), iGfun(lenG), jGvar(lenG),
     &     iw(leniw)
      double precision
     &     ObjAdd,
     &     A(lenA), xlow(maxn), xupp(maxn), Flow(maxF), Fupp(maxF),
     &     x(maxn), F(maxF), xmul(maxn), Fmul(maxF), rw(lenrw)
      character
     &     Prob*8, cw(lencw)*8

!     ==================================================================
!     Test3 defines input data for the toy problem discussed in the
!     snoptA Users Guide.
!
!     Minimize      3*x(1) + 5*x(2) + (x(1) + x(3) + x(4))**2
!
!     subject to      x(1)         +   x(3)**2 +   x(4)**2     = 2
!                                    2*x(3)    + 4*x(4)       >= 0
!                             x(2)             +   x(4)**4     = 4
!                     x(1) >= 0,                         x(4) >= 0.
!
!     On exit:
!        nF  is the number of objective and constraint functions
!               (including linear and nonlinear)
!        n    is the number of variables.
!
!        (iGfun(k),jGvar(k)), k = 1,2,...,neG, define the coordinates
!             of the nonzero problem derivatives.
!             If (iGfun(k),jGvar(k)) = (i,j), G(k) is the ijth element
!             of the problem vector F(i), i = 0,1,2,...,nF,  with
!             objective function in position 0 and constraint functions
!             in positions  1  through  m.
!
!        (iAfun(k),jAvar(k),a(k)), k = 1,2,...,neA, are the coordinates
!             of the nonzero constant problem derivatives.
!
!     ==================================================================
      integer
     &     i, Obj
!     ------------------------------------------------------------------
      double precision     zero,         one ,         two
      parameter           (zero =0.0d+0, one  =1.0d+0, two    =2.0d+0 )
      double precision     four,         five,         plInfy
      parameter           (four =4.0d+0, five =5.0d+0, plInfy =1.0d+20)
!     ------------------------------------------------------------------
!     Give the problem a name.

      Prob   = 'testerA '

!     Assign the dimensions of the constraint Jacobian.

      nF     = 6                ! two empty rows
      Obj    = 1                ! Toy problem objective row
      ObjRow = 1                ! Could be 0
      n      = 5                ! One empty variable

!     Check that there is enough storage.

      Errors = 0
      if (nF     .gt. maxF ) Errors = 1
      if (n      .gt. maxn ) Errors = 1
      if (Errors .gt.     0) return

      neG        =  0

      neG        =  neG + 1
!     G(neG)     =  two*sum + three
      iGfun(neG) =  Obj
      jGvar(neG) =  1

      neG        =  neG + 1
!     G(neG)     =  two*sum
      iGfun(neG) =  Obj
      jGvar(neG) =  3

      neG        =  neG + 1
!     G(neG)     =  two*sum
      iGfun(neG) =  Obj
      jGvar(neG) =  4

!     Nonlinear constraints (derivatives by row)

      neG        =  neG + 1
!     G(neG)     =  two*x(3)
      iGfun(neG) =  2
      jGvar(neG) =  3

      neG        =  neG + 1
!     G(neG)     =  two*x(4)
      iGfun(neG) =  2
      jGvar(neG) =  4

      neG        =  neG + 1
!     G(neG)     =  four*x(4)**3
      iGfun(neG) =  4
      jGvar(neG) =  4

!     neG        = 7 derivatives in all

!     -------------------------------------------------------
!     Next we assign the list of constant derivative entries.
!     -------------------------------------------------------
      neA        =  0

      neA        =  neA + 1
      iAfun(neA) =  Obj
      jAvar(neA) =  2
      A(neA)     =  five

      neA        =  neA + 1
      iAfun(neA) =  2
      jAvar(neA) =  1
      A(neA)     =  one

      neA        =  neA + 1
      iAfun(neA) =  3
      jAvar(neA) =  3
      A(neA)     =  two

      neA        =  neA + 1
      iAfun(neA) =  3
      jAvar(neA) =  4
      A(neA)     =  four

      neA        =  neA + 1
      iAfun(neA) =  4
      jAvar(neA) =  2
      A(neA)     =  one

!     neA        =  5  derivatives in all

!     ----------------
!     Initial x.
!     ----------------
      ObjAdd = zero

      do i = 1, n
         x(i)      =  one
         xlow(i)   = -plInfy
         xupp(i)   =  plInfy
         xstate(i) =  0
      end do

      xlow(1) = zero
      xlow(2) = zero

!     The objective row is a free row.

      Flow(Obj) = -plInfy
      Fupp(Obj) =  plInfy

      Flow(2)   = two           ! Equality constraint
      Fupp(2)   = two

      Flow(3)   = zero
      Fupp(3)   = plInfy

      Flow(4)   = four          ! Equality constraint
      Fupp(4)   = four

      Flow(5)   = zero
      Fupp(5)   = plInfy

      Flow(6)   = zero
      Fupp(6)   = plInfy

      do i = 1, nF
         Fstate(i) = 0
         F(i)      = zero
         Fmul(i)   = zero
      end do

      end ! subroutine test3

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine test4
     &   ( Errors, Prob, maxF, maxn, nF, n, neG,
     &     iAfun, jAvar, lenA, neA, A,
     &     ObjAdd, ObjRow, xlow, xupp, Flow, Fupp, x, xstate, Fmul,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     Errors, maxF, maxn, neA, lenA, nF, n, neG,
     &     ObjRow, lencw, leniw, lenrw, xstate(maxn),
     &     iAfun(lenA), jAvar(lenA), iw(leniw)
      double precision
     &     ObjAdd,
     &     A(lenA), xlow(maxn), xupp(maxn), Flow(maxF), Fupp(maxF),
     &     x(maxn), Fmul(maxF), rw(lenrw)
      character*8
     &     Prob, cw(lencw)

!     ==================================================================
!     test4   defines input data for the Diet problem of Chvatal, 1983.
!
!
!           ( 110  205  160  160  420  260 )
!     A  =  (   4   32   13    8    4   14 )
!           (   2   12   54  285   22   80 )
!           (   3   24   13    9   20   19 ) ( = objective row c')

!     Errors      is 0 if there is enough storage, 1 otherwise.
!     nF          is the number of problem functions
!                 (objective and constraints).
!     n           is the number of variables.
!     xlow        holds the lower bounds on x.
!     xupp        holds the upper bounds on x.
!     Flow        holds the lower bounds on F = Ax.
!     Fupp        holds the upper bounds on F = Ax.

!     xstate(1:n) is a set of initial states for each x  (0,1,2,3,4,5).
!     x (1:n)     is a set of initial values for x.
!     Fmul(1:nF)  is a set of initial values for the dual variables.
!
!     ==================================================================
      integer
     &     i, j
!     ------------------------------------------------------------------
      double precision     plInfy
      parameter           (plInfy = 1.0d+20)
!     ------------------------------------------------------------------
!     Give the problem a name.

      Prob   = 'testerA '

!     Assign the dimensions of the constraint Jacobian.

      nF     = 4
      ObjRow = 4
      n      = 6

!     Check that there is enough storage.

      Errors = 0
      if (nF     .gt. maxF ) Errors = 1
      if (n      .gt. maxn ) Errors = 1
      if (Errors .gt.     0) return

      ObjAdd = 0.0d+0

      neG    = 0
      neA    = 0

!     Column 1

      neA        = neA + 1
      iAfun(neA) = 1
      jAvar(neA) = 1
      A(neA)     = 110.0d+0

      neA        = neA + 1
      iAfun(neA) = 2
      jAvar(neA) = 1
      A(neA)     = 4.0d+0

      neA        = neA + 1
      iAfun(neA) = 3
      jAvar(neA) = 1
      A(neA)     = 2.0d+0

      neA        = neA + 1
      iAfun(neA) = 4
      jAvar(neA) = 1
      A(neA)     = 3.0d+0

!     Column 2.

      neA        = neA + 1
      iAfun(neA) = 1
      jAvar(neA) = 2
      A(neA)     = 205.0d+0

      neA        = neA + 1
      iAfun(neA) = 2
      jAvar(neA) = 2
      A(neA)     = 32.0d+0

      neA        = neA + 1
      iAfun(neA) = 3
      jAvar(neA) = 2
      A(neA)     = 12.0d+0

      neA        = neA + 1
      iAfun(neA) = 4
      jAvar(neA) = 2
      A(neA)     = 24.0d+0

!     Column 3.

      neA        = neA + 1
      iAfun(neA) = 1
      jAvar(neA) = 3
      A(neA)     = 160.0d+0

      neA        = neA + 1
      iAfun(neA) = 2
      jAvar(neA) = 3
      A(neA)     = 13.0d+0

      neA        = neA + 1
      iAfun(neA) = 3
      jAvar(neA) = 3
      A(neA)     = 54.0d+0

      neA        = neA + 1
      iAfun(neA) = 4
      jAvar(neA) = 3
      A(neA)     = 13.0d+0

!     Column 4.

      neA        = neA + 1
      iAfun(neA) = 1
      jAvar(neA) = 4
      A(neA)     = 160.0d+0

      neA        = neA + 1
      iAfun(neA) = 2
      jAvar(neA) = 4
      A(neA)     = 8.0d+0

      neA        = neA + 1
      iAfun(neA) = 3
      jAvar(neA) = 4
      A(neA)     = 285.0d+0

      neA        = neA + 1
      iAfun(neA) = 4
      jAvar(neA) = 4
      A(neA)     = 9.0d+0

!     Column 5.

      neA        = neA + 1
      iAfun(neA) = 1
      jAvar(neA) = 5
      A(neA)     = 420.0d+0

      neA        = neA + 1
      iAfun(neA) = 2
      jAvar(neA) = 5
      A(neA)     = 4.0d+0

      neA        = neA + 1
      iAfun(neA) = 3
      jAvar(neA) = 5
      A(neA)     = 22.0d+0

      neA        = neA + 1
      iAfun(neA) = 4
      jAvar(neA) = 5
      A(neA)     = 20.0d+0

!     Column 6.

      neA        = neA + 1
      iAfun(neA) = 1
      jAvar(neA) = 6
      A(neA)     = 260.0d+0

      neA        = neA + 1
      iAfun(neA) = 2
      jAvar(neA) = 6
      A(neA)     = 14.0d+0

      neA        = neA + 1
      iAfun(neA) = 3
      jAvar(neA) = 6
      A(neA)     = 80.0d+0

      neA        = neA + 1
      iAfun(neA) = 4
      jAvar(neA) = 6
      A(neA)     = 19.0d+0

!     ------------------------------------------------------------------
!     Set the upper and lower bounds on the variables
!     ------------------------------------------------------------------
      do j = 1, n
         xlow(j) = 0.0d+0
      end do

      xupp(1) = 4.0d+0
      xupp(2) = 3.0d+0
      xupp(3) = 2.0d+0
      xupp(4) = 8.0d+0
      xupp(5) = 2.0d+0
      xupp(6) = 2.0d+0

!     ------------------------------------------------------------------
!     Set the upper and lower bounds on  Ax.
!     The objective row is free (i.e., infinite upper and lower bounds).
!     ------------------------------------------------------------------
      Flow( 1) =  2000.0d+0
      Flow( 2) =    55.0d+0
      Flow( 3) =   800.0d+0
      Flow( 4) = - plInfy

      do i = 1, nF
         Fupp(i) =  plInfy
         Fmul(i) =  0.0d+0
      end do

!     ----------------
!     Initialize  x.
!     ----------------
      do j = 1, n
         x(j) = 1.0d+0
      end do

      do j = 1, n
         xstate(j) = 0
      end do

      end ! subroutine test4

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine usrfun1
     &   ( Status, n, x,
     &     needF, nF, f,
     &     needG, lenG, G,
     &     cu, lencu, iu, leniu, ru, lenru )

      integer
     &     lencu, lenG, leniu, lenru, n, needF, needG, nF, Status,
     &     iu(leniu)
      double precision
     &     f(nF), G(lenG), x(n), ru(lenru)
      character
     &     cu(lencu)*8

!     ==================================================================
!     Computes the nonlinear objective and constraint terms for problem
!     featured in the SNOPT user's guide.
!     ==================================================================
      integer
     &     neG, Obj, Out
      double precision
     &     sum, x1, x3, x4
!     ------------------------------------------------------------------
      double precision     zero,         one ,         two
      parameter           (zero =0.0d+0, one  =1.0d+0, two  =2.0d+0)
      double precision     three,        four,         five
      parameter           (three=4.0d+0, four =4.0d+0, five =5.0d+0)
!     ------------------------------------------------------------------
      Out = 6                   ! Output unit number
      Obj = 1                   ! Objective row of F

!     --------------------------------------------
!     Print something on the first and last entry.
!     --------------------------------------------
      if (Status .eq. 1) then       ! First
         if (Out .gt. 0) write(Out, '(/a)') ' This is problem  testerA1'
      else  if (Status .ge. 2) then ! Last
         if (Out .gt. 0) write(Out, '(/a)') ' Finished problem testerA1'
         return
      end if

      x1  = x(1)
      x2  = x(2)
      x3  = x(3)
      x4  = x(4)
      sum = x1 + x3 + x4

      if (needF .gt. 0) then
         f(Obj) = three*x1 + five*x2 + sum**2
         f(2)   =       x1 +            x3**2 +      x4**2
         f(3)   =                   two*x3    + four*x4
         f(4)   =                 x2          +      x4**4
      end if

      neG = 0
      if (needG .gt. 0) then
         neG        =  neG + 1
         G(neG)     =  two*sum + three
!        iGfun(neG) =  Obj
!        jGvar(neG) =  1

         neG        =  neG + 1
         G(neG)     =  five
!        iGfun(neG) =  Obj
!        jGvar(neG) =  2

         neG        =  neG + 1
         G(neG)     =  two*sum
!        iGfun(neG) =  Obj
!        jGvar(neG) =  3

         neG        =  neG + 1
         G(neG)     =  two*sum
!        iGfun(neG) =  Obj
!        jGvar(neG) =  4

!        Nonlinear constraints (derivatives by row)

         neG        =  neG + 1
         G(neG)     =  one
!        iGfun(neG) =  2
!        jGvar(neG) =  1

         neG        =  neG + 1
         G(neG)     =  two*x(3)
!        iGfun(neG) =  2
!        jGvar(neG) =  3

         neG        =  neG + 1
         G(neG)     =  two*x(4)
!        iGfun(neG) =  2
!        jGvar(neG) =  4

         neG        =  neG + 1
         G(neG)     =  two
!        iGfun(neG) =  3
!        jGvar(neG) =  3

         neG        =  neG + 1
         G(neG)     =  four
!        iGfun(neG) =  3
!        jGvar(neG) =  4

         neG        =  neG + 1
         G(neG)     =  one
!        iGfun(neG) =  4
!        jGvar(neG) =  2

         neG        =  neG + 1
         G(neG)     =  four*x(4)**3
!        iGfun(neG) =  4
!        jGvar(neG) =  4

      end if

      end ! subroutine usrfun1

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine usrfun2
     &   ( Status, n, x,
     &     needF, nF, f,
     &     needG, lenG, G,
     &     cu, lencu, iu, leniu, ru, lenru )

      integer
     &     lencu, lenG, leniu, lenru, n, needF, needG, nF, Status,
     &     iu(leniu)
      double precision
     &     f(nF), G(lenG), x(n), ru(lenru)
      character
     &     cu(lencu)*8

!     ==================================================================
!     Computes the nonlinear objective and constraint terms for the toy
!     problem featured in the SNOPT user's guide.
!     ==================================================================
      integer
     &     neG, Obj, Out
      double precision
     &     sum, x1, x3, x4
!     ------------------------------------------------------------------
      Out = 6                   ! Output unit number
      Obj = 1                   ! Objective row of F

!     --------------------------------------------
!     Print something on the first and last entry.
!     --------------------------------------------
      if (Status .eq. 1) then       ! First
         if (Out .gt. 0) write(Out, '(/a)') ' This is problem  testerA'
      else  if (Status .ge. 2) then ! Last
         if (Out .gt. 0) write(Out, '(/a)') ' Finished problem testerA'
         return
      end if

      x1  = x(1)
      x3  = x(3)
      x4  = x(4)
      sum = x1 + x3 + x4

      if (needF .gt. 0) then
         f(Obj) =   3.0d0*x1 + sum**2
         f(2)   =   x3**2 + x4**2
!!!      f(3)   =                    ! Linear constraint omitted!
         f(4)   =   x4**4
      end if

      neG = 0
      if (needG .gt. 0) then
         neG        =  neG + 1
         G(neG)     =  2.0d0*sum + 3.0d0
!!!      iGfun(neG) =  Obj           ! Not used, but included for clarity!
!!!      jGvar(neG) =  1             ! Not used

         neG        =  neG + 1
         G(neG)     =  2.0d0*sum
!        iGfun(neG) =  Obj
!        jGvar(neG) =  3

         neG        =  neG + 1
         G(neG)     =  2.0d0*sum
!        iGfun(neG) =  Obj
!        jGvar(neG) =  4

         neG        =  neG + 1
         G(neG)     =  2.0d0*x3
!        iGfun(neG) =  2
!        jGvar(neG) =  3

         neG        =  neG + 1
         G(neG)     =  2.0d0*x4
!        iGfun(neG) =  2
!        jGvar(neG) =  4

         neG        =  neG + 1
         G(neG)     =  4.0d0*x4**3
!        iGfun(neG) =  4
!        jGvar(neG) =  4
      end if

      end ! subroutine usrfun2

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine usrfun3
     &   ( Status, n, x,
     &     needF, nF, f,
     &     needG, lenG, G,
     &     cu, lencu, iu, leniu, ru, lenru )

      integer
     &     lencu, lenG, leniu, lenru, n, needF, needG, nF, Status,
     &     iu(leniu)
      double precision
     &     f(nF), G(lenG), x(n), ru(lenru)
      character
     &     cu(lencu)*8

!     ==================================================================
!     Computes the nonlinear objective and constraint terms for the toy
!     problem featured in the SNOPT user's guide.
!     ==================================================================
      integer
     &     neG, Obj, Out
      double precision
     &     sum, x1, x3, x4
!     ------------------------------------------------------------------
      Out = 6                   ! Output unit number
      Obj = 1                   ! Objective row of F

!     --------------------------------------------
!     Print something on the first and last entry.
!     --------------------------------------------
      if (Status .eq. 1) then       ! First
         if (Out .gt. 0) write(Out, '(/a)') ' This is problem  testerA'
      else  if (Status .ge. 2) then ! Last
         if (Out .gt. 0) write(Out, '(/a)') ' Finished problem testerA'
         return
      end if

      x1  = x(1)
      x3  = x(3)
      x4  = x(4)
      sum = x1 + x3 + x4

      if (needF .gt. 0) then
         f(Obj) =   3.0d0*x1 + sum**2
         f(2)   =   x3**2 + x4**2
!!!      f(3)   =                    ! Linear constraint omitted!
         f(4)   =   x4**4
         f(5)   =   0.0d0            ! Keep valgrind happy
         f(6)   =   0.0d0
      end if

      neG = 0
      if (needG .gt. 0) then
         neG        =  neG + 1
         G(neG)     =  2.0d0*sum + 3.0d0
!!!      iGfun(neG) =  Obj           ! Not used, but included for clarity!
!!!      jGvar(neG) =  1             ! Not used

         neG        =  neG + 1
         G(neG)     =  2.0d0*sum
!        iGfun(neG) =  Obj
!        jGvar(neG) =  3

         neG        =  neG + 1
         G(neG)     =  2.0d0*sum
!        iGfun(neG) =  Obj
!        jGvar(neG) =  4

         neG        =  neG + 1
         G(neG)     =  2.0d0*x3
!        iGfun(neG) =  2
!        jGvar(neG) =  3

         neG        =  neG + 1
         G(neG)     =  2.0d0*x4
!        iGfun(neG) =  2
!        jGvar(neG) =  4

         neG        =  neG + 1
         G(neG)     =  4.0d0*x4**3
!        iGfun(neG) =  4
!        jGvar(neG) =  4
      end if

      end ! subroutine usrfun2


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine usrfun4
     &   ( Status, n, x,
     &     needF, nF, F,
     &     needG, lenG, G,
     &     cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none
      integer
     &     Status, needF, needG, nF, n, lenG, lencu, leniu, lenru,
     &     iu(leniu)
      double precision
     &     F(nF), G(lenG), x(n), ru(lenru)
      character*8
     &     cu(lencu)

!     ==================================================================
!     Dummy  objective and constraint function for the Diet problem
!     of Chvatal, 1983.
!
!           ( 110  205  160  160  420  260 )
!     A  =  (   4   32   13    8    4   14 )
!           (   2   12   54  285   22   80 )
!           (   3   24   13    9   20   19 ) ( = objective row c')
!
!     ==================================================================
!     Relax, A*x  is computed by snOptA from A.

      end ! subroutine usrfun4
