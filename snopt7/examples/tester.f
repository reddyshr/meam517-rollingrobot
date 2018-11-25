*     ------------------------------------------------------------------
*     File tester.f
*     Illustrates using snOptB on a linear program with MPS data
*
*     16 May 1998: First   version.
*     02 Jan 2003: Current version.
*     ------------------------------------------------------------------
      program
     &     main
      implicit
     &     none
      integer
     &     maxm,        maxn,          maxne
      parameter
     &   ( maxm = 1000, maxn   = 1000, maxne  = 3000 )
      character*8
     &     PrbNms(5), Names(maxm+maxn)
      integer
     &     indA(maxne), locA(maxn+1), hs(maxn+maxm)
      double precision
     &     Acol(maxne), bl(maxn+maxm), bu(maxn+maxm), x(maxn+maxm),
     &     pi(maxm), rc(maxn+maxm)

*     ------------------------------------------------------------------
*     SNOPT workspace

      integer
     &     lencw,         leniw,         lenrw
      parameter
     &   ( lencw =   500, leniw = 50000, lenrw = 50000)
      double precision
     &     rw(lenrw)
      integer
     &     iw(leniw)
      character*8
     &     cw(lencw)
*     ------------------------------------------------------------------
      character*20
     &     lfile
      external
     &     dummyf, dummyc
      integer
     &     Errors, iMPS, INFO, iObj, iPrint, iSpecs, iSumm, m, mincw,
     &     miniw, minrw, n, ne, nInf, nName, nnCon, nnJac, nnObj, nOut,
     &     nS
      double precision
     &     Obj, ObjAdd, sInf

*     ------------------------------------------------------------------
*     Specify some of the SNOPT files.
*     iSpecs  is the Specs file   (0 if none).
*     iPrint  is the Print file   (0 if none).
*     iSumm   is the Summary file (0 if none).
*
*     nOut    is an output file used here by t7etamacro.

      iSpecs =  4
      iPrint =  9
      iSumm  =  6
      nOut   =  6

*     Open the Specs and print files.

      lfile = 'tester.spc'
      open( iSpecs, file=lfile, status='OLD',     err=800 )

      lfile = 'tester.out'
      open( iPrint, file=lfile, status='UNKNOWN', err=800 )

*     ------------------------------------------------------------------
*     First,  snInit MUST be called to initialize optional parameters
*     to their default values.
*     ------------------------------------------------------------------
      call snInit
     &   ( iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw )

*     ------------------------------------------------------------------
*     Read a Specs file (Optional).
*     ------------------------------------------------------------------
      call snSpec
     &     ( iSpecs, INFO, cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .ne. 101  .and.  INFO .ne. 107) then
         go to 990
      end if

*     ------------------------------------------------------------------
*     Set up the data structure for the linear constraints.
*     MPSinp needs to know the number of nonlinear variables, etc.
*     The following calls fetch values set in the SPECS file.
*     Optionally, these values can be set in-line.
*     ------------------------------------------------------------------
      Errors = 0

      call sngeti
     &   ( 'Nonlinear constraints',         nnCon, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call sngeti
     &   ( 'Nonlinear Jacobian  variables', nnJac, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call sngeti
     &   ( 'Nonlinear Objective variables', nnObj, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call sngeti
     &   ( 'MPS file',                       iMPS, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )


*     The problem name is not needed---it is set by MPSinp.
*     PrbNms(1) = '        '
*     Specify the OBJECTIVE, RHS, RANGES and BOUNDS to be selected
*     from the MPS file.  Blank names mean "select the first one".

      PrbNms(2) = '        '    ! OBJECTIVE name
      PrbNms(3) = '        '    ! RHS       name
      PrbNms(4) = '        '    ! RANGES    name
      PrbNms(5) = '        '    ! BOUNDS    name

      lfile = 'tester.mps'
      open( iMPS, file=lfile, status='OLD', err=800 )

      call MPSinp
     &   ( iMPS, maxm, maxn, maxne,
     &     nnCon, nnJac, nnObj,
     &     m, n, ne,
     &     iObj, ObjAdd, PrbNms,
     &     Acol, indA, locA, bl, bu, Names,
     &     hs, x, pi,
     &     INFO, mincw, miniw, minrw, nS,
     &     cw, lencw, iw, leniw, rw, lenrw )
      close( iMPS )
      if (INFO .ne. 103) go to 990
      nName = m + n

*     ------------------------------------------------------------------
*     Go for it, using a Cold start.
*     hs     need not be set if a basis file is to be input.
*            Otherwise, each hs(1:n) should be 0, 1, 2, 3, 4, or 5.
*            The values are used by the Crash procedure s2crsh
*            to choose an initial basis B.
*            If hs(j) = 0 or 1, column j is eligible for B.
*            If hs(j) = 2, column j is initially superbasic (not in B).
*            If hs(j) = 3, column j is eligible for B and is given
*                          preference over columns with hs(j) = 0 or 1.
*            If hs(j) = 4 or 5, column j is initially nonbasic.
*     ------------------------------------------------------------------
      call snOpt
     &   ( 'Cold', m, n, ne, nName,
     &     nnCon, nnObj, nnJac,
     &     iObj, ObjAdd, PrbNms(1),
     &     dummyf, dummyc,
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
      write(nOut, *) 'tester finished.'
      write(nOut, *) 'Input  errors =', Errors
      write(nOut, *) 'snOptB INFO   =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      if (iObj .gt. 0) then
         write(nOut, *) 'Obj        =', ObjAdd + x(n+iObj) + Obj
      else
         write(nOut, *) 'Obj        =', ObjAdd + Obj
      end if
      if (INFO .ge. 30) go to 910
      stop

*     ------------------------------------------------------------------
*     Error exit.
*     ------------------------------------------------------------------
  800 write(nOut, 4000) 'Error while opening file', lfile
      stop

  910 write(nOut, *) ' '
      write(nOut, *) 'STOPPING because of error condition'

  990 stop

 4000 format(/  a, 2x, a  )

      end ! program main

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine dummyf
     &   ( mode, n, x, f, g, nState, cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none
      integer
     &     mode, n, nState, lencu, leniu, lenru, iu(leniu)
      double precision
     &     f, x(n), g(n), ru(lenru)
      character*8
     &     cu(lencu)

*     ==================================================================
*     No nonlinear objective.
*     ==================================================================
      write(*,*) '***** dummyf WAS CALLED ****'

      end ! subroutine dummyf

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine dummyc
     &   ( mode, nnCon, nnJac, neJac, x, fCon, gCon,
     &     nState, cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none
      integer
     &     mode, nnCon, nnJac, neJac, nState, lencu, leniu, lenru,
     &     iu(leniu)
      double precision
     &     x(nnJac), fCon(nnCon), gCon(neJac), ru(lenru)
      character*8
     &     cu(lencu)

*     ==================================================================
*     No nonlinear constraints.
*     ==================================================================
      write(*,*) '**** dummyc WAS CALLED ****'

      end ! subroutine dummyc

