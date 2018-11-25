!------------------------------------------------------------------
! File sqmain3.f90
! Examples of calls to subroutine SQOPT.
!
! 14 Jul 2005: Date of sqmain2.f.
! 23 Apr 2009: sqmain3.f90 is an f90 version of sqmain2.f.
!              The associated subroutines have to be in a
!              separate file, sqmain3Module.f90.
!------------------------------------------------------------------
program sqmain3

  use sqmain3Module, only : sqdata, userHx, sqdata2

  implicit none
  integer, parameter   :: maxm   = 1000, &
                          maxn   = 1500, &
                          maxne  = 3000, &
                          maxnb  = maxm + maxn
  character            :: Prob*8, Names(maxnb)*8
  integer              :: indA(maxne), hEtype(maxnb), hs(maxnb), locA(maxn+1)
  real(8)              :: Acol(maxne), bl(maxnb), bu(maxnb), cObj(maxn), &
                          x(maxnb), pi(maxm), rc(maxnb)

  !-------------SQOPT workspace-------------------------------------
  integer, parameter   :: lencw =     500, &
                          leniw =  100000, &
                          lenrw =  100000
  character            :: cw(lencw)*8
  integer              :: iw(leniw)
  real(8)              :: rw(lenrw)
  !-----------------------------------------------------------------

  logical              :: byname
  integer              :: Errors, INFO, iObj, iPrint, iSumm, iSpecs, &
                          itnlim, iP, iS, j, lencObj, m, maxS,       &
                          mincw, miniw, minrw, n, ncolH, ne, nInf,   &
                          nName, nOut, nS
  real(8)              :: Obj, ObjAdd, sInf
  character            :: lfile*20

  !-----------------------------------------------------------------
  ! Specify some of the SQOPT files.
  ! iSpecs  is the Specs file   (0 if none).
  ! iPrint  is the Print file   (0 if none).
  ! iSumm   is the Summary file (0 if none).
  !
  ! nOut    is an output file used here by main.
  !-----------------------------------------------------------------
  iSpecs =  4
  iPrint =  9
  iSumm  =  6
  nOut   =  6

  byname = .true.

  if ( byname ) then ! Unix and DOS systems.
     lfile = 'sqmain3.spc'
     open( iSpecs, file=lfile, status='OLD',     err=800 )
     lfile = 'sqmain3.out'
     open( iPrint, file=lfile, status='UNKNOWN', err=800 )
  end if

  !-----------------------------------------------------------------
  ! First,  sqInit MUST be called to initialize optional parameters
  ! to their default values.
  !-----------------------------------------------------------------
  call sqInit( iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw )

  !-----------------------------------------------------------------
  ! Read a Specs file (Optional).
  !-----------------------------------------------------------------
  call sqSpec( iSpecs, INFO, cw, lencw, iw, leniw, rw, lenrw )

  if (INFO .ne. 101  .and.  INFO .ne. 107) then
     go to 910
  end if

  !-----------------------------------------------------------------
  ! 1. Solve a QP problem with an explicit linear objective.
  ! (1) Compute l, u, and A so that the constraints are ranges of the
  !     form  l <= Ax <= u.
  !     Store l and u in bl(n+1:n+m) and bu(n+1:n+m).
  !
  ! (2) Set up the constants ObjAdd and c so that the explicit
  !     objective is
  !         ObjAdd + c'*x + half*x'*H*x
  !-----------------------------------------------------------------
  call sqdata( maxm, maxn, maxne,               &
               m, n, ne, nName, lencObj, ncolH, &
               iObj, ObjAdd, Prob,              &
               Acol, indA, locA, bl, bu, cObj,  &
               Names, hEtype, hs, x )

  !-----------------------------------------------------------------
  ! Specify options not set in the Specs file.
  ! iP and iS refer to the Print and Summary file respectively.
  ! Setting them to 0 suppresses printing.
  !-----------------------------------------------------------------
  Errors = 0

  maxS   = ncolH + 1
  itnlim = 200
  iP     =  0
  iS     =  0
  call sqseti( 'Superbasics Limit', maxS  , iP, iS, Errors, &
               cw, lencw, iw, leniw, rw, lenrw )
  call sqseti( 'Iterations',        itnlim, iP, iS, Errors, &
               cw, lencw, iw, leniw, rw, lenrw )

  !-----------------------------------------------------------------
  ! Solve the QP using a Cold start.
  ! hs     need not be set if a basis file is to be input.
  !        Otherwise, each hs(1:n) should be 0, 1, 2, 3, 4, or 5.
  !        The values are used by the Crash procedure s2crsh
  !        to choose an initial basis B.
  !        If hs(j) = 0 or 1, column j is eligible for B.
  !        If hs(j) = 2, column j is initially superbasic (not in B).
  !        If hs(j) = 3, column j is eligible for B and is given
  !                      preference over columns with hs(j) = 0 or 1.
  !        If hs(j) = 4 or 5, column j is initially nonbasic.
  !-----------------------------------------------------------------
  call sqopt ( 'Cold', userHx, m, n, ne, nName, lencObj, ncolH, &
               iObj, ObjAdd, Prob,                              &
               Acol, indA, locA, bl, bu, cObj, Names,           &
               hEtype, hs, x, pi, rc,                           &
               INFO, mincw, miniw, minrw,                       &
               nS, nInf, sInf, Obj,                             &
               cw, lencw, iw, leniw, rw, lenrw,                 &
               cw, lencw, iw, leniw, rw, lenrw )

  if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
     go to 910
  end if

  write(nOut, *) ' '
  write(nOut, *) 'sqOpt finished.'
  write(nOut, *) 'Input errors  =', Errors
  write(nOut, *) 'sqOpt INFO    =', INFO
  write(nOut, *) 'nInf          =', nInf
  write(nOut, *) 'sInf          =', sInf
  write(nOut, *) 'Obj           =', ObjAdd + Obj
  if (INFO .gt. 30) go to 910

  !-----------------------------------------------------------------
  ! 2. Solve the same problem but with the objective row as part of
  ! the constraint matrix A.
  ! Use a Cold Start because the dimensions are different.
  !-----------------------------------------------------------------
  Errors = 0                ! Reset count of input errors.

  call sqset ( 'Solution    No ',        iP,    iS, Errors, &
               cw, lencw, iw, leniw, rw, lenrw )

  call sqdata2( maxm, maxn, maxne,               &
                m, n, ne, nName, lencObj, ncolH, &
                iObj, ObjAdd, Prob,              &
                Acol, indA, locA, bl, bu, cObj,  &
                Names, hEtype, hs, x )

  call sqopt  ( 'Cold', userHx, m,                     &
                n, ne, nName, lencObj, ncolH,          &
                iObj, ObjAdd, Prob,                    &
                Acol, indA, locA, bl, bu, cObj, Names, &
                hEtype, hs, x, pi, rc,                 &
                INFO, mincw, miniw, minrw,             &
                nS, nInf, sInf, Obj,                   &
                cw, lencw, iw, leniw, rw, lenrw,       &
                cw, lencw, iw, leniw, rw, lenrw )

  if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
     go to 910
  end if

  write(nOut, *) ' '
  write(nOut, *) ' '
  write(nOut, *) 'sqOpt finished.'
  write(nOut, *) 'Input errors  =', Errors
  write(nOut, *) 'sqOpt INFO    =', INFO
  write(nOut, *) 'nInf          =', nInf
  write(nOut, *) 'sInf          =', sInf
  write(nOut, *) 'Obj           =', ObjAdd + Obj + x(n+iObj)
  if (INFO .gt. 30) go to 910

  !-----------------------------------------------------------------
  ! 3. Alter some options and call sqopt again, with a Warm start.
  ! The following illustrates the use of sqset, sqseti and sqsetr
  ! to set specific options.  We could ensure that all unspecified
  ! options take default values by first calling
  ! sqset ( 'Defaults', ... ).
  ! Beware that certain parameters would then need to be redefined.
  !-----------------------------------------------------------------
  write(nOut, *) ' '
  write(nOut, *) 'Alter options and test Warm start:'

  Errors = 0
  itnlim = 500
  call sqset ( 'Defaults',           iPrint, iSumm, Errors, &
               cw, lencw, iw, leniw, rw, lenrw )
  call sqseti( 'Print frequency', 1, iPrint, iSumm, Errors, &
               cw, lencw, iw, leniw, rw, lenrw )
  call sqset ( ' ',                  iPrint, iSumm, Errors, &
               cw, lencw, iw, leniw, rw, lenrw )
  call sqset ( 'Scale option 0',     iPrint, iSumm, Errors, &
               cw, lencw, iw, leniw, rw, lenrw )
  call sqseti( 'Iterations', itnlim, iPrint, iSumm, Errors, &
               cw, lencw, iw, leniw, rw, lenrw )

  if (Errors .gt. 0) then
     write(nOut, *) 'NOTE: Some of the options were not recognized'
  end if


  !-----------------------------------------------------------------
  ! 3. Alter some options and call sqopt again, with a Warm start.
  ! The following illustrates the use of sqset, sqseti and sqsetr
  ! to set specific options.  We could ensure that all unspecified
  ! options take default values by first calling
  ! sqset ( 'Defaults', ... ).
  ! Beware that certain parameters would then need to be redefined.
  !-----------------------------------------------------------------
  write(nOut, *) ' '
  write(nOut, *) 'Alter options and test Warm start:'

  Errors = 0
  itnlim = 500
  call sqset ( 'Defaults',           iPrint, iSumm, Errors, &
               cw, lencw, iw, leniw, rw, lenrw )
  call sqseti( 'Print frequency', 1, iPrint, iSumm, Errors, &
               cw, lencw, iw, leniw, rw, lenrw )
  call sqset ( ' ',                  iPrint, iSumm, Errors, &
               cw, lencw, iw, leniw, rw, lenrw )
  call sqset ( 'Scale option 0',     iPrint, iSumm, Errors, &
               cw, lencw, iw, leniw, rw, lenrw )
  call sqseti( 'Iterations', itnlim, iPrint, iSumm, Errors, &
               cw, lencw, iw, leniw, rw, lenrw )

  if (Errors .gt. 0) then
     write(nOut, *) 'NOTE: Some of the options were not recognized'
  end if

  !-----------------------------------------------------------------
  ! Test the Warm start.
  ! hs(*) specifies a complete basis from the previous call.
  ! A Warm start uses hs(*) directly, without calling Crash.
  !
  ! Warm starts are normally used after sqopt has solved a
  ! problem with the SAME DIMENSIONS but perhaps altered data.
  ! Here we have not altered the data, so very few iterations
  ! should be required.
  !-----------------------------------------------------------------
  call sqopt  ( 'Warm', userHx, m,                     &
                n, ne, nName, lencObj, ncolH,          &
                iObj, ObjAdd, Prob,                    &
                Acol, indA, locA, bl, bu, cObj, Names, &
                hEtype, hs, x, pi, rc,                 &
                INFO, mincw, miniw, minrw,             &
                nS, nInf, sInf, Obj,                   &
                cw, lencw, iw, leniw, rw, lenrw,       &
                cw, lencw, iw, leniw, rw, lenrw )

  if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
     go to 910
  end if

  write(nOut, *) ' '
  write(nOut, *) 'sqOpt finished.'
  write(nOut, *) 'Input errors  =', Errors
  write(nOut, *) 'sqOpt INFO    =', INFO
  write(nOut, *) 'nInf          =', nInf
  write(nOut, *) 'sInf          =', sInf
  write(nOut, *) 'Obj           =', ObjAdd + Obj + x(n+iObj)
  if (INFO .gt. 30) go to 910

  !------------------------------------------------------------------
  ! 4. Call sqopt with the quasi-Newton solver and a hot start.
  !
  ! The following illustrates the use of sqset, sqseti and sqsetr
  ! to set specific options.  We could ensure that all unspecified
  ! options take default values by first calling
  ! sqset ( 'Defaults', ... ).
  ! Beware that certain parameters would then need to be redefined.
  !------------------------------------------------------------------
  write(nOut, *) ' '
  write(nOut, *) 'Calling sqopt with a hot start:'
  INFO = 0

  call sqset( 'QPsolver QN',     iPrint, iSumm, INFO, &
              cw, lencw, iw, leniw, rw, lenrw )

  call sqopt( 'Hot H', userHx, m,                    &
              n, ne, nName, lencObj, ncolH,          &
              iObj, ObjAdd, Prob,                    &
              Acol, indA, locA, bl, bu, cObj, Names, &
              hEtype, hs, x, pi, rc,                 &
              INFO, mincw, miniw, minrw,             &
              nS, nInf, sInf, Obj,                   &
              cw, lencw, iw, leniw, rw, lenrw,       &
              cw, lencw, iw, leniw, rw, lenrw )

  if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
     go to 910
  end if

  write(nOut, *) ' '
  write(nOut, *) 'sqOpt finished.'
  write(nOut, *) 'Input errors  =', Errors
  write(nOut, *) 'sqOpt INFO    =', INFO
  write(nOut, *) 'nInf          =', nInf
  write(nOut, *) 'sInf          =', sInf
  write(nOut, *) 'Obj           =', ObjAdd + Obj + x(n+iObj)

  stop

  !-----------------------------------------------------------------
  ! Error exit.
  !-----------------------------------------------------------------
  800 write(nOut, 4000) 'Error while opening file', lfile
      stop

  910 write(nOut, *) ' '
      write(nOut, *) 'STOPPING because of error condition'
      stop

 4000 format(/ a, 2x, a )

end program sqmain3
