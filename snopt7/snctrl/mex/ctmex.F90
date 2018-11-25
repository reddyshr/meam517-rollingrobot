!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File: ctmex.f90
!
! 10 Feb 2010: First version.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module ctmex
  use snctrl
  implicit none

  ! SNOPT/SNCTRL Workspace
  integer           :: leniw = 10000, lenrw = 50000, lencw = 500
  integer,          allocatable :: iw(:), iw0(:)
  double precision, allocatable :: rw(:), rw0(:)
  character*8,      allocatable :: cw(:), cw0(:)

  ! Problem structure
  type(ctProb) :: prob

  ! SNOPT mex variables
  logical              :: firstCall = .true.,  &
                          memCall   = .false., &
                          printOpen = .false., &
                          summOpen  = .false., &
                          screenON  = .false.

  integer              :: callType = 0

  mwPointer            :: objHandle, conHandle

  integer, parameter   :: iPrint     = 9, iSpecs   = 4, iSumm    = 55, &
                          systemCall = 0, userCall = 1
  integer, parameter   :: sncSolve    = 1,  &
                          sncSetXX    = 2,  &
                          sncSetIX    = 3,  &
                          sncSetRX    = 4,  &
                          sncGetXX    = 5,  &
                          sncGetCX    = 6,  &
                          sncGetIX    = 7,  &
                          sncGetRX    = 8,  &
                          sncSpecs    = 9,  &
                          sncOpenP    = 10, &
                          sncOpenS    = 11, &
                          sncClosP    = 12, &
                          sncClosS    = 13, &
                          sncSetWork  = 14, &
                          sncscrnON   = 15, &
                          sncscrnOff  = 16, &
                          sncEnd      = 999

contains

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine checkCol ( pm, n, name )
    mwPointer     :: pm
    integer       :: n
    character*(*) :: name
    !---------------------------------------------------------------------------
    ! Check column dimension of pm is equal to n.
    !---------------------------------------------------------------------------
    character*80 :: str
    mwSize       :: m, mxGetN

    m = mxGetN(pm)
    if ( m /= n ) then
       write(str,100) name, m, n
       call mexErrMsgTxt ( str )
    end if

    return

100 format ( a, ' has incorrect column dimension ', i5, &
                '.  Should be length ', i5 )

  end subroutine checkCol

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine checkRow ( pm, n, name )
    character*(*) :: name
    mwPointer     :: pm
    integer       :: n
    !---------------------------------------------------------------------------
    ! Check row dimension of pm is equal to n.
    !---------------------------------------------------------------------------
    character*80 :: str
    mwSize       :: m, mxGetM

    m = mxGetM(pm)
    if ( m /= n ) then
       write(str,100) name, m, n
       call mexErrMsgTxt ( str )
    end if

    return

100 format ( a, ' has incorrect row dimension ', i5, &
                '.  Should be length ', i5 )

  end subroutine checkRow

end module ctmex

