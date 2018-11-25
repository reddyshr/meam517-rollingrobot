!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File: ctmexModule.f90
! Module containing variables for the Matlab gateway routine.
!
! 10 Feb 2010: First version.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module ctmexModule
  use precision,   only : ip, rp
  use ctusrModule, only : ctProb
  implicit none

  ! SNOPT/SNCTRL Workspace
  integer(ip), parameter :: lencu = 1, &
                            leniu = 1, &
                            lenru = 1, &
                            lencw = 500, &
                            leniw = 10000, &
                            lenrw = 50000

  integer(ip)  :: iu(leniu), iw(leniw)
  real(rp)     :: ru(lenru), rw(lenrw)
  character(8) :: cu(lencu), cw(lencw)

  ! M-file names
  character(20) :: varbds, odecon, algcon

  ! Problem structure
  type(ctProb) :: prob

end module ctmexModule

