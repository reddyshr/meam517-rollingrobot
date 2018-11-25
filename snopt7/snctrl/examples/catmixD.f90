!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File catmixD.f90
!
!  This problem has an algebraic objective.
!  The original problem:
!      min            y1 + y2 - 1
!      subject to:
!                     y1' = u1 * (10*y2 - y1)
!                     y2' = u1 * (y1 - 10*y2) - (1-u1)*y2
!
!                     y1(0) = 1
!                     y2(0) = 0
!                     0 <= u1 <= 1
!
!  The problem input to SNCTRL with an added control to represent
!  the algebraic objective:
!      min            u2
!      subject to:
!                     y1' = u1 * (10*y2 - y1)
!                     y2' = u1 * (y1 - 10*y2) - (1-u1)*y2
!
!                     y1(0) = 1
!                     y2(0) = 0
!                     0 <= u1 <= 1
!
!                     u2 = y1 + y2 - 1
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

program catmixD
  use SNCTRL
  implicit none

  !=====================================================================

  ! SNOPT variables and workspaces
  character(20) :: lfile
  real(rp)      :: sInf
  integer(ip)   :: Start, INFO, iPrint, iSpecs, nOut, &
                   iSumm, mincw, miniw, minrw, nS, nInf, &
                   lencu, leniu, lenru, &
                   lencw, leniw, lenrw
  parameter ( lenrw = 50000, leniw = 50000, lencw = 500, &
              lenru = 1,     leniu = 1,     lencu = 1)
  integer(ip)  :: iu(leniu), iw(leniw)
  real(rp)     :: ru(lenru), rw(lenrw)
  character(8) :: cu(lencu), cw(lencw)

  external     :: odecon, algcon, varbds
  type(ctProb) :: prob

  !=====================================================================

  iSpecs =  4  ! equivalenced to .spc
  iPrint =  9  ! equivalenced to .out
  iSumm  =  6  ! summary file goes to standard output...
  nOut   =  6  ! ... as do messages from this program.

  ! Open the Specs and print files.
  lfile = 'catmixD.spc'
  open( iSpecs, file=lfile, status='OLD',     err=800 )
  lfile = 'catmixD.out'
  open( iPrint, file=lfile, status='UNKNOWN', err=800 )


  ! Set options to default values.
  call sncInit( iPrint, iSumm, prob, cw, lencw, iw, leniw, rw, lenrw )


  ! Read a Specs file.
  call sncSpec( iSpecs, INFO, cw, lencw, iw, leniw, rw, lenrw )

  if (INFO /= 101  .AND.  INFO /= 107) go to 910


  ! Setup the problem.
  prob%probname(1:8) = ' catmixD'
  prob%nY   = 2     ! number of states
  prob%nU   = 2     ! number of controls
  prob%nP   = 0     ! number of parameters
  prob%nC   = 1     ! number of algebraic constraints
  prob%nPhs = 1     ! number of phases

  ! Objective
  allocate ( prob%objL(1) )
  prob%objL(1) = prob%nY+2

  ! Interval
  allocate ( prob%npInt(1), prob%phsPt(2) )
  prob%npInt(1) = 10
  prob%phsPt(1) = 0_rp
  prob%phsPt(2) = 1_rp

  ! Linearity
  allocate ( prob%ytype(2,1), prob%ctype(1,1) )
  prob%ytype(1,1) = 0
  prob%ytype(2,1) = 0

  prob%ctype(1,1) = 0

  !---------------------------------------------------------------------
  ! Solve the problem
  !---------------------------------------------------------------------
  Start = 0
  call snctrlD ( Start, prob, odecon, algcon, varbds, &
                 INFO, mincw, miniw, minrw, nS, nInf, sInf, &
                 cu, lencu, iu, leniu, ru, lenru, &
                 cw, lencw, iw, leniw, rw, lenrw )

  if (INFO == 82 .OR. INFO == 83 .OR. INFO == 84) then
     go to 910
  end if

  write(nOut, *) ' '
  write(nOut, *) 'snctrl finished.'
  write(nOut, *) 'INFO          =', INFO
  write(nOut, *) 'nInf          =', nInf
  write(nOut, *) 'sInf          =', sInf

  call sncEnd( prob )

  if (INFO >= 30) go to 900

  stop

  !---------------------------------------------------------------------
  ! Error exit.
  !---------------------------------------------------------------------
800 write(nOut, 4000) 'Error while opening file', lfile
  stop

900 write(nOut, *) ' '
  write(nOut, *) 'STOPPING because of error condition'
910 stop

4000 format(/  a, 2x, a  )

end program catmixD

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine varbds ( curPhs, nPhs, nY, nU, nP, nC, nNodes, lbds, ubds, &
                    x, plbds, pubds, p, clbds, cubds )
  use precision, only : ip, rp
  implicit none
  integer(ip), intent(in)  :: curPhs, nPhs, nY, nU, nP, nC, nNodes
  real(rp),    intent(out) :: lbds(nY+nU,nNodes), ubds(nY+nU,nNodes), &
                              x(nY+nU,nNodes), plbds(nP), pubds(nP), p(nP), &
                              clbds(nC), cubds(nC)

  !=============================================================================
  integer(ip) :: j
  real(rp)    :: bplus, bminus
  parameter    ( bplus  = 1.0d+20, bminus = -bplus )

  ! Default set everything to free bounds.
  lbds = bminus
  ubds = bplus
  x    = 0.0

  ! At t0
  lbds(1,1) = 1.0
  lbds(2,1) = 0.0

  ubds(1,1) = 1.0
  ubds(2,1) = 0.0

  ! 0 <= u1 <= 1
  do j = 1, nNodes
     lbds(1+nY,j) = 0.0
     ubds(1+nY,j) = 1.0
  end do

  ! algebraic constraint
  clbds(1) = 1.0
  cubds(1) = 1.0

end subroutine varbds

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine odecon ( snStat, curPhs, nPhs, nY, nU, nP, nNodes, F, J, &
                    dvar, pvar, needF, needJ, cu, lencu, iu, leniu, ru, lenru )
  use precision, only : ip, rp
  implicit none
  integer(ip), intent(in) :: snStat, curPhs, nPhs, nY, nU, nP, nNodes, &
                             needF, needJ, lencu, leniu, lenru
  real(rp),    intent(in) :: dvar(nY+nU,nNodes), pvar(nP)

  integer(ip), intent(inout) :: iu(leniu)
  real(rp),    intent(inout) :: ru(lenru)
  character(8),intent(inout) :: cu(lencu)

  real(rp),    intent(out) :: F(nY,nNodes), J(nY,nY+nU+nP,nNodes)

  !=============================================================================
  integer(ip) :: jt

  if ( needF > 0 ) then
     do jt = 1, nNodes
        F(1,jt) = dvar(1+nY,jt)*(10.0*dvar(2,jt) - dvar(1,jt))
        F(2,jt) = dvar(1+nY,jt)*(dvar(1,jt) - 10*dvar(2,jt)) - &
             (1-dvar(1+nY,jt))*dvar(2,jt)
     end do
  end if

  if ( needJ > 0 ) then
     do jt = 1, nNodes
        J(1,1,jt) = -dvar(1+nY,jt)
        J(1,2,jt) = 10.0*dvar(1+nY,jt)
        J(1,3,jt) = 10.0*dvar(2,jt) - dvar(1,jt)
        J(1,4,jt) = 0.0

        J(2,1,jt) = dvar(1+nY,jt)
        J(2,2,jt) = -10.0*dvar(1+nY,jt) - (1-dvar(1+nY,jt))
        J(2,3,jt) = dvar(1,jt) - 10*dvar(2,jt) + dvar(2,jt)
        J(2,4,jt) = 0.0
     end do
  end if

end subroutine odecon

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine algcon ( snStat, curPhs, nPhs, nC, nY, nU, nP, nNodes, C, G, &
                    dvar, pvar, needC, needG, cu, lencu, iu, leniu, ru, lenru )
  use precision, only : ip, rp
  implicit none
  integer(ip), intent(in) :: snStat, curPhs, nPhs, nC, nY, nU, nP, nNodes, &
                             needC, needG, lencu, leniu, lenru
  real(rp),    intent(in) :: dvar(nY+nU,nNodes), pvar(nP)

  integer(ip), intent(inout) :: iu(leniu)
  real(rp),    intent(inout) :: ru(lenru)
  character(8),intent(inout) :: cu(lencu)

  real(rp),    intent(out) :: C(nC,nNodes), G(nC,nY+nU+nP,nNodes)

  !=============================================================================
  integer(ip) :: jt

  if ( needC > 0 ) then
     do jt = 1, nNodes
        C(1,jt) = dvar(1,jt) + dvar(2,jt) - dvar(2+nY,jt)
     end do
  end if

  if ( needG > 0 ) then
     do jt = 1, nNodes
        G(1,1,jt) =  1.0
        G(1,2,jt) =  1.0
        G(1,3,jt) =  0.0
        G(1,4,jt) = -1.0
     end do
  end if

end subroutine algcon

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
