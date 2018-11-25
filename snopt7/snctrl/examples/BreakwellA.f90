!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File BreakwellA.f90
!
!  The problem in Bolza form:
!      min            (1/2) * integral{0, 1}  u^2
!      subject to:
!                     y1' = y2
!                     y2' = u
!                     y1(0) = 0
!                     y2(0) = 1
!
!                     y1(1) = 0
!                     y2(1) = -1
!
!                     y1 <= 0.1
!
!  The problem input to SNCTRL:
!        min          y1(1)
!        subject to:
!                     y1' = (1/2)*u^2
!                     y2' = y3
!                     y3' = u
!
!                     y1(0) = 0
!                     y2(0) = 0
!                     y3(0) = 1
!
!                     y2(1) = 0
!                     y3(1) = -1
!
!                     y2 <= 0.1
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

program BreakwellA
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
  lfile = 'BreakwellA.spc'
  open( iSpecs, file=lfile, status='OLD',     err=800 )
  lfile = 'BreakwellA.out'
  open( iPrint, file=lfile, status='UNKNOWN', err=800 )


  ! Set options to default values.
  call sncInit ( iPrint, iSumm, prob, cw, lencw, iw, leniw, rw, lenrw )


  ! Read a Specs file.
  call sncSpec( iSpecs, INFO, cw, lencw, iw, leniw, rw, lenrw )

  if (INFO /= 101  .AND.  INFO /= 107) go to 910


  ! Setup the problem.
  prob%probname = 'BrkwellA'
  prob%nY   = 3     ! number of states
  prob%nU   = 1     ! number of controls
  prob%nP   = 0     ! number of parameters
  prob%nC   = 0     ! number of algebraic constraints
  prob%nPhs = 1     ! number of phases

  ! Objective
  allocate ( prob%objL(1) )
  prob%objL(1) = 1

  ! Interval
  allocate ( prob%npInt(1), prob%phsPt(2) )
  prob%npInt(1) = 10
  prob%phsPt(1) = 0
  prob%phsPt(2) = 1.0

  ! Linearity
  allocate ( prob%ytype(3,1) )
  prob%ytype(1,1) = 0
  prob%ytype(2,1) = 1
  prob%ytype(3,1) = 1


  !---------------------------------------------------------------------
  ! Solve the problem.
  !---------------------------------------------------------------------
  Start = 0
  call snctrlA ( Start, prob, odecon, algcon, varbds, &
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


end program BreakwellA

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine varbds ( curPhs, nPhs, nY, nU, nP, nC, y0low, y0upp, yflow, yfupp, &
                    ylow, yupp, ulow, uupp, plow, pupp, p, clow, cupp )
  use precision, only : ip, rp
  implicit none
  integer(ip), intent(in)  :: nPhs, curPhs, nY, nU, nP, nC
  real(rp),    intent(out) :: y0low(nY), y0upp(nY), yflow(nY), yfupp(nY), &
                             ylow(nY),  yupp(nY),  ulow(nU),  uupp(nU), &
                             plow(nP), pupp(nP), p(nP), clow(nC), cupp(nC)

  !=============================================================================
  ! Upper/lower bounds on states, controls, parameters, and
  ! algebraic constraints.
  !=============================================================================
  real(rp) :: bplus, bminus
  parameter ( bplus  = 1.0d+20, bminus = -bplus )


  ! At t0
  y0low(1) = 0.0
  y0upp(1) = 0.0

  y0low(2) = 0.0
  y0upp(2) = 0.0

  y0low(3) = 1.0
  y0upp(3) = 1.0

  ! At tf
  yflow(1) = bminus
  yfupp(1) = bplus

  yflow(2) = 0.0
  yfupp(2) = 0.0

  yflow(3) = -1.0
  yfupp(3) = -1.0

  ! For t in (t0, tf):
  ! y1 free
  ylow(1)  = bminus
  yupp(1)  = bplus

  ! y2 <= 0.1
  ylow(2)  = bminus
  yupp(2)  = 0.1

  ! y3 free
  ylow(3)  = bminus
  yupp(3)  = bplus

  ! u1 free
  ulow(1) = bminus
  uupp(1) = bplus

end subroutine varbds

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine odecon ( snStat, curPhs, nPhs, nY, nU, nP, F, J, y, u, p, &
                    needF, needJ, cu, lencu, iu, leniu, ru, lenru )
  use precision, only : ip, rp
  implicit none
  integer(ip), intent(in)    :: snStat, curPhs, nPhs, nY, nU, nP, &
                               needF, needJ, lencu, leniu, lenru
  real(rp),    intent(in)    :: y(nY), u(nU), p(nP)

  integer(ip), intent(inout) :: iu(leniu)
  real(rp),    intent(inout) :: ru(lenru)
  character(8),intent(inout) :: cu(lencu)

  real(rp),    intent(out)   :: F(nY), J(nY, nY+nU+nP)

  !=============================================================================

  if ( needF > 0 ) then
     ! nonlinear
     F(1) = (u(1)**2)/2.0

     ! linear
     F(2) = y(3)
     F(3) = u(1)
  end if

  if ( needJ > 0 ) then
     J(1,1) = 0.0
     J(1,2) = 0.0
     J(1,3) = 0.0
     J(1,4) = u(1)

     J(2,1) = 0.0
     J(2,2) = 0.0
     J(2,3) = 1.0
     J(2,4) = 0.0

     J(3,1) = 0.0
     J(3,2) = 0.0
     J(3,3) = 0.0
     J(3,4) = 1.0
  end if

end subroutine odecon

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine algcon ( snStat, curPhs, nPhs, nC, nY, nU, nP, C, G, y, u, p, &
                    needC, needG, cu, lencu, iu, leniu, ru, lenru)
  use precision, only : ip, rp
  implicit none
  integer(ip), intent(in)    :: snStat, curPhs, nPhs, nC, nY, nU, nP, &
                               needC, needG, lencu, leniu, lenru
  real(rp),    intent(in)    :: y(nY), u(nU), p(nP)

  integer(ip), intent(inout) :: iu(leniu)
  real(rp),    intent(inout) :: ru(lenru)
  character(8),intent(inout) :: cu(lencu)

  real(rp),    intent(out)   :: C(nC), G(nC, nY+nU+nP)

end subroutine algcon

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
