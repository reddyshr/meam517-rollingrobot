module sqmain3Module

  implicit none
  public   :: sqdata, userHx, sqdata2

contains

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine sqdata( maxm, maxn, maxne,               &
                   m, n, ne, nName, lencObj, ncolH, &
                   iObj, ObjAdd, Prob,              &
                   Acol, indA, locA, bl, bu, cObj,  &
                   Names, hEtype, hs, x )

    integer            :: maxm, maxn, maxne, m, n, ne,      &
                          nName, lencObj, ncolH, iObj,      &
                          hEtype(maxn+maxm), hs(maxn+maxm), &
                          indA(maxne), locA(maxn+1)
    character          :: Prob*8, Names(maxn+maxm)*8
    real(8)            :: ObjAdd, Acol(maxne),              &
                          bl(maxn+maxm), bu(maxn+maxm),     &
                          cObj(maxn), x(maxn+maxm)

    !=================================================================
    ! sqdat0   defines the problem discusses in the SQOPT Users Guide.
    !
    ! (1) Compute l, u, and A so that the constraints are ranges of the
    !     form  l <= Ax <= u.
    !     Store l and u in bl(n+1:n+m) and bu(n+1:n+m).
    !
    ! (2) Set up the constants ObjAdd and c so that the explicit
    !     objective is
    !           ObjAdd + cObj'*x + half*x'*H*x
    !     If necessary, include an additional linear objective terms
    !     as row iObj of A.
    !=================================================================
    integer              :: i, j, neA
    real(8)              :: infBnd
    real(8), parameter   :: zero = 0.0d+0, one = 1.0d+0


    Prob   = 'sqProb  '  ! Give the problem a name

    !-----------------------------------------------------------------
    ! Assign the constraint nonzeros to Acol, column by column.
    ! indA(i) gives the row index of element Acol(i).
    ! locA(j) gives the index in a of the start of column j.
    !-----------------------------------------------------------------
    !
    !   -inf     ( 1  -1                           )    0
    !   -inf     (     1  -1                       )    0
    !   -inf     (         1  -1                   )    0
    !   -inf le  (             1  -1               ) le 0
    !   -inf     (                 1  -1           )    0
    !   -inf     (                     1  -1       )    0
    !   -inf     (                         1  -1   )    0
    !      1 =   ( 1   1   1   1   1   1   1   1   ) =  1
    !-----------------------------------------------------------------

    n       = 30
    m       = n               ! Does not include an objective row
    ne      = n + 2*(n-1)

    ObjAdd  = one
    lencObj = n
    iObj    = 0

    ncolH   = n
    nName   = 1

    infBnd  =  1.0d+20

    neA = 0

    do j = 1, n                   ! Set the elements of column j

       locA( j) =  neA + 1

       if (j .gt. 1) then
          neA       =  neA + 1
          indA(neA) =  j   - 1
          Acol(neA) = -one
       endif

       if (j .lt. n) then
          neA       =  neA + 1
          indA(neA) =  j
          Acol(neA) =  one
       end if

       neA       =  neA + 1
       indA(neA) =  m
       Acol(neA) =  one

    end do

    locA(n+1) =  neA + 1     ! Don't forget to finish off  locA.

    !-----------------------------------------------------------------
    ! Set the upper and lower bounds on the variables
    !-----------------------------------------------------------------
    do i = 1, n-1
       bl(i) = zero          ! Bounds on  x
       bu(i) = infBnd
       j     = n + i         ! Bounds on Ax
       bl(j) = -infBnd
       bu(j) =  zero
    end do

    bl(n)   = zero
    bu(n)   = infBnd
    bl(n+m) =  one
    bu(n+m) =  one

    !-----------------------------------------------------------------
    ! Set the objective terms.
    ! The objective linear term is explicit.
    !-----------------------------------------------------------------
    cObj(1) = one
    do j = 2, n
       cObj(j) = - cObj(j-1)
    end do

    !-----------------------------------------------------------------
    ! Set the initial value and status of each variable.
    ! For want of something better to do, make the variables x(1:n)
    ! temporarily fixed at their current values.
    ! Crash can set the rest.
    !-----------------------------------------------------------------
    hs(1:n) = 0

    !-----------------------------------------------------------------
    ! Fix the column variables to be non-elastic
    ! and the row variables to be elastic.
    !-----------------------------------------------------------------
    hEtype(1:n)     = 0
    hEtype(n+1:n+m) = 3

    !-----------------------------------------------------------------
    ! Set the initial estimate of the solution.
    !-----------------------------------------------------------------
    x(1) = -one
    do j = 2, n
       x(j) = -x(j-1)
    end do

  end subroutine sqdata

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine userHx( ncolH, x, Hx, State, cu, lencu, iu, leniu, ru, lenru )

    integer              :: ncolH, State, lencu, leniu, lenru, iu(leniu)
    character            :: cu(lencu)*8
    real(8)              :: x(ncolH), Hx(ncolH), ru(lenru)

    !=================================================================
    ! This is the user-defined Hessian-vector product for the
    ! sqopt example program.
    !
    !          (  2  -1   0   0   0  0  0 . 0  0 )
    !          ( -1   2  -1   0   0  0  0 . 0  0 )
    !          (  0  -1   2  -1   0  0  0 . 0  0 )
    !          (  0   0   .   .   .  0  0 . 0  0 )
    !          (  0   0   0   .   2 -1  0 . 0  0 )
    !     H =  (  0   0   0   0  -1  2  0 . 0  0 )   rank(H) = ncolH
    !          (  0   0   0   0   0  0  0 . 0  0 )
    !          (  .   .   .   .   .  .  . . .  . )
    !          (  0   0   0   0   0  0  0 . 0  0 )
    !          (  0   0   0   0   0  0  0 . 0  0 )
    !          (  0   0   0   0   0  0  0 . 0  0 )
    !==================================================================

    integer              :: j
    real(8)              :: rho
    real(8), parameter   :: two   = 2.0d+0, five  = 5.0d+0
    integer, parameter   :: nOut  = 6


    if (State .eq. 1) then    ! First entry.  Print on standard output.
       if (nOut .gt. 0) write(nOut, 1000) ncolH
    end if

    rho = five                ! Smoothing parameter

    Hx(1) = - x(2) + two*x(1)
    do j = 2, ncolH-1
       Hx(j) = - x(j+1) + two*x(j) - x(j-1)
    end do
    Hx(ncolH) =  two*x(ncolH) - x(ncolH-1)

    do j = 1, ncolH
       Hx(j) = rho*Hx(j)
    end do

    if (State .ge. 2) then    ! Final entry.
       if (nOut .gt. 0) write(nOut, 2000)
    end if
    return

1000 format(/ ' This is problem  sqmain.   ncolH =', i4)
2000 format(/ ' Finished         sqmain.')

  end subroutine userHx

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine sqdata2( maxm, maxn, maxne,               &
                      m, n, ne, nName, lencObj, ncolH, &
                      iObj, ObjAdd, Prob,              &
                      Acol, indA, locA, bl, bu, cObj,  &
                      Names, hEtype, hs, x )

    integer              :: maxm, maxn, maxne, m, n, ne,      &
                            nName, lencObj, ncolH, iObj,      &
                            hEtype(maxn+maxm), hs(maxn+maxm), &
                            indA(maxne), locA(maxn+1)
    character            :: Prob*8, Names(maxn+maxm)*8
    real(8)              :: ObjAdd, Acol(maxne),              &
                            bl(maxn+maxm), bu(maxn+maxm),     &
                            cObj(maxn), x(maxn+maxm)

    !================================================================
    ! sqdat0   defines the problem discusses in the SQOPT Users Guide.
    !
    ! (1) Compute l, u, and A so that the constraints are ranges of the
    !     form  l <= Ax <= u.
    !     Store l and u in bl(n+1:n+m) and bu(n+1:n+m).
    !
    ! (2) Set up the constants ObjAdd and cObj so that the explicit
    !     objective is
    !         ObjAdd + (row iObj)*x + half*x'*H*x
    !================================================================
    integer              :: i, j, neA
    real(8)              :: infBnd, cj
    real(8), parameter   :: zero = 0.0d+0, one = 1.0d+0


    Prob   = 'sqProb 2'       ! Give the problem a name

    !----------------------------------------------------------------
    ! Assign the constraint nonzeros to Acol, column by column.
    ! indA(i) gives the row index of element Acol(i).
    ! locA(j) gives the index in a of the start of column j.
    !----------------------------------------------------------------
    !
    !   -inf     ( 1  -1                           )    0
    !   -inf     (     1  -1                       )    0
    !   -inf     (         1  -1                   )    0
    !   -inf le  (             1  -1               ) le 0
    !   -inf     (                 1  -1           )    0
    !   -inf     (                     1  -1       )    0
    !   -inf     (                         1  -1   )    0
    !   -inf le  (-1   1  -1   1  -1   1  -1   1   ) le inf
    !      1 =   ( 1   1   1   1   1   1   1   1   ) =  1
    !----------------------------------------------------------------

    n       = 30
    m       = n + 1           ! Includes the objective row
    ne      = 2*n + 2*(n-1)

    ObjAdd  = one
    iObj    = n
    lencObj = 0

    ncolH   = n
    nName   = 1

    infBnd  =  1.0d+20

    neA     = 0
    cj      = one

    do j = 1, n     ! Set the elements of column j

       locA( j) =  neA + 1

       if (j .gt. 1) then
          neA       =  neA + 1
          indA(neA) =  j   - 1
          Acol(neA) = -one
       endif

       if (j .lt. n) then
          neA       =  neA + 1
          indA(neA) =  j
          Acol(neA) =  one
       end if

       neA       =  neA + 1
       indA(neA) =  iObj
       Acol(neA) =   cj
       cj        =  -cj

       neA       =  neA + 1
       indA(neA) =  m
       Acol(neA) =  one
    end do

    locA(n+1) =  neA + 1     ! Don't forget to finish off  locA.

    !---------------------------------------------------------------
    ! Set the upper and lower bounds on the variables
    !---------------------------------------------------------------
    do i = 1, n-1
       bl(i) = zero        ! Bounds on  x
       bu(i) = infBnd
       j     = n + i       ! Bounds on Ax
       bl(j) = -infBnd
       bu(j) =  zero
    end do

    bl(n)   = zero
    bu(n)   = infBnd
    bl(n+m) =  one
    bu(n+m) =  one

    bl(n+iObj) = -infBnd
    bu(n+iObj) =  infBnd

    !---------------------------------------------------------------
    ! Set the initial value and status of each variable.
    ! For want of something better to do, make the variables x(1:n)
    ! temporarily fixed at their current values.
    ! The crash can set the rest.
    !---------------------------------------------------------------
    hs(1:n) = 0

    !---------------------------------------------------------------
    ! Fix the column variables to be non-elastic
    ! and the row variables to be elastic.
    !---------------------------------------------------------------
    hEtype(1:n)     = 0
    hEtype(n+1:n+m) = 3

    !---------------------------------------------------------------
    ! Set the initial estimate of the solution.
    !---------------------------------------------------------------
    x(1) = -one
    do j = 2, n
       x(j) = -x(j-1)
    end do

  end subroutine sqdata2

end module sqmain3Module
