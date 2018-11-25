      program
     &     stats
      implicit
     &     double precision (a-h,o-z)

      integer
     &     INFO
      character
     &     name*35, aexit*35
      character
     &     msg*43
      character*3
     &     dit
      double precision
     &     time_limit

      parameter       ( infl  = 4, iscr = 6, nTeX = 7, nrun = 8 )
      parameter       ( isumm = 9 )
      parameter       ( zero  = 0.0d+0 )

*     open  ( unit=infl, name='data.txt',status='OLD' )
*     rewind( unit=infl, err=999 )

      timtot = zero
      nline  = 1
      numcpu = 0
      numGnf = 0
      numLnf = 0
      numopt = 0
      numfp  = 0
      numinf = 0
      nummax = 0
      numacc = 0
      numSB  = 0
      numcbi = 0
      numUnb = 0
      numvlm = 0
      mintot = 0
      majtot = 0
      nfntot = 0
      timtot = zero
      time_limit = 3600         ! one hour
      time_limit = 6000000

      write(nrun, 5000)

*     Note: "Minor print level" > 0 for stats.f to work for an LP.

      do 100 i = 1, 2000
         nInf   = 0
         read (infl, 1004, end=900) name ! 1. BEGIN SNOPT Problem
         numtot = i

         read (infl, 1111) nnCon, linCon ! 2. Nonlinear constraints
         read (infl, 1111) nnL  , linVar ! 3. Nonlinear variables

         nconst = nnCon + linCon
         nVars  = nnL   + linVar
         majits = 0
         nfObj  = 0
         nfCon  = 0
         nfmax  = 0
         vmax   = zero

         read (infl, 1050) aexit         ! 4. SNOPTB EXIT
         read (infl, 1055) INFO          ! 5. SNOPTB INFO
         read (infl, 1008) name          ! 6. Problem name
         read (infl, 1010) itn, obj      ! 7. No. of iterations

         if (nnL .gt. 0) then ! NLP or QP (SNOPT)
            read (infl, 1022) msg, itmp  ! 8.

            if (INFO .eq. 31 .or. INFO .eq. 32
     &     .or. INFO .eq. 41 .or. INFO .eq. 42 .or. INFO .eq. 43) then
               if (msg .eq. ' No. of infeasibilities       ') then
*                            123456789012345678901234567890
                  read (infl, 1022) msg, itmp ! 8.
               end if
            end if

            if (INFO .eq. 11 .or. INFO .eq. 12 .or. INFO .eq. 13) then ! infeasible
               nInf   = itmp
               if      (INFO .eq. 11) then ! infeasible linear constraints
*                 Relax
               else if (INFO .eq. 12) then ! infeasible linear equalities
*                 Relax
               else if (INFO .eq. 13) then ! infeasible nonlinear constraints
                  read (infl, 1020) majits
               end if
            else
               majits = itmp
            end if

            if (INFO .ne. 11  .and. INFO .ne. 12) then
               read (infl, 1021) nfObj, nfCon           !  9
               if (nnCon .gt. 0) read (infl, 1025) vmax ! 10
            end if
         else                ! LP        (SQOPT)
            if (INFO .eq. 11 .or. INFO .eq. 12 .or. INFO .eq. 14) then
               read (infl, 1020) nInf
            end if
         end if

         read (infl, 1030) time
         if (time .gt. time_limit) then
            time = time_limit
            INFO = -1
         end if

*        Write the Header for the files.

         if ( mod (nline, 75) .eq. 1 ) write(isumm,2000)
         if ( mod (nline, 50) .eq. 1 ) then
                   if (nline  .gt. 1 ) write(nTeX, 3000)
                                       write(nTeX, 3010)
         end if
         if ( mod (nline, 10) .eq. 1 ) write(iscr, 3050)

*        Write the summary of the problem statistics

         if (nnL .gt. 0) then
            nfmax = max( nfObj, nfCon )
         end if

         write(isumm, 2100) i, name, nvars, linCon, nnCon, obj

         mintot = mintot + itn
         majtot = majtot + majits
         nfntot = nfntot + nfmax
         timtot = timtot + time

*        1 hour = 3600

         if      (     INFO .eq. -1) then
            dit    = 'Cpu'
            msg    = 'Cpu time limit exceeded'
            numcpu = numcpu  + 1
         else if (     INFO .eq.  1) then
            dit    = '   '
            msg    = 'Optimal solution found'
            numopt = numopt + 1
         else if (INFO .eq.  2) then
            dit    = 'fp '
            msg    = 'Feasible point found'
            numfp  = numfp  + 1
         else if (INFO .eq.  3) then
            dit    = 'acc'
            msg    = 'Optimal, but accuracy not achvd'
            numacc = numacc + 1
         else if (INFO .eq. 11 .or. INFO .eq. 12
     &       .or. INFO .eq. 14 .or. INFO .eq. 15) then
            dit    = 'Lnf'
            msg    = 'Infeasible linear constraints'
            numLnf = numLnf + 1
         else if (INFO .eq. 13) then
            dit    = 'Inf'
            msg    = 'Infeasible nonlinear constraints'
            numinf = numinf + 1
         else if (INFO .eq. 31) then
            dit    = 'Itr'
            msg    = 'Iteration limit exceeded'
            nummax = nummax + 1
         else if (INFO .eq. 32) then
            dit    = 'Itr'
            msg    = 'Major iterations limit'
            nummax = nummax + 1
         else if (INFO .eq. 33) then
            dit    = 'SB '
            msg    = 'Superbasics limit too small'
            numSB  = numSB + 1
         else if (INFO .eq. 41) then
            dit    = 'Cbi'
            msg    = 'Current point cannot be improved'
*                     12345678901234567890123456789012
            numcbi = numcbi + 1
         else if (INFO .eq. 42) then
            dit    = 'Gnf'
            msg    = 'Singular basis'
         else if (INFO .eq. 43) then
            dit    = 'Gnf'
            msg    = 'Constraints cannot be satisfied'
            numGnf = numGnf + 1
         else if (INFO .eq. 44) then
            dit    = 'Gnf'
            msg    = 'Ill-conditioned null=space basis'
            numGnf = numGnf + 1
         else if (INFO .eq. 21) then
            dit    = 'Unb'
            msg    = 'Unbounded problem'
            numUnb = numUnb + 1
         else if (INFO .eq. 22) then
            dit    = 'vlm'
            msg    = 'Violation limit exceeded'
            numvlm = numvlm + 1
         else if (INFO .gt. 0) then
            dit    = '???'
            msg    = ' '
         else
            write(iscr, '(a)') ' Some problem died before EXIT'
            stop
         end if

         write(nTeX, 3100) i, name, itn, majits, nfmax, obj, vmax,
     &                     time, dit
         write(iscr, 3200) i, name, itn, majits, nfmax, obj, vmax,
     &                     time, dit
         write(nrun, 5100) name, nvars, nconst, majits, nfObj, nfCon,
     &                     time, obj, nInf, msg

         nline = nline + 1
  100 continue

  900 write(nTeX, 4000)
      write(nTeX, 4100) numtot, numopt, numfp , numinf, numacc, nummax,
     &                  numSB , numLnf, numcbi, numGnf, numUnb, numvlm,
     &                  numcpu, majtot, mintot, nfntot, timtot
      write(iscr, 4100) numtot, numopt, numfp , numinf, numacc, nummax,
     &                  numSB , numLnf, numcbi, numGnf, numUnb, numvlm,
     &                  numcpu, majtot, mintot, nfntot, timtot
      write(nTeX, 4200)

  999 stop

 1004 format(29x,a)
 1005 format(a40)
 1008 format(30x,a)
 1010 format(30x,i8,23x,e17.10)
 1020 format(30x,i8)
 1022 format(a30,i8)
 1021 format(30x,i8,32x,i8)
 1025 format(29x,e9.2)
 1030 format(43x,f9.2)
 1050 format(a)
 1055 format(12x, i4)
 1111 format(23x,i7,24x,i7)

 2000 format(
     & //    52x,  'Linear   Nonlinear         Optimal'
     &  / '  No.    Problem Name', 20x, 'Variables',
     &           ' constrnts constrnts        objective'
     &  /)
 2100 format(i4, 1x, a35, 3i10, 1p, e18.6)

 3000 format(' \\end{tabular}' //)
 3010 format(' \\begin{tabular}{rlrrrrrrl}' /
     &   '  No.&          Problem Name&       Min',
     &   ' &     Maj &     Fun &',
     &   '    Objective   &  norm(c)   &    Time& Status \\\\[2ex]')
 3050 format(/ '   No.', 10x, 'Problem Name', 11x, '  Min   Maj   Fun',
     &       '    Objective    norm(c)     Time Status')
 3100 format( 1x, i4, '&', a22, '&', i10, 2(' &',i8), ' &$',
     &        1p, e13.6, ' $&$', e9.2 , ' $&', 0p, f8.2, '&{\\em ',
     &        a3, '} \\\\' )
 3200 format( 1x, i4, 1x, a32, 1x, i5, 2(1x,i5), 1x,
     &        1p, e13.6, 1x, e9.2 ,1x, 0p, f8.2, 1x, a3)

 4000 format(' \\end{tabular}'
     &    // ' \\begin{verbatim}')
 4100 format(/ 5x, '  Summary ' /
     &         5x, '  ================================= ' /
     &         5x, '  No. problems attempted :', i9 /
     &         5x, '  No. optimal            :', i9 /
     &         5x, '  No. feasible point only:', i9 /
     &         5x, '  No. infeasible         :', i9 /
     &         5x, '  No. accuracy not achvd :', i9 /
     &         5x, '  No. iterations limit   :', i9 /
     &         5x, '  No. superbasics limit  :', i9 /
     &         5x, '  No. infeasible lin con :', i9 /
     &         5x, '  No. cannot be improved :', i9 /
     &         5x, '  No. infeasible gen con :', i9 /
     &         5x, '  No. Unbounded          :', i9 /
     &         5x, '  No. violation limit    :', i9 /
     &         5x, '  No. cpu limit          :', i9 /
     &         5x, '  ---------------------------------', /
     &         5x, '  No. Major iterations   :', i9 /
     &         5x, '  No. Minor iterations   :', i9 /
     &         5x, '  No. Function evals.    :', i9 /
     &         5x, '  Total Time             :', f9.2 )
 4200 format(' \\end{verbatim}'
     &    // ' \\end{document}')

 5000 format( 'Name            n      m   Majr   nFun   Ncon',
     &        '    cpu(s)',
     &        '      Objective      inf  Result'/)
 5100 format( a8, 2x, 2i7, i7, 2i7, 1x, f9.2, 1p, 2x, e16.9, i6, 2x, a )

      end ! stats
