*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     File  sncmexlog.f
*
*     snPRNT  snAbort    snPROB
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snREAD
     &   ( unitno, string, nchar, endfile )

      implicit
     &     none
      character*(*)
     &     string
      integer
     &     endfile, nchar, unitno
*     ==================================================================
*     snREAD reads a string of length nchar from file  unitno.
*
*     30 Apr 2006: First version of snREAD.
*     30 Apr 2006: Matlab version.
*     ==================================================================
      call sncmxread ( unitno, string, nchar, endfile )

      end ! subroutine snREAD

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snPRNT
     &   ( mode, string, iw, leniw )

      implicit
     &     none
      character*(*)
     &     string
      integer
     &     mode, leniw, iw(leniw)

*     ==================================================================
*     snPRNT  prints a trimmed form of "string" on various files.
*     If mode = 0,      nothing is output.
*     If mode = 1,      string is output to iPrint.
*     If mode = 2,      string is output to iSumm.
*     If mode = 3 or 4, string is output to iPrint and iSumm.
*     If mode = 4,      string is output to the screen.
*                       This mode is intended for error messages.
*     If mode = 5,      string is output to iStdo (standard output)
*                       This mode is to be used when the elements of
*                       the integer work array iw cannot be trusted.
*
*     mode 11-15 are the same as mode 1-5 with blank line before output.
*
*     If mode > 15 then nothing is printed unless  lvlSys > 0.
*     mode 21-25 are the same as mode 1-5
*     mode 31-35 are the same as mode 11-15
*
*     25 Sep 2002: First version of snPRNT.
*     31 Jul 2003: mode 11-14 added.  form introduced.
*     27 Dec 2003: mode 5 added to allow printing before iw is set.
*     12 Mar 2004: s1trim called to trim the string.
*     22 Jun 2004: System printing option added.
*     14 Oct 2004: Matlab version of snPRNT.
*     30 Apr 2006: Files opened and closed in C.
*     ==================================================================
      integer
     &     iPrint, iSumm, length, lvlSys, m, s1outpt, newline,
     &     screenOK, summaryOK, printOK
      character
     &     Buff*140
*     ------------------------------------------------------------------
      lvlSys    = iw( 71) ! > 0   => print system info

      newline = 0
      m       = 0
      if (mode .le.  0) then
!        Relax
      else if (mode   .lt. 10) then
         m       = mode
      else if (mode   .lt. 20) then ! Blank line first
         m       = mode - 10
         newline = 1
      else if (lvlSys .gt.  0) then ! Print system Info
         if (mode .lt. 30) then
            m       = mode - 20
         else
            m       = mode - 30
            newline = 1
         end if
      end if

      if (m .gt. 0) then

         call sncmxfilestatus( screenOK, summaryOK, printOK  )

!        length = len_trim(string)     ! An F90 intrinsic
         call s1trim( string, length ) ! The F77 equivalent
         Buff = string

         if (m .eq. 5) then
            call sncmxwritescreen( Buff, length)
         else
            iPrint = iw( 12) ! Print file
            iSumm  = iw( 13) ! Summary file

            if (m .eq. 1  .or.  m .ge. 3) then
               if (printOK .gt. 0) then
                  call sncmxwritefile( newline, iPrint, Buff, length )
               end if
            end if

            if (m .eq. 2  .or.  m .ge. 3) then
               if (screenOK  .gt. 0) then
                  call sncmxwritescreen( Buff, length)
               end if
               if (summaryOK .gt. 0) then
                  call sncmxwritefile( newline, iSumm , Buff, length )
               end if
            end if

            if (m .eq. 4) then
               call sncmxwritescreen( Buff, length)
            end if
         end if
      end if

      end ! subroutine snPRNT

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snAbort
     &   ( iAbort, info, Htype, KTcond, MjrPrt, minimz,
     &     m, maxS, n, nb, nnCon0, nnCon, nnObj0, nnObj, nS,
     &     itn, nMajor, nMinor, nSwap,
     &     condHz, iObj, sclObj, ObjAdd, fMrt, PenNrm, step,
     &     prInf, duInf, vimax, virel, hs,
     &     ne, nlocJ, locJ, indJ, Jcol, negCon,
     &     Ascale, bl, bu, fCon, gCon, gObj,
     &     yCon, pi, rc, rg, x,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      logical
     &     KTcond(2)
      integer
     &     Htype, iAbort, info(6), iObj, itn,
     &     lencu, lencw, leniu, leniw, lenru, lenrw,
     &     MjrPrt, minimz, m, maxS, n, nb, ne, negCon, nlocJ,
     &     nnCon0, nnCon, nnObj0, nnObj, nMajor, nMinor, nS, nSwap,
     &     hs(nb), locJ(nlocJ), indJ(ne), iu(leniu), iw(leniw)
      double precision
     &     condHz, sclObj, ObjAdd, fMrt, PenNrm, virel, vimax, step,
     &     prInf, duInf, Ascale(nb), bl(nb), bu(nb), fCon(nnCon0),
     &     gCon(negCon), gObj(nnObj0), Jcol(ne), pi(m),
     &     rc(nb), rg(maxS), yCon(nnCon0), x(nb), ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

*     ==================================================================
*     snAbort  is called every major iteration.
*     If iAbort > 0 on exit, the run is terminated.
*     By specifying a custom version of snSTOP, the user can arrange for
*     snopt to be terminated at any given major iteration.
*
*     14 Oct 2004: First version of   snAbort.
*     01 Sep 2007: Parameter list expanded.
*     ==================================================================
      iAbort    = 0
      call sncmxabort( iAbort )

      end ! subroutine snAbort

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snPROB
     &     ( Prob )

      implicit
     &     none
      character
     &     Prob*8

*     ==================================================================
*     Assigns an empty problem name.
*
*     31 Dec 2002: First version of snPROB
*     ==================================================================
      character          Blank*8
      parameter         (Blank ='        ')
*     ------------------------------------------------------------------
      Prob = Blank

      end ! subroutine snPROB

