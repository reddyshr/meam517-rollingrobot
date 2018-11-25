*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     snFileOpenRead   snFileOpenAppend   snFileClose
*  This files is here for the C/C++ interface.
*
*     snFilewrapper   snOpenappend   snClose   snOpen
*  This files is only used in the F2C'd version of SNOPT
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snFileOpenRead
     &     ( unitno, filename )
      implicit
     &     none
      integer
     &     unitno
      character*(*)
     &     filename
!     ==================================================================
!     Intel compiler doesn't like opening a file in one library and
!     passing it to another.
!
!     Wrapper to open files for reading inside SNOPT library.
!
!     10 Jan 2017: First version.
!     ==================================================================
      integer len

      call s1trim(filename,len)
      open(unit=unitno, file=filename, status='old')

      end ! subroutine snFileOpenRead

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snFileOpenAppend
     &     ( unitno, filename )
      implicit
     &     none
      integer
     &     unitno
      character*(*)
     &     filename
!     ==================================================================
!     Intel compiler doesn't like opening a file in one library and
!     passing it to another.
!
!     Wrapper to open files for append inside SNOPT library.
!
!     10 Jan 2017: First version.
!     ==================================================================
      integer len

      call s1trim(filename,len)
      open(unit=unitno, file=filename,
     &     status='unknown', position='append')

      end ! subroutine snFileOpenAppend

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snFileClose
     &     ( unitno )
      implicit
     &     none
      integer
     &     unitno
!     ==================================================================
!     Close the given unit file.
!
!     10 Jan 2017: First version.
!     ==================================================================

      close(unitno)

      end ! subroutine snFileClose

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c$$$
c$$$      subroutine snFilewrapper
c$$$     &   ( name, ispec, inform, cw, lencw, iw, leniw, rw, lenrw )
c$$$
c$$$      implicit
c$$$     &     none
c$$$      external
c$$$     &     snSpec
c$$$      character*(*)
c$$$     &     name
c$$$      integer
c$$$     &     inform, lencw, leniw, lenrw, iw(leniw), ispec
c$$$      double precision
c$$$     &     rw(lenrw)
c$$$      character
c$$$     &     cw(lencw)*8
c$$$
c$$$*     ==================================================================
c$$$*     Read options for snopt from the file named name. inform .eq 0 if
c$$$*     successful.
c$$$*
c$$$*     09 Jan 2000: First version of snFilewrapper.
c$$$*     ==================================================================
c$$$      integer
c$$$     &     iostat
c$$$
c$$$      open( ispec, iostat=iostat, file=name, status='old' )
c$$$      if (iostat .ne. 0) then
c$$$         inform = 2 + iostat
c$$$      else
c$$$         call snSpec( ispec, inform, cw, lencw, iw, leniw, rw, lenrw )
c$$$         close( ispec )
c$$$      end if
c$$$
c$$$      end ! subroutine snFilewrapper
c$$$
c$$$*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c$$$
c$$$      subroutine snOpenappend
c$$$     &   ( iunit, name, inform )
c$$$
c$$$      implicit
c$$$     &     none
c$$$      integer
c$$$     &     iunit
c$$$      character*(*)
c$$$     &     name
c$$$      integer
c$$$     &     inform
c$$$
c$$$*     ==================================================================
c$$$*     Open file named name to Fortran unit iunit. inform .eq. 0 if
c$$$*     sucessful.  Although opening for appending is not in the f77
c$$$*     standard, it is understood by f2c.
c$$$*
c$$$*     09 Jan 2000: First version of snOpenappend
c$$$*     ==================================================================
c$$$      open( iunit, iostat=inform, file=name, access='append' )
c$$$
c$$$      end ! subroutine snOpenappend
c$$$
c$$$*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c$$$
c$$$      subroutine snClose
c$$$     &   ( iunit )
c$$$
c$$$      implicit
c$$$     &     none
c$$$      integer
c$$$     &     iunit
c$$$
c$$$*     ==================================================================
c$$$*     Close unit iunit.
c$$$*
c$$$*     09 Jan 2000: First version of snClose
c$$$*     ==================================================================
c$$$      close( iunit )
c$$$
c$$$      end ! subroutine snClose
c$$$
c$$$*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c$$$
c$$$      subroutine snOpen
c$$$     &   ( iunit )
c$$$
c$$$      implicit
c$$$     &     none
c$$$      integer
c$$$     &     iunit
c$$$
c$$$*     =================================================================
c$$$*     Open file named name to Fortran unit iunit.  inform .eq. 0
c$$$*     if successful.
c$$$*     =================================================================
c$$$      character
c$$$     &     lfile*20
c$$$
c$$$      lfile='testing.out'
c$$$      open( iunit, file=lfile, status='UNKNOWN')
c$$$
c$$$      end ! subroutine snOpen
