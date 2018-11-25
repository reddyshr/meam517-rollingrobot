*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     File  nlp.f
*
*     newnlp   setnlp   s2tnlp   getnlp   g2tnlp   ininlp
*
*     Written by Mike Gertz - 07-Mar-2000
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine newnlp
     &   ( ia, lenia, lusdi, morei )

      implicit
     &     none
      integer
     &     lenia, lusdi, morei, ia(lenia)

*     ==================================================================
*     newnlp
*
*     07 Mar 2000: First version of newnlp.
*     29 Jan 2001: Reduced the size of the nlp data-structure to 10
*     07 Feb 2001: Current version
*     ==================================================================

      if ( lenia .lt. 10 ) then
         morei = 10 - lenia
         lusdi = 0
      else
         morei = 0
         lusdi = 10
         ia(1) = 0
         ia(2) = 0
      end if

      end ! subroutine newnlp

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine setnlp
     &   ( lenia, lenra,
     &     iiadmx, liadmx, iradmx, lradmx,
     &     iiwork, liwork, irwork, lrwork,
     &     inlp, linlp )

      implicit
     &     none
      integer
     &     lenia, lenra, iiadmx, liadmx, iradmx, lradmx, linlp,
     &     iiwork, liwork, irwork, lrwork, inlp(linlp)

*     ==================================================================
*     setnlp
*
*     07 Mar 2000: First version of setnlp.
*     29 Jan 2001: Eliminated copy of the Jacobian from the structure
*     07 Feb 2001: Current version
*     ==================================================================

      inlp(1)   = lenia
      inlp(2)   = lenra
      inlp(3)   = iiadmx
      inlp(4)   = liadmx
      inlp(5)   = iradmx
      inlp(6)   = lradmx
      inlp(7)   = iiwork
      inlp(8)   = liwork
      inlp(9)   = irwork
      inlp(10)  = lrwork

      end ! subroutine setnlp

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine getnlp
     &   ( lenia, lenra,
     &     iiadmx, liadmx, iradmx, lradmx,
     &     iiwork, liwork, irwork, lrwork,
     &     inlp, linlp )

      implicit
     &     none
      integer
     &     lenia, lenra, iiadmx, liadmx, iradmx, lradmx, linlp,
     &     iiwork, liwork, irwork, lrwork, inlp(linlp)

*     ==================================================================
*     getnlp
*
*     07 Mar 2000: First version of getnlp.
*     29 Jan 2001: Eliminated copy of the Jacobian from the structure
*     07 Feb 2001: Current version
*     ==================================================================

      lenia  = inlp(1)
      lenra  = inlp(2)
      iiadmx = inlp(3)
      liadmx = inlp(4)
      iradmx = inlp(5)
      lradmx = inlp(6)
      iiwork = inlp(7)
      liwork = inlp(8)
      irwork = inlp(9)
      lrwork = inlp(10)

      end ! subroutine getnlp

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine ininlp
     &   ( x, n, blx, bux,
     &     bl, mCon, bu,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     inlp, linlp, rnlp, lrnlp,
     &     lusdi, lusdr, morei, morer )

      implicit
     &     none
      integer
     &     n, mCon, lusdi, lusdr, morei, morer, linlp, lrnlp,
     &     lencu, leniu, lenru, lencw, leniw, lenrw, inlp(linlp),
     &     iu(leniu), iw(leniw)
      double precision
     &     x(n), blx(n), bux(n), bl(mCon), bu(mCon), rnlp(lrnlp),
     &     ru(lenru), rw(lenrw)
      character*8
     &     cu(lencu), cw(lencw)

*     ==================================================================
*     ininlp
*
*
*     Initialize the integer and real parts of the nlp data structure.
*
*     07 Mar 2000: First version of ininlp.
*     07 Feb 2001: Eliminated the copy of the Jacobian
*     07 Feb 2001: Current version
*     ==================================================================
      integer
     &     iiadmx, liadmx, iradmx, lradmx, iiwork, liwork, irwork,
     &     lrwork, k
*     ------------------------------------------------------------------

      morei   = 0
      morer   = 0

      iiadmx  = lusdi + 1
      liadmx  = 0
      iradmx  = 1
      lradmx  = 0

      call nwadmx
     &   ( inlp(iiadmx), linlp - lusdi, liadmx, morei )
      if ( morei .ne. 0 ) goto 600

      call iniad
     &   ( x, n, blx, bux,
     &     bl, mCon, bu,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     inlp(iiadmx), linlp - lusdi,
     &     rnlp(iradmx), lrnlp - lusdr,
     &     liadmx, lradmx, morei, morer )

      if ( morei .ne. 0  .or.  morer .ne. 0 ) goto 600
      lusdi = lusdi + liadmx
      lusdr = lusdr + lradmx

*     Allocate some workspace. First,  make sure space is available.

      if ( lusdi + n .gt. linlp  .or.  lusdr + n .gt. lrnlp ) then
         morei = max( 0, lusdi + n - linlp )
         morer = max( 0, lusdr + n - lrnlp )
         goto 600
      end if

*     Now do it.

      iiwork = lusdi + 1
      liwork = n
      lusdi  = lusdi + n

      irwork = lusdr + 1
      lrwork = n
      lusdr  = lusdr + n

      call setnlp
     &   ( lusdi, lusdr,
     &     iiadmx, liadmx, iradmx, lradmx,
     &     iiwork, liwork, irwork, lrwork,
     &     inlp, linlp )

      goto 900

*     Memory allocation errors

  600 lusdi = 0
      lusdr = 0
  900 continue

      end ! subroutine ininlp
