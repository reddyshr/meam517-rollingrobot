      program            snmain

      implicit           none

      external           userfg

      integer            n, neF, nName, nFname
      integer            lencw, leniw, lenrw, lencu, leniu, lenru

*
      parameter( lencw = 501, leniw = 10000, lenrw = 20000 )
      parameter( lencu = 1, leniu = 1, lenru = 1 )
      parameter( n = 5, neF = 6 )
      parameter( nName = 1, nFname = 1 )

*

      character*8        Names(nName), Fnames( nFname )

      integer            nsnr
      parameter        ( nsnr = 500 )

      character*8        cu(lencu),          cw(lencw)
      integer            iu(leniu),          iw(leniw)
      double precision   ru(lenru),          rw(lenrw)

      integer            nbnd
      parameter        ( nbnd = n + neF )

      double precision   bl(nbnd), bu(nbnd)
      integer            xstate(n)
      double precision   x(n), xmul(n)
      integer            Fstate(neF)
      double precision   F(neF), Fmul(neF)

      double precision   nInfty,           pInfty
      parameter         (nInfty = -1.0d20, pInfty = 1.0d20)
      double precision   zero
      parameter         (zero = 0.0d0)

      character*8        Prob

      integer            INFO, mincw, miniw, minrw
      integer            nS, nInf
      double precision   sInf, ObjAdd

      integer            iSpecs, iPrint, iSumm, iErr

      integer            iiAfun, ijAvar, lenA, neA, iA
      integer            iiGfun, ijGvar, lenG, neG, mm

      integer            lPrint, lSumm
      integer            lusdi, lusdr, moreiw, morerw
      integer            ObjRow, iCold

      integer            iiadmx, liadmx, iradmx, lradmx
      integer            linlp,  lrnlp
      integer            iiwork, liwork, irwork, lrwork
      integer            lenia,  lenra
*
*     Data statements
*
      data bl          / nbnd * nInfty /
      data bu          / nbnd * pInfty /
      data xstate      / n * 0 /
      data x           / n * zero /
      data xmul        / n * zero /
      data Fstate      / neF * 0 /
      data F           / neF * zero/
      data Fmul        / neF * zero/

      data ObjAdd      / zero /
      data ObjRow      / 1 /
      data neA, neG    / 0, 0 /
      data Prob        /'hs47    '/

      data iSpecs, iPrint, iSumm, iErr /0, 0, 6, 0/
      data         lPrint, lSumm       /   0, 0   /
*
*     End Data statements
*

*     ------------------------------------------------------------------
*     First,  snInit MUST be called to initialize optional parameters
*     to their default values.
*     ------------------------------------------------------------------
      call snInit
     &  ( iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw )

      call usrini( ObjAdd, ObjRow, Prob,
     &     x, bl(1), bu(1), xstate, Names,
     &     Fmul, bl(n+1), bu(n+1), Fstate, Fnames,
     &     iSpecs, iPrint, iSumm, iErr,
     &     cu, iu, ru, cw, iw, rw )


*     ------------------------------------------------------------------
*     Reset iPrint and iSumm from usrini.
*     ------------------------------------------------------------------
      call snseti
     &   ( 'iw 12 = ', iPrint, lPrint, lSumm, INFO,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call snseti
     &   ( 'iw 13 = ', iSumm , lPrint, lSumm, INFO,
     &     cw, lencw, iw, leniw, rw, lenrw )

      lenia          = leniw - nsnr
      lenra          = lenrw - nsnr
      lusdi          = 0
      lusdr          = 0

      call newnlp( iw(nsnr + 1), leniw - nsnr, lusdi, moreiw )
      if ( moreiw .ne. 0 ) goto 600

      call ininlp( x, n, bl(1), bu(1),
     &     bl(n+1), neF, bu(n+1),
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     iw(nsnr + 1), lenia, rw(nsnr + 1), lenra,
     &     lusdi, lusdr, moreiw, morerw )

      if (moreiw .ne. 0  .or.  morerw .ne. 0) goto 600
      lenia = lusdi
      lenra = lusdr

*     Get the copy of the Jacobian (it is a copy because snopt
*     scrambles it.

      call getnlp( linlp, lrnlp,
     &     iiadmx, liadmx, iradmx, lradmx,
     &     iiwork, liwork, irwork, lrwork,
     &     iw(nsnr + 1), lenia )

      call getmtx( iiAfun, neA, ijAvar, iA,
     &     iiGfun, neG, ijGvar, mm,
     &     iw(nsnr + iiadmx), liadmx, rw(nsnr + iradmx), lradmx )

*     Adjust iiAfun, etc. to be indices into iw and rw
      iiAfun = iiAfun + nsnr + iiadmx - 1
      ijAvar = ijAvar + nsnr + iiadmx - 1
      iiGfun = iiGfun + nsnr + iiadmx - 1
      ijGvar = ijGvar + nsnr + iiadmx - 1

      iA     = iA     + nsnr + iradmx - 1

      call snseti
     &   ( 'User real    workspace', lusdr + nsnr,
     &     lPrint, lSumm, INFO,
     &     cw, lencw, iw, leniw, rw, lenrw )

      call snseti
     &   ( 'User integer workspace', lusdi + nsnr,
     &     lPrint, lSumm, INFO,
     &     cw, lencw, iw, leniw, rw, lenrw )

*     ------------------------------------------------------------------
*     Read a Specs file (Optional).
*     ------------------------------------------------------------------
      if (iSpecs .gt. 0) then
         call snSpec
     &     ( iSpecs, INFO, cw, lencw, iw, leniw, rw, lenrw )

         if (INFO .ne. 101) then
             write(iErr, *) 'Trouble reading the Specs file'
             stop 1
         end if
      end if

      iCold = 0
      lenA  = max( neA, 1 )
      lenG  = max( neG, 1 )

      call snOptA
     &   ( iCold, neF, n, nName, nFname,
     &     ObjAdd, ObjRow, Prob, userfg,
     &     iw(iiAfun), iw(ijAvar), lenA, neA, rw(iA),
     &     iw(iiGfun), iw(ijGvar), lenG, neG,
     &     bl(1), bu(1), Names, bl(n+1), bu(n+1), FNames,
     &     x, xstate, xmul, F, Fstate, Fmul,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .ge. 10) goto 900

      stop

*     Begin error handling

 600  INFO  = 1
      miniw = leniw + moreiw
      minrw = lenrw + morerw

 900  if ( INFO .eq. 1 ) then

*        Memory allocation errors

         write (iErr, *) 'There is not enough storage ',
     &        'to start solving the problem'

         write (iErr, *) 'Total integer workspace should be ',
     &        'significantly more than ', miniw

         write (iErr, *) 'Total real    workspace should be ',
     &        'significantly more than ', minrw

         stop 1
      else
         write (iErr, *) Prob, ' terminated with INFO = ',
     &        INFO
         stop 2
      end if

      end

