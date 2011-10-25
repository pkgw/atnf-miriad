c***********************************************************************
c* Title -  Write title in standard format into LogFile.
c& mchw
c: image analysis,log-file
c+
c $Id$

      subroutine Title(lIn,naxis,blc,trc,cbof)

      integer lIn,naxis,blc(naxis),trc(naxis)
c  ---------------------------------------------------------------------
c  Write title in standard format into LogFile.
c
c  Inputs:
c    lIn        The handle of the Image.
c    naxis      Number of axes of the Image.
c    blc,trc    Corners of region of interest.
c  Output:
c    cbof       Beam oversampling factor.
c-----------------------------------------------------------------------
      include 'mirconst.h'

      logical   more
      integer   lblc, ltrc
      real      bmaj, bmin, cbof, ifrq, KperJy, rfreq, omega
      character algo*3, bunit*20, ctype1*9, ctype2*9, line*80,
     *          txtblc*20, txttrc*20

      external  itoaf, len1
      integer   len1
      character itoaf*2
c-----------------------------------------------------------------------
c     Get beam size etc.
      call GetBeam(lIn,naxis,bunit,bmaj,bmin,omega,cbof)

c     Get the rest frequency.
      call rdhdr(lIn, 'restfreq', rfreq, 0.0)
      if (rfreq.eq.0.0) then
c       Rest frequency not recorded, use the frequency reference value.
        call coInit(lIn)
        call coSpcSet(lIn, 'FREQ', ' ', ifrq, algo)
        if (ifrq.gt.0) then
          call rdhdr(lIn, 'crval'//itoaf(ifrq), rfreq, 0.0)
        endif
        call coFin(lIn)
      endif

      KperJy = 1.0
      if (rfreq.ne.0.0) then
        KperJy = (0.3/rfreq)**2 / (2.0*1.38e3*omega)
      endif

      write(line,'(a,a,a,f10.6,a,f15.4)')
     *      '  bunit: ', bunit(1:len1(bunit)),
     *      '  frequency: ', rfreq,
     *      '    K/Jy: ', KperJy
      call LogWrite(line,more)

      call rdhda(lIn,'ctype1',ctype1,' ')
      call rdhda(lIn,'ctype2',ctype2,' ')
      write(line,'(a,a,a,a,a,f6.2,a,f6.2,a)')
     *  '  Axes: ',ctype1,' x ',ctype2,
     *  '  Beam: ',bmaj*DR2AS,' x ',bmin*DR2AS,' arcsec'
      call LogWrite(line,more)

      call mitoaf(blc,3,txtblc,lblc)
      call mitoaf(trc,3,txttrc,ltrc)
      line = '  Bounding region is Blc=('//txtblc(1:lblc)//
     *                          '),Trc=('//txttrc(1:ltrc)//')'
      call LogWrite(line,more)

      end
