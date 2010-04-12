c***********************************************************************
        program Moment
c
c= MOMENT -  Calculate moments of a Miriad image.
c& mchw
c: image analysis
c+
c       MOMENT calculates the nth moment of a Miriad image.  Currently
c       only the moment of axis 1 or axis 3 can be computed for three-
c       dimensional images.  To obtain the moments for other axes use
c       task REORDER to reorder the axes.
c@ in
c       The input image.  No default.
c@ region
c       The region of the input image to be used. See documentation for
c       region on how to specify this. Only the bounding box is
c       supported.
c@ out
c       The output image. No default.
c@ mom
c       -3      Velocity at peak, using a 3 point polynomial fit around
c               the peak
c       -2      Peak temperature (units are same as individual channels)
c       -1      Average intensity. (units are same as individual
c               channels)
c        0      Integrated intensity. (e.g. units I x km/s)
c        1      1st moment = sum(I*v)/sum(I)
c        2      dispersion = sqrt(sum(I*(v-M1)**2)/sum(I))
c               The v(i) are the axis values (e.g. velocities).  M1 is
c               the first moment.  The moment is calculated
c               independently for each pixel within the specified clip
c               range.  Default = 0.
c               Note that mom=2 produces a map of sigma ('gaussian'
c               dispersions), with FWHM = 2*sqrt(2*ln(2)) = 2.355*sigma
c               if the profile was a true gaussian.
c@ axis
c       The axis for which the moment is calculated. Default = 3.
c@ clip
c       Two values.  Exclude pixels with values in the range clip(1) to
c       clip(2) inclusive.  If only one value is given, then exclude
c       -abs(clip) to +abs(clip).
c@ rngmsk
c       Mask pixels as bad when the 1st moment gives a 1st moment out of
c       range of that of the axis. This option can only be activated for
c       mom=1. Default: false
c--
c   Other things to be improved in this code:
c
c       - when items such as cdeltX, crvalX are missing, assume
c         reasonable defaults.  I recently fixed the ctypeX problem.
c         (X = the moment axis)
c
c  History:
c    26may89 Robert Loushin  original version - v-axis of vxy image
c                  only.
c    25jun90 mchw  Major update and rework to allow for other image
c                  axes.
c                       Fixed a bug for 2nd and higher moments.
c    13nov90 mchw  Added suport for xyv images, and pixel blanking.
c                       Made 2nd moment work and removed higher moments.
c    08feb91 mchw  Two values for clip range, and version in history.
c    25feb91 mjs   Changed references of itoa to itoaf.
c    01aug91 rjs   Use boxSet(...'s')
c    05aug91 alr/pjt Changed clip from include to exclude range.
c                       Write out dummy 3rd. axis.
c    06aug91 mchw       Merge in above changes.
c    11jul92 pjt   Fixed rather serious flagging bug when axis=1
c    26nov92 nebk  Add btype
c     4jan93 pjt   Fixed rather serious flagging bug when axis=3
c     5feb93 rjs   Use memalloc.
c     5apr93 pjt   Since rjs is walkabouting, I had to fix that amazing
c                  bug: indexor's (i0,j0) in moment3 wrong if blc,trc
c                  used
c    04jan96 nebk  Write header descriptors as double precision
c    26nov96 rjs   Increase length of input and output file names.
c    07mar97 rjs   Improve a message.
c    30jun97 rjs   Change units of 0th moment image to be jy/beam.km/s
c    23jul97 rjs   Added pbtype.
c    27feb98 mwp   Added mom=-2 for peak temperature
c    13jul00 rjs   Copy across llrot keyword.
c    12feb01 pjt   Mask pixels if their velocity is out of range
c    15feb01 dpr   Truncate rangemask key to rngmsk - doh!
c    16feb01 pjt   Added mom=-3 for velocity of peak fit to poly=2
c    18jan02 pjt   Turned rngmask typo into rngmsk (duh)
c     4mar02 pjt   documented FWHM/sigma, fixed units of mom=2 map
c
c $Id$
c-----------------------------------------------------------------------
      include 'maxdim.h'
      integer maxnax,maxboxes,maxruns,naxis
      parameter(maxnax=3,maxboxes=2048)
      parameter(maxruns=3*maxdim)

      logical rngmsk
      integer axis, blc(maxnax), boxes(maxboxes), i, j, lIn, lOut,
     :          mom, nsize(maxnax), size(maxnax), trc(maxnax)
      real    blo, bhi, clip(2)
      character in*64, line*72, out*64, version*80
c
c  External
c
      logical keyprsnt
      character versan*80
c-----------------------------------------------------------------------
      version = versan ('moment',
     :                  '$Revision$',
     :                  '$Date$')

c
c Get inputs.
c
      call keyini
      call keya('in',in,' ')
      call BoxInput('region',in,boxes,maxboxes)
      call keya('out',out,' ')
      call keyi('mom',mom,0)
      call keyi('axis',axis,3)
      call keyr('clip',clip(1),0.0)
      if(keyprsnt('clip'))then
        call keyr('clip',clip(2),0.0)
      else
        clip(2) = abs(clip(1))
        clip(1) = -clip(2)
      endif
      call keyl('rngmsk',rngmsk,.FALSE.)
      call keyfin
c
c Check inputs.
c
      if(in .eq. ' ') call bug('f','No input specified. (in=)')
      if(out .eq. ' ') call bug('f','No output specified. (out=)')
      if(mom.lt.-3 .or. mom.gt.2)
     :   call bug('f','moment must be between -3 and 2')
      if(clip(2).lt.clip(1)) call bug('f','clip range out of order')
      call xyopen(lIn,in,'old',maxnax,nsize)
      call rdhdi(lIn,'naxis',naxis,0)
      naxis = min(naxis,maxnax)
      if(nsize(1).gt.maxdim)call bug('f','Input file too big for me')
      if(axis.ne.1.and.axis.ne.3) then
        call output('Reorder axes to take moment of axis 2')
        call bug('f','axis 2 is not implemented.')
      endif
c
c  Determine the min and max value.
c
      call ImMinMax(lIn,naxis,nsize,blo,bhi)
      if(blo.eq.bhi)then
        call xyclose(lIn)
        write(line,'(''All pixels are '',1pg10.3)')blo
        call output(line)
        stop
      endif
c
c  Set up the region of interest. Warn if not rectangle.
c
c     call BoxMask(lIn,boxes,maxboxes)          - not yet used
      call BoxSet(boxes,maxnax,nsize,'s')
      call BoxInfo(boxes,maxnax,blc,trc)
c
c  Open output image and fill in its header.
c
      j = 0
      do i = 1,naxis
        if(blc(i).lt.1) blc(i) = 1
        if(trc(i).gt.nsize(i)) trc(i) = nsize(i)
        if(i.ne.axis) then
          j = j + 1
          size(j) = trc(i) - blc(i) + 1
        endif
      enddo
      size(naxis) = 1
      call xyopen(lOut,out,'new',naxis,size)
      call header(lIn,lOut,naxis,blc,trc,mom,axis)
c
c  Calculate moment.
c
      call makemom(lIn,lOut,naxis,blc,trc,mom,axis,clip,rngmsk)
c
c  Update history and close files.
c
      call Hisopen(lOut,'append')
      call HisWrite(lOut,'MOMENT: '//version)
      call HisInput(lOut,'MOMENT')
      call HisClose(lOut)
      call xyclose(lIn)
      call xyclose(lOut)
      end
c***********************************************************************
      subroutine header(lIn,lOut,naxis,blc,trc,mom,axis)
      integer lIn,lOut,naxis,blc(naxis),trc(naxis),mom,axis
c
c  Copy keywords to output file.
c
c  Inputs:
c    lIn,lOut   Handle of input and output files.
c    naxis      The number of input axes.
c    blc,trc    The corners of the input image.
c    mom        The moment to be calculated.
c    axis       The axis for which the moment is calculated.
c-----------------------------------------------------------------------
      character atemp*16,ctype*10,cin*2,cout*2
      double precision rtemp,cdelt,crval,crpix,idx
      integer i,j,k,l
c
c  Externals.
c
      character itoaf*2
      integer len1
      logical hdprsnt
c
c  Be careful that nkeys and nckeys match the number of keywords.
c
      integer nkeys,nckeys
      parameter (nkeys=23, nckeys=4)
      character keyw(nkeys)*8, ckeyw(nckeys)*5

      data keyw/   'bmaj    ','bmin    ','bpa     ',
     :    'obstime ','epoch   ','history ','llrot   ',
     :    'ltype   ','lstart  ','lstep   ','lwidth  ','pbfwhm  ',
     :    'instrume','niters  ','object  ','telescop','pbtype  ',
     :    'restfreq','vobs    ','observer','obsra   ',
     :    'obsdec  ','mostable'/
c
c  Keyword values which must be changed as they are passed from in to
c  out.
c
      data ckeyw/'ctype','cdelt','crval','crota'/
c
c  Copy across unchanged header keywords.
c
      do i = 1,nkeys
        call hdcopy(lIn,lOut,keyw(i))
      enddo
c
c  Handle the keywords which must be moved to another axis.
c
      j = 0
      do i = 1,naxis
        if(i.ne.axis) then
          j = j + 1
          cin = itoaf(i)
          cout = itoaf(j)
          atemp = ckeyw(1)//cin
          call rdhda(lIn,ckeyw(1)//cin,atemp,' ')
          if(atemp.ne.' ')call wrhda(lOut,ckeyw(1)//cout,atemp)
          do k = 2,nckeys
            atemp = ckeyw(k)//cin
            if(hdprsnt(lIn,atemp)) then
              call rdhdd(lIn,ckeyw(k)//cin,rtemp,0.0d0)
              call wrhdd(lOut,ckeyw(k)//cout,rtemp)
            endif
          enddo
c
c  Special cases: the crpixes will change if the user uses a subcube.
c
          if(hdprsnt(lIn,'crpix'//cin)) then
            call rdhdd(lIn,'crpix'//cin,rtemp,0.0d0)
            rtemp = rtemp - dble(blc(i)) + 1
            call wrhdd(lOut,'crpix'//cout,rtemp)
          endif
        endif
      enddo
c
c  bunit is usually a pain.  For this program it is relatively simple.
c
      if(mom.le.-1)then
        call hdcopy(lIn,lOut,'bunit')
      else if(mom.eq.0)then
        call rdhda(lIn,'bunit',atemp,' ')
        l = len1(atemp)
        if(l.gt.0)then
          atemp(l+1:) = '.KM/S'
          call wrhda(lOut,'bunit',atemp)
        endif
      else if(mom.eq.1)then
        call wrhda(lOut,'bunit','KM/S')
        call wrbtype(lOut,'velocity')
      else
        call wrhda(lOut,'bunit','KM/S')
        call wrbtype(lOut,'velocity_dispersion')
      endif
c
c  Write out additional information about the ``third'' dummy axis.
c  also the third dimension is 1 (naxis3=1), the axes are labeled
c  as much as possible from the input cube
c
      cin = itoaf(axis)
      idx = (blc(axis)+trc(axis))/2.0
      call rdhdd(lIn,'crpix'//cin, crpix, 1d0)
      call rdhdd(lIn,'cdelt'//cin, cdelt, 1d0)
      call rdhdd(lIn,'crval'//cin, crval, 0d0)
      call rdhda(lIn,'ctype'//cin, ctype, ' ')
      crval = crval + (idx-crpix)*cdelt
      cdelt = cdelt * (trc(axis)-blc(axis)+1)
      call wrhdd(lOut,'crpix3', 1d0)
      call wrhdd(lOut,'cdelt3', cdelt)
      call wrhdd(lOut,'crval3', crval)
      call wrhda(lOut,'ctype3', ctype)
      end
c********1*********2*********3*********4*********5*********6*********7*c
      subroutine makemom(lIn,lOut,naxis,blc,trc,mom,axis,clip,rngmsk)
      integer lIn,lOut,naxis,blc(naxis),trc(naxis),mom,axis
      real clip(2)
      logical rngmsk
c
c  Calculate the scale factor and channel offset for moments.
c
c  Inputs:
c    lIn,lOut   Handle of input and output files.
c    naxis      The number of input axes.
c    blc,trc    The corners of the input image.
c    mom        The moment to be calculated. Default = 0.
c    axis       The axis for which the moment is calculated.
c    clip       Pixel values in range clip(1) to clip(2) are excluded.
c    rngmsk     Also set pixels flagged if mom=1 and outside axis range
c-----------------------------------------------------------------------
      include 'maxdim.h'
      include 'mirconst.h'
      include 'mem.h'
      integer n1, n2, pFlags, pSum, pTemps
      real    cdelt, crpix, crval, offset, restfreq, scale
      character cin*1, ctype*9
c
c  Externals.
c
      logical hdprsnt
      character itoaf*1
c-----------------------------------------------------------------------
c     Calculate offset and scale to convert from channels to km/sec.
      cin = itoaf(axis)
      if (hdprsnt(lIn, 'crpix'//cin) .and.
     :    hdprsnt(lIn, 'crval'//cin)) then
        call rdhdr(lIn,'crpix'//cin, crpix, 0.0)
        call rdhdr(lIn,'crval'//cin, crval, 0.0)
      else
        call bug('f', 'crpix and/or crval not in header.')
      endif

      call rdhdr(lIn, 'cdelt'//cin, cdelt, 0.0)
      if (cdelt.eq.0.0) call bug('f','cdelt is 0 or not present.')

      call rdhda(lIn, 'ctype'//cin, ctype, ' ')
      if (ctype(1:4).ne.'VELO' .and.
     :    ctype(1:4).ne.'FELO' .and.
     :    ctype(1:4).ne.'FREQ')
     :      call bug('w','Axis is not VELO, FELO, or FREQ.')

      offset = crpix - crval/cdelt
      if (mom.eq.-2) then
c       Peak channel value.
        scale  = 1.0
      else if (mom.eq.-1) then
c       Average channel value.
        scale  = 1.0 / real(trc(axis)-blc(axis)+1)
      else
c       Velocity needed.
        if (ctype(1:4).eq.'FREQ') then
          call rdhdr(lIn, 'restfreq', restfreq, 0.0)
          if (restfreq.eq.0.0) then
            call bug('f','restfreq not present in header.')
          endif

          offset = crpix - (crval - restfreq) / cdelt
          scale  = -CMKS * (cdelt / restfreq) * 1e-3
        else
          scale  = cdelt
        endif

        if (mom.eq.0) scale = abs(scale)
      endif


c     Compute the moment.
      if (axis.eq.1) then
        call moment1(lIn,lOut,naxis,blc,trc,mom,scale,offset,clip,
     :               rngmsk)
      else if (axis.eq.2) then
        call bug('f','axis 2 is not implemented.')
      else if (axis.eq.3) then
        n1 = trc(1) - blc(1) + 1
        n2 = trc(2) - blc(2) + 1
        call memalloc(pSum,3*n1*n2,'r')
        call memalloc(pFlags,n1*n2,'l')
        call memalloc(pTemps,n1*n2,'r')
        call moment3(lIn,lOut,naxis,blc,trc,mom,scale,offset,clip,
     :    memr(pSum),meml(pFlags),n1,n2,memr(pTemps),rngmsk)
        call memfree(pFlags,n1*n2,'l')
        call memfree(pSum,3*n1*n2,'r')
        call memfree(pTemps,n1*n2,'r')
      endif

      end
c********1*********2*********3*********4*********5*********6*********7*c
      subroutine moment1(lIn,lOut,naxis,blc,trc,mom,scale,offset,clip,
     :                   rngmsk)
      integer lIn,lOut,naxis,blc(naxis),trc(naxis),mom
      real scale,offset,clip(2)
      logical rngmsk
c
c  Calculate the moment for axis 1.
c
c  Inputs:
c    lIn,lOut   Handle of input and output files.
c    naxis      The number of input axes.
c    blc,trc    The corners of the input image.
c    mom        The moment to be calculated. Default = 0.
c    scale      Scale factor to convert from channels to km/s
c    offset     Offset in channels.
c    clip       Pixels with values in range clip(1:2) are excluded.
c    rngmsk     Also set pixels flagged if mom=1 and outside axis range
c-----------------------------------------------------------------------
      include 'maxdim.h'
      real buf(maxdim),sum(maxdim),sum1(maxdim),sum2(maxdim)
      real chan,flux,sigsq,vmin,vmax,vpeak,dpeak
      integer i,j,k,ipeak
      logical flags(maxdim), outflags(maxdim)
c-----------------------------------------------------------------------
      vmin = (blc(1) - offset)*scale
      vmax = (trc(1) - offset)*scale
      write(*,*) 'Velocity range: ',vmin,vmax,scale

      do k = blc(3), trc(3)
        call xysetpl(lIn,1,k)
        do j = blc(2), trc(2)
          call xyread(lIn,j,buf)
          call xyflgrd(lIn,j,flags)

          sum(j)  = 0.0
          sum1(j) = 0.0
          sum2(j) = 0.0

c
c  Loop through the velocity channels and accumulate the moments.
c  (note tricky dpeak initialization)
c
          do i = blc(1), trc(1)
            if (flags(i).and.
     :           (buf(i).lt.clip(1) .or. clip(2).lt.buf(i))) then
              chan = (i - offset)*scale
              if (mom.eq.-3 .and.
     :           (sum(j).eq.0.0 .or. buf(i).gt.dpeak)) then
                 dpeak = buf(i)
                 vpeak = chan
                 ipeak = i
              endif

              sum(j) = sum(j) + buf(i)
              if (mom.ge.1) then
                sum1(j) = sum1(j) + buf(i)*chan
                if (mom.eq.2) then
                  sum2(j) = sum2(j) + buf(i)*chan*chan
                endif
              endif
            endif
          enddo
c
c  Normalize and scale the moments.
c
          flux = sum(j)
          if (flux.ne.0.0) then
            sum(j) = sum(j) * scale
            if (mom.ge.1) sum1(j) = sum1(j)/flux
            outflags(j) = .true.
            if (mom.eq.2) then
              sigsq = sum2(j)/flux - sum1(j)*sum1(j)
              if (sigsq.gt.0.0) then
                sum2(j) = sqrt(sigsq)
              else
                sum2(j) = 0.0
                flags(j) = .false.
              endif
            endif
            if (rngmsk) then
               if (vmin.lt.vmax. and.
     :              (sum1(j).lt.vmin .or. sum1(j).gt.vmax) .or.
     :             vmin.gt.vmax. and.
     :              (sum1(j).lt.vmax .or. sum1(j).gt.vmin)) then
                  flags(j) = .false.
                  sum(j)  = 0.
                  sum1(j) = 0.
                  sum2(j) = 0.
               endif
            endif
            if (mom.eq.-3) then
               sum(j) = vpeak + 0.5*(buf(ipeak-1)-buf(ipeak+1))*
     :                  scale/(buf(ipeak-1)+buf(ipeak+1)-2*buf(ipeak))
            endif
          else
            sum1(j) = 0.0
            sum2(j) = 0.0
            outflags(j) = .false.
          endif
        enddo
c
c  Now write out the moment map.
c
        if (mom.le.0) then
          call xywrite(lOut,k-blc(3)+1,sum(blc(2)))
        else if (mom.eq.1) then
          call xywrite(lOut,k-blc(3)+1,sum1(blc(2)))
        else if (mom.eq.2) then
          call xywrite(lOut,k-blc(3)+1,sum2(blc(2)))
        endif
        call xyflgwr(lOut,k-blc(3)+1,outflags(blc(2)))
      enddo

      end
c********1*********2*********3*********4*********5*********6*********7*c
      subroutine moment3(lIn,lOut,naxis,blc,trc,mom,scale,offset,clip,
     :  sum,outflag,n1,n2,pTemps,rngmsk)

      integer lIn,lOut,naxis,blc(naxis),trc(naxis),mom,n1,n2
      real scale,offset,clip(2),sum(n1,n2,3)
      logical outflag(n1,n2)
      real pTemps(n1,n2)
      logical rngmsk
c
c  Calculate the moments for axis 3.
c
c  Inputs:
c    lIn,lOut   Handle of input and output files.
c    naxis      The number of input axes.
c    blc,trc    The corners of the input image.
c    mom        The moment to be calculated. Default = 0.
c    scale      Scale factor to convert from channels to km/s
c    offset     Offset in channels.
c    clip       Pixels with value in range clip(1),clip(2) are excluded.
c    rngmsk     Also set pixels flagged if mom=1 and outside axis range
c  Scratch:
c    sum        Used to accumulate the moments.
c    outflag    Flagging info for the output array.
c    pTemps     use to calculate max temperature (flux) map if mom=-2
c-----------------------------------------------------------------------
      include 'maxdim.h'

      logical flags(maxdim)
      integer i, i0, j, j0, k
      real    buf(maxdim), chan, flux, sigsq, vmax, vmin, vtmp
c-----------------------------------------------------------------------
c
c  Check consistency.
c
      if (trc(1)-blc(1)+1.ne.n1 .or. trc(2)-blc(2)+1.ne.n2) then
        call bug('f', 'Dimension inconsistency in MOMENT3')
      endif
      if (mom.eq.-3) then
        call bug('f', 'mom=-3 not yet available for axis=3; ' //
     :                'try: reorder in= out= mode=312')
      endif
c
c  Intialize the max temperature array
c    (should really use some kind of POSIX-type MINFLOAT here)
      do j = 1, n2
        do i = 1, n1
          pTemps(i,j) = -1e38
        enddo
      enddo
c
c  Zero the array to accumulate the moments.
c
      do k = 1, 3
        do j = 1, n2
          do i = 1, n1
            sum(i,j,k) = 0.0
          enddo
        enddo
      enddo
c
c  Loop through the velocity channels.
c
      if (rngmsk) then
        vmin = (blc(3) - offset)*scale
        vmax = (trc(3) - offset)*scale

        if (vmax.lt.vmin) then
          vtmp = vmin
          vmin = vmax
          vmax = vtmp
        endif

        write(*,*) 'moment3: rangemask ', vmin, vmax
      endif

      do k = blc(3), trc(3)
        call xysetpl(lIn,1,k)
        chan = (k - offset)*scale

c       Accumulate the moments, one row at a time.
        j0 = 1
        do j = blc(2), trc(2)
          call xyread(lIn,j,buf)
          call xyflgrd(lIn,j,flags)

          i0 = 1
          do i = blc(1), trc(1)
            if (flags(i) .and.
      :          (buf(i).lt.clip(1) .or. buf(i).gt.clip(2))) then
              if (mom.eq.-2) then
                if (buf(i).gt.pTemps(i,j)) then
                  pTemps(i,j)  = buf(i)
                  sum(i0,j0,1) = pTemps(i,j)
                endif
              else
                sum(i0,j0,1) = sum(i0,j0,1) + buf(i)
              endif

              if (mom.ge.1) then
                sum(i0,j0,2) = sum(i0,j0,2) + buf(i)*chan

                if (mom.ge.2) then
                  sum(i0,j0,3) = sum(i0,j0,3) + buf(i)*chan*chan
                endif
              endif
            endif

            i0 = i0 + 1
          enddo

          j0 = j0 + 1
        enddo
      enddo
c
c  Normalize and scale the moments.
c
      do j = 1,n2
        do i = 1,n1
          flux = sum(i,j,1)
          if (flux.ne.0.0) then
            sum(i,j,1) = sum(i,j,1) * scale
            if (mom.ge.1) sum(i,j,2) = sum(i,j,2)/flux

            outflag(i,j) = .true.
            if (mom.eq.2) then
              sigsq = sum(i,j,3)/flux - sum(i,j,2)*sum(i,j,2)
              if (sigsq.gt.0.) then
                sum(i,j,3) = sqrt(sigsq)
              else
                sum(i,j,3) = 0.0
                outflag(i,j) = .false.
              endif
            endif

            if (rngmsk) then
              if ((sum(i,j,2).lt.vmin .or. vmax.lt.sum(i,j,2))) then
                outflag(i,j) = .false.
                sum(i,j,1) = 0.0
                sum(i,j,2) = 0.0
                sum(i,j,3) = 0.0
              endif
            endif
          else
            sum(i,j,2) = 0.
            sum(i,j,3) = 0.
            outflag(i,j) = .false.
          endif
        enddo
      enddo
c
c  Now write out the moment map.
c
      k = max(mom+1,1)
      do j = 1,n2
        do i = 1,n1
          buf(i) = sum(i,j,k)
          flags(i) = outflag(i,j)
        enddo

        call xywrite(lOut,j,buf)
        call xyflgwr(lOut,j,flags)
      enddo

      end
