      program imom

c= IMOM - Compute intensity-weighted moment of an image
c& vjm
c: image analysis
c+
c       IMOM finds the 'centroid' of a Miriad image, and the 'spread'
c       about that centroid, by evaluating the expressions
c         mom1 = SUM{I(x)*x} / SUM{I(x)}
c         mom2 = sqrt(SUM{I(x)*(x-mom1)**2} / SUM{I(x)})
c       along the first two axes of an image.
c
c       Notes:
c       * The first moment should always be in the range 1..naxis(n).
c       * mom2 is equivalent to the standard deviation, not the FWHM,
c         for a gaussian intensity distribution.
c       * mom2 shows the spread along the directions parallel to the
c         *image axes*. So if your object of interest is elongated along
c         some random angle, the mom2 values returned will be a
c         combination of the major and minor-axis sizes, projected onto
c         the x- and y-axes.
c       * If all pixels in the region are blank, moments of -1.0 are
c         returned and a warning is printed.
c       * In case of a rounding error etc, moments of 0.0 are returned
c         and a warning is printed.
c
c@ in
c       Input image name. No default.
c@ options
c       The following options can be given (minimum match applies):
c         skew
c           Compute the 'third moment'. Since the formal third moment
c           has the rather useless dimensions of pixel**3, what
c           is actually returned is the dimensionless quantity
c                SUM{ I(x)*(x - mom1)**3 } / SUM{ I(x) } / mom2**3.
c           A negative value means the distribution is skewed towards
c           values less than the first moment; i.e. the pixels at
c           lower pixel number than mom1 are generally brighter than
c           those at higher pixel number. A positive value
c           means the pixels at higher pixel number than the first
c           moment are brighter.
c         clipmean
c           When you use this program to estimate the peak of a source,
c           better accuracy can be had by ignoring low-intensity values.
c           Setting this option forms the moments using only pixels that
c           are brighter than the (local) mean.
c           If a min or max value is specified (see below), the mean is
c           computed using only the values within the selected range.
c         clip1sigma
c           Forms moments using just the pixels brighter than
c           local_mean + local_rms.
c           If a min or max value is specified, the mean and rms are
c           computed using only the values within the selected range.
c@ region
c       The region of interest. See the help on "region" for more
c       information. This gives the region of pixels to calculate the
c       moments within. Only rectangular regions are supported.
c       The default is the entire image.
c       If a cube is given as the input image, the program computes
c       moments for each plane along axis 3 (within the selected range).
c@ min
c       The minimum pixel value for inclusion in the moment calculation.
c       The default is the image minimum.
c@ max
c       The maximum pixel value for inclusion in the moment calculation
c       The default is the image maximum.
c@ log
c       The output log file. The default is the terminal. Output to
c       the terminal is paged.
c
c$Id$
c--
c  History:
c    Refer to the RCS log, v1.1 includes prior revision information.
c-----------------------------------------------------------------------
      include 'maxdim.h'

      integer   MAXBOXES, MAXNAX
      parameter (MAXBOXES=2048, MAXNAX=3)

      logical   clip1sig, clipmean, doskew, getmax, getmin
      integer   blc(MAXNAX), boxes(MAXBOXES), lin, naxis, nsize(MAXNAX),
     *          trc(MAXNAX)
      real      imhi, imlo, rhi, rlo
      character in*64, line*80, out*64, version*72

      external  keyprsnt, versan
      logical   keyprsnt
      character versan*72
c-----------------------------------------------------------------------
c     Get input parameters.
      version = versan('imom',
     *                 '$Revision$',
     *                 '$Date$')

      call keyini
      call keya('in', In, ' ')
      if (in.eq.' ') call bug('f', 'Image name not specified')

      call GetOpt(doskew,clipmean,clip1sig)

      call BoxInput('region',In,boxes,MAXBOXES)

      getmin = .not.keyprsnt('min')
      getmax = .not.keyprsnt('max')
      call keyr('min',rlo,0.0)
      call keyr('max',rhi,0.0)
      call keya('log',out,' ')
      call keyfin

c     Open the output text file.
      call LogOpen(out,'q')

c     Open input image and check dimensions.
      call xyopen(lIn,In,'old',MAXNAX,nsize)
      if (nsize(1).gt.maxbuf) call bug('f','Image too big for buffer')
      call rdhdi(lIn,'naxis',naxis,0)
      naxis = min(naxis,MAXNAX)

c     Determine portion of image to list.
      call BoxSet(boxes,MAXNAX,nsize,'s')
      call BoxInfo(boxes,MAXNAX,blc,trc)

c     Determine the min and max value, if needed.
      call ImMinMax(lIn,naxis,nsize,imlo,imhi)
      if (getmin .or. getmax) then
        if (imlo.eq.imhi) then
          call xyclose(lIn)
          write(line,'(''All pixels are '',1pg10.3)') rlo
          call output(line)
        endif

c       rlo,rhi are set to user's value if one is given, else 0.0.
c       Reset them to image low, high if user gave no value.
        if (getmin) rlo=imlo
        if (getmax) rhi=imhi

        if (rlo.gt.rhi) then
          call xyclose(lIn)
          write(line,
     *        '(''Minimum '',1pg10.3,'' greater than max '',g10.3)')
     *         rlo,rhi
          call bug('f',line)
        endif
      endif

c     Title line.
      call LogWrit('Moments for image '//In)

      if (clip1sig) then
        call LogWrit('Will clip pixels below (RoI mean + 1 sigma)')
      else
        if (clipmean) then
          call LogWrit('Will clip pixels below the RoI mean')
        endif
      endif

c     Compute
      call ListMom(lIn,naxis,imlo,imhi,blc,trc,rlo,rhi,
     *             getmin,getmax,doskew,clipmean,clip1sig)

c     All done.
      call xyclose(lIn)
      call LogClose

      end

c***********************************************************************

      subroutine GetOpt(doskew,clipmean,clip1sig)

      logical   doskew, clipmean, clip1sig
c-----------------------------------------------------------------------
c  Determine which of the options is to be done. Default is dodata.
c
c  Outputs:
c    doskew     Compute the skew, or not?
c-----------------------------------------------------------------------
      integer   NOPTS
      parameter (NOPTS=3)

      logical   present(NOPTS)
      character opts(NOPTS)*10

      data opts/'skew      ','clipmean  ','clip1sigma'/
c-----------------------------------------------------------------------
      call options('options',opts,present,NOPTS)

      doskew   = present(1)
      clipmean = present(2)
      clip1sig = present(3)

      end

c***********************************************************************

      subroutine ListMom(lIn,naxis,imlo,imhi,blc,trc,rlo,rhi,
     *                   getmin,getmax,doskew,clipmean,clip1sig)
      logical getmin,getmax,doskew,clipmean,clip1sig
      integer lIn,naxis,blc(naxis),trc(naxis)
      real    rlo,rhi,imlo,imhi
c-----------------------------------------------------------------------
c  List moments and write in standard format into LogFile.
c
c      Inputs:
c        lIn        The handle of the Image.
c        naxis      Number of axes of the Image.
c        imlo,imhi  Min, max values in the entire image
c        blc,trc    Corners of region of interest.
c        rlo,rhi    Intensity range to use
c        getmin     Find the minimum in each region
c        getmax     Find the maximum in each region
c                   (getmin or getmax=false means user has set it)
c        doskew     Compute the skew statistic
c        clipmean   Set lower clip limit to region mean
c        clip1sig   Set lower clip limit to region mean+rms
c-----------------------------------------------------------------------
      integer   axis, length, plane
      real      cbof, momx(3), momy(3), stats(3)
      character cax*2, ctype*9, line*80

      external  itoaf
      character itoaf*2
c-----------------------------------------------------------------------
c     Write title lines.
      call LogWrit(' ')
      call LogWrit('Image Moments')
      call Title(lIn,naxis,blc,trc,cbof)

c     List moments for each plane in each axis in selected region.
      if (naxis.ge.3) then
        axis = 3
        do while (axis.le.naxis)
          cax = itoaf(axis)
          call rdhda(lIn, 'ctype'//cax, ctype, ' ')
          write(line,'(a,i5,5x,a)') 'Axis: ',axis,ctype
          length=6+5+5+10
          call LogWrit(' ')
          call LogWrit(line(1:length))

c         Do a separate estimate for each hyperplane.
          plane = blc(axis)
          do while (plane.le.trc(axis))
            call xysetpl(lIn,1,plane)

            call momit(lin,naxis,blc,trc,imlo,imhi,rlo,rhi,
     *                 getmin,getmax,doskew,clipmean,clip1sig,
     *                 momx,momy,stats)

            write(line,'(a,x,i4,2x,a10,x,1pe12.5,1p,2(x,a,x,e12.5))')
     *        'Plane:',plane,'Npix:',stats(1),'good, range',stats(2),
     *        'to',stats(3)
            call LogWrit(line(1:80))

            write(line,'(a,x,i4,2x,a10,1p,2(x,e12.5),x,a)')
     *        'Plane:',plane,'Centroid:',momx(1),momy(1),'pixels'
            call LogWrit(line(1:80))

            write(line,'(a,x,i4,2x,a10,1p,2(x,e12.5),x,a)')
     *        'Plane:',plane,'Spread:',momx(2),momy(2),'pixels'
            call LogWrit(line(1:80))

            if (doskew) then
              write(line,'(a,x,i4,2x,a10,1p,2(x,e12.5),x,a)')
     *          'Plane:',plane,'Skew:',momx(3),momy(3),' '
              call LogWrit(line(1:80))
            endif
            plane = plane + 1
          enddo

          axis = axis + 1
        enddo

      else if (naxis.eq.2) then
        call momit(lin,naxis,blc,trc,imlo,imhi,rlo,rhi,getmin,getmax,
     *             doskew,clipmean,clip1sig,momx,momy,stats)

        write(line,'(x,a10,x,1pe12.5,2x,a,x,1pe12.5,x,a,x,1pe12.5))')
     *    'Npix:',stats(1),'good, range',stats(2),'to',stats(3)
        call LogWrit(line(1:80))

        write(line,'(x,a10,1p,2(x,e12.5),x,a)')
     *    'Centroid:',momx(1),momy(1),' pixels'
        call LogWrit(line(1:80))

        write(line,'(x,a10,1p,2(x,e12.5),x,a)')
     *    'Spread:',momx(2),momy(2),' pixels'
        call LogWrit(line(1:80))

        if (doskew) then
          write(line,'(x,a10,1p,2(x,e12.5),x,a)')
     *      'Skew:',momx(3),momy(3),' '
          call LogWrit(line(1:80))
        endif
      endif

      end

c***********************************************************************

      subroutine momit(lin,naxis,blc,trc,imlo,imhi,rlo,rhi,
     *                 getmin,getmax,doskew,clipmean,clip1sig,
     *                 momx,momy,stats)

      integer lIn, naxis, blc(naxis), trc(naxis)
      real    imlo, imhi, rlo, rhi
      logical getmin, getmax, doskew, clipmean, clip1sig
      real    momx(3), momy(3), stats(3)
c-----------------------------------------------------------------------
c  Compute image moments within a bounding box.
c
c  Inputs:
c    lIn        The handle of the Image.
c    naxis      Number of axes of the Image.
c    imhi,imlo  Max,min values in image (from image header)
c    blc,trc    Corners of region of interest.
c    rlo,rhi    Intensity range of interest.
c    getmin     Find the minumum value in the region.
c    getmax     Find the maximum value in the region.
c    doskew     Compute the skew statistic
c    clipmean   Set low end of clip range = mean of subregion
c    clip1sig   Set low end of clip range = mean+rms of subregion
c
c  Output:
c    momx       moments along "x" axis, of unflagged pixels within box.
c    momy       moments along "y" axis, of unflagged pixels within box.
c    stats      (no. pixels, lowest value, highest value) used in the
c               calculation of the moments.
c
c  Note:
c    a moment < 0 means no unflagged pixels
c    a moment = 0 means the sum of intensities (used as weights) = 0,
c                 or there was a rounding error (mom(2)=0).
c
c Doesn't yet handle a user-fixed min or max; it recomputes rlo,rhi.
c-----------------------------------------------------------------------
      include 'maxdim.h'

      logical   apply, flags(MAXBUF)
      integer   i, j, num
      real      avg, epsilon, nrmi, nrmj, offset, ri, rj, rms, sum,
     *          sumcbi, sumcbj, sumi, sumj, sumrow, sumsq, sumsqi,
     *          sumsqj, val, vmax, vmin, wt(MAXBUF)
c-----------------------------------------------------------------------
      epsilon = 1e-12

c     Sentinel values; these mean "no unflagged pixels"
      do i = 1, 3
        momx(i)  = -1.0
        momy(i)  = -1.0
        stats(i) = -1.0
      enddo

c     We need the local minimum value to get the moments right - the
c     weights must be +ve, and the minimum weight should be zero, to
c     avoid biasing the result.  So get the min,max, mean & rms of all
c     unmasked pixels.  The mean and rms are needed for clipping so we
c     get the right answer when finding peaks at low s/n.

c     Get statistics for the current region.
      num   = 0
      sum   = 0.0
      sumsq = 0.0
      avg   = 0.0
      rms   = 0.0
      vmin  = imhi
      vmax  = imlo

      if (getmin .and. getmax) then
c       No user-specified limits.
        do j = trc(2),blc(2),-1
          call xyread(lIn,j,wt)
          call xyflgrd(lIn,j,flags)
          do i = blc(1), trc(1)
            val = wt(i)
            if (flags(i)) then
              vmin  = min(vmin,val)
              vmax  = max(vmax,val)
              num   = num   + 1
              sum   = sum   + val
              sumsq = sumsq + val*val
            endif
          enddo
        enddo

      else
        if ((.not.getmin) .and. (.not.getmax)) then
c         User has fixed range; clip pixels outside it.
          do j = trc(2),blc(2),-1
            call xyread(lIn,j,wt)
            call xyflgrd(lIn,j,flags)
            do i = blc(1), trc(1)
              val = wt(i)
              if (flags(i) .and. val.gt.rlo .and. val.lt.rhi) then
                vmin  = min(vmin,val)
                vmax  = max(vmax,val)
                num   = num   + 1
                sum   = sum   + val
                sumsq = sumsq + val*val
              endif
            enddo
          enddo

        else
          if (.not.getmin) then
c           User has fixed minimum only: clip pixels below rlo.
            do j = trc(2),blc(2),-1
              call xyread(lIn,j,wt)
              call xyflgrd(lIn,j,flags)
              do i = blc(1), trc(1)
                val = wt(i)
                if (flags(i) .and. val.gt.rlo) then
                  vmin  = min(vmin,val)
                  vmax  = max(vmax,val)
                  num   = num   + 1
                  sum   = sum   + val
                  sumsq = sumsq + val*val
                endif
              enddo
            enddo
          else
c           User has fixed maximum; clip pixels above rhi.
            do j = trc(2),blc(2),-1
              call xyread(lIn,j,wt)
              call xyflgrd(lIn,j,flags)
              do i = blc(1), trc(1)
                val = wt(i)
                if (flags(i) .and. val.lt.rhi) then
                  vmin  = min(vmin,val)
                  vmax  = max(vmax,val)
                  num   = num   + 1
                  sum   = sum   + val
                  sumsq = sumsq + val*val
                endif
              enddo
            enddo
          endif
        endif
      endif

c     Don't waste time if no data
      if (num.eq.0) then
        call bug('w','no valid pixels selected')
        return
      endif

      avg = sum/num
      rms = sumsq/num - avg*avg
      if (rms.gt.0.0) then
        rms = sqrt(rms)
      else
        rms = 0.0
        call bug('w','rms was root of a negative number')
      endif

c     Set up intensity range of interest.
      offset = -1.0*vmin
      if (clipmean) vmin = avg
      if (clip1sig) vmin = avg + rms

      apply = (.not.getmin) .or. (.not.getmax)
      apply = apply .or. (clipmean .or. clip1sig)

      num    = 0
      sum    = 0.0
      sumi   = 0.0
      sumsqi = 0.0
      sumcbi = 0.0
      sumj   = 0.0
      sumsqj = 0.0
      sumcbj = 0.0

c     Keep position values near 1, to minimise overflow errors.
      nrmi = 2.0/(blc(1)+trc(1))
      nrmj = 2.0/(blc(2)+trc(2))

c     Accumulate moment statistics for unflagged data.
      if (doskew) then
        if (apply) then
          do j = trc(2),blc(2),-1
            call xyread(lIn,j,wt)
            call xyflgrd(lIn,j,flags)
            sumrow = 0.0
            do i = blc(1), trc(1)
              if (flags(i)) then
                if (wt(i).gt.vmin .and. wt(i).le.vmax) then
                  val    = wt(i) + offset
                  ri     = real(i)*nrmi
                  num    = num    + 1
                  sumrow = sumrow + val
                  sumi   = sumi   + val*ri
                  sumsqi = sumsqi + val*ri*ri
                  sumcbi = sumcbi + val*ri*ri*ri
                endif
              endif
            enddo

c           Since rj constant for each row, we can do this for y-moment:
            rj     = real(j)*nrmj
            sum    = sum    + sumrow
            sumj   = sumj   + sumrow*rj
            sumsqj = sumsqj + sumrow*rj*rj
            sumcbj = sumcbj + sumrow*rj*rj*rj
          enddo

        else
          do j = trc(2),blc(2),-1
            call xyread(lIn,j,wt)
            call xyflgrd(lIn,j,flags)
            sumrow = 0.0
            do i = blc(1), trc(1)
              if (flags(i)) then
                ri     = real(i)*nrmi
                num    = num    + 1
                val    = wt(i)  + offset
                sum    = sum    + val
                sumrow = sumrow + val
                sumi   = sumi   + val*ri
                sumsqi = sumsqi + val*ri*ri
                sumcbi = sumcbi + val*ri*ri*ri
              endif
            enddo

            rj     = real(j)*nrmj
            sumj   = sumj   + sumrow*rj
            sumsqj = sumsqj + sumrow*rj*rj
            sumcbj = sumcbj + sumrow*rj*rj*rj
          enddo
        endif

      else
c       Only first 2 moments.
        if (apply) then
          do j = trc(2),blc(2),-1
            call xyread(lIn,j,wt)
            call xyflgrd(lIn,j,flags)
            sumrow = 0.0
            do i = blc(1), trc(1)
              if (flags(i)) then
                if (wt(i).gt.vmin .and. wt(i).le.vmax) then
                  val    = wt(i) + offset
                  ri     = real(i)*nrmi
                  num    = num    + 1
                  sumrow = sumrow + val
                  sumi   = sumi   + val*ri
                  sumsqi = sumsqi + val*ri*ri
                endif
              endif
            enddo

            rj     = real(j)*nrmj
            sum    = sum    + sumrow
            sumj   = sumj   + sumrow*rj
            sumsqj = sumsqj + sumrow*rj*rj
          enddo

        else
          do j = trc(2),blc(2),-1
            call xyread(lIn,j,wt)
            call xyflgrd(lIn,j,flags)
            sumrow = 0.0
            do i = blc(1), trc(1)
              if (flags(i)) then
                ri      = real(i)*nrmi
                num     = num    + 1
                val     = wt(i)  + offset
                sumrow  = sumrow + val
                sumi    = sumi   + val*ri
                sumsqi  = sumsqi + val*ri*ri
              endif
            enddo

            rj     = real(j)*nrmj
            sum    = sum    + sumrow
            sumj   = sumj   + sumrow*rj
            sumsqj = sumsqj + sumrow*rj*rj
          enddo
        endif
      endif

      if (num.eq.0) then
        call bug('w','No valid data in chosen range')
        return
      endif

c     Calculate moments:
c
c       mom1 = sum[ x(i)*wt(i) ] / sum[wt(i)]
c       mom2 = sqrt(  sum[ wt(i) * (x(i)-mom1)**2 ] / sum[wt(i)] )
c
c            = sqrt(   sum[          wt(i)*x(i)*x(i)] / sum[wt(i)]
c                    - sum[   2*mom1*wt(i)*x(i)]      / sum[wt(i)]
c                    + sum[mom1*mom1*wt(i)]           / sum[wt(i)] )
c
c       mom3 = cbrt ( sum [ wt(i) * (x(i)-mom1)**3 ] / sum[ wt(i) ] )
c       skew =        sum [ wt(i) * (x(i)-mom1)**3 ] / mom2**3
c
c     Numerical Details:
c       You don't get the right answer for the first or third moments if
c       the data cross zero.
c
c       For the first moment, this can be done after the summation at
c       the cost of summimg the pixel numbers (ri, rj).
c       However for higher moments, sums of rj**3 etc have to be
c       accumulated, so it's cheaper to offset every pixel; only 1 extra
c       addition is needed, once the offset is known.

      stats(1) = real(num)
      stats(2) = vmin
      stats(3) = vmax

      if (abs(sum).gt.epsilon) then
        momx(1) = sumi / sum
        momx(2) = (sumsqi - 2.0*momx(1)*sumi)/sum + momx(1)*momx(1)
        if (momx(2).lt.epsilon) then
          momx(2) = 0.0
          momx(3) = 0.0
          call bug('w','mom2 was root of a negative number')
        else
          momx(2) = sqrt(momx(2))
          if (doskew) then
            momx(3) = sumcbi - 3*sumsqi*momx(1) + 3*sumi*(momx(1)**2)
            momx(3) = momx(3)/sum - (momx(1)**3)
            momx(3) = momx(3)/(momx(2)**3)
          endif
        endif
        momx(1) = momx(1) / nrmi
        momx(2) = momx(2) / nrmi

        momy(1) = sumj/sum
        momy(2) = (sumsqj - 2*momy(1)*sumj)/sum + momy(1)*momy(1)
        if (momy(2).lt.epsilon) then
          momy(2) = 0.0
          momy(3) = 0.0
          call bug('w','mom2 was root of a negative number')
        else
          momy(2) = sqrt(momy(2))
          if (doskew) then
            momy(3) = sumcbj - 3*sumsqj*momy(1) + 3*sumj*(momy(1)**2)
            momy(3) = momy(3)/sum - (momy(1)**3)
            momy(3) = momy(3)/(momy(2)**3)
          endif
        endif

        momy(1) = momy(1) / nrmj
        momy(2) = momy(2) / nrmj
      endif

      end
