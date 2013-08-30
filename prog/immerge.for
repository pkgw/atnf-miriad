      program immerge

c= immerge -- Linear merging of images.
c& rjs
c: map combination
c+
c       IMMERGE is a Miriad task to merge linearly two images with
c       different resolutions.  The two images must be of the same field
c       and use the same coordinate system.
c
c       In combining the data, it is assumed that the low resolution
c       image better represents the short spacing data, whereas the
c       high resolution best represents the fine structure. Commonly,
c       the low resolution image will be a single-dish observation,
c       and the high resolution will be a mosaiced interferometric
c       observation.
c@ in
c       This gives the two images to be merged. The first input must be
c       the high resolution image, and the second the low. There is
c       no default. The two images must be on the same coordinate grid.
c       If necessary, use REGRID to achieve this.
c@ out
c       The output image. No default.
c@ factor
c       The flux calibration factor. Ideally the two inputs will have
c       correctly calibrated flux scales, in units like Jy/beam. In this
c       case, the flux calibration factor would be 1. In practise, the
c       calibration may not be perfect, and an extra calibration factor
c       will be required to align the flux scales.  "factor" gives the
c       factor to scale the low resolution image to put it on the same
c       flux scale as the high resolution image.  Note that this factor
c       is in addition to the scaling needed to convert Jy/beam at a low
c       resolution to Jy/beam at the high resolution.
c
c       If no value is given for "factor", IMMERGE determines this by
c       comparing the values in the Fourier plane, of data within an
c       particular annulus (after accounting for the differing
c       resolutions of the two inputs).  IMMERGE finds the scale factor
c       that minimises differences between the data of the two (in a
c       robust/L1 sense) in the annulus.
c@ uvrange
c       This specifies an annulus, in the Fourier domain, where the
c       high and low resolution images should agree (after allowing
c       for resolution differences). This is the annulus of data used to
c       deduce the flux calibration factor, and with options=feather.
c
c       Two or three values can be given.  The first two give the inner
c       and outer radius of the annulus in the Fourier domain.  The
c       third argument gives the units. Possible units are "klambda"
c       (kilo-wavelengths), "meters", "feet" and "nanoseconds".  The
c       default is "klambda".
c
c       Values for "uvrange" must be given either if "options=feather"
c       is used or if the flux calibration factor is being deduced.
c@ region
c       Region-of-interest parameter.  See the help on "region"
c       for more information.  NOTE: This parameter is ONLY used for
c       determining the flux calibration factor.  Only plane selection
c       (e.g. via the "image" command) is allowed.  Typically you
c       would want to select a range of planes which contains
c       significant signal in the overlap region.
c@ device
c       PGPLOT device for a plot.  When determining the flux calibration
c       factor, IMMERGE can produce a plot showing the correspondence
c       between the high and low resolution data points in the annulus
c       (after correcting for resolution effects and the deduced flux
c       calibration factor).  Ideally it will show a line with "y=x".
c       The default is not to produce a plot.  It also plots the
c       difference from this "y=x" line as a function of spatial
c       frequency.
c@ guard
c       Before Fourier transforming, the images are padded with a guard
c       band.  "guard" gives one or two values, being the minimum width
c       of this guard band in pixels, in the x and y directions.  The
c       actual guard band used is such that the size of the image plus
c       guard band is a power of 2.
c@ options
c       Task enrichment parameters.  Several can be given, separated by
c       commas.  Minimum match is used.
c         normalize  Rather than the output being the merged images, the
c                    output is the low resolution image corrected by
c                    the flux calibration factor.
c         zero       By default, IMMERGE pads the guard band with data
c                    that minimizes FFT edge effects.  If the input
c                    images are really zero beyond the edges of the two
c                    input images, then padding with zeros might be
c                    preferable.  This is particularly so if IMMERGE is
c                    deducing the flux calibration scale factor.
c         feather    This merges the two images in a fashion similar to
c                    AIPS IMERG.  This method is generally less
c                    desirable than the default scheme used by IMMERGE.
c         shift      Determine the optimum shift to apply to the low
c                    resolution image to make it align with the high
c                    resolution one.
c         notaper    Normally the low-resolution image is tapered to
c                    match any residual primary beam response in the
c                    high-resolution image.  This option causes this
c                    step to be skipped.
c
c$Id$
c--
c  History:
c    Refer to the RCS log, v1.1 includes prior revision information.
c-----------------------------------------------------------------------
      include 'maxdim.h'
      include 'maxnax.h'
      include 'mirconst.h'
      include 'mem.h'

      integer MAXBOX, NUNITS
      parameter (MAXBOX = 2048, NUNITS = 4)

      logical   dofac, dofeath, domerge, doout, doshift, dotaper,
     *          dozero, notaper
      integer   box(MAXBOX), i, iax, ifail, k, lIn1, lIn2, lOut, n,
     *          naxis, ngx, ngy, nin(3), nout(MAXNAX), npnt, ntemp(3),
     *          pIn1, pIn2, pWt1, pWt2, xoff, yoff, zoff
      real      bmaj, bmaj1, bmaj2, bmajt, bmin, bmin1, bmin2, bmint,
     *          bpa, bpa1, bpa2, bpat, du, dv, fac, freq, lambda, norm,
     *          pfac, sfac, sxx, sxy, syy, temp, uvhi, uvlo, xs, ys
      double precision cdelt1, cdelt2, freq1, freq2, x1(2), x2(2)
      character algo*3, device*64, in1*80, in2*80, line*80, mess1*64,
     *          mess2*64, out*80, unit*12, units(NUNITS)*12, version*72

      logical   keyprsnt
      integer   nextpow2
      character versan*72
      external keyprsnt, nextpow2, versan

      data units/'klambda     ','meters      ','nanoseconds ',
     *           'feet        '/
c-----------------------------------------------------------------------
      version = versan('immerge',
     *                 '$Revision$',
     *                 '$Date$')

      call keyini
      call GetOpt(domerge,dozero,dofeath,doshift,notaper)
      call keyf('in',in1,' ')
      call keyf('in',in2,' ')
      if (in1.eq.' ' .or. in2.eq.' ')
     *  call bug('f','Two input files must be given')
      call keyi('guard',ngx,0)
      call keyi('guard',ngy,ngx)
      if (ngx.lt.0 .or. ngy.lt.0)
     *  call bug('f','Invalid values for guard')
      dofac = .not.keyprsnt('factor')
      call keyr('factor',fac,1.0)
      if (dofac) then
        call keya('device',device,' ')
        call BoxInput('region',in1,box,MAXBOX)
      endif
      if (dofac .or. dofeath) then
        call keyr('uvrange',uvlo,0.0)
        call keyr('uvrange',uvhi,-1.0)
        call keymatch('uvrange',NUNITS,units,1,unit,n)
        if (n.eq.0) unit = units(1)
        if (uvlo.ge.uvhi)
     *    call bug('f','Invalid uvrange value')
      else
        uvlo = 0
        uvhi = -1
      endif
      call keya('out',out,' ')
      doOut = out.ne.' '
      if (.not.dofac .and. .not.doout)
     *  call bug('f','Inputs do not require any work')
      call keyfin
      dotaper = (domerge .or. dofac) .and. .not.notaper

c     Open the input datasets, and get their beam parameters.
      call xyopen(lIn1,in1,'old',3,nIn)
      call coInit(lIn1)
      call GetBeam(lIn1,in1,bmaj1,bmin1,bpa1,mess1)
      call coGetD(lIn1,'cdelt1',cdelt1)
      call coGetD(lIn1,'cdelt2',cdelt2)

      call xyopen(lIn2,in2,'old',3,ntemp)
      call coInit(lIn2)
      call GetBeam(lIn2,in2,bmaj2,bmin2,bpa2,mess2)

c     Check that the two images are compatible.
      do i = 1, 3
        if (nIn(i).ne.nTemp(i)) call bug('f','Incompatible image sizes')
      enddo
      call alignIni(lIn1,lIn2,nIn(1),nIn(2),nIn(3),xoff,yoff,zoff)
      if (xoff.ne.0 .or. yoff.ne.0 .or. zoff.ne.0)
     *  call bug('f','Inputs are not aligned')

      call rdhdi(lIn1,'naxis',naxis,1)
      naxis = min(naxis,MAXNAX)

c     Round up the total sizes to the next power of two.
      ngx = nextpow2(nIn(1)+ngx)
      ngy = nextpow2(nIn(2)+ngy)
      if (max(ngx,ngy).gt.MAXDIM)
     *  call bug('f','Image+guard size too big for me')

c     Determine which of the two images has the coarser resolution, and
c     determine the gaussian to convolve with to bring one to the same
c     resolution as the other.
      if (bmaj1*bmin1.gt.bmaj2*bmin2) call bug('f',
     *  'The first input image must be the high resolution')
      call GauDFac(bmaj2,bmin2,bpa2,bmaj1,bmin1,bpa1,
     *                        norm,bmaj,bmin,bpa,ifail)
      if (ifail.eq.0) call Gaufac(bmaj1,bmin1,bpa1,bmaj,bmin,bpa,
     *                        norm,bmajt,bmint,bpat,ifail)
      if (ifail.ne.0)
     *  call bug('f','Could not determine convolving parameters')
c       write(line,'(a,1pe10.3)')'Flux unit conversion factor:',
c     *                                 abs(cdelt1*cdelt2/norm)
c       call trimout(line)

c     Convert the gaussian parameters to the Fourier domain.
      call CvtParam(norm,bmaj,bmin,bpa,cdelt1,cdelt2,ngx,ngy,
     *                                        sfac,sxx,sxy,syy)

c     Determine the scaling factor to convert uvlo and uvhi to lambda.
      if (dofeath .or. dofac) then
        if (unit.eq.'nanoseconds') then
          temp = 1e-9*CMKS
        else if (unit.eq.'feet') then
          temp = 12*0.0254
        else if (unit.eq.'klambda') then
          temp = 1000
        else if (unit.eq.'meters') then
          temp = 1
        else
          call bug('f','Unrecognised units')
        endif

        if (unit.ne.'klambda') then
          call coSpcSet(lIn1, 'FREQ', ' ', iax, algo)
          if (iax.eq.0) call bug('f','No spectral axis in image 1')
          call coCvt1(lIn1, iax, 'op', 0d0, 'aw', freq1)

          call coSpcSet(lIn2, 'FREQ', ' ', iax, algo)
          if (iax.eq.0) call bug('f','No spectral axis in image 2')
          call coCvt1(lIn2, iax, 'op', 0d0, 'aw', freq2)

          if (min(freq1,freq2).le.0d0)
     *      call bug('f','Could not determine observing frequency')

          freq = sqrt(freq1*freq2)
          write(line,'(a,f8.3,a)')'Using a frequency of ',freq,
     *        ' GHz, to convert annulus to wavelengths'
          call trimout(line)
          lambda = CMKS/(freq*1e9)
          temp = temp/lambda
        endif
        uvlo = uvlo * temp
        uvhi = uvhi * temp
      endif
      du = 1/(ngx*cdelt1)
      dv = 1/(ngy*cdelt2)

c     Allocate memory if needed.
      if (dofac .or. domerge) then
        call memAlloc(pIn1,(ngx+2)*ngy,'r')
        call memAlloc(pIn2,(ngx+2)*ngy,'r')
        if (dotaper) then
          call mosLoad(lIn1,npnt)
          call memAlloc(pWt1,nIn(1)*nIn(2),'r')
          call memAlloc(pWt2,nIn(1)*nIn(2),'r')
        else
          pWt1 = pIn1
          pWt2 = pIn2
        endif
      endif

c     Determine the scale factor.
      if (dofac) then
        call boxSet(box,3,nIn,' ')
        pfac = (4.0*log(2.0))/PI*abs(cdelt1*cdelt2)/(bmaj2*bmin2)
        call GetFac(lIn1,lIn2,box,memr(pIn1),memr(pIn2),
     *    nIn(1),nIn(2),ngx,ngy,dozero,device,
     *      sfac,sxx,sxy,syy,uvlo,uvhi,du,dv,fac,pfac,
     *      dotaper,memr(pWt1),memr(pWt2),xs,ys,doshift)
        write(line,'(a,1pe11.3)')'Flux calibration factor:',fac
        call trimout(line)
        if (doshift) then
          x1(1) = xs
          x1(2) = ys
          call coCvt(lIn1,'op/op',x1,'ow/ow',x2)
          write(line,'(a,f10.2,a,f10.2,a)')
     *      'Applying a shift of',
     *      180*3600/PI*x2(1),', ',180*3600/PI*x2(2),
     *      ' arcseconds to low the resolution data'
          call trimout(line)
        endif
      endif

c     Open the output.
      if (doOut) then
        nOut(1) = nIn(1)
        nOut(2) = nIn(2)
        nOut(3) = nIn(3)
        do i = 4, naxis
          nOut(i) = 1
        enddo
        call xyopen(lOut,out,'new',naxis,nOut)
        if (domerge) then
          call MkHead(lIn1,lOut,version,fac,mess1,mess2)
        else
          call MkHead(lIn2,lOut,version,fac,mess1,mess2)
          call hdcopy(lIn2,lOut,'mask')
        endif

        do k = 1, nOut(3)
          if (k.ne.1) then
            call xySetpl(lIn1,1,k)
            call xySetpl(lIn2,1,k)
            call xySetpl(lOut,1,k)
          endif

          if (domerge) then
            call GetDat(lIn1,memr(pIn1),nIn(1),nIn(2),
     *          ngx,ngy,dozero,.false.,memr(pWt1),memr(pWt2))
            if (dotaper) call mosMIni(lIn1,real(k))
            call GetDat(lIn2,memr(pIn2),nIn(1),nIn(2),
     *          ngx,ngy,dozero,dotaper,memr(pWt1),memr(pWt2))
            if (abs(xs)+abs(ys).gt.0)
     *        call shiftit1(memr(pIn2),ngx,ngy,xs,ys)
            if (dotaper) call mosMFin
            call Mergeit(dofeath,memr(pIn1),memr(pIn2),ngx,ngy,
     *        uvlo,uvhi,du,dv,fac/sfac,sxx,sxy,syy)
            call WriteOut(lOut,memr(pIn1),ngx,nOut(2))
          else
            call Normalis(lIn2,lOut,nOut(1),nOut(2),fac)
          endif
        enddo
        call xyClose(lOut)
      endif

      call xyClose(lIn1)
      call xyClose(lIn2)
      if (dofac .or. domerge) then
        call memFree(pIn1,(ngx+2)*ngy,'r')
        call memFree(pIn2,(ngx+2)*ngy,'r')
        if (dotaper) then
          call memFree(pWt1,nIn(1)*nIn(2),'r')
          call memFree(pWt2,nIn(1)*nIn(2),'r')
        endif
      endif

      end

c***********************************************************************

      subroutine CvtParam(norm,bmaj,bmin,bpa,cdelt1,cdelt2,ngx,ngy,
     *                                        sfac,sxx,sxy,syy)

      real norm,bmaj,bmin,bpa
      double precision cdelt1,cdelt2
      integer ngx,ngy
      real sfac,sxx,syy,sxy
c-----------------------------------------------------------------------
      include 'mirconst.h'
      real theta,bmajf,bminf,c2,s2,a,b,dx,dy
c-----------------------------------------------------------------------
      sfac = PI/(4*log(2.0)) * bmaj * bmin / norm
      bmajf = 4*log(2.0)/PI/bmaj
      bminf = 4*log(2.0)/PI/bmin
      dx = 1/(cdelt1*ngx)
      dy = 1/(cdelt2*ngy)

      theta = pi/180.0 * bpa
      s2 = -sin(2*theta)
      c2 = -cos(2*theta)
      a = 4*log(2.0) / (bmajf*bmajf)
      b = 4*log(2.0) / (bminf*bminf)
      sxx = -0.5*(a*(c2+1) + b*(1-c2)) * dx*dx
      syy = -0.5*(b*(c2+1) + a*(1-c2)) * dy*dy
      sxy = -(b-a)*s2 * dx*dy

      end

c***********************************************************************

      subroutine Normalis(lIn,lOut,nx,ny,fac)

      integer lIn,lOut,nx,ny
      real fac
c-----------------------------------------------------------------------
      include 'maxdim.h'
      real Data(MAXDIM)
      integer i,j
c-----------------------------------------------------------------------
      do j = 1, ny
        call xyread(lIn,j,Data)
        do i = 1, nx
          Data(i) = fac * Data(i)
        enddo
        call xywrite(lOut,j,Data)
      enddo

      end

c***********************************************************************

      subroutine WriteOut(lOut,Data,ngx,ny)

      integer ngx,ny,lOut
      real Data(ngx+2,ny)
c-----------------------------------------------------------------------
      integer j
c-----------------------------------------------------------------------
      do j = 1, ny
        call xywrite(lOut,j,Data(1,j))
      enddo

      end

c***********************************************************************

      subroutine GetOpt(domerge,dozero,dofeath,doshift,notaper)

      logical domerge,dozero,dofeath,doshift,notaper
c-----------------------------------------------------------------------
c  Get extra processing parameters.
c-----------------------------------------------------------------------
      integer NOPTS
      parameter (NOPTS=5)
      character opts(NOPTS)*12
      logical   present(NOPTS)

      data opts/'normalize   ','zero        ','feather     ',
     *          'shift       ','notaper     '/
c-----------------------------------------------------------------------
      call options('options',opts,present,NOPTS)
      domerge = .not.present(1)
      dozero  =      present(2)
      dofeath =      present(3)
      if (dofeath .and. .not.domerge) call bug('f',
     *  'The feather and normalize options make no sense together')
      doshift =      present(4)
      notaper =      present(5)

      end

c***********************************************************************

      subroutine MkHead(lIn,lOut,version,fac,mess1,mess2)

      integer   lIn, lOut
      character version*(*), mess1*(*), mess2*(*)
      real      fac
c-----------------------------------------------------------------------
c  Make the output dataset header.
c-----------------------------------------------------------------------
      character line*64
c-----------------------------------------------------------------------
c     Copy the header verbatim.
      call headcp(lIn, lOut, 0, 0, 0, 0)

c     Create the output history.
      call hisopen(lOut,'append')
      line = 'IMMERGE: Miriad '//version
      call hiswrite(lOut,line)
      call hisinput(lOut,'IMMERGE')

      call hiswrite(lOut,
     *  'IMMERGE: The inputs are assumed to have Gaussian beams')

      line = 'IMMERGE: '//mess1
      call hiswrite(lOut,line)
      line = 'IMMERGE: '//mess2
      call hiswrite(lOut,line)

      write(line,'(a,1pe13.5)')
     *  'IMMERGE: Using a flux calibration factor of',fac
      call hiswrite(lOut,line)
      call hisclose(lOut)

      end

c***********************************************************************

      subroutine GetDat(lIn,In,nx,ny,ngx,ngy,dozero,dotaper,Wt1,Wt2)

      integer lIn,nx,ny,ngx,ngy
      logical dozero,dotaper
      real In(ngx+2,ngy),Wt1(nx,ny),Wt2(nx,ny)
c-----------------------------------------------------------------------
      include 'maxdim.h'
      integer i,j,length,ntaper
      real scale,x,fac,staper
      logical flags(MAXDIM)
      character ctaper*4
c-----------------------------------------------------------------------
c     Get the weight image
      ntaper = 0
      staper = 0
      if (dotaper) call mosWts(Wt1,Wt2,nx,ny,0,0)

      do j = 1, ny
c       Get pixels and flags and set bad pixels to zero.
        call xyRead(lIn,j,In(1,j))
        call xyflgrd(lIn,j,flags)
        do i = 1, nx
          if (.not.flags(i)) In(i,j) = 0
        enddo

        if (dotaper) then
          do i = 1, nx
            ntaper = ntaper + 1
            if (Wt1(i,j).gt.0) then
              In(i,j) = In(i,j)/Wt1(i,j)
              staper = staper + 1/Wt1(i,j)
            else
              In(i,j) = 0
            endif
          enddo
        endif

c       Zero extend the image.
        do i = nx+1, ngx
          In(i,j) = 0
        enddo

c       Reflect the image in a bizzar fashion, if needed.
        if (.not.dozero) then
          if (nx+1.eq.ngx) then
            In(nx+1,j) = 0.5*(In(nx,j) + In(1,j))
          else
            length = min(ngx-nx-1,nx-1)
            scale = 2.0/real(length)
            do i = 1, length
              x = scale*(i-1) - 1
              fac = 0.5 - x*(0.75 - 0.25*x*x)
              In(nx+i,j) = In(nx+i,j) + fac*In(nx-i,j)
              In(ngx-i+1,j) = In(ngx-i+1,j) + fac*In(i+1,j)
            enddo
          endif
        endif
      enddo

c     Zero the rest of the guard band.
      do j = ny+1, ngy
        do i = 1, ngx
          In(i,j) = 0
        enddo
      enddo

c     Reflect the image into the guard band if needed.
      if (.not.dozero) then
        if (ny+1.eq.ngy) then
          do i = 1, ngx
            In(i,ny+1) = 0.5*(In(i,ny) + In(i,1))
          enddo
        else
          length = min(ngy-ny-1,ny-1)
          scale = 2.0/real(length)
          do j = 1, length
            x = scale*(j-1) - 1
            fac = 0.5 - x*(0.75 - 0.25*x*x)
            do i = 1, ngx
              In(i,ny+j) = In(i,ny+j) + fac*In(i,ny-j)
              In(i,ngy-j+1) = In(i,ngy-j+1) + fac*In(i,j+1)
            enddo
          enddo
        endif
      endif

c     Report on the tapering.
      if (ntaper.gt.0) then
        staper = staper/ntaper
        write(ctaper,'(f4.2)') staper
        if (staper.lt.0.5) call bug('w','Taper correction for low '//
     *    'resolution image is small (mean='//ctaper//')')
        if (staper.lt.0.1) call bug('f','This is too small')
      endif

c     Fourier transform the image.
      call FFFT(In,ngx,ngy)

      end

c***********************************************************************

      subroutine GetBeam(lIn,in,bmaj,bmin,bpa,message)

      integer lIn
      character in*(*),message*(*)
      real bmaj,bmin,bpa
c-----------------------------------------------------------------------
c  Get the gaussian beam parameters corresponding to this observation.
c-----------------------------------------------------------------------
      include 'mirconst.h'
      character line*80,pbtype*16,label*18
      double precision ra,dec
      real rms,pbfwhm,cutoff,maxrad
      integer npnt,pbObj

      character itoaf*6
      external  itoaf
c-----------------------------------------------------------------------
      call rdhdr(lIn,'bmaj',bmaj,0.0)
      call rdhdr(lIn,'bmin',bmin,0.0)
      call rdhdr(lIn,'bpa',bpa,0.0)

c     If the Gaussian parameters were not in the header, assume its
c     a single dish, and try to get the primary beam size.
      if (bmaj*bmin.le.0) then
        call mosLoad(lIn,npnt)
        line = 'Cannot determine gaussian beam '//
     *         'parameters for mosaic in '//in
        if (npnt.gt.1) call bug('f',line)
        call mosGet(1,ra,dec,rms,pbtype)
        call pbInit(pbObj,pbtype,lIn)
        call pbInfo(pbObj,pbfwhm,cutoff,maxrad)
        call pbFin(pbObj)
        bmaj = pbfwhm
        bmin = pbfwhm
        bpa = 0
        line = 'Assuming single-dish observation in '//in
      else
        line = 'Found gaussian beam parameters for '//in
      endif

      call trimout(line)
      label = ' arcsec; pa='//itoaf(nint(bpa))
      write(line,'(a,f8.2,a,f8.2,a)')' ... with beam fwhm:',
     *  (3600*180/PI)*bmaj,' by',(3600*180/PI)*bmin,
     *  label
      call trimout(line)
      message = line

      end

c***********************************************************************

      subroutine GetFac(lIn1,lIn2,box,In1,In2,nx,ny,ngx,ngy,
     *  dozero,device,sfac,sxx,sxy,syy,uvlo,uvhi,du,dv,fac,pfac,
     *  dotaper,Wt1,Wt2,xs,ys,doshift)

      integer nx,ny,ngx,ngy,lIn1,lIn2
      integer box(*)
      logical dozero,doshift,dotaper
      complex In1(ngx/2+1,ngy),In2(ngx/2+1,ngy)
      real sfac,sxx,sxy,syy,uvlo,uvhi,fac,du,dv,pfac,xs,ys
      character device*(*)
      real Wt1(nx,ny),Wt2(nx,ny)
c-----------------------------------------------------------------------
c  Determine the scale factor to normalise everything.
c
c  Inputs:
c    In1        Fine-resolution input data.
c    In2        Coarse-resolution input data.
c    ngx,ngy    Fourier transform size.
c    sfac,sxx,sxy,syy Gaussian parameters.
c    uvlo,uvhi  Annulus of interest, in lambda.
c    du,dv      Cell increments in u and v.
c    device     PGPLOT device to plot the fit on.
c    pfac       Scale factor used in plotting.
c  Output:
c    fac        Scale factor to apply to the coarser resolution data.
c-----------------------------------------------------------------------
      include 'mirconst.h'
      include 'maxdim.h'
      include 'mem.h'
      integer MAXRUNS
      parameter (MAXRUNS=MAXDIM+1)
      integer nul,nuh,nvl,nvh,n,pX,pY,pUU,pVV,np,npd
      integer imin,imax,jmin,jmax,kmin,kmax,k
      integer nruns,runs(3,MAXRUNS)
      integer blc(3),trc(3)

      character itoaf*8
      external  itoaf
c-----------------------------------------------------------------------
c     Determine the range of planes to process.
      call boxInfo(box,3,blc,trc)
      kmin = blc(3)
      kmax = trc(3)

c     Determine the number of pixels which are in the overlap annulus.
      nul = nint(abs(uvlo/du)-0.5)
      nuh = nint(abs(uvhi/du)+0.5)
      nvl = nint(abs(uvlo/dv)-0.5)
      nvh = nint(abs(uvhi/dv)+0.5)
      n = (kmax-kmin+1)*
     *        nint(0.5*PI*(nuh*nvh - nul*nvl) + nuh - nul + 1.5)
      call memAlloc(pX,n,'c')
      call memAlloc(pY,n,'c')
      call memAlloc(pUU,n,'i')
      call memAlloc(pVV,n,'i')

c     Get all the data in the annuli.
      np = 0
      do k = kmin, kmax
        call boxRuns(1,k,' ',box,runs,MAXRUNS,nruns,
     *                                        imin,imax,jmin,jmax)
        if (imin.ne.1 .or. imax.ne.nx .or.
     *      jmin.ne.1 .or. jmax.ne.ny .or.
     *    (nruns.ne.ny .and. nruns.ne.0)) call bug('f',
     *      'Only plane selection supported in region keyword')

        if (nruns.gt.0) then
          if (k.ne.1) then
            call xysetpl(lIn1,1,k)
            call xysetpl(lIn2,1,k)
          endif
          call GetDat(lIn1,In1,nx,ny,ngx,ngy,dozero,
     *                                        .false.,Wt1,Wt2)
          if (dotaper) call mosMIni(lIn1,real(k))
          call GetDat(lIn2,In2,nx,ny,ngx,ngy,dozero,
     *                                        dotaper,Wt1,Wt2)
          if (dotaper) call mosMFin
          call AnnExt(In1,In2,ngx,ngy,memc(pX+np),memc(pY+np),
     *        memi(pUU+np),memi(pVV+np),
     *        n-np,npd,sfac,sxx,sxy,syy,uvlo,uvhi,du,dv)
          np = np + npd
        endif
      enddo

      call trimout('Number of data points in the annulus: '
     *                                        //itoaf(np))
      if (np.eq.0) call bug('f','No Data points found in annulus')

c     Determine the offset.
      if (doshift) then
        call GetShift(np,memi(pUU),memi(pVV),memc(pX),memc(pY),
     *                                        In1,ngx,ngy,xs,ys)
        call Shiftit2(np,memi(pUU),memi(pVV),memc(pY),xs,ys,ngx,ngy)
      else
        xs = 0
        ys = 0
      endif

c     Fit for the scale factor.
      call DetFac(memc(pX),memc(pY),2*np,fac)

c     Plot the scale factor, if the user wanted this.
      if (device.ne.' ') call PlotFac(device,memc(pX),memc(pY),
     *                                2*np,fac,pfac,doshift)

c     Free the allocated memory.
      call memFree(pX,n,'c')
      call memFree(pY,n,'c')
      call memFree(pUU,n,'i')
      call memFree(pVV,n,'i')

c     Reset to select the first plane.
      if (kmax.ne.1) then
        call xysetpl(lIn1,1,1)
        call xysetpl(lIn2,1,1)
      endif

      end

c***********************************************************************

      subroutine GetShift(np,uu,vv,x,y,data,ngx,ngy,xs,ys)

      integer np,ngx,ngy
      integer uu(np),vv(np)
      complex x(np),y(np),data(ngx/2+1,ngy)
      real xs,ys
c-----------------------------------------------------------------------
c  Determine the required shift to optimally align the X and Y data.
c-----------------------------------------------------------------------
      include 'maxdim.h'
      integer i,j,k,l
c-----------------------------------------------------------------------
      do j = 1, ngy
        do i = 1, ngx/2+1
          data(i,j) = (0.0,0.0)
        enddo
      enddo

c     Accumulate the cross correlation between X and Y.
      do k = 1, np
        l = uu(k) + vv(k)
        if (2*(l/2).eq.l) then
          data(uu(k),vv(k)) = data(uu(k),vv(k)) + x(k)*conjg(y(k))
        else
          data(uu(k),vv(k)) = data(uu(k),vv(k)) - x(k)*conjg(y(k))
        endif
      enddo

c     Fourier transform it now.
      call IFFT(data,ngx,ngy)

c     Find the peak of the data.
      call findpk(data,ngx,ngy,xs,ys)

      end

c***********************************************************************

      subroutine Shiftit2(np,uu,vv,y,xs,ys,ngx,ngy)

      integer np,ngx,ngy
      integer uu(np),vv(np)
      real xs,ys
      complex y(np)
c-----------------------------------------------------------------------
      include 'mirconst.h'
      integer i
      complex W
      real theta
c-----------------------------------------------------------------------
      do i = 1, np
        theta = 2*PI*(xs*(uu(i)-1)/real(ngx)+ys*(vv(i)-1)/real(ngy))
        W = cmplx(cos(theta),sin(theta))
        y(i) = W*y(i)
      enddo

      end

c***********************************************************************

      subroutine shiftit1(data,ngx,ngy,xs,ys)

      integer ngx,ngy
      real xs,ys
      complex data(ngx/2+1,ngy)
c-----------------------------------------------------------------------
      include 'mirconst.h'
      integer i,j
      complex W
      real theta
c-----------------------------------------------------------------------
      do j = 1, ngy
        do i = 1, ngx/2+1
          theta = 2*PI*(xs*(i-1)/real(ngx)+ys*(j-1)/real(ngy))
          W = cmplx(cos(theta),sin(theta))
          data(i,j) = W*data(i,j)
        enddo
      enddo

      end

c***********************************************************************

      subroutine findpk(data,ngx,ngy,xs,ys)

      integer ngx,ngy
      real data(ngx+2,ngy),xs,ys
c-----------------------------------------------------------------------
c  Find the position, to a fraction of a pixel of the maximum pixel.
c-----------------------------------------------------------------------
      integer i,j,k,imax,jmax
      real pkval,fit(9),coeffs(6),fmax
      double precision pixmax(2)
c-----------------------------------------------------------------------
      imax = 1
      jmax = 1
      pkval = data(1,1)
      do j = 1, ngy
        do i = 1, ngx
          if (data(i,j).gt.pkval) then
            pkval = data(i,j)
            imax = i
            jmax = j
          endif
        enddo
      enddo

      if (imax.eq.1 .or. imax.eq.ngx .or. jmax.eq.1 .or. jmax.eq.ngy)
     *  call bug('f','Cannot fit to peak next to image edge')
      k = 0
      do j = jmax-1, jmax+1
        do i = imax-1, imax+1
          k = k + 1
          fit(k) = data(i,j)
        enddo
      enddo

c     Locate the peak.
      call pkfit(fit,3,fmax,pixmax,coeffs)
      xs = imax + pixmax(1) - (ngx/2+1)
      ys = jmax + pixmax(2) - (ngy/2+1)

      end

c***********************************************************************

      subroutine PlotFac(device,X,Y,n,a,pfac,doshift)

      integer n
      character device*(*)
      real X(n),Y(n),a,pfac
      logical doshift
c-----------------------------------------------------------------------
      real xmin,xmax,xlo,xhi,xp(2)
      integer i

      integer  pgbeg
      external pgbeg
c-----------------------------------------------------------------------
c     Find the min and max.
      xmin = pfac*x(1)
      xmax = xmin
      do i = 1, n
        x(i) = pfac*x(i)
        y(i) = pfac*a*y(i)
        xmin = min(xmin,x(i),y(i))
        xmax = max(xmax,x(i),y(i))
      enddo

c     Create the plot of the normal data.
      if (pgbeg(0,device,1,1).ne.1) then
        call pgldev
        call bug('f','Error opening graphics device')
      endif
      call pgscf(2)
      call pgpage
      call pgvstd
      call pgrnge(xmin,xmax,xlo,xhi)
      call pgwnad(xlo,xhi,xlo,xhi)
      call pgbox('BCNST',0.0,0,'BCNST',0.0,0)
      if (n.lt.100) then
        call pgpt(n,x,y,17)
      else
        call pgpt(n,x,y,1)
      endif

c     Plot the y = x line.
      call pgsci(2)
      xp(1) = xlo
      xp(2) = xhi
      call pgline(2,xp,xp)

c     Label it.
      call pgsci(1)
      if (doshift) then
        call pglab('High Resolution Data (Jy)',
     *             'Scaled/Shifted Low Resolution Data (Jy)',
     *             'Plot of Data in Fourier Annulus')
      else
        call pglab('High Resolution Data (Jy)',
     *             'Scaled Low Resolution Data (Jy)',
     *             'Plot of Data in Fourier Annulus')
      endif
      call pgend

      end

c***********************************************************************

      subroutine DetFac(X,Y,n,a)

      integer n
      real X(n),Y(n),a
c-----------------------------------------------------------------------
c  Determine "a" such that
c    X \approx a * Y
c
c  Input:
c    X,Y        Data.
c    n          Number of data points.
c  Output:
c    a          Scale factor.
c-----------------------------------------------------------------------
      double precision SumXX,SumYY,SumXY
      real rms,a1,a2,f1,f2,f
      integer i

      real     medfunc
      external medfunc
c-----------------------------------------------------------------------
c     Determine the least squares fit first.
      SumXX = 0d0
      SumYY = 0d0
      SumXY = 0d0
      do i = 1, n
        SumXX = SumXX + X(i)*X(i)
        SumYY = SumYY + Y(i)*Y(i)
        SumXY = SumXY + X(i)*Y(i)
      enddo

      if (SumYY.eq.0) call bug('f','Low resolution image is zero')
      if (SumXX.eq.0) call bug('f','High resolution image is zero')
      a = SumXY/SumYY
      rms = (SumXX + a*a*SumYY - 2*a*SumXY)/(n*SumYY)
      rms = max(0.01*a,sqrt(max(0.0,rms)))

      a1 = a
      f1 = medfunc(a1,X,Y,n)
      a2 = a1 + sign(3*rms,f1)
      f2 = medfunc(a2,X,Y,n)
      do while (f1*f2.gt.0)
        a = 2*a2 - a1
        a1 = a2
        f1 = f2
        a2 = a
        f2 = medfunc(a2,X,Y,n)
      enddo

      do while (abs(a1-a2).gt.0.001*rms)
        a = 0.5*(a1+a2)
        if (a.eq.a1 .or. a.eq.a2) return
        f = medfunc(a,X,Y,n)
        if (f*f1.ge.0) then
          f1 = f
          a1 = a
        else
          f2 = f
          a2 = a
        endif
      enddo

      end

c***********************************************************************

      real function medfunc(a,X,Y,n)

      integer n
      real a,X(n),Y(n)
c-----------------------------------------------------------------------
c  Determine the function needed to solve for the median.
c-----------------------------------------------------------------------
      integer i
c-----------------------------------------------------------------------
      medfunc = 0
      do i = 1, n
        medfunc = medfunc + y(i)*sign(1.0,x(i)-a*y(i))
     *                    - x(i)/a/a*sign(1.0,y(i)-x(i)/a)
      enddo

      end

c***********************************************************************

      subroutine AnnExt(In1,In2,ngx,ngy,X,Y,UU,VV,nmax,np,
     *  sfac,sxx,sxy,syy,uvlo,uvhi,du,dv)

      integer ngx,ngy,nmax,np
      complex X(nmax),Y(nmax)
      integer UU(nmax),VV(nmax)
      complex In1(ngx/2+1,ngy),In2(ngx/2+1,ngy)
      real sfac,sxx,sxy,syy,uvlo,uvhi,du,dv
c-----------------------------------------------------------------------
c  Extract the pixels from the annulus of overlap of the coarse and
c  fine resolution images.
c
c  Output:
c    X,Y        Pixels from the low and high resolution data, that we
c               are scaling to equalize.
c    R          Radius in wavelengths
c    np         Number of pixels extracted.
c-----------------------------------------------------------------------
      integer i,j,ic,jc
      real uv2,t
c-----------------------------------------------------------------------
      np = 0

c     Do the top half,
      ic = 1
      jc = 1
      do j = 1, ngy/2+1
        do i = 1, ngx/2+1
          uv2 = ((i-ic)*du)**2 + ((j-jc)*dv)**2
          if (uv2.ge.uvlo*uvlo .and. uv2.le.uvhi*uvhi) then
            t = sxx*(i-ic)*(i-ic) + sxy*(i-ic)*(j-jc) +
     *          syy*(j-jc)*(j-jc)
            t = sfac*exp(t)
            np = np + 1
            if (np.gt.nmax) call bug('f','Too many points in Annext')
            X(np) = t*In1(i,j)
            Y(np) =   In2(i,j)
            UU(np) = i
            VV(np) = j
          endif
        enddo
      enddo

c     Do the bottom half.
      ic = 1
      jc = ngy+1
      do j = ngy/2+2, ngy
        do i = 1, ngx/2+1
          uv2 = ((i-ic)*du)**2 + ((j-jc)*dv)**2
          if (uv2.ge.uvlo*uvlo .and. uv2.le.uvhi*uvhi) then
            t = sxx*(i-ic)*(i-ic) + sxy*(i-ic)*(j-jc) +
     *          syy*(j-jc)*(j-jc)
            t = sfac*exp(t)
            np = np + 1
            if (np.gt.nmax) call bug('f','Too many points in Annext')
            X(np) = t*In1(i,j)
            Y(np) =   In2(i,j)
            UU(np) = i
            VV(np) = j
          endif
        enddo
      enddo

      end

c***********************************************************************

      subroutine Mergeit(dofeath,In1,In2,ngx,ngy,uvlo,uvhi,du,dv,
     *                                        fac,sxx,sxy,syy)

      integer ngx,ngy
      complex In1(ngx/2+1,ngy),In2(ngx/2+1,ngy)
      real fac,sxx,sxy,syy,uvlo,uvhi,du,dv
      logical dofeath
c-----------------------------------------------------------------------
c  Merge two Fourier transforms.
c
c  Input/Output:
c    In1        On input,  the transform of the finer-resolution data.
c               On output, the merged data in the image domain.
c  Input:
c    In2        This is the transform of the coarser-resolution data.
c    ngx,ngy
c    uvlo,uvhi
c    fac,sxx,sxy,syy Gaussian parameters.
c    du,dv      Cell increments in u and v.
c-----------------------------------------------------------------------
      integer i,j,ic,jc
      real uv2,t,x,y
c-----------------------------------------------------------------------
c     Do the top half.
      ic = 1
      jc = 1
      do j = 1, ngy/2+1
        do i = 1, ngx/2+1
          uv2 = ((i-ic)*du)**2 + ((j-jc)*dv)**2
          t = sxx*(i-ic)*(i-ic) + sxy*(i-ic)*(j-jc) +
     *        syy*(j-jc)*(j-jc)
          if (dofeath) then
            if (uv2.lt.uvlo*uvlo) then
              t = fac*exp(-t)
              In1(i,j) = t*In2(i,j)
            else if (uv2.le.uvhi*uvhi) then
              t = fac*exp(-t)
              x = 2*(sqrt(uv2) - uvlo)/(uvhi-uvlo) - 1
              y = 0.5 - x*(0.75 - 0.25*x*x)
              In1(i,j) = (1-y)*In1(i,j) + y*t*In2(i,j)
            endif
          else if (t.gt.-20) then
            In1(i,j) = fac*In2(i,j) + (1-exp(t))*In1(i,j)
          endif
        enddo
      enddo

c     Do the bottom half.
      ic = 1
      jc = ngy+1
      do j = ngy/2+2, ngy
        do i = 1, ngx/2+1
          uv2 = ((i-ic)*du)**2 + ((j-jc)*dv)**2
          t = sxx*(i-ic)*(i-ic) + sxy*(i-ic)*(j-jc) +
     *        syy*(j-jc)*(j-jc)
          if (dofeath) then
            if (uv2.lt.uvlo*uvlo) then
              t = fac*exp(-t)
              In1(i,j) = t*In2(i,j)
            else if (uv2.le.uvhi*uvhi) then
              t = fac*exp(-t)
              x = 2*(sqrt(uv2) - uvlo)/(uvhi-uvlo) - 1
              y = 0.5 - x*(0.75 - 0.25*x*x)
              In1(i,j) = (1-y)*In1(i,j) + y*t*In2(i,j)
            endif
          else if (t.gt.-20) then
            In1(i,j) = fac*In2(i,j) + (1-exp(t))*In1(i,j)
          endif
        enddo
      enddo

c     Fourier transform back to the image domain.
      call IFFT(In1,ngx,ngy)

      end

c***********************************************************************

      subroutine IFFT(In,ngx,ngy)

      integer ngx,ngy
      complex In(ngx/2+1,ngy)
c-----------------------------------------------------------------------
      include 'maxdim.h'
      integer i,j
      real scale
      complex Temp1(MAXDIM),temp2(MAXDIM)
c-----------------------------------------------------------------------
      scale = 1/real(ngx*ngy)
      do i = 1, ngx/2+1
        do j = 1, ngy
          Temp1(j) = scale*In(i,j)
        enddo
        call fftcc(Temp1,Temp2,-1,ngy)
        do j = 1, ngy
          In(i,j) = Temp2(j)
        enddo
      enddo

      do j = 1, ngy
        do i = 1, ngx/2+1
          Temp1(i) = In(i,j)
        enddo
        call fftcr(Temp1,In(1,j),-1,ngx)
      enddo

      end

c***********************************************************************

      subroutine FFFT(In,ngx,ngy)

      integer ngx,ngy
      complex In(ngx/2+1,ngy)
c-----------------------------------------------------------------------
      include 'maxdim.h'
      integer i,j
      complex Temp1(MAXDIM),Temp2(MAXDIM)
c-----------------------------------------------------------------------
c     First pass of the FFT.
      do j = 1, ngy
        call fftrc(In(1,j),Temp1,1,ngx)
        do i = 1, ngx/2+1
          In(i,j) = Temp1(i)
        enddo
      enddo

c     Second pass of the FFT.
      do i = 1, ngx/2+1
        do j = 1, ngy
          Temp1(j) = In(i,j)
        enddo
        call fftcc(Temp1,Temp2,1,ngy)
        do j = 1, ngy
          In(i,j) = Temp2(j)
        enddo
      enddo

      end

c***********************************************************************

      subroutine trimout(line)

      character line*(*)
c-----------------------------------------------------------------------
      integer i,l
      logical blank,add
      character message*128
c-----------------------------------------------------------------------
      blank = .false.
      l = 0
      do i = 1, len(line)
        add = line(i:i).ne.' ' .or. .not.blank
        if (add) then
          l = l + 1
          message(l:l) = line(i:i)
        endif
        blank = line(i:i).eq.' '
      enddo
      if (l.gt.0) call output(message(1:l))

      end
