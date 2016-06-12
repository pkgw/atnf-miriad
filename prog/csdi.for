      program csdi

c= csdi - Complex Steer Clean
c& rjs mchw
c: deconvolution
c+
c       CSDI is a MIRIAD task that performs the complex Steer Clean 
c       algorithm, which takes a Stokes Q and U dirty map and beam, 
c       and produces two output maps that consist of the Stokes Q and  
c       Stokes U Clean components.
c       These outputs can be input to RESTOR to produce a "clean" image.  
c       The model could be from a previous CSDI run, or from other
c       deconvolution tasks.
c
c       Full details of the algorithm are found in: 
c       Pratley & Johnston-Hollitt, “An improved method for polarimetric 
c       image restoration in interferometry", MNRAS, 2016. ArXiv: 1606.01482. 
c	Please acknowledge this work in publications using the code.
c       
c@ map
c       The input Stokes Q and U dirty maps, which should have units of
c       Jy/beam. No default.
c@ model
c       Initial models of the Stokes Q and U deconvolved images. This  
c       could be the output from a previous run of CSDI, or the output of
c       other deconvolution tasks. It must have flux units of
c       Jy/pixel. The default is no model (i.e. a zero map).
c@ beam
c       The input of Stokes I, Q or U dirty beam. No default
c@ out
c       The names of the Stokes Q and U output map. The units of the  
c       output will be Jy/pixel.  The files will contain the contribution
c       of the input models.  They should have a different name to the 
c       input model (if any).  They can be input to RESTOR or CSDI (as a
c       model,to do more cleaning)   
c@ gain
c       The minor iteration loop gain. Default is 0.1.
c@ options
c       Extra processing options. Several can be given, separated
c       by commas. Minimum match is used. Possible values are:
c         asym      The beam is asymmetric.  By default CLEAN assumes
c                   the beam has a 180 degree rotation symmetry, which
c                   is the norm for beams in radio-astronomy.
c         pad       Double the beam size by padding it with zeros. This
c                   will give you better stability if you are daring enough 
c                   to CLEAN an area more than half the size 
c                   (in each dimension) of the dirty beam.
c@ cutoff
c       CLEAN finishes either when the absolute maximum residual falls
c       below CUTOFF, or when the criteria described below is
c       satisfied. The default CUTOFF is 0. It is recommended that the
c       cutoff is 3 times the rms noise of Stokes Q or U.
c@ niters
c       The maximum number of minor iterations.  The default is 250,
c       which is too small for all but the simplest of images.  CLEAN
c       will stop when either the maximum number of iterations is
c       performed, or the cutoff (see above) is reached.
c@ region
c       This specifies the region to be Cleaned.  See the help on
c       "region" for more information.  The default is the largest
c       region that can be deconvolved safely.
c@ clip
c       This sets the relative clip level in Steer mode. Values are
c       typically 0.75 to 0.9. The default is image dependent.
c
c$Id$
c--
c  CSDI History:
c    lp  02aug15 - Original version modified from clean.for, v1.13
c    lp & mjh Jun 16 - updated task explanations and added paper reference
c    mhw 10jun16 - Some tidying up for inclusion in Miriad distribution      
c
c  Important Constants:
c    MaxDim     The max linear dimension of an input (or output) image.
c
c  Runs:
c    Run(3,nrun) We may want to clean several boxes (which may overlap).
c               Run is an 3 x N array, with each entry being of the
c               form:
c                 j, imin, imax
c               This gives a range of i in line j to process.  There
c               could be several entries for a given line, describing
c               disconnected ranges of i to process (there is no overlap
c               between runs).  Specifying boxes in this way is
c               reasonably concise for the programmer, and yet makes
c               vectorisable code straightforward to write.
c-----------------------------------------------------------------------
      include 'maxdim.h'
      integer MAXBEAM,MAXCMP1,MAXCMP2,MAXBOX,MAXRUN,MAXP
      parameter (MAXCMP1 = 66000, MAXCMP2 = 32000, MAXP = 257,
     *           MAXBEAM = MAXP*MAXP, MAXBOX = 3*MAXDIM,
     *           MAXRUN = 3*MAXDIM)

      real Data(MaxBuf)
      integer Boxes(MAXBOX),Run(3,MAXRUN),MaxMap
      ptrdiff pBem,pMap,pEst,pRes,upMap,upEst,upRes
      real Histo(MAXP/2+1)

      character Mode*8,Moded*8,Text*7,flags*8,version*72
      real Cutoff,Gain,Phat,Speed,Clip,defClip,Limit
      logical NegStop,Positive,Pad,Asym,NegFound,More,FFTIni,steermsg
      integer MaxNiter,oNiter,Niter,totNiter,minPatch,maxPatch
      integer totNiter1,totNiter2
      integer naxis,n1,n2,icentre,jcentre,nx,ny
      integer blc(3),trc(3),xmin,xmax,ymin,ymax
      integer k,nRun,nPoint,xoff,yoff,zoff
      character MapNam*256,uMapNam*256,BeamNam*256,ModelNam*256
      character uModelNam*256,OutNam*256,uOutNam*256,line*72
      integer lMap,ulMap,lBeam,lModel,ulModel,lOut,ulOut
      integer nMap(3),nBeam(3),nModel(3),nOut(4)
      real EstASum,Flux,uEstASum,uFlux
      real ResMin,ResMax,ResAMax,ResRms
      real uResMin,uResMax,uResAMax,uResRms
      real maxPres
      common Data

c     Externals.
      character itoaf*8, versan*72
c-----------------------------------------------------------------------
      version = versan ('csdi',
     *                  '$Revision$',
     *                  '$Date$')
c
c  Get the input parameters.
c
      call inputs(MapNam,uMapNam,BeamNam,ModelNam,uModelNam,
     *     OutNam,uOutNam,MaxNiter,
     *     NegStop,positive,pad,asym,Cutoff,Boxes,MAXBOX,MinPatch,Gain
     *     ,PHat,Speed,Clip,mode)
c
c  Open the beam, get some characteristics about it, then read in the
c  beam patch. The beam patch is not required if we are performing Steer
c  iterations only. However, some of the statistics returned by BeamChar
c  are required.
c
      call xyopen(lBeam,BeamNam,'old',3,nBeam)
      n1 = nBeam(1)
      n2 = nBeam(2)
      if (nBeam(3).ne.1) call bug('w',
     *  'Beam contains more than one plane')
      FFTIni = .false.
c
c  Fiddle the min and max patch sizes.
c
      maxPatch = min(MAXP,2*((min(n1,n2)-1)/2) + 1)
      if (maxPatch.le.0) call bug('f','Bad patch size')
      call BeamChar(lBeam,n1,n2,icentre,jcentre,Histo,maxPatch)
      maxPatch = min(maxPatch,
     *            2*min(icentre-1,n1-icentre,jcentre-1,n2-jcentre)+1)
      if (minPatch.gt.maxPatch) then
        call bug('w','Setting min patch size to '//itoaf(maxPatch))
        minPatch = maxPatch
      endif

        defClip = 0.2*Histo(1) + 0.8*Histo(2)
        if (Clip.eq.0) then
          Clip = defClip
        else if (Clip.lt.defClip) then
          call bug('w','Clip level seems small')
        endif

c  Open the map, and determine the area being cleaned.
c

      call xyopen(lMap,MapNam,'old',3,nMap)
      call xyopen(ulMap,uMapNam,'old',3,nMap)
      call rdhdi(lMap,'naxis',naxis,3)
      call rdhdi(ulMap,'naxis',naxis,3)
      naxis = min(naxis,4)
      call defregio(boxes,nMap,nBeam,icentre,jcentre)
      call BoxMask(lMap,boxes,MAXBOX)
      call BoxMask(ulMap,boxes,MAXBOX)
      call BoxSet(boxes,3,nMap,' ')
      call BoxInfo(boxes,3,blc,trc)
      nOut(1) = trc(1) - blc(1) + 1
      nOut(2) = trc(2) - blc(2) + 1
      nOut(3) = trc(3) - blc(3) + 1
      nOut(4) = 1
      if (nOut(1).gt.n1 .or. nOut(2).gt.n2)
     *  call bug('f','Region of map to deconvolve is too big')
      if (2*nOut(1)-1.gt.n1 .or. 2*nOut(2)-1.gt.n2)
     *  call bug('w','Size of region of map to deconvolve is unsafe')
c
c  Allocate space for the Map,Estimate and Residuals.
c
      MaxMap = nOut(1)*nOut(2)
      call MemAllop(pMap,MaxMap,'r')
      call MemAllop(pEst,MaxMap,'r')
      call MemAllop(pRes,MaxMap,'r')
      call MemAllop(upMap,MaxMap,'r')
      call MemAllop(upEst,MaxMap,'r')
      call MemAllop(upRes,MaxMap,'r')
c
c  Open the model if there is one.  Note that currently the model
c  must agree exactly in size with the output map (an unfortunate
c  restriction).
c
      if (ModelNam.ne.' ') then
        totNiter1 = 0
        totNiter2 = 0
        call xyopen(lModel,ModelNam,'old',3,nModel)
        call rdhdi(lModel,'niters',totNiter1,0)
        call AlignIni(lModel,lMap,nMap(1),nMap(2),nMap(3),
     *                                      xoff,yoff,zoff)
        call xyopen(ulModel,uModelNam,'old',3,nModel)
        call rdhdi(ulModel,'niters',totNiter2,0)
        call AlignIni(ulModel,ulMap,nMap(1),nMap(2),nMap(3),
     *                                      xoff,yoff,zoff)
        totNiter = NINT((totNiter1+totNiter2)*0.5)        
      else
        totNiter = 0
      endif
      steermsg = .true.
c
c  Open the output.
c
      call xyopen(lOut,OutNam,'new',naxis,nOut)
      call xyopen(ulOut,uOutNam,'new',naxis,nOut)
c
c  Loop over all the planes of interest.
c
      do k = blc(3), trc(3)
        if (blc(3).ne.trc(3)) call output('Plane: '//itoaf(k))
c
c  Get the Map, Estimate and Residual.
c
        call BoxRuns(1,k,'r',boxes,Run,MAXRUN,nRun,
     *                                xmin,xmax,ymin,ymax)
        nx = xmax - xmin + 1
        ny = ymax - ymin + 1

        call xysetpl(lMap,1,k)
        call xysetpl(ulMap,1,k)
        call GetPlane(lMap,Run,nRun,xmin-1,ymin-1,nMap(1),nMap(2),
     *                        Data(pMap),MaxMap,nPoint)
        call GetPlane(ulMap,Run,nRun,xmin-1,ymin-1,nMap(1),nMap(2),
     *                        Data(upMap),MaxMap,nPoint)
c
c  Determine the CLEAN algorithm that is to be used.
c
        if (nPoint.gt.0) then
          moded = mode
c
c  Initialise the FFT of the beam if needed.
c
          if ((moded.ne.'hogbom' .or. ModelNam.ne.' ')
     *                                .and. .not.FFTIni) then
            FFTIni = .true.
            flags = 'p'
            if (.not.asym) flags(2:2) = 's'
            if (pad)       flags(3:3) = 'e'
            call CnvlIniF(pBem,lBeam,n1,n2,icentre,jcentre,PHat,flags)
          endif
c
c  Initialise the estimate, and determine the residuals if the the user
c  gave an estimate. Determine statistics about the estimate and the
c  residuals.
c
          if (ModelNam.eq.' ') then
            EstASum = 0
            call NoModel(Data(pMap),Data(pEst),Data(pRes),nPoint)
            call NoModel(Data(upMap),Data(upEst),Data(upRes),nPoint)
          else
              call output('Subtracting initial model ...')
              call AlignGet(lModel,Run,nRun,k,xmin+xoff-1,ymin+yoff-1,
       *        zoff,nModel(1),nModel(2),nModel(3),
       *        Data(pEst),MaxMap,nPoint)
              call Diff(pBem,Data(pEst),Data(pMap),Data(pRes),
       *        nPoint,nx,ny,Run,nRun)
              call SumAbs(EstASum,Data(pEst),nPoint)
              call AlignGet(ulModel,Run,nRun,k,xmin+xoff-1,ymin+yoff-1,
       *        zoff,nModel(1),nModel(2),nModel(3),
       *        Data(upEst),MaxMap,nPoint)
              call Diff(pBem,Data(upEst),Data(upMap),Data(upRes),
       *        nPoint,nx,ny,Run,nRun)
              call SumAbs(uEstASum,Data(upEst),nPoint)
          endif
          call Stats(Data(pRes),nPoint,ResMin,ResMax,ResAMax,ResRms)
          call Stats(Data(upRes),nPoint,uResMin,uResMax,uResAMax,
     *         uResRms)
          call output('Begin iterating')
        else
          ResMin = 0
          ResMax = 1
        endif
c
c  Perform the appropriate iteration until no more.
c
        call maxPol(Data(pRes),Data(upRes),nPoint,maxPres)
        Niter = 0
        oNiter = 0
        negFound = .false.
        More = (nPoint.gt.0 .and. ResMin.ne.ResMax .and.
        * uResMin.ne.uResMax)
        Limit = 0
        do while (More)
          oNiter = Niter
c
c  Give some information about steer mode clip level.
c
          if (steermsg .and. moded.eq.'steer') then
            write(line,'(a,f6.3)')'Steer Clip Level:',Clip
            call output(line)
            steermsg = .false.
          endif
          call ComplexSteer(pBem,Data(pRes),Data(upRes),Data(pEst),
     *        Data(upEst),Data(pMap),Data(upMap),
     *        nPoint,nx,ny,Clip*maxPres,Gain,Niter,Run,nRun)
          text = ' Steer'

c
c  Output some messages to assure the user that the computer has not
c  crashed.
c
            
            call Stats(Data(pRes),nPoint,ResMin,ResMax,ResAMax,ResRms)
            call Stats(Data(upRes),nPoint,uResMin,uResMax,uResAMax
     *         ,uResRms)
            call maxPol(Data(pRes),Data(upRes),nPoint,maxPres)
            call SumFlux(Flux,Data(pEst),nPoint)
            call SumFlux(uFlux,Data(upEst),nPoint)
            if ( Clip.le.0.98 .or. (MOD(Niter,100).eq.0 )) Then
            call output(Text//' Iterations: '//itoaf(Niter))
            write(line,'(a,1p3e12.3)')' qResidual min,max,rms: ',
     *         ResMin,ResMax,ResRms
            call output(line)
            write(line,'(a,1p3e12.3)')' uResidual min,max,rms: ',
     *                              uResMin,uResMax,uResRms
            call output(line)
            write(line,'(a,1p3e12.3)')' maxPres, clip*maxPres: ',
     *         maxPres, Clip*maxPres
            call output(line)
            write(line,'(a,1pe12.3)')' Total CLEANed qflux: ',Flux
            call output(line)
            write(line,'(a,1pe12.3)')' Total CLEANed uflux: ',uFlux
            call output(line)
          endif
c
c  If we are doing Steer iterations, see if the next iteration would run
c  into negative components.
c
          if (moded.eq.'steer')
     *      negFound = negFound .or. (-Clip*ResAMax.gt.ResMin)
c
c  Check for convergence.
c
          more = .not.((negFound .and. negStop) .or. (Niter.eq.oNiter)
     *              .or. (ResAMax.le.Cutoff) .or. (Niter.ge.MaxNiter))
        enddo
c
c  Give a message about what terminated the iterations.
c
        if ((ResMin.eq.ResMax) .or. (uResMin.eq.uResMax)) then
          call bug('w','All pixels for this plane are identical')
          call bug('w','No cleaning performed')
          nPoint = 0
        else if (nPoint.eq.0) then
          call output('No region selected in this plane')
       else if (maxPres.le.Cutoff) then
          call output(' Stopping -- Clean cutoff limit reached')
        else if (Niter.ge.MaxNiter) then
          call output(' Stopping -- Maximum iterations performed')
        else if (NegStop .and. NegFound) then
          call output(' Stopping -- Negative components encountered')
        else if (oNiter.eq.Niter) then
          call output(' Stopping -- Could not find more components')
        endif
c
c  Write out this plane.
c
        totniter = totniter + Niter
        call xysetpl(lOut,1,k-blc(3)+1)
        call xysetpl(ulOut,1,k-blc(3)+1)
        call PutPlane(lOut,Run,nRun,xmin-blc(1),ymin-blc(2),
     *                        nOut(1),nOut(2),Data(pEst),nPoint)
        call PutPlane(ulOut,Run,nRun,xmin-blc(1),ymin-blc(2),
     *                        nOut(1),nOut(2),Data(upEst),nPoint)
      enddo
c
c  Construct a header for the output file, and give some history
c  information.
c
      call Header(lMap,lOut,blc,trc,totNiter,minpatch,clip,mode,
     *                                                version)
      call Header(ulMap,ulOut,blc,trc,totNiter,minpatch,clip,mode,
     *                                                version)
c
c  Close up the files. Ready to go home.
c
      call xyclose(lMap)
      call xyclose(ulMap)
      call xyclose(lBeam)
      if (ModelNam.ne.' ') call xyclose(lModel)
      if (uModelNam.ne.' ') call xyclose(ulModel)
      call xyclose(lOut)
      call xyclose(ulOut)
c
c  Thats all folks.
c
      end

c***********************************************************************
      subroutine Stats(Data,n,Dmin,Dmax,DAmax,Drms)

      integer n
      real Data(n)
      real Dmin,Dmax,DAmax,Drms
c-----------------------------------------------------------------------
c  Calculate every conceivably wanted statistic.
c
c  Input:
c    n          Number of points.
c    Data       Input data array.
c
c  Output:
c    Dmin       Data minima.
c    Dmax       Data maxima.
c    DAmax      Data absolute maxima.
c    Drms       Rms value of the data.
c-----------------------------------------------------------------------
      integer i

c     Externals.
      integer ismax,ismin
c-----------------------------------------------------------------------
c
c  Calculate the minima and maxima.
c
      i = ismax(n,Data,1)
      Dmax = Data(i)
      i = ismin(n,Data,1)
      Dmin = Data(i)
      DAmax = max(abs(Dmax),abs(Dmin))
c
c  Calculate the sums.
c
      Drms = 0
      do i = 1, n
        Drms = Drms + Data(i)*Data(i)
      enddo
      Drms = sqrt(Drms/n)

      end

c***********************************************************************
      subroutine maxPol(qData,uData,n,Dmax)

      integer n
      real qData(n),uData(n)
      real Dmax
c-----------------------------------------------------------------------
c  Calculate every conceivably wanted statistic.
c
c  Input:
c    n          Number of points.
c    Data       Input data array.
c
c  Output:
c    Dmax       Data maxima.
c-----------------------------------------------------------------------
      integer i
      real temp
c-----------------------------------------------------------------------
c
c  Calculate the maxima.
c

      Dmax = 0
      do i = 1, n
         temp = qData(i)*qData(i)+uData(i)*uData(i)
         if (temp.gt.Dmax) then
            Dmax = temp
         endif
      enddo
      Dmax = sqrt(Dmax)
      end

c***********************************************************************

      subroutine inputs(qmap,umap,beam,qestimate,uestimate,qout,uout,
     *     Niter,negStop,positive,pad,asym,cutoff,box,maxbox,minpatch,
     *     gain,phat,speed,clip,mode)

      integer Niter, minpatch, maxbox
      integer box(maxbox)
      real cutoff,gain,phat,speed,clip
      logical negStop,positive,pad,asym
      character qmap*(*),umap*(*),beam*(*),qestimate*(*),uestimate*(*)
      character qout*(*),uout*(*),mode*(*)
c-----------------------------------------------------------------------
c       Get user supplied inputs
c
c    Input:
c      maxbox      The maximum number of boxes.
c    Output:
c      Map         Name of the input map. No default.
c      Beam        Name of the input beam. No default.
c      Estimate    Name of the input estimate of the deconvolved image.
c                  Default is a estimate which is zero.
c      Out         Name of the output deconvolved map. No default.
c      Gain        Clean loop gain.  Default 0.5 (a bit high).
c      Cutoff)     The iterations stop when either the absolute max
c                  residual remaining is less than Cutoff, of Niter
c      Niter)     minor iterations have been performed (whichever comes
c                  first).  The default Cutoff is 0 (i.e. iterate
c                  forever), and the default number of minor iterations
c                  is 250.  This is the total number of iterations to
c                  do.
c      negStop     Stop on first negative component.
c      positive    Constrain the output components to be positive
c                  valued.
c      pad         Zero pad the beam to a bigger size.
c      asym        The beam is asymmetric.
c      minpatch    The minimum beam patch width.
c      phat        Prussian hat. Default is 0.
c      speed       Speedup factor. Default is 0.
c      box         The boxes specification.
c      clip        The Steer clip level.
c      mode        Either "clark" (default), "steer" or "any".
c-----------------------------------------------------------------------
      include 'maxdim.h'
c-----------------------------------------------------------------------
      call keyini
      call keya('map', qmap, ' ')
      call keya('map', umap, ' ')
      call keya('beam', beam, ' ')
      call keya('model', qestimate,' ')
      call keya('model', uestimate,' ')
      call keya('out', qout, ' ')
      call keya('out', uout, ' ')
      if (qmap.eq.' ' .or. umap.eq.' ' .or. beam.eq.' ' 
     *     .or. qout.eq.' ' .or. uout.eq.' ' .or.
     *    (qestimate.ne.' '.and.uestimate.eq.' ')) then
         call bug('f', 'A file name was missing from the parameters')
      endif
      call GetOpt(negstop,positive,pad,asym)
      call keyi('niters', Niter, 250)
      if (Niter.le.0) call bug('f', 'NITERS must be positive')
      call keyr('cutoff', cutoff,0.0)

      call BoxInput('region',qmap,box,maxbox)

      call keyi('minpatch', minpatch, 51)
      minpatch = 51
      call keyr('gain', gain, 0.1)
      if (gain.le.0 .or. gain.gt.1)
     *  call bug('f','Bad gain value, it must be in the range (0,1]')
      call keyr('phat', phat, 0.0)
      phat = 0.0
      call keyr('speed', speed, 0.0)
      speed = 0.0
      call keyr('clip', clip, 0.0)
      if (clip.lt.0 .or. clip.gt.1)
     *  call bug('f','Bad clip value, it must be in the range [0,1]')
c      call keya('mode',mode,'steer')
c      if (mode.ne.'clark' .and.
c     *    mode.ne.'steer' .and.
c     *    mode.ne.'any'   .and.
c     *    mode.ne.'hogbom')
      mode='steer'
      call keyfin

      end

c***********************************************************************

      subroutine GetOpt(negstop,positive,pad,asym)

      logical negstop,positive,pad,asym
c-----------------------------------------------------------------------
c  Get extra processing options.
c
c  Output:
c    negstop
c    positive
c-----------------------------------------------------------------------
      integer NOPTS
      parameter (NOPTS = 4)
      character opts(NOPTS)*8
      logical present(NOPTS)
      data opts/'negstop ','positive','pad     ','asym    '/
c-----------------------------------------------------------------------
      call options('options',opts,present,NOPTS)
      negstop  = .false.
      positive = .false.
      pad      = present(3)
      asym     = present(4)

      end

c***********************************************************************

      subroutine Header (lIn, lOut, blc, trc, niters, minpatch, clip,
     *  mode, version)

      integer   lIn, lOut, blc(3), trc(3), niters, minpatch
      real      clip
      character mode*(*), version*(*)
c-----------------------------------------------------------------------
c Copy the header to the model.
c-----------------------------------------------------------------------
      integer   lblc, ltrc
      double precision crpix1, crpix2, crpix3
      character line*72, txtblc*32, txttrc*32

      character itoaf*8
      external  itoaf
c-----------------------------------------------------------------------
c     Copy all header keywords.
      call headcp(lIn, lOut, 0, 0, 0, 0)

c     Update keywords that have changed.
      call rdhdd(lIn,  'crpix1', crpix1, 1d0)
      call rdhdd(lIn,  'crpix2', crpix2, 1d0)
      call rdhdd(lIn,  'crpix3', crpix3, 1d0)
      crpix1 = crpix1 - blc(1) + 1d0
      crpix2 = crpix2 - blc(2) + 1d0
      crpix3 = crpix3 - blc(3) + 1d0
      call wrhdd(lOut, 'crpix1', crpix1)
      call wrhdd(lOut, 'crpix2', crpix2)
      call wrhdd(lOut, 'crpix3', crpix3)
      call wrhda(lOut, 'bunit', 'JY/PIXEL')
      call wrhdi(lOut, 'niters', niters)

c     Write crap to the history file, to attempt (ha!) to appease Neil.
      call hisopen(lOut,'append')
      line = 'CLEAN: Miriad ' // version
      call hiswrite(lOut,line)
      call hisinput(lOut,'CLEAN')

      call mitoaf(blc,3,txtblc,lblc)
      call mitoaf(trc,3,txttrc,ltrc)
      line = 'CLEAN: Bounding region is Blc = (' // txtblc(1:lblc) //
     *       '), Trc = (' // txttrc(1:ltrc) // ')'
      call hiswrite(lOut,line)

      if (mode.eq.'steer' .or. mode.eq.'any') then
        write(line,'(''CLEAN: Steer Clip Level = '',f6.3)') Clip
        call hiswrite(lOut,line)
      endif
      call hiswrite(lOut,'CLEAN: Minpatch = ' // itoaf(minpatch))
      call hiswrite(lOut,'CLEAN: Total Iterations = ' // itoaf(niters))
      call hisclose(lOut)

      end

c***********************************************************************

      subroutine BeamChar(lBeam,n1,n2,ic,jc,Histo,maxPatch)

      integer lBeam,n1,n2,ic,jc,maxPatch
      real Histo(maxPatch/2+1)
c-----------------------------------------------------------------------
c  Determine the location of the max value in the beam, and an array
c  giving the max abs value of the beam outside a particular area.
c
c  Input:
c    lBeam      Handle of the file containing the beam.
c    n1,n2      Beam size.
c    maxPatch   Maximum patch size.
c
c  Output:
c    ic,jc      Pixel coordinate of the max pixel.
c    Histo      Histo(i) is the absolute maximum outside the patch
c               with width 2*i-1 around the beam maximum.
c-----------------------------------------------------------------------
      include 'maxdim.h'
      integer i,j,k,imin,imax,jmin,jmax,nHisto
      real Data(MAXDIM),bmax

c     External.
      integer isamax,ismax
c-----------------------------------------------------------------------
c
c  Initialise.
c
      nHisto = maxPatch/2 + 1
c
c  Determine the value and location of the beam maximum. Just search
c  around the centre of the beam.
c
      bmax = 0
      jmin = max(n2/2-3,1)
      jmax = min(n2/2+3,n2)
      imin = max(n1/2-3,1)
      imax = min(n1/2+3,n1)
      do j = jmin, jmax
        call xyread(lBeam,j,Data)
        i = ismax(imax-imin+1,Data(imin),1) + imin - 1
        if (j.eq.jmin .or. Data(i).gt.bmax) then
           bmax = Data(i)
           ic = i
           jc = j
        endif
      enddo

      if (abs(bmax-1.0).gt.0.001 .and. abs(bmax-1.0).le.0.01)
     *  call bug('w','Beam peak value is not 1')
      if (abs(bmax-1.0).gt.0.01)
     *  call bug('f','Beam peak value differs from 1 by more than 1%')
c
c  Initialise the "histo" array.
c
      do k = 1, nHisto
        Histo(k) = 0
      enddo
c
c  Determine some limits.
c
      jmin = max(1, jc - (nHisto-2))
      jmax = min(n2,jc + (nHisto-2))
      imin = max(1, ic - (nHisto-2))
      imax = min(n1,ic + (nHisto-2))

      do j = 1, n2
        call xyread(lBeam,j,Data)
        if (j.lt.jmin .or. j.gt.jmax) then
          i = isamax(n1,Data,1)
          Histo(nHisto) = max(abs(Data(i)),Histo(nHisto))
        else
          if (imin.gt.1) then
            i = isamax(imin-1,Data,1)
            Histo(nHisto) = max(abs(Data(i)),Histo(nHisto))
          endif
          if (imax.lt.n1) then
            i = isamax(n1-imax,Data(imax+1),1) + imax
            Histo(nHisto) = max(abs(Data(i)),Histo(nHisto))
          endif

          do i = imin, imax
            k = max(abs(i-ic),abs(j-jc)) + 1
            Histo(k) = max(Histo(k),abs(Data(i)))
          enddo
        endif
      enddo
c
c  Now Histo(k) contains the max abs value occurring at distance k.
c  Collapse this do so that Histo(k) contains the max abs value
c  occurring at a distance greater of equal to k.
c
      do k = nHisto-1, 1, -1
        Histo(k) = max(Histo(k),Histo(k+1))
      enddo
c
c  If Histo(1) is greater than bmax, then the sidelobes are greater than
c  the CLEAN peak -- and CLEAN might as well crap out.
c
      if (Histo(1)-bmax.gt.0.001) then
        call bug('w','Beam sidelobes appear bigger than beam peak?')
        call bug('f','Try imaging a larger field of view')
      endif

      end

c***********************************************************************

      subroutine SumFlux(Flux,Estimate,nPoint)

      integer nPoint
      real Estimate(nPoint),Flux
c-----------------------------------------------------------------------
c  Find the sum of the estimate.
c
c  Input:
c    nPoint
c    Estimate
c  Output:
c    Flux
c-----------------------------------------------------------------------
      integer i
c-----------------------------------------------------------------------
      Flux = 0
      do i = 1, nPoint
        Flux = Flux + Estimate(i)
      enddo

      end

c***********************************************************************

      subroutine SumAbs(EstASum,Estimate,nPoint)

      integer nPoint
      real Estimate(nPoint),EstASum
c-----------------------------------------------------------------------
c  Find the sum of the absolute value of the estimate.
c
c  Input:
c    nPoint
c    Estimate
c  Output:
c    EstASum
c-----------------------------------------------------------------------
      integer i
c-----------------------------------------------------------------------
      EstASum = 0
      do i = 1, nPoint
        EstASum = EstASum + abs(Estimate(i))
      enddo
      end

c***********************************************************************

      subroutine NoModel(Map,Estimate,Residual,nPoint)

      integer nPoint
      real Map(nPoint),Estimate(nPoint),Residual(nPoint)
c-----------------------------------------------------------------------
c  This initialises the estimate and the residuals, for the case where
c  there is no model.
c
c  Input:
c    nPoint     Number of points.
c    Map        The original dirty map.
c  Output:
c    Residual   The residuals, which are the same as the dirty map.
c    Estimate   The estimate, which are initially zero.
c-----------------------------------------------------------------------
      integer i
c-----------------------------------------------------------------------
      do i = 1, nPoint
        Estimate(i) = 0.0
        Residual(i) = Map(i)
      enddo

      end


c***********************************************************************

      subroutine ComplexSteer(pBem,qResidual,uResidual,qEstimate,
     *     uEstimate,qTemp,uTemp,nPoint,nx,ny,
     *     Limit,Gain,Niter,Run,nrun)

      integer Niter,nrun,Run(3,nrun),nPoint,nx,ny
      ptrdiff pBem
      real Gain,Limit
      real qResidual(nPoint),qTemp(nPoint),qEstimate(nPoint)
      real uResidual(nPoint),uTemp(nPoint),uEstimate(nPoint)
c-----------------------------------------------------------------------
c  Perform an iteration of the Complex Steer "Clean".
c
c  Input:
c    Run        The Run specifying the area of interest.
c    nrun       The number of runs. A null run follows run number nrun.
c    limit      The limit to cut at.
c    clip       Steer clip level.
c
c  Input/Output:
c    Niter      Number of Niter iterations.
c
c-----------------------------------------------------------------------
      real MinOptGain
      parameter (MinOptGain = 0.02)
      integer i
      real qSumRE,SumEE,qg,uSumRE,ug,zr,zi,mag
c-----------------------------------------------------------------------
c
c  Form the new Steer estimate.
c
      do i = 1, nPoint
         if ((qResidual(i)*qResidual(i)
         *   +uResidual(i)*uResidual(i)).ge.
         *      Limit*Limit) then
          qTemp(i) = qResidual(i)
          uTemp(i) = uResidual(i)
          Niter = Niter + 1
        else
          qTemp(i) = 0.0
          uTemp(i) = 0.0
        endif
      enddo
c
c  Convolve it with the dirty beam.
c
      call CnvlR(pBem,qTemp,nx,ny,Run,nRun,qTemp,'c')
      call CnvlR(pBem,uTemp,nx,ny,Run,nRun,uTemp,'c')      
c
c  "Temp" now contains the new estimate convolved with the dirty beam.
c  Determine the optimum gain to apply to minimise the residuals.
c  However apply the users damping factor to this gain.
c
      qSumRE = 0
      SumEE = 0
      uSumRE = 0
      zr = 0
      zi = 0
      do i = 1, nPoint
        Call ComplexMultiply(qResidual(i),uResidual(i),qTemp(i),
        *   -uTemp(i),zr,zi)
        qSumRE = qSumRE + zr
        uSumRE = uSumRE + zi
        SumEE = SumEE + qTemp(i)*qTemp(i) + uTemp(i)*uTemp(i)
      enddo
c
c       SumRE can be negative, so it is better to take the
c       absolute value of it when determining the optimum
c       gain (gmx - 07mar04)
c
c       abs(SumRE)/SumEE may be close to zero, in which case
c       a semi-infinite loop can be the result. We apply a
c       lower limit to abs(SumRE)/SumEE. A good value for it
c       is empirically determined to be 0.02 (MinOptGain),
c       which may however not be the best choice in all cases.
c       In case of problems, you can try a lower value for the
c       task option Gain before changing MinOptGain (gmx - 07mar04).
c
c       Choosing to scale the complex damping by magnitude rather than 
c       individual components. This will keep the method consistent under
c       choice of Q and U axes.
      mag = (qSumRE*qSumRE + uSumRE*uSumRE)/(SumEE*SumEE)
      qg = Gain*max(MinOptGain,mag)*(qSumRE/SumEE)/mag
      ug = Gain*max(MinOptGain,mag)*(uSumRE/SumEE)/mag
c
c  Now go through and update the estimate, using this gain. Also
c  determine the new residuals.
c
      do i = 1, nPoint
         if ((qResidual(i)*qResidual(i)
         *   +uResidual(i)*uResidual(i)).ge.
         *      Limit*Limit) then
            Call ComplexMultiply(qg, ug, qResidual(i), uResidual(i), zr
              * , zi)
            qEstimate(i) = qEstimate(i) + zr
            uEstimate(i) = uEstimate(i) + zi
         endif
         Call ComplexMultiply(qg, ug, qTemp(i), uTemp(i), zr, zi)        
         qResidual(i) = qResidual(i) - zr
         uResidual(i) = uResidual(i) - zi 
      enddo

      end

c***********************************************************************
      subroutine ComplexMultiply(pr, pi, qr, qi, zr, zi)

      real pr, pi, qr, qi, zr, zi
c-----------------------------------------------------------------------
c  Multiply two complex numbers p and q, to give an output of z
c
c  Input:
c    pr, qr     Real parts of each input
c    pi, q      Imaginary parts of each input
c  Output:
c    zr, zi     Real and imaginary part of output
c-----------------------------------------------------------------------
      zr = pr*qr - pi*qi
      zi = pr*qi + pi*qr
      end

c***********************************************************************

      subroutine Diff(pBem,Estimate,Map,Residual,nPoint,nx,ny,
     *  Run,nRun)

      integer nPoint,nx,ny,nRun,Run(3,nRun)
      ptrdiff pBem
      real Estimate(nPoint),Map(nPoint),Residual(nPoint)
c-----------------------------------------------------------------------
      integer i
c-----------------------------------------------------------------------
      call CnvlR(pBem,Estimate,nx,ny,Run,nRun,Residual,'c')

      do i = 1, nPoint
        Residual(i) = Map(i) - Residual(i)
      enddo

      end

c***********************************************************************

      subroutine defregio(boxes,nMap,nBeam,icentre,jcentre)

      integer boxes(*),nMap(3),nBeam(2),icentre,jcentre
c-----------------------------------------------------------------------
c  Set the region of interest to the lastest area that can be safely
c  deconvolved.
c-----------------------------------------------------------------------
      integer blc(3),trc(3),width
c-----------------------------------------------------------------------
      width  = min(icentre-1,nBeam(1)-icentre) + 1
      blc(1) = max(1,(nMap(1)-width)/2)
      trc(1) = min(nMap(1),blc(1)+width-1)

      width  = min(jcentre-1,nBeam(2)-jcentre) + 1
      blc(2) = max(1,(nMap(2)-width)/2)
      trc(2) = min(nMap(2),blc(2)+width-1)

      blc(3) = 1
      trc(3) = nMap(3)

      call BoxDef(boxes,3,blc,trc)

      end
