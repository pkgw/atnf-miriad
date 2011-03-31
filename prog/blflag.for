      program blflag

c= blflag -- Interactive flagging task.
c& rjs
c: uv analysis
c+
c       BLFLAG is a Miriad task for flagging visibilities interactively.
c       It plots visibilities (e.g. amplitude vs time), either one
c       baseline at a time or all together, and allows discrepant points
c       to be flagged with a cursor.
c
c       Commands are entered as a single character at the keyboard:
c         Left-Button  Left mouse button flags the nearest visibility.
c         Right-Button Right mouse button causes BLFLAG to precede to
c                      the next baseline.
c         <CR>         Carriage-return gives help.
c         ?            Help.
c         a            Flag nearest visibility (same as left mouse
c                      button).
c         c            Clear the flagging for this baseline and redraw
c                      plot.
c         h            Give help (same as carriage return).
c         p            Define a polygonal region, and flag visibilities
c                      within this region.  Define the vertices of the
c                      polygon by moving the cursor and then hitting the
c                      left mouse button (or a).  Finish defining the
c                      polygon by hitting the right mouse button (or x).
c                      You can delete vertices with the middle mouse
c                      button (or d).
c         q            Abort completely.  This does not apply flagging.
c         r            Redraw plot.
c         u            Unzoom.
c         x            Move to next baseline (same as right mouse
c                      button).
c         z            Zoom in.  You follow this by clicking the mouse
c                      on the left and right limits to zoom.
c
c@ vis
c       Input visibility dataset to be flagged.  No default.
c
c@ line
c       The normal Miriad linetype specification.  BLFLAG averages all
c       channels together before displaying them, flags being applied
c       to all channels.  The default is all channels.
c
c@ device
c       Normal PGPLOT plot device.  An interactive device, e.g. /xserve,
c       must be selected.  No default.
c
c@ stokes
c       Normal Stokes/polarisation parameter selection.  The default is
c       'ii' (i.e. Stokes-I assuming the source is unpolarised).  NOTE
c       BLFLAG plots the average of all the selected Stokes/polarisation
c       quantities.  Also it flags ALL quantities, regardless of whether
c       they were selected or not.
c
c@ select
c       Normal visibility data selection.  Only selected data can be
c       flagged.  The default is to select all data.
c
c@ axis
c       Two character strings, giving the X and Y axes of the plot.
c       Possible axis values are:
c         time         (the default for the X axis)
c         lst          Local apparent sidereal time.
c         uvdistance   sqrt(u**2+v**2)
c         hangle       (hour angle)
c         channel      (implies nofqav)
c         amplitude    (the default for the Y axis)
c         phase
c         real
c         imaginary
c         rms          Theoretical rms noise.
c
c@ xrange
c       Plot range in the x-direction
c         If axis = uvdistance            [kilo-lambda;   2 values]
c         If axis = time                  [dd,hh,mm,ss.s; 8 values]
c         If axis = amplitude, real, imag [natural units; 2 values]
c         If axis = phase                 [degrees;       2 values]
c         If axis = hangle                [hh,mm,ss.s;    6 values]
c         If axis = rms                   [flux units;    2 values]
c         If axis = lst                   [decimal hours; 2 values]
c         If axis = freq                  [GHz;           2 values]
c
c       For axis types other than 'time' or 'hangle', one or other of
c       the limits may be set with the other self-scaled by specifying
c       the lower limit as 'min' and the upper as 'max' (or simply
c       omitting it).  For example,
c
c         xrange=min,0
c
c       self-scales the lower limit while pinning the upper limit to
c       zero, whereas either of the following
c
c         xrange=0,max
c         xrange=0
c
c       set the lower limit to zero while self-scaling the upper limit.
c
c       Default is to self-scale both limits.
c
c@ yrange
c       Plot range for the y-axis as for the x-axis.  The default is to
c       self-scale.  For amplitude type plots you can greatly reduce the
c       number of points to plot by using something like yrange=0.3 to
c       cut out noise.
c
c@ options
c       Task enrichment parameters.  Several can be given, separated by
c       commas.  Minimum match is used.  Possible values are:
c         nobase  Normally BLFLAG plots a single baseline at a time.
c                 This option causes all baselines to be plotted on
c                 a single plot.
c         selgen  Generate a file appropriate for selecting the bad
c                 data (via a "select" keyword).  The output is a text
c                 file called "blflag.select".
c         noapply Do not apply the flagging.
c         rms     When processing spectra, blflag normally plots the
c                 mean value of the spectra.  Using options=rms causes
c                 it to plot the rms value instead.
c         scalar  When processing spectra, blflag normally forms an
c                 average value by vector averaging.  This option
c                 causes it to generate the scalar average.  It should
c                 be used with significant caution.
c         nofqaver Do not average spectra - the resulting number of
c                 points may be too large to handle.  Use select to
c                 break up the data in time ranges or use yrange to
c                 exclude noise.
c       The following options can be used to disable calibration.
c         nocal   Do not apply antenna gain calibration.
c         nopass  Do not apply bandpass correction.
c         nopol   Do not apply polarisation leakage correction.
c
c$Id$
c--
c  History:
c    Refer to the RCS log, v1.1 includes prior revision information.
c-----------------------------------------------------------------------
      include 'maxdim.h'
      include 'mem.h'

c     Before increasing any of these, particularly MAXDAT, be sure that
c     it does not cause linking to fail with truncated relocations.
      integer    MEBI, MAXDAT, MAXPLT, MAXEDIT
      parameter (MEBI    = 1024*1024,
     *           MAXDAT  = 32*MEBI,
     *           MAXPLT  =  4*MEBI,
     *           MAXEDIT =    MEBI)

      logical   noapply, nobase, nofqaver, present(MAXBASE), rms,
     *          scalar, selgen
      integer   i, j, k, length, npol, tno
      real      xmax, xmin, ymax, ymin
      character device*64, title*32, uvflags*12, val*16, version*72,
     *          xaxis*12, yaxis*12

c     Data store 768MiB ((4*4+8) * 32*MEBI).
      logical   ltemp(MAXDAT)
      integer   bldat(MAXDAT), chdat(MAXDAT), ndat
      real      xdat(MAXDAT), ydat(MAXDAT)
      double precision timedat(MAXDAT)

c     Plot buffer 96MiB ((4*4+8) * 4*MEBI).
      integer   blplt(MAXPLT), chplt(MAXPLT), nplt
      real      xplt(MAXPLT), yplt(MAXPLT)
      double precision timeplt(MAXPLT)

c     Editing buffer 16MiB ((2*4+8) * MEBI).
      integer   bledit(MAXEDIT), chedit(MAXEDIT), nedit
      double precision timeedit(MAXEDIT)

      external  itoaf, len1, pgbeg, uvDatOpn, versan
      logical   uvDatOpn
      integer   len1, pgbeg
      character itoaf*10, versan*72
c-----------------------------------------------------------------------
      version = versan('blflag',
     *                 '$Revision$',
     *                 '$Date$')

c     Get the input parameters.
      call keyini
      call keya('device',device,' ')
      if (device.eq.' ') call bug('f','A PGPLOT device must be given')
      call GetAxis(xaxis,yaxis)
      call GetOpt(nobase,selgen,noapply,rms,scalar,nofqaver,
     *  uvflags)
      if (xaxis.eq.'channel' .or. yaxis.eq.'channel') nofqaver=.true.

c     Get axis ranges
      call getrng('xrange', xaxis, xmin, xmax)
      call getrng('yrange', yaxis, ymin, ymax)
      call uvDatInp('vis',uvflags)
      call keyfin

c     Set the default polarisation type if needed.
      call uvDatGti('npol',npol)
      if (npol.eq.0) call uvDatSet('stokes',0)

c     Open the input data.
      if (.not.uvDatOpn(tno)) call bug('f','Error opening input')

c     Open the plot device.
      if (pgbeg(0,device,1,1).ne.1) then
        call pgldev
        call bug('f','Unable to open PGPLOT device')
      endif
      call pgqinf('CURSOR',val,length)
      if (val.eq.'NO') call bug('f','PGPLOT device is not interactive')
      call pgask(.false.)

c     Get the data.
      call GetDat(tno,rms,scalar,nofqaver,xaxis,yaxis,xmin,xmax,
     *  ymin,ymax,MAXBASE,present,MAXDAT,xdat,ydat,bldat,chdat,timedat,
     *  ndat)
      call uvDatCls
      if (ndat.eq.0) call bug('f','No points to flag')

c     Loop over the baselines.
      call output('Entering interactive mode ...')
      nedit = 0
      if (nobase) then
        call output('Processing '//itoaf(ndat)//' points')
        call Edit(ndat,xdat,ydat,bldat,chdat,timedat,ltemp,xaxis,yaxis,
     *    'All baselines',MAXEDIT,timeedit,bledit,chedit,nedit)
      else
        k = 0
        do j = 1, MAXANT
          do i = 1, j
            k = k + 1
            if (present(k)) then
              title = 'Baseline '//itoaf(i)
              length = len1(title)
              title(length+1:) = '-'//itoaf(j)
              call Extract(k,ndat,xdat,ydat,bldat,chdat,timedat,
     *          MAXPLT,xplt,yplt,blplt,chplt,timeplt,nplt)
              if (nplt.gt.0) then
                call Edit(nplt,xplt,yplt,blplt,chplt,timeplt,ltemp,
     *            xaxis,yaxis,title,MAXEDIT,timeedit,bledit,chedit,
     *            nedit)
              endif
            endif
          enddo
        enddo
      endif

      call pgend

c     Generate the "blflag.select" file, if needed.
      if (selgen) then
        if (nedit.eq.0) then
          call bug('w','No edit commands to write out!')
        else
          call doSelGen(nedit,timeedit,bledit,chedit)
        endif
      endif

c     Apply the changes.
      if (nedit.gt.0 .and. .not.noapply) then
        call output('Applying the flagging ...')
        call uvDatRew
        call uvDatSet('disable',0)
        if (.not.uvDatOpn(tno)) call bug('f','Error reopening input')
        call FlagApp(tno,nedit,timeedit,bledit,chedit,version)
        call uvDatCls
      endif

      end

c***********************************************************************

      subroutine doSelGen(n,time,bl,chn)

      integer   N
      double precision time(N)
      integer   bl(N),chn(N)
c-----------------------------------------------------------------------
c  Generate a file of select commands.
c-----------------------------------------------------------------------
      double precision TTOL
      parameter (TTOL=1d0/86400d0)
      include 'maxdim.h'
      integer i,j,k,lu,iostat,length
      integer i1(MAXBASE),i2(MAXBASE)
      character line*80,time1*24,time2*24
      logical warn

      external  itoaf, len1
      integer   len1
      character itoaf*4
c-----------------------------------------------------------------------
c     Generate a table to allow translation from baseline number to
c     antenna pairs.
      k = 0
      do j = 1, MAXANT
        do i = 1, j
          k = k + 1
          i1(k) = i
          i2(k) = j
        enddo
      enddo

      call txtopen(lu,'blflag.select','new',iostat)
      if (iostat.ne.0) then
        call bug('w','Error opening output text file blflag.select')
        call bugno('f',iostat)
      endif

      warn=.true.
      do i = 1, N
        if (chn(i).eq.0) then
          call julday(time(i)-TTOL,'H',time1)
          call julday(time(i)+TTOL,'H',time2)
          line = 'ant('//itoaf(i1(bl(i)))
          length = len1(line)
          line(length+1:) = ')('//itoaf(i2(bl(i)))
          length = len1(line)
          line(length+1:) = '),time('//time1
          length = len1(line)
          line(length+1:) = ','//time2
          length = len1(line)
          line(length+1:) = ')'
          length = length + 1
          if (i.ne.n) then
            line(length+1:length+3) = ',or'
            length = length + 3
          endif
          call txtwrite(lu,line,length,iostat)
          if (iostat.ne.0) then
            call bug('w','Error writing to text file blflag.select')
            call bugno('f',iostat)
          endif
        else
          if (warn) then
            call bug('w',
     *       'Omitting channel specific flags in text file')
            warn=.false.
          endif
        endif
      enddo

      call txtclose(lu)

      end

c***********************************************************************

      subroutine FlagApp(tno,NEDIT,timeedit,bledit,chedit,version)

      integer   tno, NEDIT
      double precision timeedit(NEDIT)
      integer   bledit(NEDIT),chedit(NEDIT)
      character version*(*)
c-----------------------------------------------------------------------
c  Apply flagging to the dataset.
c-----------------------------------------------------------------------
      include 'maxdim.h'
      double precision TTOL
      parameter (TTOL=1d0/86400d0)

      integer nchan,bl,i1,i2,i,k
      double precision preamble(4),time
      complex data(MAXCHAN)
      logical flags(MAXCHAN),match,chflag
      integer nflag,ncorr
      character line*64

      external  itoaf
      character itoaf*10
c-----------------------------------------------------------------------
      nflag = 0
      ncorr = 0

      call uvDatRd(preamble,data,flags,MAXCHAN,nchan)
      do while (nchan.gt.0)
        ncorr = ncorr + nchan
        time = preamble(3)
        call basant(preamble(4),i1,i2)
        bl = (i2*(i2-1))/2 + i1
        chflag=.false.

c       Search for this integration.
        do i = 1, NEDIT
          match = bledit(i).eq.bl .and.
     *               abs(timeedit(i)-time).le.TTOL

c         Flag bad channels if found.
          if (match) then
            if (chedit(i).eq.0) then
              do k = 1, nchan
                if (flags(k)) then
                  flags(k) = .false.
                  nflag = nflag + 1
                  chflag=.true.
                endif
              enddo
            else
              if (flags(chedit(i))) then
                flags(chedit(i))=.false.
                nflag = nflag +1
                chflag=.true.
              endif
            endif
          endif
        enddo
        if (chflag) call uvflgwr(tno,flags)

c       Go back for more.
        call uvDatRd(preamble,data,flags,MAXCHAN,nchan)
      enddo

c     Write history.
      call hisOpen(tno,'append')
      line = 'BLFLAG: Miriad '//version
      call hisWrite(tno,line)
      call hisInput(tno,'BLFLAG')
      call hisWrite(tno,'BLFLAG: Number of correlations flagged: '//
     *                        itoaf(nflag))
      call hisClose(tno)

c     Report how much we have done.
      call output('Total number of correlations:   '//itoaf(ncorr))
      call output('Number of correlations flagged: '//itoaf(nflag))

      end

c***********************************************************************

      subroutine Edit(NPLT,xplt,yplt,blplt,chplt,timeplt,flag,
     *  xaxis,yaxis,title,MAXEDIT,timeedit,bledit,chedit,nedit)

      integer   NPLT
      real      xplt(NPLT), yplt(NPLT)
      integer   blplt(NPLT), chplt(NPLT)
      double precision timeplt(NPLT)
      logical   flag(NPLT)
      character xaxis*(*), yaxis*(*), title*(*)
      integer   MAXEDIT
      double precision timeedit(MAXEDIT)
      integer   bledit(MAXEDIT), chedit(MAXEDIT)
      integer   nedit
c-----------------------------------------------------------------------
      character mode*1
      integer nedit0,i
      logical more
      real xv,yv,xs,ys,xmin,xmax
c-----------------------------------------------------------------------
      nedit0 = nedit
      mode = 'c'
      more = .true.
      do while (more)
        call lcase(mode)
        if (mode.eq.'q') then
          call pgend
          call bug('f','Aborting at users request')
        else if (mode.eq.'c') then
          nedit = nedit0
          xmin = 0
          xmax = 0
          do i = 1, nplt
            flag(i) = .true.
          enddo
          call Draw(xmin,xmax,nplt,xplt,yplt,flag,xaxis,yaxis,title,
     *      xs,ys)

        else if (mode.eq.'r') then
          call Draw(xmin,xmax,nplt,xplt,yplt,flag,xaxis,yaxis,title,
     *      xs,ys)

        else if (mode.eq.'h' .or. mode.le.' ' .or. mode.eq.'?') then
          call output('-------------------------------------')
          call output('Single key commands are')
          call output(' Left-button  Delete nearest point')
          call output(' Right-button Next baseline')
          call output(' <CR>      Help')
          call output(' ?         Help')
          call output(' a         Delete nearest point')
          call output(' c         Clear flagging of this baseline')
          call output(' h         Help -- these messages')
          call output(' p         Define and delete polygonal region')
          call output(' q         Quit -- Abort completely')
          call output(' r         Redraw')
          call output(' u         Unzoom')
          call output(' x         Next baseline')
          call output(' z         Zoom in')
          call output('-------------------------------------')
        else if (mode.eq.'u') then
          xmin = 0
          xmax = 0
          call Draw(xmin,xmax,nplt,xplt,yplt,flag,xaxis,yaxis,title,
     *      xs,ys)

        else if (mode.eq.'x') then
          more = .false.

        else if (mode.eq.'z') then
          call output('Click on left-hand edge of the zoomed region')
          call pgcurs(xmin,yv,mode)
          call output('Click on right-hand edge of the zoomed region')
          call pgcurs(xmax,yv,mode)
          call Draw(xmin,xmax,nplt,xplt,yplt,flag,xaxis,yaxis,title,
     *      xs,ys)
        else if (mode.eq.'a') then
          call Nearest(xv,yv,xs,ys,nplt,xplt,yplt,blplt,chplt,timeplt,
     *      flag,MAXEDIT,bledit,chedit,timeedit,nedit)
        else if (mode.eq.'p') then
          call Region(nplt,xplt,yplt,blplt,chplt,timeplt,flag,
     *      MAXEDIT,bledit,chedit,timeedit,nedit)
        else
          call bug('w','Unrecognised keystroke - use h for help')
        endif
        if (more) call pgcurs(xv,yv,mode)
      enddo

      end

c***********************************************************************

      subroutine Region(NPLT,xplt,yplt,blplt,chplt,timeplt,flag,
     *  MAXEDIT,bledit,chedit,timeedit,nedit)

      integer   NPLT
      real      xplt(NPLT), yplt(NPLT)
      integer   blplt(NPLT),chplt(NPLT)
      double precision timeplt(nplt)
      logical   flag(NPLT)
      integer   MAXEDIT
      integer   bledit(MAXEDIT),chedit(MAXEDIT)
      double precision timeedit(MAXEDIT)
      integer   nedit
c-----------------------------------------------------------------------
      integer MAXV
      parameter (MAXV=100)
      real xv(MAXV+1),yv(MAXV+1)
      integer nv,i,j
      logical within
c-----------------------------------------------------------------------
      call output('Define a region - exit with x')
      call pgsci(3)
      nv = 0
      call pgolin(MAXV,nv,xv,yv,17)
      if (nv.lt.3) then
        call bug('w','Too few vertices')
      else if (nv.gt.MAXV) then
        call bug('w','Too many vertices for me!')
      else
        call pgsfs(2)
        call pgslw(2)
        call pgpoly(nv,xv,yv)
        call pgslw(1)
        call pgsci(2)

        xv(nv+1) = xv(1)
        yv(nv+1) = yv(1)

c       Find all points that are within the poly.
        call pgbbuf
        do i = 1, nplt
          if (flag(i)) then
            within = .false.
            do j = 1, nv
              if ((xplt(i)-xv(j))*(xplt(i)-xv(j+1)).le.0 .and.
     *           abs(xplt(i)-xv(j))*(yv(j+1)-yv(j)).lt.
     *           abs(xv(j+1)-xv(j))*(yplt(i)-yv(j)))
     *           within = .not.within
            enddo

            if (within) then
              nedit = nedit + 1
              if (nedit.gt.MAXEDIT)
     *          call bug('f','Too many editing ops')
              bledit(nedit)   = blplt(i)
              chedit(nedit)   = chplt(i)
              timeedit(nedit) = timeplt(i)
              flag(i) = .false.
              call pgpt(1,xplt(i),yplt(i),1)
            endif
          endif
        enddo
        call pgebuf
      endif

      end

c***********************************************************************

      subroutine Draw(xmin,xmax,nplt,xplt,yplt,flag,xaxis,yaxis,title,
     *  xs,ys)

      real      xmin, xmax
      integer   NPLT
      real      xplt(NPLT), yplt(NPLT)
      logical   flag(NPLT)
      character xaxis*(*), yaxis*(*), title*(*)
      real      xs, ys
c-----------------------------------------------------------------------
      real xlo,xhi,ylo,yhi
      integer i,n
      logical good
      character xtitle*32,ytitle*32,xflags*12,yflags*12
c-----------------------------------------------------------------------
c     Determine the min and max values.
      call pgbbuf
      call SetUp(NPLT,xplt,flag,xaxis,xlo,xhi,xtitle,xflags)
      if (xmin.lt.xmax) then
        xlo = xmin
        xhi = xmax
      endif
      call SetUp(NPLT,yplt,flag,yaxis,ylo,yhi,ytitle,yflags)

      xs = 1/(xhi-xlo)**2
      ys = 1/(yhi-ylo)**2

c     Draw the plot.
      call pgsci(1)
      call pgpage
      if ((xaxis.eq.'real' .or. xaxis.eq.'imaginary') .and.
     *   (yaxis.eq.'real' .or. yaxis.eq.'imaginary')) then
        call pgvstd
        call pgwnad(xlo,xhi,ylo,yhi)
      else
        call pgvstd
        call pgswin(xlo,xhi,ylo,yhi)
      endif
      call pgtbox(xflags,0,0.0,yflags,0,0.0)
      call pglab(xtitle,ytitle,title)

c     Plot all the good data.
      n = 1
      good = flag(1)
      do i = 2, NPLT
        if (good.neqv.flag(i)) then
          if (good) then
            call pgpt(i-n,xplt(n),yplt(n),1)
          else
            n = i
          endif
          good = flag(i)
        endif
      enddo
      if (good) call pgpt(NPLT-n+1,xplt(n),yplt(n),1)

c     Change the colour to red.
      call pgsci(2)
      call pgebuf

      end

c***********************************************************************

      subroutine SetUp(NPLT,plt,flag,axis,lo,hi,title,flags)

      integer   NPLT
      real      plt(NPLT)
      logical   flag(NPLT)
      character axis*(*)
      real      lo, hi
      character title*(*), flags*(*)
c-----------------------------------------------------------------------
      integer i
      logical first
      real x1,x2,delta,absmax
c-----------------------------------------------------------------------
      first = .true.
      do i = 1, NPLT
        if (flag(i)) then
          if (first) then
            x1 = plt(i)
            x2 = x1
            first = .false.
          else
            x1 = min(x1,plt(i))
            x2 = max(x2,plt(i))
          endif
        endif
      enddo

      delta = 0.05*(x2-x1)
      absmax = max(abs(x2),abs(x1))
      if (delta.le.1e-4*absmax) delta = 0.01*absmax
      if (delta.eq.0) delta = 1
      lo = x1 - delta
      hi = x2 + delta

      if (axis.eq.'time' .or. axis.eq.'lst' .or. axis.eq.'hangle') then
        flags = 'BCNSTHZ0'
      else
        flags = 'BCNST'
      endif

      if (axis.eq.'uvdistance') then
        title = '(u\u2\d+v\u2\d)\u1/2\d (k\gl)'
      else if (axis.eq.'phase') then
        title = 'Phase (degrees)'
      else if (axis.eq.'rms') then
        title = 'Theoretical rms noise'
      else
        title = axis
        call ucase(title(1:1))
      endif

      end

c***********************************************************************

      subroutine Nearest(xv,yv,xs,ys,NPLT,xplt,yplt,blplt,chplt,timeplt,
     *  flag,MAXEDIT,bledit,chedit,timeedit,nedit)

      real    xv, yv, xs, ys
      integer NPLT
      real    xplt(NPLT), yplt(NPLT)
      integer blplt(NPLT),chplt(NPLT)
      double precision timeplt(NPLT)
      logical flag(NPLT)
      integer MAXEDIT
      integer bledit(MAXEDIT), chedit(MAXEDIT)
      double precision timeedit(MAXEDIT)
      integer nedit
c-----------------------------------------------------------------------
      integer i,k
      logical first
      real r2,r2d
c-----------------------------------------------------------------------
      first = .true.
      r2 = 0

      do i = 1, NPLT
        if (flag(i)) then
          r2d = xs*(xplt(i)-xv)**2 + ys*(yplt(i)-yv)**2
          if (first .or. r2d.lt.r2) then
            r2 = r2d
            k = i
            first = .false.
          endif
        endif
      enddo

      if (first) then
        call bug('w','No points left to edit')
      else
        nedit = nedit + 1
        if (nedit.gt.MAXEDIT) call bug('f','Too many ops')
        bledit(nedit)   = blplt(k)
        timeedit(nedit) = timeplt(k)
        chedit(nedit) = chplt(k)
        flag(k) = .false.
        call pgpt(1,xplt(k),yplt(k),1)
      endif

      end

c***********************************************************************

      subroutine Extract(k,NDAT,xdat,ydat,bldat,chdat,timedat,
     *  MAXPLT,xplt,yplt,blplt,chplt,timeplt,nplt)

      integer   k, NDAT
      real      xdat(NDAT), ydat(NDAT)
      integer   bldat(NDAT), chdat(NDAT)
      double precision timedat(NDAT)

      integer   MAXPLT
      real      xplt(MAXPLT),yplt(MAXPLT)
      integer   blplt(MAXPLT), chplt(MAXPLT)
      double precision timeplt(MAXPLT)

      integer   nplt
c-----------------------------------------------------------------------
      integer i
c-----------------------------------------------------------------------
      nplt = 0
      do i = 1, ndat
        if (bldat(i).eq.k) then
          nplt = nplt + 1
          if (nplt.gt.MAXPLT) call bug('f','Too many points')
          blplt(nplt) = k
          timeplt(nplt) = timedat(i)
          chplt(nplt)   = chdat(i)
          xplt(nplt)    = xdat(i)
          yplt(nplt)    = ydat(i)
        endif
      enddo

      end

c***********************************************************************

      subroutine GetDat(tno,rms,scalar,nofqaver,xaxis,yaxis,xmin,xmax,
     *  ymin,ymax,MAXBASE1,present,MAXDAT,xdat,ydat,bldat,chdat,
     *  timedat,ndat)

      integer   tno
      logical   rms, scalar, nofqaver
      character xaxis*(*), yaxis*(*)
      real      xmin, xmax, ymin, ymax
      integer   MAXBASE1
      logical   present(MAXBASE1)
      integer   MAXDAT
      real      xdat(MAXDAT), ydat(MAXDAT)
      integer   bldat(MAXDAT), chdat(MAXDAT)
      double precision timedat(MAXDAT)
      integer   ndat
c-----------------------------------------------------------------------
      include 'maxdim.h'

      double precision TTOL
      parameter (TTOL=1d0/86400d0)

c     The default MAXCHAN makes blflag fail to link, use one that will
c     fit a single spectrum
      integer    MAXCHAN1
      parameter (MAXCHAN1=15000)

      logical flags(MAXCHAN1),ok
      complex data(MAXCHAN1)
      complex corr(MAXBASE),corr1(MAXBASE),corr2(MAXBASE)
      complex fcorr(MAXBASE,MAXCHAN1),fcorr1(MAXBASE,MAXCHAN1),
     *  fcorr2(MAXBASE,MAXCHAN1)
      double precision preamble(4),time,time0,tprev,lst,ra
      real uvdist2(MAXBASE),var(MAXBASE),temp
      integer i,j,n,bl,i1,i2,nants,npnt(MAXBASE),
     *  fnpnt(MAXBASE,MAXCHAN1),mbase,nchan,nchanprev
c-----------------------------------------------------------------------
c     Miscellaneous initialisation.
      mbase = min(MAXBASE,maxbase1)
      do i = 1, MAXBASE
        present(i) = .false.
      enddo

      do i = 1, mbase
        npnt(i)    = 0
        uvdist2(i) = 0
        corr(i)    = 0
        corr1(i)   = 0
        corr2(i)   = 0
        var(i) = 0
        do j = 1, maxchan1
          fnpnt(i,j)  = 0
          fcorr(i,j)  = 0
          fcorr1(i,j) = 0
          fcorr2(i,j) = 0
        enddo
      enddo
      ndat = 0

c     Let's get going.
      call output('Reading the data ...')
      call uvDatRd(preamble,data,flags,MAXCHAN1,nchan)
      if (nchan.eq.0) call bug('f','No visibility data found')
      if (nchan.eq.MAXCHAN1) call bug('f','Too many channels for me')
      call flagchk(tno)
      nants = 0
      nchanprev = 0
      tprev = preamble(3)
      time0 = int(tprev - 0.5d0) + 0.5d0
      call uvrdvrd(tno,'lst',lst,0d0)
      call uvrdvrd(tno,'ra',ra,0d0)
      do while (nchan.gt.0)
        call BasAnt(preamble(4),i1,i2)
        bl = (i2*(i2-1))/2 + i1
        ok = bl.lt.mbase
        if (ok) then
          time = preamble(3)
          if (nofqaver) then
            if (abs(time-tprev).gt.TTOL) then
              if (nants.gt.0) call IntFlushF(nants,rms,scalar,ra,lst,
     *          tprev,uvdist2,var,fcorr,fcorr1,fcorr2,xaxis,yaxis,
     *          xmin,xmax,ymin,ymax,fnpnt,time0,present,mbase,
     *          xdat,ydat,timedat,bldat,chdat,ndat,MAXDAT,nchanprev)
              nants = 0
              tprev = time
              call uvrdvrd(tno,'lst',lst,0d0)
              call uvrdvrd(tno,'ra',ra,0d0)
            endif

            n = 0
            do i = 1, nchan
              if (flags(i)) then
                n = n + 1
                fnpnt(bl,i) = fnpnt(bl,i) + 1
                fcorr(bl,i) = fcorr(bl,i) + data(i)
                fcorr1(bl,i) = fcorr1(bl,i) + abs(data(i))
                fcorr2(bl,i) = fcorr2(bl,i) +
     *                    cmplx(real(data(i))**2,aimag(data(i))**2)
              endif
            enddo

          else
            if (abs(time-tprev).gt.TTOL) then
              if (nants.gt.0) call IntFlush(nants,rms,scalar,ra,lst,
     *          tprev,uvdist2,var,corr,corr1,corr2,xaxis,yaxis,
     *          xmin,xmax,ymin,ymax,npnt,time0,present,mbase,
     *          xdat,ydat,timedat,bldat,chdat,ndat,MAXDAT)
              nants = 0
              tprev = time
              call uvrdvrd(tno,'lst',lst,0d0)
              call uvrdvrd(tno,'ra',ra,0d0)
            endif

            n = 0
            do i = 1, nchan
              if (flags(i)) then
                n = n + 1
                npnt(bl) = npnt(bl) + 1
                corr(bl) = corr(bl) + data(i)
                corr1(bl) = corr1(bl) + abs(data(i))
                corr2(bl) = corr2(bl) +
     *                    cmplx(real(data(i))**2,aimag(data(i))**2)
              endif
            enddo
          endif

          if (n.gt.0) then
            call uvDatGtr('variance',temp)
            var(bl) = var(bl) + n*temp
            uvdist2(bl) = uvdist2(bl) +
     *       n * (preamble(1)*preamble(1)+preamble(2)*preamble(2))
            nants = max(nants,i1,i2)
          endif
        endif

        nchanprev = nchan
        call uvDatRd(preamble,data,flags,MAXCHAN,nchan)
      enddo

      if (nants.gt.0) then
        if (nofqaver) then
          call IntFlushF(nants,rms,scalar,ra,lst,time,uvdist2,var,
     *        fcorr,fcorr1,fcorr2,xaxis,yaxis,xmin,xmax,ymin,ymax,
     *        fnpnt,time0,present,mbase,xdat,ydat,timedat,bldat,
     *        chdat,ndat,MAXDAT,nchanprev)
        else
          call IntFlush(nants,rms,scalar,ra,lst,time,uvdist2,var,
     *        corr,corr1,corr2,xaxis,yaxis,xmin,xmax,ymin,ymax,npnt,
     *        time0,present,mbase,xdat,ydat,timedat,bldat,chdat,
     *        ndat,MAXDAT)
        endif
      endif

      end

c***********************************************************************

      subroutine IntFlush(nants,rms,scalar,ra,lst,time,uvdist2,var,
     *  corr,corr1,corr2,xaxis,yaxis,xmin,xmax,ymin,ymax,npnt,
     *  time0,present,MAXBASE,xdat,ydat,timedat,bldat,chdat,
     *  ndat,MAXDAT)

      integer MAXBASE,MAXDAT,nants,npnt(MAXBASE),bldat(MAXDAT),
     *  chdat(MAXDAT),ndat
      double precision ra,lst,time,time0,timedat(MAXDAT)
      real uvdist2(MAXBASE),var(MAXBASE),xdat(MAXDAT),ydat(MAXDAT)
      real xmin,xmax,ymin,ymax
      complex corr(MAXBASE),corr1(MAXBASE),corr2(MAXBASE)
      logical present(MAXBASE),rms,scalar
      character xaxis*(*),yaxis*(*)
c-----------------------------------------------------------------------
      integer i,j,k,ic
      real x,y

      external getVal
      real     getVal
c-----------------------------------------------------------------------

      ic=0
      k = 0
      do j = 1, nants
        do i = 1, j
          k = k + 1
          if (npnt(k).gt.0) then
            x = GetVal(xaxis,uvdist2(k),var(k),corr(k),
     *        corr1(k),corr2(k),npnt(k),lst,time,ic,ra,time0,
     *        rms,scalar)
            y = GetVal(yaxis,uvdist2(k),var(k),corr(k),
     *        corr1(k),corr2(k),npnt(k),lst,time,ic,ra,time0,
     *        rms,scalar)
            if (x.ge.xmin .and. x.le.xmax .and.
     *          y.ge.ymin .and. y.le.ymax) then
              ndat = ndat + 1
              if (ndat.gt.MAXDAT) call bug('f','Too many points.')
              xdat(ndat) = x
              ydat(ndat) = y
              bldat(ndat) = k
              timedat(ndat) = time

c             Use 0 to indicate whole spectrum.
              chdat(ndat) = 0
            endif

            present(k) = .true.
            npnt(k) = 0
            uvdist2(k) = 0
            var(k) = 0
            corr(k) = 0
            corr1(k) = 0
            corr2(k) = 0
          endif
        enddo
      enddo

      end

c***********************************************************************

      subroutine IntFlushF(nants,rms,scalar,ra,lst,time,uvdist2,var,
     *  corr,corr1,corr2,xaxis,yaxis,xmin,xmax,ymin,ymax,npnt,
     *  time0,present,MAXBASE,xdat,ydat,timedat,bldat,chdat,
     *  ndat,MAXDAT,nchan)

      integer MAXBASE,MAXDAT,nants,nchan,npnt(MAXBASE,nchan),
     *  bldat(MAXDAT),chdat(MAXDAT),ndat
      double precision ra,lst,time,time0,timedat(MAXDAT)
      real uvdist2(MAXBASE),var(MAXBASE),xdat(MAXDAT),ydat(MAXDAT)
      real xmin,xmax,ymin,ymax
      complex corr(MAXBASE,nchan),corr1(MAXBASE,nchan),
     *  corr2(MAXBASE,nchan)
      logical present(MAXBASE),rms,scalar
      character xaxis*(*),yaxis*(*)
c-----------------------------------------------------------------------
      integer i,j,k,ic
      real x,y

      external getVal
      real     getVal
c-----------------------------------------------------------------------
      k = 0
      do j = 1, nants
        do i = 1, j
          k = k + 1
          do ic = 1, nchan
            if (npnt(k,ic).gt.0) then
              x = GetVal(xaxis,uvdist2(k),var(k),corr(k,ic),
     *          corr1(k,ic),corr2(k,ic),npnt(k,ic),lst,time,ic,ra,
     *          time0,rms,scalar)
              y = GetVal(yaxis,uvdist2(k),var(k),corr(k,ic),
     *          corr1(k,ic),corr2(k,ic),npnt(k,ic),lst,time,ic,ra,
     *          time0,rms,scalar)
              if (x.ge.xmin .and. x.le.xmax .and. y.ge.ymin .and.
     *            y.le.ymax) then
                ndat = ndat + 1
                if (ndat.gt.MAXDAT) call bug('f','Too many points!')
                xdat(ndat) = x
                ydat(ndat) = y
                bldat(ndat) = k
                chdat(ndat) = ic
                timedat(ndat) = time
              endif
              npnt(k,ic) = 0
              uvdist2(k) = 0
              corr(k,ic) = 0
              corr1(k,ic) = 0
              corr2(k,ic) = 0
            endif
          enddo
          present(k) = .true.
          var(k) = 0
        enddo
      enddo

      end

c***********************************************************************

      real function GetVal(axis,uvdist2,var,corr,corr1,corr2,npnt,
     *  lst,time,chn,ra,time0,rms,scalar)

c-----------------------------------------------------------------------
      character axis*(*)
      real uvdist2,var
      complex corr,corr1,corr2
      double precision time,time0,lst,ra
      integer npnt,chn
      logical rms,scalar
c-----------------------------------------------------------------------
      include 'mirconst.h'
      complex data
      double precision dtemp
c-----------------------------------------------------------------------
      if (rms) then
        data = cmplx(sqrt(real(corr2)/npnt - real(corr/npnt)**2),
     *               sqrt(aimag(corr2)/npnt- aimag(corr/npnt)**2))
      else if (scalar) then
        data = corr1/npnt
      else
        data = corr/npnt
      endif

      if (axis.eq.'real') then
        GetVal = real(data)
      else if (axis.eq.'imaginary') then
        GetVal = aimag(data)
      else if (axis.eq.'amplitude') then
        GetVal = abs(data)
      else if (axis.eq.'phase') then
        GetVal = 180/pi * atan2(aimag(data),real(data))
      else if (axis.eq.'uvdistance') then
        GetVal = 0.001 * sqrt(uvdist2/npnt)
      else if (axis.eq.'rms') then
        GetVal = sqrt(var/npnt)
      else if (axis.eq.'time') then
        GetVal = 86400*(time - time0)
      else if (axis.eq.'lst') then
        GetVal = 86400*lst/(2*pi)
      else if (axis.eq.'hangle') then
        dtemp = lst - ra
        if (dtemp.gt.DPI) then
          dtemp = dtemp - 2*DPI
        else if (dtemp.lt.-DPI) then
          dtemp = dtemp + 2*DPI
        endif
        GetVal = 86400d0*dtemp/(2*DPI)
      else if (axis.eq.'channel') then
        GetVal = chn
      else
        call bug('f','I should never get here')
      endif

      end

c***********************************************************************

      subroutine GetAxis(xaxis,yaxis)

      character xaxis*(*), yaxis*(*)
c-----------------------------------------------------------------------
      integer NAX
      parameter (NAX=10)
      integer n
      character axes(NAX)*12
      data axes/'amplitude   ','phase       ',
     *          'real        ','imaginary   ',
     *          'time        ','uvdistance  ',
     *          'lst         ','hangle      ',
     *          'rms         ','channel     '/
c-----------------------------------------------------------------------
      call keymatch('axis',NAX,axes,1,xaxis,n)
      if (n.eq.0) xaxis = 'time'
      call keymatch('axis',NAX,axes,1,yaxis,n)
      if (n.eq.0) yaxis = 'amplitude'

      end

c***********************************************************************

      subroutine GetOpt(nobase,selgen,noapply,rms,scalar,nofqaver,
     *  uvflags)

      logical   nobase, selgen, noapply, rms, scalar, nofqaver
      character uvflags*(*)
c-----------------------------------------------------------------------
c  Get extra processing options.
c-----------------------------------------------------------------------
      integer NOPTS
      parameter (NOPTS=9)
      logical present(NOPTS)
      character opts(NOPTS)*8
      data opts/'nobase  ','nocal   ','nopass  ','nopol   ',
     *          'selgen  ','noapply ','rms     ','scalar  ',
     *          'nofqaver'/
c-----------------------------------------------------------------------
      call options('options',opts,present,NOPTS)

      nobase = present(1)
      selgen = present(5)
      noapply= present(6)
      rms    = present(7)
      scalar = present(8)
      nofqaver=present(9)
      if (scalar .and. rms)
     *  call bug('f','Options scalar and rms cannot be used together')
      uvflags = 'sdlwb'
      if (.not.present(2)) uvflags(6:6) = 'c'
      if (.not.present(3)) uvflags(7:7) = 'f'
      if (.not.present(4)) uvflags(8:8) = 'e'

      end

c***********************************************************************

      subroutine flagchk(tno)

      integer tno
c-----------------------------------------------------------------------
c  Check that the user's linetype is not going to cause the flagging
c  routine to vomit when the flagging is applied.
c-----------------------------------------------------------------------
      integer CHANNEL,WIDE
      parameter (CHANNEL=1,WIDE=2)
      double precision line(6)
c-----------------------------------------------------------------------
      call uvinfo(tno,'line',line)
      if (nint(line(1)).ne.CHANNEL .and. nint(line(1)).ne.WIDE)
     *  call bug('f','Can only flag "channel" or "wide" linetypes')
      if (nint(line(4)).ne.1)
     *  call bug('f','Cannot flag when the linetype width is not 1')

      end

c***********************************************************************

      subroutine getrng(keyw, axis, rmin, rmax)

      character keyw*(*), axis*(*)
      real      rmin, rmax
c-----------------------------------------------------------------------
c  Get the axis ranges given by the user
c
c  Input
c    keyw     Keyword to get from user
c    axis     Axis type
c  Output
c    rmin,max Range in appropriate units
c-----------------------------------------------------------------------
      logical ok
      integer il, len1, nt, s
      real trange(8)
      character cval*64
c-----------------------------------------------------------------------
      il = len1(keyw)
      if (axis.eq.'time') then
        call mkeyr(keyw(1:il), trange, 8, nt)
        if (nt.gt.0) then
          if (nt.ne.8) then
            call bug('f',
     *        'You must specify 8 numbers for the time range')
          else
c           Convert to seconds.
            rmin = 24.0*3600.0*trange(1) + 3600.0*trange(2) +
     *                    60.0*trange(3) + trange(4)
            rmax = 24.0*3600.0*trange(5) + 3600.0*trange(6) +
     *                    60.0*trange(7) + trange(8)
          endif
        else
          rmin = -1e32
          rmax =  1e32
        endif

      else if (axis.eq.'hangle') then
        call mkeyr(keyw(1:il), trange, 6, nt)
        if (nt.gt.0) then
          if (nt.ne.6) then
            call bug('f',
     *        'You must specify 6 numbers for the hangle range')
          else
c           Convert to seconds.
            s = 1
            if (trange(1).lt.0.0) s = -1
            rmin = 3600.0*abs(trange(1)) + 60.0*trange(2) + trange(3)
            rmin = s * rmin

            s = 1
            if (trange(4).lt.0.0) s = -1
            rmax = 3600.0*abs(trange(4)) + 60.0*trange(5) + trange(6)
            rmax = s * rmax
          endif
        else
          rmin = -1e32
          rmax =  1e32
        endif

      else
        call keya(keyw(:il), cval, 'min')
        if (cval.eq.'min') then
          rmin = -1e32
        else
          call atorf(cval, rmin, ok)
          if (.not.ok) then
            cval = 'Conversion error decoding parameter ' // keyw(:il)
            call bug('f', cval)
          endif
        endif

        call keya(keyw(:il), cval, 'max')
        if (cval.eq.'max') then
          rmax = 1e32
        else
          call atorf(cval, rmax, ok)
          if (.not.ok) then
            cval = 'Conversion error decoding parameter ' // keyw(:il)
            call bug('f',cval)
          endif
        endif

c       Because atorf actually uses atodf and conversion between
c       double and real may introduce rounding errors.
        if (-1.000001e32.lt.rmin .and. rmin.lt.-0.999999e32) then
          rmin = -1e32
        endif

        if (0.999999e32.lt.rmax .and. rmax.lt.1.000001e32) then
          rmax =  1e32
        endif
      endif

      end
