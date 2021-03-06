c************************************************************************
	program sdip
c
	implicit none
c
c= sdip - Compute receiver, sky and diode temperatures.
c& rjs
c: analysis
c	SDIP is a Miriad task to reduce sky dip and hot load measurement
c	data from the ATCA. The input is the so-called DSD log.
c
c	The absolute calibration is either determined by a hot load, or
c	by specifying the zenith sky brightness.
c+
c@ dsdlog
c	DSD log file.
c@ tsky
c	A value, giving the zenith sky brightness temperature (including the
c	contribution of the CMB), in Kelvin. Zero, one of four values can be
c	given. If no values are given, then "thot" and "htime" (see below)
c	must be given.
c@ thot
c	One or two values, giving the hot load temperature and the
c	antenna temperature, in Kelvin. The default is 300,0.
c@ htime
c	Time range of hot loads.
c@ stime
c	Time range of sky dips.
c@ freq
c	Observing frequency, in GHz. This value is used to determine
c	sky temperatures based on meteorology data. A number of values
c	can be given.
c@ ant
c	The antenna to process. Use either 3 or 4. The default is 3.
c@ options
c	Extra processing options.
c	  onoff  Use this option if the antenna has a noise diode.
c	  thin   Use the optically thin atmosphere approximation in
c	         the solution. The default is to not assume this (which
c	         mades is a non-linear problem). Note the plots never use
c	         the optically thin approximation.
c@ device
c	PGPLOT device to plot the result. The data are plotted in black, and
c	the "fit" in red. The default is to not plot the results.
c--
c  History:
c    27apr01 rjs   Original version.
c    08aug01 rjs   Added model sky calculations and met info.
c    17nov01 rjs   Handle "dsdn" format.
c    02dec01 rjs   Get rid of the optically thin approximation.
c------------------------------------------------------------------------
	include 'mirconst.h'
	character version*(*)
	parameter(version='Sdip: version 02-Dec-01')
	integer MAXT,MAXFREQ
	parameter(MAXT=10,MAXFREQ=3)
c
	integer ITMAX
	real EPS1,EPS2
	parameter(ITMAX=100,EPS1=0.05,EPS2=0.005)
c
c  Commons.
c
	integer MAXEQ
	real Tcmb
	parameter(MAXEQ=2000,Tcmb=2.72)
	real t0,tant,eqs(3,MAXEQ,4),bs(MAXEQ,4)
	integer neq,j0
	common/sdipcom/T0,Tant,bs,eqs,neq,j0
c
	logical onoff,ht,st,dotsky,thin
	character dsdlog*64,device*64,line*64,chan*5
	double precision htime(2,MAXT),stime(2,MAXT),time,time0
	integer nht,nst,nhts,nsts,i,j,ifail,ant,offsets(4,6)
	real dsdons(4),dsdoffs(4)
	real times(MAXEQ),y1(MAXEQ),y2(MAXEQ)
	real Thot,Trec,Tcal,Tsky,Tmod,Tatm
	real press,humid,temp0,press0,humid0,fac,tau
	real xlo,xhi,ylo,yhi,el,dsdon,dsdoff
	integer nmet,nfreq,ntskys
	real ymin,ymax,freqs(MAXFREQ)
	real eqs2(2,MAXEQ,4),bs2(MAXEQ,4)
	real temps(3),tskys(4)
c
c  Scratch arrays.
c
	real f(MAXEQ),fp(MAXEQ),dx(3),dfdx(3*MAXEQ),aa(9)
	real rtmp(9)
	integer itmp(3)
c
c  Externals.
c
	character itoaf*4,stcat*64,streal*64
	integer pgBeg,tinNext
	logical gotData
	external func,deriv
c
	data offsets/
     *       0,    0,    0,    0,
     *       0,    0,    0,    0,
     *	  -286,  155, -259,  437,
     *	   -38,  -35,  608, -324,
     *       0,    0,    0,    0,
     *       0,    0,    0,    0/
c
	call output(version)
	call keyini
	call keyi('ant',ant,3)
	call keya('dsdlog',dsdlog,' ')
	call mkeyr('tsky',tskys,4,ntskys)
	dotsky = ntskys.gt.0
	if(ntskys.gt.0)then
	  do i=ntskys+1,4
	    tskys(i) = tskys(i-1)
	  enddo
	endif
c
	call keyr('thot',thot,300.0)
	call keyr('thot',tant,0.0)
	call keya('device',device,' ')
	call mkeyt('htime',htime,2*MAXT,nht,'atime')
	if(2*(nht/2).ne.nht)call bug('f','Odd no. of hot load times')
	if(nht.lt.2.and..not.dotsky)
     *		call bug('f','Hot load times must be given')
	call mkeyt('stime',stime,2*MAXT,nst,'atime')
	if(2*(nst/2).ne.nst)call bug('f','Odd no. of sky dip times')
	if(nst.lt.2)call bug('f','Sky dip times must be given')
	call mkeyr('freq',freqs,MAXFREQ,nfreq)
	call getopt(onoff,thin)
	call keyfin
c
	nmet = 0
	Tatm = 0
	press = 0
	humid = 0
c
	neq = 0
	nhts = 0
	nsts = 0
	call tinOpen(dsdlog,'n')
	dowhile(tinNext().gt.0)
	  if(gotData(ant,time,el,dsdons,dsdoffs,
     *				temp0,press0,humid0))then
	    ht = .false.
	    st = .false.
	    do j=1,nht
	      ht = ht.or.(htime(1,j).le.time.and.time.le.htime(2,j))
	    enddo
	    ht = ht.and..not.dotsky
	    do j=1,nst
	      st = st.or.(stime(1,j).le.time.and.time.le.stime(2,j))
	    enddo
	    if(ht.or.st)then
	      neq = neq + 1
	      if(neq.gt.MAXEQ)call bug('f','Too many equations')
	      if(neq.eq.1)time0 = int(time-0.5) + 0.5
	      times(neq) = 86400*(time-time0)
	      if(ht)nhts = nhts + 1
	      if(st)nsts = nsts + 1
c
	      do j=1,4
		dsdon = dsdons(j)
		dsdoff = dsdoffs(j)
		if(dsdon.gt.61000.or.dsdoff.gt.61000.or.
     *		  (dsdon.le.dsdoff.and.onoff))then
		  dsdon = 0
		else if(onoff)then
		  dsdon = (dsdoff-offsets(j,ant))/(dsdon-dsdoff)
		else
		  dsdon = 0.5*(dsdon+dsdoff)-offsets(j,ant)
		endif
		if(dsdon.eq.0)then
		  eqs(1,neq,j) = 0
		  eqs(2,neq,j) = 0
		  eqs(3,neq,j) = 0
		  bs(neq,j) = 0
		else if(ht)then
		  eqs(1,neq,j) = 1
		  eqs(2,neq,j) = -dsdon
		  eqs(3,neq,j) = 0
		  bs(neq,j) = -thot
		else
		  eqs(1,neq,j) = 1
		  eqs(2,neq,j) = -dsdon
		  eqs(3,neq,j) = 1./sin(el)
		  bs(neq,j) = - tant - tcmb
		endif
	      enddo
	      Tatm = Tatm + temp0
	      press = press + press0
	      humid = humid + humid0
	      nmet = nmet + 1
	    endif
	  endif
	enddo
	call tinClose
	call output('Hot load samples: '//itoaf(nhts))
	call output('Sky dip samples:  '//itoaf(nsts))
	if((nhts.eq.0..and..not.dotsky).or.nsts.eq.0)
     *				call bug('f','Ill conditioned')
c
c  Check that we have some met data.
c
	if(nmet.le.0)call bug('f','No met data!')
	Tatm = Tatm/nmet
	press = press/nmet
	humid = humid/nmet
	call output('-------------------------------------------')
	call output('Mean meteorological parameters:')
	line = stcat('  T =',
     *	       stcat(' '//streal(Tatm-273.15,'(f10.1)'),' Celsius'))
	call output(line)
	line = stcat('  P =',
     *	       stcat(' '//streal(press/100/0.975,'(f10.1)'),' hPa'))
	call output(line)
	line = stcat('  h =',
     *		stcat(' '//streal(humid*100,'(f10.1)'),'%'))
	call output(line)
	call output('-------------------------------------------')
c
	T0 = Tatm
	do i=1,nfreq
	  if(i.eq.1)call output('Model sky parameters:')
	  if(freqs(i).lt.0.1.or.freqs(i).gt.200)
     *	    call bug('f','Invalid frequency')
	  call opacGet(1,freqs(i)*1e9,0.5*PI,Tatm,press,humid,
     *							fac,Tmod)
	  line = stcat('Observing frequency:   nu   =',
     *		stcat(' '//streal(freqs(i),'(f10.4)'),' GHz'))
	  call output(line)
	  line = stcat('Model sky brightness:  Tsky =',
     *		stcat(' '//streal(Tmod,'(f10.1)'),' Kelvin'))
	  call output(line)
	  line = stcat('Model sky opacity:     tau  =',
     *		stcat(' '//streal(-log(fac),'(f10.3)'),' nepers'))
	  call output(line)
	  call output(' ')
	  if(i.eq.1)then
	    t0 = (tmod-tcmb*fac)/(1-fac)
	  endif
	enddo
	line = stcat('Effective temperature: T0   =',
     *	       stcat(' '//streal(T0,'(f10.1)'),' Kelvin'))
	call output(line)
	call output('-------------------------------------------')
c
c  Fix up the equations.
c
	if(dotsky)then
	  do j=1,4
	    fac = (T0-Tskys(j))/(T0-Tcmb)
	    tau = -log(fac)
	    do i=1,neq
	      eqs2(1,i,j) = eqs(1,i,j)
	      eqs2(2,i,j) = eqs(2,i,j)
	      if(eqs(3,i,j).eq.0)then
		bs2(i,j) = 0
	      else
	        bs2(i,j) = -tant - tcmb - 
     *				(t0-tcmb)*(1-exp(-tau*eqs(3,i,j)))
	      endif
	    enddo
	  enddo
	endif
c
c  Now solve for all the parameters.
c
	if(device.ne.' ')then
	  if(pgbeg(0,device,2,2).ne.1)
     *	    call bug('f','Failed to open plot device')
	  call pgscf(2)
	  call pgsch(2.0)
	endif
	call output('Channel Trec    Tcal   Tsky')
	call output('------- ----    ----   ----')
c
	do j=1,4
	  j0 = j
	  chan = char(ichar('A')+j-1)//char(ichar('0')+ant)
	  if(dotsky)then
	    call llsqu(bs2(1,j),eqs2(1,1,j),2,neq,temps,
     *						ifail,rtmp,itmp)
	    temps(3) = tskys(j)
	  else
	    call llsqu(bs(1,j),eqs(1,1,j),3,neq,temps,
     *						ifail,rtmp,itmp)
	    temps(3) = temps(3) + tcmb
	    if(.not.thin.and.ifail.eq.0)then
	      Tsky = temps(3)
	      fac = (T0-Tsky)/(T0-Tcmb)
	      temps(3) = -log(fac)
	      call nllsqu(3,neq,temps,temps,ITMAX,EPS1,EPS2,.true.,
     *		ifail,func,deriv,f,fp,dx,dfdx,aa)
	      fac = exp(-temps(3))
	      temps(3) = (T0-Tcmb)*(1-fac) +Tcmb
	    endif
	  endif
	  if(ifail.ne.0)then
	    call bug('w','Ifail = '//itoaf(ifail))
	    trec = 0
	    tsky = 0
	    tcal = 0
	  else
	    trec = temps(1)
	    tcal = temps(2)
	    tsky = temps(3)
	  endif
	  write(line,'(a,f7.1,f8.2,f7.1)')chan,Trec,Tcal,Tsky
	  call output(line)
c
	  if(device.ne.' ')then
	    fac = (T0-Tsky)/(T0-Tcmb)
	    tau = -log(fac)
	    do i=1,neq
	      y1(i) = -eqs(2,i,j)*Tcal
	      if(eqs(1,i,j).eq.0)then
		y2(i) = 0
	      else if(eqs(3,i,j).gt.0)then
	        y2(i) = trec+tant+T0-(T0-tcmb)*exp(-tau*eqs(3,i,j))
	      else
	        y2(i) = (trec+thot)
	      endif
	    enddo
	    xlo = times(1)
	    xhi = times(neq)
	    xlo = xlo - 0.1*(xhi-xlo)
	    xhi = xhi + 0.1*(xhi-xlo)
	    ymin = min(y1(1),y2(1))
	    ymax = max(y1(1),y2(1))
	    do i=1,neq
	      ymin = min(ymin,y1(i),y2(i))
	      ymax = max(ymax,y1(i),y2(i))
	    enddo
	    call pgrnge(ymin,ymax,ylo,yhi)
	    yhi = max(600.0,yhi+0.1*(yhi-ylo))
	    ylo = min(0.0,ylo-0.1*(yhi-ylo))
c
	    call pgpage
	    call pgvstd
	    call pgswin(xlo,xhi,ylo,yhi)
	    call pgtbox('BCSTNHZO',0.,0,'BCNST',0.,0)
	    call pglab('Time (UTC)','DSD Value (K)',
     *				'Channel '//chan)
	    call pgpt(neq,times,y1,1)
	    call pgsci(2)
	    call pgpt(neq,times,y2,1)
	    call pgsci(1)
	  endif
	enddo
	if(device.ne.' ')call pgend
c
	end
c************************************************************************
	subroutine getopt(onoff,thin)
c
	implicit none
	logical onoff,thin
c------------------------------------------------------------------------
	integer NOPTS
	parameter(NOPTS=2)
	character opts(NOPTS)*8
	logical present(NOPTS)
	data opts/'onoff   ','thin    '/
	call options('options',opts,present,NOPTS)
	onoff = present(1)
	thin  = present(2)
	end
c************************************************************************
	logical function gotData(ant,time,el,dsdons,dsdoffs,
     *						temp,press,humid)
c
	implicit none
	integer ant
	double precision time
	real el,dsdons(4),dsdoffs(4),temp,press,humid
c------------------------------------------------------------------------
	include 'mirconst.h'
	character ftype*16,state*2
	integer i,iant,nants
	logical ok
c
c  Externals.
c
	integer tinNext
c
c  Determine the file type.
c
	call tinGeta(ftype,' ')
c
c  The old format.
c
	if(ftype.eq.'dsd34')then
	  call tinGett(time,0.d0,'atime')
	  call tinGeta(state,' ')
	  ok = state.eq.'i'.or.state.eq.'t'.and.
     *				(ant.eq.3.or.ant.eq.4)
	  if(ok)then
	    call tinSkip(1)
	    call tinGetr(el,0.)
	    if(ant.eq.4)call tinSkip(8)
	    do i=1,4
	      call tinGetr(dsdons(i),0.)
	      call tinGetr(dsdoffs(i),0.)
	    enddo
	    if(ant.eq.3)call tinSkip(8)
	    call tinSkip(2)
	    call tinGetr(temp,0.)
	    call tinGetr(press,0.)
	    call tinGetr(humid,0.)
	  endif
c
c  The newer format.
c
	else if(ftype.eq.'dsdn')then
	  call tinGeti(nants,0)
	  call tinGett(time,0.d0,'atime')
	  call tinGeta(state,' ')
	  ok = state.eq.'i'.or.state.eq.'t'
	  if(ok)then
	    call tinSkip(1)
	    call tinGetr(el,0.)
	    call tinSkip(2)
	    call tinGetr(temp,0.)
	    call tinGetr(press,0.)
	    call tinGetr(humid,0.)
	    ok = .false.
	    dowhile(nants.gt.0.and..not.ok)
	      i = tinNext()
	      nants = nants - 1
	      call tinGeti(iant,0)
	      if(iant.eq.ant)then
		ok = .true.
		do i=1,4
		  call tinGetr(dsdons(i),0.)
		  call tinGetr(dsdoffs(i),0.)
		enddo
	      endif
	    enddo
	  endif
	  dowhile(nants.gt.0)
	    i = tinNext()
	    nants = nants - 1
	  enddo
	else
	  call bug('f','Unrecognised file type')
	endif
	if(ok)then
	  temp = temp + 273.15
	  press = 0.975*100*press
	  humid = 0.01*humid
	  el = PI/180.0 * el
	endif
c
	gotData = ok
c
	end
c************************************************************************
	subroutine func(temps,f,n,m)
c
	implicit none
	integer n,m
	real temps(n),f(m)
c------------------------------------------------------------------------
	real trec,tcal,tau
	integer i
c
c  Commons.
c
	integer MAXEQ
	real Tcmb
	parameter(MAXEQ=2000,Tcmb=2.72)
	real t0,tant,eqs(3,MAXEQ,4),bs(MAXEQ,4)
	integer neq,j0
	common/sdipcom/T0,Tant,bs,eqs,neq,j0
c
	if(neq.ne.m)call bug('f','Inconsistent number of equations')
	if(n.ne.3)call bug('f','Inconsistent number of unknowns')
c
	trec = temps(1)
	tcal = temps(2)
	tau  = temps(3)
c
	do i=1,neq
	  f(i) = trec + tcal*eqs(2,i,j0) - bs(i,j0)
	  if(eqs(3,i,j0).ne.0)
     *	    f(i) = f(i) + (t0-tcmb)*(1-exp(-tau*eqs(3,i,j0)))
	enddo
c
	end
c************************************************************************
	subroutine deriv(temps,dfdx,n,m)
c
	implicit none
	integer n,m
	real temps(n),dfdx(n,m)
c
c------------------------------------------------------------------------
	real tau
	integer i
c
c  Commons.
c
	integer MAXEQ
	real Tcmb
	parameter(MAXEQ=2000,Tcmb=2.72)
	real t0,tant,eqs(3,MAXEQ,4),bs(MAXEQ,4)
	integer neq,j0
	common/sdipcom/T0,Tant,bs,eqs,neq,j0
c
	if(neq.ne.m)call bug('f','Inconsistent number of equations')
	if(n.ne.3)call bug('f','Inconsistent number of unknowns')
c
	tau  = temps(3)
c
	do i=1,neq
	  dfdx(1,i) = 1
	  dfdx(2,i) = eqs(2,i,j0)
	  if(eqs(3,i,j0).eq.0)then
	    dfdx(3,i) = 0
	  else
	    dfdx(3,i) = (t0-tcmb)*eqs(3,i,j0)*exp(-tau*eqs(3,i,j0))
	  endif
	enddo
c
	end
