c************************************************************************
c
c  A set of routines to do extended precision integer arithmetic
c  in standard FORTRAN. These routines approximately double the
c  range of integers.
c
c  History:
c    27mar09 rjs  Original version
c    21apr09 rjs  Complete rewrite to use a simpler, more limited,
c		  more efficient, representation.
c    10may09 rjs  Added simple division
c    17jul09 rjs  Added sign routine and corrected a documentation error.
c
c  The set of routines manipulate these triples.
c
c  subroutine mpCvtim(v,int)             v = int
c  integer function mpCvtmi(v)           mpCvtmi = v
c  double precision function mpCvtmd(v)  mpCvtmd = v
c  subroutine mpFmt(string,v)            
c  subroutine mpAddmi(v,int)             v = v + int
c  subroutine mpAddmm(v1,v2)             v1 = v1 + v2
c  subroutine mpSubmi(v,int)             v  = v - int
c  subroutine mpSubmm(v1,v2)             v1 = v1 - v2
c  subroutine mpMulmi(v,int)             v  = v * int
c  subroutine mpMulmm(v1,v2)             v1 = v1 * v2
c  subroutine mpDivmi(v,int,rem)         v  = v/int, rem = v%int
c  integer function mpCmp(v1,v2)         v1 cmp v2
c  integer function mpSign(v)		 sign(v) [-1,0,+1]
c  subroutine mpSet(v1,v2)               v1 = v2
c  subroutine mpNeg(v)                   v = -v
c  subroutine mpAbs(v)                   v = abs(v)
c
c  A multi-precision integer is stored as a triple
c
c    value = v(1) + v(2)*v(3)
c
c  The "standard" representation is to have v(3) be a number approximately
c  a half the size of the maximum allowed integer and also a square.
c
c
c************************************************************************
c* mpCvtim -- Convert integer to multi-precision format.
c& rjs
c+
	subroutine mpCvtim(v,k)
c
	implicit none
	integer v(3),k
c
c  Convert an integer to multi-precision format.
c
c  Input:
c    k		The integer to be converted.
c  Output:
c    v		The multi-precision representation of the integer.
c--
c------------------------------------------------------------------------
	include 'mp.h'
	logical first
	save first
c
	data first/.true./
c
	if(first)call mpInit
	first = .false.
	call mpStd2(k,0,mpBase2,v(1),v(2))
	v(3) = mpBase2
	end
c************************************************************************
	subroutine mpStd2(v1,v2,v3,d1,d2)
c
	implicit none
	integer v1,v2,v3,d1,d2
c
c  Convert an arbitrary triple, v1+v2*v3, into a "standard triple",
c  d1+d2*mpBase2, where d1 and d2 have the same sign.
c------------------------------------------------------------------------
	include 'mp.h'
	integer d1a,a(3),b(3),d(6),carry,i,j,k
c
	d1 = v1
	d1a = 0
	if(v3.eq.mpBase2)then
	  d2 = v2
	else if(v2.ne.0.and.v3.ne.0)then
	  a(1) = v2
	  a(2) = a(1)/mpBase
	  a(1) = a(1) - a(2)*mpBase
	  a(3) = a(2)/mpBase
	  a(2) = a(2) - a(3)*mpBase
	  b(1) = v3
	  b(2) = b(1)/mpBase
	  b(1) = b(1) - b(2)*mpBase
	  b(3) = b(2)/mpBase
	  b(2) = b(2) - b(3)*mpBase
	  do k=1,6
	    d(k) = 0
	  enddo
	  do j=1,3
	    carry = 0
	    do i=1,3
	      k = i + j - 1
	      d(k) = d(k) + a(i)*b(j) + carry
	      carry = d(k)/mpBase
	      d(k) = d(k) - carry*mpBase
	    enddo
	    d(j+3) = d(j+3) + carry
	  enddo
	  d1a = d(2)*mpBase + d(1)
	  d2  = d(4)*mpBase + d(3)
	  if(abs(d(5))+abs(d(6)).ne.0)
     *		call bug('f','Integer overflow in mpStd2')
	else
	  d2 = 0
	endif
c
c  Add in the different components with care to avoid overflow.
c
	carry = d1/mpBase2
	d1 = d1 - carry*mpBase2
	d2 = d2 + carry
c
	d1 = d1 + d1a
	carry = d1/mpBase2
	d1 = d1 - carry*mpBase2
	d2 = d2 + carry
c
c  Make d1 and d2 have the same signs.
c
	if(d1.lt.0.and.d2.gt.0)then
	  d1 = d1 + mpBase2
	  d2 = d2 - 1
	else if(d1.gt.0.and.d2.lt.0)then
	  d1 = d1 - mpBase2
	  d2 = d2 + 1
	endif
c
	end
c************************************************************************
c* mpCvtmd -- Convert multi-precision integer to double precision.
c& rjs
c+
	double precision function mpCvtmd(v)
c
	implicit none
	integer v(3)
c
c  Convert multi-precision integer to double precision.
c
c  Input:
c    v		The multi-precision integer.
c  Output:
c    mpCvtmd	The double precision representation.
c--
c------------------------------------------------------------------------
	double precision temp
c
c  Order the operations to avoid the possibility of overflow.
c
	temp = v(3)
	temp = temp*v(2) + v(1)
	mpCvtmd = temp
c
	end
c************************************************************************
c* mpCvtmi -- Convert multi-precision integer to a standard integer.
c& rjs
c+
	integer function mpCvtmi(v)
c
	implicit none
	integer v(3)
c
c  Convert multi-precision integer to a standard integer.
c  If the number exceeds the maximum allowed integer, the program
c  is aborted.
c
c  Input:
c    v		The multi-precision integer.
c  Output:
c    mpCvtmi	The standard integer representation.
c--
c------------------------------------------------------------------------
	include 'mp.h'
	integer t(3)
c
c  Externals.
c
	integer mpCmp
c
c  Check for overflow first.
c
	t(1) = v(1)
	t(2) = v(2)
	t(3) = v(3)
	call mpAbs(t)
	if(mpCmp(t,mpMAXInt).gt.0)
     *		call bug('f','Integer overflow in mpCvtmi')
	mpCvtmi = v(1) + v(2)*v(3)
	end
c************************************************************************
c* mpSet -- Copy one multi-precision integer to another.
c& rjs
c+
	subroutine mpSet(v1,v2)
c
	implicit none
	integer v1(3),v2(3)
c
c  Copy one multi-precision integer to another.
c
c  Input:
c    v1		The input multi-precision integer.
c  Output:
c    v2		The copy of the input multi-precision integer.
c--
c------------------------------------------------------------------------
	include 'mp.h'
	call mpStd2(v2(1),v2(2),v2(3),v1(1),v1(2))
	v1(3) = mpBase2
	end
c************************************************************************
c* mpAddmi -- Add an integer to a multi-precision integer.
c& rjs
c+
	subroutine mpAddmi(v,d)
c
	implicit none
	integer v(3),d
c
c  Add an integer to a multi-precision integer.
c    v = v + d
c
c  Input:
c    d		The integer to be added.
c  Input/Output:
c    v		The multi-precision integer to be added to.
c--
c------------------------------------------------------------------------
	include 'mp.h'
	integer t(3),carry

c
c  If the triple is in standard format and the input integer is
c  small, then do the operation quickly.
c
	if(v(3).eq.mpBase2.and.abs(d).lt.mpBase2
     *			  .and.abs(v(1)).lt.mpBase2)then
	  v(1) = v(1) + d
	  carry = v(1)/mpBase2
	  v(1) = v(1) - carry*mpBase2
	  v(2) = v(2) + carry
	else
	  call mpCvtim(t,d)
	  call mpAddmm(v,t)
	endif
c
	end
c************************************************************************
c* mpAddmm -- Add two multi-precision integers.
c& rjs
c+
	subroutine mpAddmm(v1,v2)
c
	implicit none
	integer v1(3),v2(3)
c
c--
c------------------------------------------------------------------------
	include 'mp.h'
	integer d11,d12,d21,d22,carry
c
c  If the triple is in the standard form, then do the operation
c  quickly. Otherwise convert to standard and do it.
c
	if(abs(v1(1)).lt.mpBase2.and.abs(v2(1)).lt.mpBase2.and.
     *	       v1(3) .eq.mpBase2.and.    v2(3) .eq.mpBase2)then
	  v1(1) = v1(1) + v2(1)
	  v1(2) = v1(2) + v2(2)
	else
	  call mpStd2(v1(1),v1(2),v1(3),d11,d12)
	  call mpStd2(v2(1),v2(2),v2(3),d21,d22)
	  v1(1) = d11 + d21
	  v1(2) = d12 + d22
	  v1(3) = mpBase2
	endif
	carry = v1(1)/mpBase2
	v1(1) = v1(1) - carry*mpBase2
	v1(2) = v1(2) + carry
	end
c************************************************************************
c* mpSubmi -- Subtract an integer from a multi-precision integer.
c& rjs
c+
	subroutine mpSubmi(v,d)
c
	implicit none
	integer v(3),d
c
c  Subtract two numbers.
c--
c------------------------------------------------------------------------
	include 'mp.h'
	integer t(3),carry

c
c  If the triple is in standard format and the integer is small, then do
c  the operation quickly.
c
	if(v(3).eq.mpBase2.and.abs(d).lt.mpBase2
     *			  .and.abs(v(1)).lt.mpBase2)then
	  v(1) = v(1) - d
	  carry = v(1)/mpBase2
	  v(1) = v(1) - carry*mpBase2
	  v(2) = v(2) + carry
	else
	  call mpCvtim(t,d)
	  call mpSubmm(v,t)
	endif
c
	end
c************************************************************************
c* mpSubmm -- Subtract two multi-precision integers.
c& rjs
c+
	subroutine mpSubmm(v1,v2)
c
	implicit none
	integer v1(3),v2(3)
c
c  Subtract two multi-precision integers.
c
c  Input:
c    v2		Number to be subtracted.
c  Input/Output:
c    v1		On output, v1 = v1 - v2.
c--
c------------------------------------------------------------------------
	include 'mp.h'
	integer carry,d11,d12,d21,d22
c
	if(abs(v1(1)).lt.mpBase2.and.abs(v2(1)).lt.mpBase2.and.
     *	       v1(3) .eq.mpBase2.and.    v2(3) .eq.mpBase2)then
	  v1(1) = v1(1) - v2(1)
	  v1(2) = v1(2) - v2(2)
	else
	  call mpStd2(v1(1),v1(2),v1(3),d11,d12)
	  call mpStd2(v2(1),v2(2),v2(3),d21,d22)
	  v1(1) = d11 - d21
	  v1(2) = d12 - d22
	  v1(3) = mpBase2
	endif
	carry = v1(1)/mpBase2
	v1(1) = v1(1) - carry*mpBase2
	v1(2) = v1(2) + carry
	end
c************************************************************************
c* mpMulmi -- Multiply a multi-precision integer by an integer.
c& rjs
c+
	subroutine mpMulmi(v,d)
c
	implicit none
	integer v(3),d
c
c  Multiply a multi-precision integer by an integer.
c
c  Input:
c    d		Integer to multiply by.
c  Input/Output:
c    v		On output, v = v * d.
c--
c------------------------------------------------------------------------
	include 'mp.h'
	integer q(4),t(3),carry,i
c
	if(abs(d).lt.mpBase)then
	  call mpStd4(v,q)
	  carry = 0
	  do i=1,4
	    q(i) = q(i)*d + carry
	    carry = q(i)/mpBase
	    q(i) = q(i) - carry*mpBase
	  enddo
	  if(carry.ne.0)call bug('f','Integer overflow in mpMulmi')
	  v(1) = q(2)*mpBase + q(1)
	  v(2) = q(4)*mpBase + q(3)
	  v(3) = mpBase2
	else
	  call mpCvtim(t,d)
	  call mpMulmm(v,t)
	endif
	end
c************************************************************************
c* mpMulmm -- Multiply two multi-precision integers.
c& rjs
c+
	subroutine mpMulmm(v1,v2)
c
	implicit none
	integer v1(3),v2(3)
c
c  Multiply two multi-precision integers.
c
c  Input:
c    v2		The multi-precision integer to multiply by.
c  Input/Output:
c    v1		On output, v1 = v1 * v2.
c--
c------------------------------------------------------------------------
	include 'mp.h'
	integer q1(4),q2(4),d(8),i,j,k,carry
c
	call mpStd4(v1,q1)
	call mpStd4(v2,q2)
c
	do k=1,8
	  d(k) = 0
	enddo
c
	do j=1,4
	  carry = 0
	  do i=1,4
	    k = i + j - 1
	    d(k) = d(k) + q1(i)*q2(j) + carry
	    carry = d(k)/mpBase
	    d(k) = d(k) - carry*mpBase
	  enddo
	  d(j+4) = d(j+4) + carry
	enddo
c
	if(d(8).ne.0.or.d(7).ne.0.or.d(6).ne.0.or.d(5).ne.0)
     *	  call bug('f','Integer overflow in mpMulmm')
	v1(1) = d(2)*mpBase + d(1)
	v1(2) = d(4)*mpBase + d(3)
	v1(3) = mpBase2
c
	end
c************************************************************************
c* mpDivmi -- Divide a multi-precision integer by an integer.
c& rjs
c+
	subroutine mpDivmi(v,k,rem)
c
	implicit none
	integer v(3),k,rem
c
c  Divide a multi-precision integer by an integer, returning the
c  quotient and remainder.
c
c  NOTE: The current version is a very simple one, which limits the
c  maximum size of the divisor to mpBase.
c
c  Input:
c    k		The divisor
c  Input/Output:
c    v		On output v = v / k.
c  Output:
c    rem	The remainder. Specifically, rem = v % k.
c    
c--
c------------------------------------------------------------------------
	include 'mp.h'
	integer q(4),t,carry
c
	if(abs(k).gt.mpBase)
     *	  call bug('f','Algorithmic failure in mpDivmi')
	call mpStd4(v,q)
	carry = q(4)*mpBase + q(3)
	v(2) = carry/k
	carry = (carry - v(2)*k)*mpBase + q(2)
	v(1) = carry/k
	carry = (carry - v(1)*k)*mpBase + q(1)
	t = carry/k
	v(1) = v(1)*mpBase + t
	rem = (carry - t*k)
	end
c************************************************************************
	subroutine mpStd4(v,q)
c
	implicit none
	integer v(3),q(4)
c
c  Convert a triple into a standard "quad" representation.
c
c  v(1)+v(2)*v(3) = q(1) + q(2)*mpBase + q(3)*mpBase**2 + q(4)*mpBase**3
c
c------------------------------------------------------------------------
	include 'mp.h'
	integer d1,d2
c
	call mpStd2(v(1),v(2),v(3),d1,d2)
	q(2) = d1/mpBase
	q(1) = d1 - q(2)*mpBase
	q(4) = d2/mpBase
	q(3) = d2 - q(4)*mpBase
	end
c************************************************************************
c* mpSign -- Return the sign of a multi-precision integer.
c& rjs
c+
	integer function mpSign(v)
c
	implicit none
	integer v(3)
c
c  Return the sign of a multi-precision integer.
c
c  Input:
c    v		Multi-precision integer to be considered.
c  Output:
c    mpSign	-1, 0 or +1 depending on whether v is negative, zero
c		or positive respectively.
c--
c------------------------------------------------------------------------
	integer d1,d2,s
c
	call mpStd2(v(1),v(2),v(3),d1,d2)
	s = d2
	if(s.eq.0)s = d1
	if(s.lt.0)then
	  s = -1
	else if(s.gt.0)then
	  s = 1
	endif
c
	mpSign = s
c
	end
c************************************************************************
c* mpCmp -- Compare two multi-precision integers.
c& rjs
c+
	integer function mpCmp(v1,v2)
c
	implicit none
	integer v1(3),v2(3)
c
c  Compute two multi-precision integers.
c  Specifically this computes v1-v2 and return -1, 0 or +1 depending
c  on whether the result is negative, zero or positive.
c
c  Input:
c    v1,v2	The multi-precision integers to be compared.
c  Output:
c    mpCmp	-1,0 or +1 depending on whether v1 is less than, equal
c		to or greater than v2 respectively.
c--
c------------------------------------------------------------------------
	integer t(3),v
c
	t(1) = v1(1)
	t(2) = v1(2)
	t(3) = v1(3)
	call mpSubmm(t,v2)
c
	if(t(3).le.0)call bug('f','Assertion in mpCmp failed')
	v = t(2)
	if(v.eq.0)v = t(1)
	if(v.lt.0)then
	  v = -1
	else if(v.gt.0)then
	  v = 1
	endif
c
	mpCmp = v
	end
c************************************************************************
c* mpAbs -- Take the absolute value of a multi-precision integer.
c& rjs
c+
	subroutine mpAbs(v)
c
	implicit none
	integer v(3)
c
c  Take the absolute value of a multi-precision integer.
c
c  Intput/Output:
c    v		The output, v = abs(v).
c--
c------------------------------------------------------------------------
	include 'mp.h'
	integer d1,d2
c
	if(((v(1).le.0.and.v(2).le.0).or.(v(1).ge.0.and.v(2).ge.0))
     *	    .and.v(3).ge.0)then
	  v(1) = abs(v(1))
	  v(2) = abs(v(2))
	else
	  call mpStd2(v(1),v(2),v(3),d1,d2)
	  v(1) = abs(d1)
	  v(2) = abs(d2)
	  v(3) = mpBase2
	endif
c
	end
c************************************************************************
c* mpNeg -- Negate a multi-precision integer.
c& rjs
c+
	subroutine mpNeg(v)
c
	implicit none
	integer v(3)
c
c  Negate a multi-precision integer.
c
c  Input/Output:
c    v		The output, v = -v
c
c--
c------------------------------------------------------------------------
	v(1) = -v(1)
	v(2) = -v(2)
	end
c************************************************************************
	subroutine mpInit
c
	implicit none
c
c------------------------------------------------------------------------
	integer k
	include 'mp.h'
c
c  Externals.
c
	integer ipmpar
	k = ipmpar(3)
	mpBase = sqrt(0.99*real(k/2))
	mpBase2 = mpBase*mpBase
	call mpStd2(k,0,mpBase2,mpMaxInt(1),mpMaxInt(2))
	mpMaxInt(3) = mpBase2
c
	end
c************************************************************************
c* mpFmt -- Format a multi-precision integer as a string.
c& rjs
c+
	subroutine mpFmt(out,v)
c
	implicit none
	integer v(3)
	character out*(*)
c
c  Format a multi-precision integer as a string.
c
c  Input:
c    v		The multi-precision integer being formatted.
c  Output:
c    out	The string containing a representation (base 10) of the 
c		multi-precision integer.
c--
c------------------------------------------------------------------------
	include 'mp.h'
	integer maxtd,maxdig
	parameter(maxtd=9,maxdig=10)
	character fmt*6
	integer dout(2*maxdig),ndo,ndi,ntd,baseo,l1,l2,i,din(4)
	logical neg
c
c  Externals.
c
	character itoaf*(maxtd)
	integer len1
c
	call mpStd4(v,din)
	ndi = 0
	do i=1,4
	  if(din(i).ne.0)then
	    ndi = i
	    neg = din(i).lt.0
	  endif
	  din(i) = abs(din(i))
	enddo
c
	if(ndi.eq.0)then
	  out = '0'
	else
	  ntd = 1
	  baseo = 10
	  dowhile(10*baseo.lt.mpBase.and.ntd.lt.maxtd)
	    ntd = ntd + 1
	    baseo = 10*baseo
	  enddo
	  write(fmt,'(a,i1,a,i1,a)') '(i',ntd,'.',ntd,')'
c
	  call mpNewBas(ndi,din,mpBase,2*maxdig,ndo,dout,baseo)
c
	  if(len(out).lt.ntd+1)call bug('f','Format overflow in mpFmt')
	  if(neg)then
	    out = '-'//itoaf(dout(ndo))
	  else
	    out = itoaf(dout(ndo))
	  endif
	  l1 = len1(out(1:ntd+1)) + 1
	  l2 = len(out)
c
	  if((ndo-1)*ntd+l1-1.gt.l2)call bug('f','Format overflow')
	  do i=ndo-1,1,-1
	    write(out(l1:l1+ntd-1),fmt)dout(i)
	    l1 = l1 + ntd
	  enddo
	endif
c
	end
c************************************************************************
	subroutine mpNewBas(ndi,din,basei,maxdo,ndo,dout,baseo)
c
	implicit none
	integer ndi,din(ndi),basei,maxdo,ndo,dout(maxdo),baseo
c
c  Internal routine used to change the radix base of the representation
c  of a multi-precision integer.
c
c  Input:
c    nd		The number of digits in the input.
c    din	An array of the digits of the input.
c    basei	Radix base of the input.
c    maxdo	The maximum possible number of digits for the output.
c    baseo	Radix base of the output.
c
c  Output:
c    dout	The digits for the output.
c    ndo	The number of digits in the output.
c
c------------------------------------------------------------------------
	integer carry,i,j
c
	dout(1) = 0
	ndo = 0
	do i=ndi,1,-1
	  carry = din(i)
	  do j=1,ndo
	    carry = carry + basei*dout(j)
	    dout(j) = mod(carry,baseo)
	    carry = carry/baseo
	  enddo
	  dowhile(carry.ne.0)
	    ndo = ndo + 1
	    if(ndo.gt.maxdo)
     *	      call bug('f','Integer overflow in mpNewBas')
	    dout(ndo) = mod(carry,baseo)
	    carry = carry/baseo
	  enddo
	enddo
c
	end
