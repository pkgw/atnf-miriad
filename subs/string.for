c************************************************************************
c History
c   rjs     long time ago
c   nebk    19Sep89 initialize RESULT in ATOI
c   jm      29Jan90 included tabs as field delimiters in GetField
c   jm      29Jan90 removed any leading spaces in ATOx
c   pjt      5may90 getfield: check if enough space
c   jm/pxt  24aug90 added dtoa, rtoa, and padleft.
c   bpw     11sep90 add rtfmt, match, matchdcd, equals, substr, isupper,
c                   islower, isdigit, isalpha, isalnum, indek
c   rjs     22feb91 Deleted getparm, both because it was obsolete,
c                   and because it conflicted with another routine.
c   bpw     25feb91 Changed atoi, etc to atoif, etc to avoid naming
c                   conflicts.
c   bpw     10apr91 Fixed removal of leading/trailing quotes in getfield.
c   bpw     07aug91 Add matchnr
c   bpw     24sep91 Add atorf, matorf, matodf
c   bpw     30sep91 extract match routines
c   mjs     29mar92 "dtoaf" converted to use miriad's "numbpg" rather
c		    than PGPLOT's "pgnumb" (relevent to shared lib
c		    implementations).
c   mchw    30nov92 added hex and octal conversions to atoif.
c   bpw     16aug93 changed an accidental len to len1 in atodf
c   rjs      3nov94 Treat [] and {} as brackets in GetField.
c   rjs     06feb95 Fixed handling of quotes in getfield. What did bpw do?
c   rjs     25jul97 Treat " and ' as quote characters.
c   rjs     03aug98 Included updated version of matodf and matorf.
c   rjs     05feb01 Added st routines.
c
c $Id$
c************************************************************************
c* stcat - Concatenate two strings together (avoiding blank pads).
c& rjs
c: string
c+
	character*(*) function stcat(a,b)
c
	character a*(*),b*(*)
c
c  Concatenate two strings together, ignoring blank padding.
c
c  Inputs:
c    a,b  The strings to be concatenated.
c  Output:
c    stcat The concatenated strings a//b.
c------------------------------------------------------------------------
	integer length
c
c  Externals.
c
	integer len1
c
	length = len1(a)
	stcat = a(1:length)//b
	end
c************************************************************************
	character*(*) function streal(a,b)
c
	real a
	character b*(*)
c------------------------------------------------------------------------
	integer i
	character line*32
c
	write(line,b)a
	i = 1
	dowhile(line(i:i).eq.' ')
	  i = i + 1
	enddo
	streal = line(i:)
	end
c************************************************************************
c* atoif -- Convert a string into an integer.
c& rjs
c: strings
c+
	subroutine atoif(string,result,ok)
c
	character string*(*)
	integer result
	logical ok
c
c  Convert a string into its corresponding integer.
c  The string can be input as a hex, octal or decimal number using a
c  prefix 0x, %x or h for hex, o or %o for octal, and +, - or nothing
c  for decimal. Either case is ok.
c
c  Input:
c    string	The ascii string, containing the integer.
c  Output:
c    result	The integer value.
c    ok		This will be true if the decoding succeeded.
c--
c------------------------------------------------------------------------
	integer k,i,i0,sign1,nbase
	character ch*1, ch2*2
c
c  External
c
	integer len1
c
c  Ignore leading blanks.
c
	do i = 1, len(string)
	  if (string(i:i).gt.' ') go to 10
	enddo

 10     result = 0
	ok = i.le.len(string)
	if (.not.ok) return

c       Find number base.
	call lcase(string)
	nbase = 10
	sign1 = 1
	i0    = i

	ch = string(i:i)
	if(ch.eq.'+')then
	  i0 = i + 1
	else if(ch.eq.'-')then
	  sign1 = -1
	  i0 = i + 1
	else if(ch.eq.'h')then
	  nbase = 16
	  i0 = i + 1
	else if(ch.eq.'o')then
	  nbase = 8
	  i0 = i + 1
	else if (i.lt.len(string)) then
	  ch2 = string(i:i+1)
	  if(ch2.eq.'0x' .or. ch2.eq.'%x')then
	    nbase = 16
	    i0 = i + 2
	  else if(ch2.eq.'%o')then
	    nbase = 8
	    i0 = i + 2
	  endif
	endif

c	Decode value.
	ok = i0.le.len1(string)
	if (.not.ok) return

	do i = i0, len1(string)
	  ch = string(i:i)
	  if(ch.ge.'0' .and. ch.le.'9')then
	    k = ichar(ch) - ichar('0')
	  else if(ch.ge.'a' .and. ch.le.'f')then
	    k = ichar(ch) - ichar('a') + 10
	  else
	    k = -1
	  endif

	  if(k.ge.0 .and. k.lt.nbase)then
	    result = nbase*result + k
	  else
	    result = 0
	    ok = .false.
	    return
	  endif
	enddo

	result = sign1 * result

	end
c************************************************************************
c* atorf -- Convert a string into a real
c& rjs
c: strings
c+
	subroutine atorf(string,result,ok)
c
	character string*(*)
	real result
	logical ok
c
c  Convert a string into its corresponding real.
c
c  Input:
c    string	The ascii string, containing the real.
c  Output:
c    result	The double precision value.
c    ok		This will be true if the decoding succeeded.
c--
c------------------------------------------------------------------------
       double precision dresult
       call atodf(string,dresult,ok)
       result= dresult
c
       end
C************************************************************************
c* atodf -- Convert a string into a double precision.
c& rjs
c: strings
c+
	subroutine atodf(string,d,ok)
c
	character string*(*)
	double precision d
	logical ok
c
c  Convert a string into its corresponding double precision.
c
c  Input:
c    string	The ascii string, containing the double precision.
c  Output:
c    result	The double precision value.
c    ok		This will be true if the decoding succeeded.
c--
c------------------------------------------------------------------------
	integer whole,frac,bad,expon
	parameter(whole=1,frac=2,bad=3,expon=4)
	integer length,l,expo,digits,state,iw0,temp,len1
	double precision wd,w0,ww0,dtemp
c
	length = len1(string)
	wd = 10
	ww0 = 1
	w0 = 1
	digits = 0
	state = whole
	dtemp = 0
c
	l = 1
	do while (string(l:l) .le. ' ')
	  l = l + 1
	enddo
	if(string(l:l).eq.'-')then
	  w0 = -1
	else if(string(1:1).ne.'+')then
	  l = l - 1
	endif
c
c  Handle the mantissa.
c
	do while(l.lt.length.and.(state.eq.whole.or.state.eq.frac))
	  l = l + 1
	  if(string(l:l).ge.'0'.and.string(l:l).le.'9')then
	    digits = digits + 1
	    temp = ichar(string(l:l)) - ichar('0')
	    dtemp = wd * dtemp + w0 * temp
	    w0 = ww0 * w0
	  else if(string(l:l).eq.'.')then
	    if(state.eq.whole)then
	      state = frac
	      wd = 1
	      w0 = 0.1 * w0
	      ww0 = 0.1
	    else
	      state = bad
	    endif
	  else if(string(l:l).eq.'e'.or.string(l:l).eq.'E'.or.
     *		  string(l:l).eq.'d'.or.string(l:l).eq.'D')then
	    state = expon
	  else
	    state = bad
	  endif
	enddo
c
c  Handle the exponent if there is one.
c
	if(state.eq.expon.and.digits.gt.0)then
	  digits = 0
	  expo = 0
	  iw0 = 1
	  if(l.lt.length)then
	    l = l + 1
	    if(string(l:l).eq.'-')then
	      iw0 = -1
	    else if(string(l:l).ne.'+')then
	      l = l - 1
	    endif
	  endif
	  do while(l.lt.length)
	    l = l + 1
	    if(string(l:l).ge.'0'.and.string(l:l).le.'9')then
	      temp = ichar(string(l:l))-ichar('0')
	      expo = 10*expo + iw0 * temp
	      digits = digits + 1
	    else
	      state = bad
	    endif
	  enddo
	  if(digits.gt.0.and.state.ne.bad)dtemp = dtemp * 10.d0**expo
	endif
c
c  Finish up.
c
	ok = state.ne.bad.and.digits.gt.0
	if ( ok ) d = dtemp
	end
c************************************************************************
c* matodf - convert a string to many double precision numbers
c& rjs
c: strings
c+
      subroutine matodf( string, array, n, ok )
c
      integer          n
      character*(*)    string
      double precision array(n)
      logical          ok
c
c Convert a string to many double precision numbers.
c
c Input:
c   string      The ascii string containing the numbers
c   n           The number of values wanted
c Output:
c   array       The double precision numbers in the string
c   ok          True if the decoding succeeded
c--
c-----------------------------------------------------------------------
      character*30     substr
      double precision dval
      integer          tlen,k,k1,k2
c
c  External
c
	integer len1
c
        k1 = 1
        k2 = len1(string)
        k  = 0
        do while (k1.le.k2.and.k.lt.n)
          call getfield(string, k1, k2, substr, tlen)
          call atodf(substr,dval,ok)
          if(ok) then
            k = k+1
            array(k) = dval
          else
            return
          endif
        enddo
	end
c************************************************************************
c* matorf - convert a string to many real numbers
c& rjs
c: strings
c+
      subroutine matorf( string, array, n, ok )

      character*(*)    string
      integer          n
      real             array(n)
      logical          ok
c
c Convert a string to many real numbers.
c
c Input:
c   string      The ascii string containing the numbers
c   n           The number of values wanted
c Output:
c   array       The real numbers in the string
c   ok          True if the decoding succeeded
c--
c-----------------------------------------------------------------------
      character*30     substr
      double precision dval
      integer          tlen,k,k1,k2
c
c  External
c
	integer len1
c
        k1 = 1
        k2 = len1(string)
        k  = 0
        do while (k1.le.k2.and.k.lt.n)
          call getfield(string, k1, k2, substr, tlen)
          call atodf(substr,dval,ok)
          if(ok) then
            k = k+1
            array(k) = dval
          else
            return
          endif
        enddo
	end
c************************************************************************
c* Itoaf -- Convert an integer into a string.
c& rjs
c: strings
c+
	character*(*) function itoaf(n)
c
	integer n
c
c  Convert an integer into its ascii representation. It is returned
c  left justified.
c
c  Input:
c    n		The integer to convert.
c  Output:
c    itoaf      The formated integer, left justified.
c--
c------------------------------------------------------------------------
	integer nono,lstring
	parameter(nono=30,lstring=12)
	character string*(lstring)
	integer i
	character numbers(nono)*2
	save numbers
	data numbers/'1 ','2 ','3 ','4 ','5 ','6 ','7 ','8 ','9 ',
     *		'10','11','12','13','14','15','16','17','18','19',
     *		'20','21','22','23','24','25','26','27','28','29','30'/
c
	if(n.gt.0.and.n.le.nono)then
	  itoaf = numbers(n)
	else
	  write(string,'(i12)')n
	  i = 1
	  dowhile(string(i:i).eq.' ')
	    i = i + 1
	  enddo
	  itoaf = string(i:lstring)
	endif
	end
c************************************************************************
c* mitoaf -- Convert many integers into a string.
c& rjs
c: strings
c+
	subroutine mitoaf(array,n,line,length)
c
	integer n,array(n),length
	character line*(*)
c
c  Convert multiple integers to a neatly formatted character string.
c
c  Input:
c    array	Array of integers to be converted to a string.
c    n		Number of integers.
c
c  Output:
c    line	The formatted string. Integers are separated by
c		commas (no blanks).
c    length	The length of the line.
c--
c------------------------------------------------------------------------
	integer lmax
	parameter(lmax=12)
	integer i,newlen,l
	character txt*(lmax)
c
	length = 0
	do i=1,n
	  if(i.ne.n)then
	    write(txt,'(i11,a)')array(i),','
	  else
	    write(txt,'(i12)')array(i)
	  endif
	  l = 1
	  dowhile(txt(l:l).eq.' ')
	    l = l + 1
	  enddo
	  newlen = min(len(line),length + lmax - l + 1)
	  if(length.lt.newlen)line(length+1:newlen) = txt(l:lmax)
	  length = length + lmax - l + 1
	enddo
	end
c***********************************************************************
c* dtoaf -- Convert a double precision value into a string.
c& rjs
c: strings
c+
      character*(*) function dtoaf(value, form, nsf)
c
      double precision value
      integer form, nsf
c
c  Format a double precision value into a string using miriad's numbpg.
c
c  Inputs:
c    value         The value (double precision)
c    form          How the number is formatted:
c                    form = 0 - uses either decimal or
c                               exponential notation
c                    form = 1 - uses decimal notation
c                    form = 2 - uses exponential notation
c    nsf           Number of significant figures for output
c                    (best range is between 1 and 6)
c--
c    pxt july 1990
c-----------------------------------------------------------------------
      integer ns, np, nv
      character tmpstr*25
c
      if (value .eq. 0.0D0) then
        np = 0
        nv = 0
      else
        np = log10(abs(value)) - nsf
        nv = nint(value * (10.0 ** (-np)))
      endif
      call numbpg(nv, np, form, tmpstr, ns)
      dtoaf = tmpstr(1:ns)
      end
c***********************************************************************
c* rtoaf -- Convert a real value into a string.
c& rjs
c: strings
c+
      character*(*) function rtoaf(value, form, nsf)
c
      real value
      integer form, nsf
c
c  Format a real value into a string using miriad's numbpg.
c
c  Inputs:
c    value         The value (must be real).
c    form          How the number is formatted:
c                    form = 0 - uses either decimal or
c                               exponential notation
c                    form = 1 - uses decimal notation
c                    form = 2 - uses exponential notation
c    nsf           Number of significant figures for output
c                    (best range is between 1 and 6)
c--
c    pxt july 1990
c-----------------------------------------------------------------------
      character dtoaf*25
c
      rtoaf = dtoaf(dble(value), form, nsf)
      end
c************************************************************************
c* GetTok -- Extract a token from a string.
c& rjs
c: strings
c+
	subroutine gettok(string,k1,k2,token,length)
c
	integer k1,k2,length
	character string*(*),token*(*)
c
c  Extract an alphanumeric token from the string. White space at the
c  start of the string is first skipped.
c
c  Input:
c    string	The string containing the tokens.
c  Input/Output:
c    k1,k2	These delimit the characters in the string that are
c		to be processed. On output, k1 is incremented to point
c		beyond the token just returned.
c  Output:
c    token	The returned token.
c    length	The length of the token.
c--
c------------------------------------------------------------------------
	integer k
	logical more
c
c  Skip leading white.
c
	more = .true.
	do while(k1.le.k2.and.more)
	  if(string(k1:k1).le.' ')then
	    k1 = k1 + 1
	  else
	    more = .false.
	  endif
	enddo
c
c  Go over ths string.
c
	k = k1
	more = .true.
	do while(k1.le.k2.and.more)
	  if(  (string(k1:k1).ge.'0'.and.string(k1:k1).le.'9').or.
     *	       (string(k1:k1).ge.'a'.and.string(k1:k1).le.'z').or.
     *         (string(k1:k1).ge.'A'.and.string(k1:k1).le.'Z').or.
     *		string(k1:k1).eq.'$')then
	    k1 = k1 + 1
	  else
	    more = .false.
	  endif
	enddo
	length = k1 - k
	if(length.gt.0)token = string(k:k1-1)
	end
c************************************************************************
c* GetField -- Extract a field from a string.
c& rjs
c: strings
c+
	subroutine getfield(string,k1,k2,token,length)
c
	integer k1,k2,length
	character string*(*),token*(*)
c
c  Extract a "field". This looks for a separator, which may be a comma
c  or white character. It also worries about bracketing with () characters
c  and strings enclosed within quotes.
c
c  Input:
c    string	The string containing the fields.
c  Input/Output:
c    k1,k2	These delimit the characters in the string that are
c		to be processed. On output, k1 is incremented to point
c		beyond the token just returned.
c  Output:
c    token	The returned token.
c    length	The length of the token.
c--
c------------------------------------------------------------------------
	integer k,l,depth,k1old
	logical more,quoted
	character c*1,line*80,quotec*1
c
c  Skip leading white.
c
	k1old = k1
	more = .true.
	do while(k1.le.k2.and.more)
	  if(string(k1:k1).le.' '.or.string(k1:k1).eq.',')then
	    k1 = k1 + 1
	  else
	    more = .false.
	  endif
	enddo
c
c  Go over this string.
c
	k = k1
	more = .true.
	quoted = .false.
	depth = 0
	do while(k1.le.k2.and.more)
	  c = string(k1:k1)
	  if(quoted)then
	    quoted = c.ne.quotec
	  else if(c.eq.'"'.or.c.eq.'''')then
	    quoted = .true.
	    quotec = c
	  else
	    if(c.eq.'('.or.c.eq.'['.or.c.eq.'{')then
	      depth = depth + 1
	    else if(c.eq.')'.or.c.eq.']'.or.c.eq.'}')then
	      depth = depth - 1
	    else if(c.le.' '.or.c.eq.',')then
	      more = depth.gt.0
	    endif
	  endif
	  if(more)k1 = k1 + 1
	enddo
c
c  Remove leading and trailing quotes.
c
	l = k1 - 1
	if( (string(k:k).eq.'"' .and. string(l:l).eq.'"' ) .or.
     *	    (string(k:k).eq.''''.and. string(l:l).eq.'''')) then
           k = k + 1
	   l = l - 1
	endif
c
c  Return the string and its length (plus stupid bug construction -- pjt)
c
	length = l - k + 1
	if(length.gt.0) then
            if (length.gt.len(token)) then
		line = string(k1old:k2)
                call bug('w','GETFIELD: string too long: '//line)
            endif
	    token = string(k:l)
	endif
	end
c************************************************************************
c* Spanchar -- Skip over a particular character.
c& rjs
c: strings
c+
	subroutine spanchar(string,k1,k2,c)
c
	character string*(*),c*1
	integer k1,k2
c
c  Skip over any occurrences of the character "c" in the string.
c
c  Input:
c    string	The string containing the characters.
c  Input/Output:
c    k1,k2	These delimit the characters in the string that are
c		to be processed. On output, k1 is incremented to point
c		beyond the character that was being skipped.
c
c------------------------------------------------------------------------
	logical more
c
	more = .true.
	do while(k1.le.k2.and.more)
	  if(string(k1:k1).eq.c)then
	    k1 = k1 + 1
	  else
	    more = .false.
	  endif
	enddo
	end
c************************************************************************
c* ScanChar -- Scan a string for a character.
c& rjs
c: strings
c+
	subroutine scanchar(string,k1,k2,c)
c
	character string*(*),c*1
	integer k1,k2
c
c  This scans forward through a string, until it finds a character "c".
c
c  Input:
c    string	The string containing the characters.
c  Input/Output:
c    k1,k2	These delimit the characters in the string that are
c		to be processed. On output, k1 is incremented to point
c		to the character that was being scanned for.
c--
c------------------------------------------------------------------------
	logical more
c
	more = .true.
	do while(k1.le.k2.and.more)
	  if(string(k1:k1).ne.c)then
	    k1 = k1 + 1
	  else
	    more = .false.
	  endif
	enddo
	end
c************************************************************************
c* Len1 -- Determine the unpadded length of a character string.
c& rjs
c: strings
c+
	integer function len1(string)
c
	character string*(*)
c
c  This determines the unblanked length of a character string.
c
c  Input:
c    string	The character string that we are interested in.
c  Output:
c    len1	The unpadded length of the character string.
c--
c------------------------------------------------------------------------
	integer l
	logical more
c
	l = len(string)
	more = .true.
	do while(l.gt.0.and.more)
	  if(string(l:l).le.' '.or.string(l:l).gt.'~')then
	    l = l - 1
	  else
	    more = .false.
	  endif
	enddo
c
	len1 = l
	end
c************************************************************************
c* Lcase -- Convert a string to lower case.
c& rjs
c: strings
c+
	subroutine lcase(string)
c
	character string*(*)
c
c  Convert a string to lower case.
c
c  Input/Output:
c    string	The string to be converted to lower case.
c--
c------------------------------------------------------------------------
	integer i
	character c*1
c
	do i=1,len(string)
	  c = string(i:i)
	  if(c.ge.'A'.and.c.le.'Z')
     *	    string(i:i) = char(ichar(c)-ichar('A')+ichar('a'))
	enddo
	end
c************************************************************************
c* Ucase -- Convert string to upper case.
c& rjs
c: strings
c+
	subroutine ucase(string)
c
	character string*(*)
c
c  Convert a string to upper case.
c
c  Input/Output:
c    string	The string to be converted to upper case.
c--
c------------------------------------------------------------------------
	integer i
	character c*1
c
	do i=1,len(string)
	  c = string(i:i)
	  if(c.ge.'a'.and.c.le.'z')
     *	    string(i:i) = char(ichar(c)-ichar('a')+ichar('A'))
	enddo
	end
c***********************************************************************
c* PadLeft -- Right justify a string to length characters.
c& rjs
c: strings
c+
      subroutine padleft(string, length)
c
      character string*(*)
      integer length
c
c  Right-justify a string to a specified length.
c
c  Inputs:
c    string        The string to be padded.
c    length        The length of the string after padding.
c  Outputs:
c    string        The string with spaces added to the left.
c--
c    pxt july 1990
c-----------------------------------------------------------------------
      character temp*132
      integer ns, i
c  Externals
      integer Len1
c
      ns = Len1(string)
      if ((len(string) .lt. length) .or. (ns .gt. len(temp))) then
        call bug('w', 'PADLEFT:  Input string overflows buffer.')
        call bug('w', 'No padding performed.')
        return
      endif
      if (ns .lt. length) then
        i = length - ns + 1
        temp = string(1:ns)
        string = ' '
        string(i:length) = temp(1:ns)
      endif
c  The string(1:length) now has the blanks and the value.
      end
************************************************************************
c* indek - get position of substring in a string, returning length if not found
c& rjs
c: strings
c+
      integer function indek ( string, substrng )
c
      character*(*) string, substrng

c Indek works basically the same as the intrinsic function index, but if the
c substring is not found, it does not return 0, but the length of the string
c plus 1.
c Example of use:
c filename = dir(:nelc(dir)) // name( : indek(name,'.')-1 ) // ' '
c Now only one statement is needed instead of three (find index(name,'.'); test'
c if ok; the above)
c
c Input:
c   string:      the string to search in
c   substrng:   the substring to search for
c--
      integer       idx, len1
      idx                  = index ( string, substrng )
      if( idx .eq. 0 ) idx = len1( string ) + 1
      indek = idx
c
      end
