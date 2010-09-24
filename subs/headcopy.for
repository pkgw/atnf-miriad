c* headcopy - Copy image header
c& bpw
c: utilities, image-i/o
c+
      subroutine headcopy(lIn, lOut, axMap, nAxMap, blc, trc)

      integer lIn, lOut, axMap(*), nAxMap, blc(*), trc(*)
c-----------------------------------------------------------------------
c Headcopy copies the small items (header keywords) and the history from
c one image to another.  It transfers all official header keywords with
c only a few exceptions:
c   - naxis and naxis# which are maintained by XYOPEN and XYZOPEN,
c   - datamin, datamax, and rms must be recalculated for the output
c     image and updated with WRHDR.
c
c The items listed in the KEYW array (below) are copied verbatim except
c for crpix, crval, cdelt, crota and ctype which may be deleted,
c exchanged, or reversed depending on axis permutations.  The axMap
c array defines the relation between new and old axes in the sense
c axMap(new) = old, e.g.
c
c                                     axMap
c                       nAxMap  (1)  (2)  (3)  (4)
c                       ------  ---  ---  ---  ---
c           direct copy:   4     1    2    3    4
c         delete z-axis:   2     1    2    0    0
c         delete x-axis:   2     2    3    0    0
c    output is zxy cube:   3     3    1    2    0
c       x-axis reversed:   4    -1    2    3    4
c         verbatim copy:  any    0    -    -    -
c
c   The last entry shows the shorthand used in the common case where
c   there are no axis permutations and the input and output images have
c   the same dimensions (no sub-imaging).
c
c   For reversed axes cdelt is multiplied by -1 and 180 is added to
c   crota.
c
c   For output images whose corners differ from the input image, crpix
c   may have to be corrected.  This is done using the arrays blc and trc
c   which give the corners of the output image relative to the pixel
c   numbers of the input image.  blc is used for non-reversed axes and
c   trc for reversed axes.  Often one can use the output of subroutine
c   boxinfo to get blc and trc.
c
c   Input:
c      lIn        Handle of input image.
c      lOut       Handle of output image.
c      axMap      Array that relates new axes with old (see above).
c                 Can be specified as 0 (scalar) meaning verbatim copy.
c      nAxMap     Dimension of axMap array - ignored if axMap(1) = 0.
c      blc        List of bottom-left-corners of input image region.
c                 Can be specified as 0 (scalar) meaning (1,1,...).
c      trc        List of top-right-corners of input image region.
c                 Can be specified as 0 (scalar) meaning
c                 (naxis1,naxis2,...) for the input image.
c
c $Id$
c-----------------------------------------------------------------------
      integer   CRPIX, CDELT, CRVAL, CROTA, CTYPE
      parameter (CRPIX = 1, CDELT = 2, CRVAL = 3, CROTA = 4, CTYPE = 5)

      integer   NKEYS
      parameter (NKEYS = 32)

      logical   verbtm
      integer   axLen, iAxIn, iAxOut, nAxOut, k
      double precision defVal(5), dvalue
      character avalue*80, keyIn*8, keyOut*8, keyw(NKEYS)*8

      logical   hdprsnt
      character itoaf*1
      external  hdprsnt, itoaf

      data defVal /0d0, 1d0, 0d0, 0d0, 0d0/
      data keyw /
     *    'crpix   ', 'cdelt   ', 'crval   ', 'crota   ', 'ctype   ',
     *    'bmaj    ', 'bmin    ', 'bpa     ', 'btype   ', 'bunit   ',
     *    'cellscal', 'date-obs', 'epoch   ', 'instrume', 'llrot   ',
     *    'ltype   ', 'lstart  ', 'lstep   ', 'lwidth  ', 'mostable',
     *    'niters  ', 'object  ', 'observer', 'obsra   ', 'obsdec  ',
     *    'obstime ', 'pbfwhm  ', 'pbtype  ', 'restfreq', 'telescop',
     *    'vobs    ', 'history '/
c-----------------------------------------------------------------------
c     All axes in the output image must have coordinate keywords.  
      call rdhdi(lOut, 'naxis', nAxOut, 0)

c     Loop for crpix, crval, cdelt, crota, and ctype.  
      verbtm = axMap(1).eq.0
      do k = 1, 5
c       Set default values.
        avalue = ' '
        dvalue = defVal(k)

        do iAxOut = 1, nAxOut
          keyOut = keyw(k)(1:5) // itoaf(iAxOut)

          if (verbtm .and. hdprsnt(lIn, keyOut)) then
c           Copy verbatim.
            call hdcopy(lIn, lOut, keyOut)

          else
c           Handle axis permutations.
            if (.not.verbtm .and. iAxOut.le.nAxMap) then
              iAxIn = axMap(iAxOut)
            else
c             Use default values.
              iAxIn = 0
            endif

            if (iAxIn.ne.0) then
              keyIn = keyw(k)(1:5) // itoaf(abs(iAxIn))

              if (hdprsnt(lIn, keyIn)) then
c               Read it from the input image.
                if (k.eq.CTYPE) then
                  call rdhda(lIn, keyIn, avalue, ' ')
                else
                  call rdhdd(lIn, keyIn, dvalue, 0d0)
                endif

                if (iAxIn.gt.0) then
                  if (k.eq.CRPIX) then
c                   Sub-imaging.
                    if (blc(1).ne.0) then
                      dvalue = dvalue - dble(blc(iAxIn)-1)
                    endif
                  endif

                else if (iAxIn.lt.0) then
c                 Axis reversal.
                  iAxIn = -iAxIn

                  if (k.eq.CRPIX) then
c                   Sub-imaging.
                    if (trc(1).ne.0) then
                      axLen = trc(iAxIn)
                    else
                      keyIn = 'naxis' // itoaf(iAxIn)
                      call rdhdi(lOut, keyIn, axLen, 0)
                    endif

                    dvalue = dble(axLen+1) - dvalue

                  else if (k.eq.CDELT) then
                    dvalue = -dvalue
                  else if (k.eq.CROTA) then
                    dvalue = dvalue - 180d0
                  endif
                endif
              endif
            endif

c           Write it to the output image.
            if (k.eq.CTYPE) then
              call wrhda(lOut, keyOut, avalue)
            else
              call wrhdd(lOut, keyOut, dvalue)
            endif
          endif
        enddo
      enddo

c     Copy the remaining items verbatim.
      do k = 6, NKEYS
         call hdcopy(lIn, lOut, keyw(k))
      enddo

      return
      end
