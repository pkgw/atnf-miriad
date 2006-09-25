#-----------------------------------------------------------------------------
# GNUmakefile used to compile Miriad.
#
# Original: 2006/08/28, Mark Calabretta, ATNF
# $Id$
#-----------------------------------------------------------------------------
# Get common makefile variables and rules.
include $(MIR)/GNUmakedefs

ifeq "$(MAKEMODE)" "system"
  # Subdirectories in which to invoke "allsys", in order.
  #------------------------------------------------------
  ALLSYSD  := $(findstring linpack,$(SUBDIRS))
  ALLSYSD  += scripts tools inc subs prog spec guides

  # The Miriad distribution is split into RCS, code, compiled documentation,
  # and platform-specific binary kits.
  DISTRCS  := .rcs ./RCS ./*/RCS ./*/*/RCS ./*/*/*/RCS
  DISTCODE := DISCLAIMER GNUmake* MIRRC* README.html cat guides inc linpack
  DISTCODE += prog scripts spec specdoc subs tests tools 
  DISTDOC  := doc html man progguide.ps userguide.ps
  DISTBINS := $(subst /bin,,$(wildcard */bin))

  show ::
	-@ echo ""
	-@ echo "Variables defined in the top-level GNUmakefile"
	-@ echo "=============================================="
	-@ echo ""
	-@ echo "ALLSYSD  = $(ALLSYSD)"
	-@ echo ""
	-@ echo "DISTRCS  = $(DISTRCS)"
	-@ echo "DISTCODE = $(DISTCODE)"
	-@ echo "DISTDOC  = $(DISTDOC)"
	-@ echo "DISTBINS = $(DISTBINS)"

  # Static and static pattern rules.
  #---------------------------------
  .PHONY : announce bookings dist updates

  allsys :: initial MIRRC MIRRC.sh sysdirs $(ALLSYSD)

  initial :: FORCE
	-@ echo "Rebuilding/updating Miriad for $(ARCH) machines."

  help ::
	-@ echo ""
	-@ echo "Targets defined in the top-level GNUmakefile"
	-@ echo "--------------------------------------------"
	-@ echo "     allsys: recursively rebuild or update Miriad."


  # The following rules are for ATNF use only.
  ifdef ATNF
    # Update the copy of the RPFITS library and include file via allsys.
    initial :: $(MIRINCD)/rpfits.inc $(MIRLIBD)/librpfits.a

    $(MIRLIBD)/librpfits.a : /usr/local/lib/librpfits.a
	-@ $(RM) $@
	   cp $< $@

    $(MIRINCD)/rpfits.inc : /usr/local/include/rpfits.inc
	-@ rcs -l $@ < /dev/null
	-@ $(RM) $@
	   cp $< $@
	 @ ci -u -m"Updated from /usr/local/include/rpfits.inc." $@

    # Regenerate the Miriad ftp distribution kits.
    dist : allsys
	-@ echo " "
	   tar cf miriad-rcs.tar $(DISTRCS)
	   gzip miriad-rcs.tar
	-@ $(RM) $(MIRFTPD)/miriad-rcs.tar.gz
	 @ mv miriad-rcs.tar.gz $(MIRFTPD)/
	-@ $(RM) .rcs-list
	 @ find . -name RCS | sort > .rcs-list
	   tar cXf .rcs-list miriad-code.tar $(DISTCODE)
	   gzip miriad-code.tar
	-@ $(RM) $(MIRFTPD)/miriad-code.tar.gz
	 @ mv miriad-code.tar.gz $(MIRFTPD)/
	 @ $(RM) .rcs-list
	   tar cf miriad-doc.tar $(DISTDOC)
	   gzip miriad-doc.tar
	-@ $(RM) $(MIRFTPD)/miriad-doc.tar.gz
	 @ mv miriad-doc.tar.gz $(MIRFTPD)/
	 @ for bin in $(DISTBINS) ; do \
	     echo "tar cf miriad-$$bin.tar $$bin" ; \
	     tar cf miriad-$$bin.tar $$bin ; \
	     echo "gzip miriad-$$bin.tar" ; \
	     gzip miriad-$$bin.tar ; \
	     $(RM) $(MIRFTPD)/miriad-$$bin.tar.gz ; \
	     mv miriad-$$bin.tar.gz $(MIRFTPD)/ ; \
	   done

    bookings :
	-@ echo ""
	-@ echo "Checking for users needing help..."
	 @ cd $(MIRROOT)/at_friends && ./CheckBookings.csh

    updates :
	-@ echo " "
	-@ echo "Creating updates..."
	 @ $(MIRBIND)/mirexport

    help ::
	-@ echo "       dist: allsys, then generate distribution kits."
	-@ echo "   bookings: check disk bookings."
	-@ echo "    updates: create ftp update tar files."
  endif

else
  # Programmer-oriented rules.
  #---------------------------
  help ::
	-@ echo ""
	-@ echo "No programmer-oriented rules are defined in the top-level"
	-@ echo "GNUmakefile."
endif
