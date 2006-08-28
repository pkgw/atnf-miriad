#-----------------------------------------------------------------------------
# GNUmakefile used to compile Miriad.
#
# Original: 2006/08/28, Mark Calabretta, ATNF
# $Id$
#-----------------------------------------------------------------------------
# Get common makefile variables and rules.
include $(MIR)/GNUmakedefs

# Subdirectories in which to invoke "allsys", in order.
#------------------------------------------------------
ALLSYSD  := $(findstring linpack,$(SUBDIRS))
ALLSYSD  += tools subs prog

# Static and static pattern rules.
#---------------------------------
allsys :: sysdirs $(ALLSYSD)

show ::
	-@ echo ""
	-@ echo "Variables defined in the top-level GNUmakefile"
	-@ echo "=============================================="
	-@ echo ""
	-@ echo "ALLSYSD  = $(ALLSYSD)"
