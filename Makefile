#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *#
#*                                                                           *#
#*                  This file is part of the class library                   *#
#*       SoPlex --- the Sequential object-oriented simPlex.                  *#
#*                                                                           *#
#*    Copyright (C) 1996-2015 Konrad-Zuse-Zentrum                            *#
#*                            fuer Informationstechnik Berlin                *#
#*                                                                           *#
#*  SoPlex is distributed under the terms of the ZIB Academic Licence.       *#
#*                                                                           *#
#*  You should have received a copy of the ZIB Academic License              *#
#*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  *#
#*                                                                           *#
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *#

#@file    Makefile
#@brief   SoPlex Makefile
#@author  Thorsten Koch
#@author  Ambros Gleixner

#-----------------------------------------------------------------------------
# paths variables
#-----------------------------------------------------------------------------

# define to be able to locate library files
ifeq ($(OSTYPE),mingw)
SPXDIR		=	./
else
SPXDIR		=	$(realpath .)
endif

INSTALLDIR	=


#-----------------------------------------------------------------------------
# detect host architecture
#-----------------------------------------------------------------------------

include make/make.detecthost


#-----------------------------------------------------------------------------
# default settings
#-----------------------------------------------------------------------------

VERSION		:=	2.2.0
SPXGITHASH	=

VERBOSE		=	false
SHARED		=	false
OPT		=	opt
STATICLIBEXT	=	a
SHAREDLIBEXT	=	so
LIBEXT		=	$(STATICLIBEXT)
EXEEXTENSION	=	
TEST		=	quick
ALGO		=  1 2 3 4
LIMIT		=  #
SETTINGS	=	default
TIME		=	3600
RESDIR		=	results
MAKESOFTLINKS	=	true
SOFTLINKS	=
LINKSINFO	=

# these variables are needed for cluster runs
MEM		=	2000
CONTINUE	=	false

# will this be compiled for PARASCIP? (disables output because it uses global variables)
PARASCIP	=	false

# will this be compiled with the 1.x interface?
LEGACY		=	false

GMP		=	true
ZLIB		=	true
EGLIB		=	false

COMP		=	gnu
CXX		=	g++
CXX_c		=	-c # the trailing space is important
CXX_o		=	-o # the trailing space is important
LINKCXX		=	$(CXX)
LINKCXX_L	=	-L
LINKCXX_l	=	-l
LINKCXX_o	=	-o # the trailing space is important
LINKLIBSUFFIX	=
DCXX		=	$(CXX)
LINT		=	flexelint
AR		=	ar
AR_o		=
RANLIB		=	ranlib
DOXY		=	doxygen
VALGRIND	=	valgrind

READ		=	read -e
LN_s		=	ln -s

LIBBUILD	=	$(AR)
LIBBUILD_o	=	$(AR_o)
LIBBUILDFLAGS	=       $(ARFLAGS)

CPPFLAGS	=	-Isrc
CXXFLAGS	=
BINOFLAGS	=	
LIBOFLAGS	=	
LDFLAGS		=	
ARFLAGS		=	cr
DFLAGS		=	-MM
VFLAGS		=	--tool=memcheck --leak-check=yes --show-reachable=yes #--gen-suppressions=yes

SOPLEXDIR	=	$(realpath .)
SRCDIR		=	src
BINDIR		=	bin
LIBDIR		=	lib
INCLUDEDIR	=	include
NAME		=	soplex
LIBHEADER	=	array.h \
			basevectors.h \
			classarray.h \
			clufactor.h \
			clufactor_rational.h \
			cring.h \
			dataarray.h \
			datahashtable.h \
			datakey.h \
			dataset.h \
			didxset.h \
			dsvectorbase.h \
			dsvector.h \
			dvectorbase.h \
			dvector.h \
			exceptions.h \
			gzstream.h \
			idlist.h \
			idxset.h \
			islist.h \
			lpcolbase.h \
			lpcol.h \
			lpcolsetbase.h \
			lpcolset.h \
			lprowbase.h \
			lprow.h \
			lprowsetbase.h \
			lprowset.h \
			mpsinput.h \
			nameset.h \
			random.h \
			rational.h \
			ratrecon.h \
			slinsolver.h \
			slinsolver_rational.h \
			slufactor.h \
			slufactor_rational.h \
			solbase.h \
			sol.h \
			soplex.h \
			soplexlegacy.h \
			sorter.h \
			spxalloc.h \
			spxautopr.h \
			spxbasis.h \
			spxboundflippingrt.h \
			spxdantzigpr.h \
			spxdefaultrt.h \
			spxdefines.h \
			spxdevexpr.h \
			spxequilisc.h \
			spxfastrt.h \
			spxfileio.h \
			spxgeometsc.h \
			spxgithash.h \
			spxharrisrt.h \
			spxhybridpr.h \
			spxid.h \
			spxlpbase.h \
			spxlp.h \
			spxmainsm.h \
			spxout.h \
			spxparmultpr.h \
			spxpricer.h \
			spxratiotester.h \
			spxscaler.h \
			spxsimplifier.h \
			spxsolver.h \
			spxstarter.h \
			spxsteepexpr.h \
			spxsteeppr.h \
			spxsumst.h \
			spxvectorst.h \
			spxweightpr.h \
			spxweightst.h \
			ssvectorbase.h \
			ssvector.h \
			statistics.h \
			svectorbase.h \
			svector.h \
			svsetbase.h \
			svset.h \
			timer.h \
			unitvectorbase.h \
			unitvector.h \
			updatevector.h \
			vectorbase.h \
			vector.h
LIBOBJ		= 	changesoplex.o \
			clufactor.o \
			clufactor_rational.o \
			didxset.o \
			enter.o \
			gzstream.o \
			idxset.o \
			leave.o \
			mpsinput.o \
			nameset.o \
			rational.o \
			ratrecon.o \
			slufactor.o \
			slufactor_rational.o \
			solverational.o \
			solvereal.o \
			soplex.o \
			soplexlegacy.o \
			spxautopr.o \
			spxbasis.o \
			spxboundflippingrt.o \
			spxbounds.o \
			spxchangebasis.o \
			spxdantzigpr.o \
			spxdefaultrt.o \
			spxdefines.o \
			spxdesc.o \
			spxdevexpr.o \
			spxequilisc.o \
			spxfastrt.o \
			spxfileio.o \
			spxgeometsc.o \
			spxgithash.o \
			spxharrisrt.o \
			spxhybridpr.o \
			spxid.o \
			spxlpbase_rational.o \
			spxlpbase_real.o \
			spxmainsm.o \
			spxout.o \
			spxparmultpr.o \
			spxquality.o \
			spxscaler.o \
			spxshift.o \
			spxsolve.o \
			spxsolver.o \
			spxstarter.o \
			spxsteeppr.o \
			spxsumst.o \
			spxvecs.o \
			spxvectorst.o \
			spxweightpr.o \
			spxweightst.o \
			spxwritestate.o \
			statistics.o \
			usertimer.o \
			wallclocktimer.o \
			updatevector.o
BINOBJ		=	soplexmain.o
EXAMPLEOBJ	=	example.o
REPOSIT		=	# template repository, explicitly empty  #spxproof.o 

BASE		=	$(OSTYPE).$(ARCH).$(COMP).$(OPT)

LINKSMARKERFILE	=	$(LIBDIR)/linkscreated.$(OSTYPE).$(ARCH).$(COMP)$(LINKLIBSUFFIX).$(EGLIB)
LASTSETTINGS	=	$(OBJDIR)/make.lastsettings

SPXGITHASHFILE	=	$(SRCDIR)/git_hash.cpp

#------------------------------------------------------------------------------
#--- NOTHING TO CHANGE FROM HERE ON -------------------------------------------
#------------------------------------------------------------------------------

GCCWARN		=	-pedantic -Wall -W -Wpointer-arith -Wcast-align -Wwrite-strings \
			-Wconversion -Wsign-compare -Wshadow \
			-Wredundant-decls -Wdisabled-optimization \
			-Wctor-dtor-privacy -Wnon-virtual-dtor -Wreorder \
			-Woverloaded-virtual -Wsign-promo -Wsynth -Wundef \
			-Wcast-qual \
			-Wmissing-declarations \
			-Wno-unused-parameter -Wno-strict-overflow -Wno-long-long
#			-Wold-style-cast
#			-Weffc++


#-----------------------------------------------------------------------------
include make/make.$(BASE)
-include make/local/make.$(HOSTNAME)
-include make/local/make.$(HOSTNAME).$(COMP)
-include make/local/make.$(HOSTNAME).$(COMP).$(OPT)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# SHARED Libaries
#-----------------------------------------------------------------------------

ifeq ($(SHARED),true)
CPPFLAGS	+=	-fPIC
LIBEXT		=	$(SHAREDLIBEXT)
LIBBUILD	=	$(LINKCXX)
LIBBUILDFLAGS	+=      -shared
LIBBUILD_o	= 	-o # the trailing space is important
ARFLAGS		=
RANLIB		=
endif

CPPFLAGS	+=	$(USRCPPFLAGS)
CXXFLAGS	+=	$(USRCXXFLAGS)
LDFLAGS		+=	$(USRLDFLAGS)
ARFLAGS		+=	$(USRARFLAGS)
DFLAGS		+=	$(USRDFLAGS)

#-----------------------------------------------------------------------------
# PARASCIP
#-----------------------------------------------------------------------------

PARASCIPDEP	:=	$(SRCDIR)/depend.parascip
PARASCIPSRC	:=	$(shell cat $(PARASCIPDEP))

ifeq ($(PARASCIP),true)
CPPFLAGS	+=	-DDISABLE_VERBOSITY
endif

#-----------------------------------------------------------------------------
# LEGACY
#-----------------------------------------------------------------------------

LEGACYDEP	:=	$(SRCDIR)/depend.legacy
LEGACYSRC	:=	$(shell cat $(LEGACYDEP))

ifeq ($(LEGACY),true)
CPPFLAGS	+=	-DSOPLEX_LEGACY
endif

#-----------------------------------------------------------------------------
# Main Program
#-----------------------------------------------------------------------------

BINNAME		=	$(NAME)-$(VERSION).$(BASE)
EXAMPLENAME	=	example.$(BASE)
LIBNAME		=	$(NAME)-$(VERSION).$(BASE)
BINFILE		=	$(BINDIR)/$(BINNAME)$(EXEEXTENSION)
EXAMPLEFILE	=	$(BINDIR)/$(EXAMPLENAME)$(EXEEXTENSION)
LIBFILE		=	$(LIBDIR)/lib$(LIBNAME).$(LIBEXT)
LIBSHORTLINK	=	$(LIBDIR)/lib$(NAME).$(LIBEXT)
LIBLINK		=	$(LIBDIR)/lib$(NAME).$(BASE).$(LIBEXT)
BINLINK		=	$(BINDIR)/$(NAME).$(BASE)$(EXEEXTENSION)
BINSHORTLINK	=	$(BINDIR)/$(NAME)$(EXEEXTENSION)
DEPEND		=	src/depend

# potential valgrind suppression file name
VSUPPNAME	= 	$(OSTYPE).$(ARCH).$(COMP).supp

OBJDIR		=	obj/O.$(BASE)
BINOBJDIR	=	$(OBJDIR)/bin
LIBOBJDIR	=	$(OBJDIR)/lib
BINOBJFILES	=	$(addprefix $(BINOBJDIR)/,$(BINOBJ))
EXAMPLEOBJFILES	=	$(addprefix $(BINOBJDIR)/,$(EXAMPLEOBJ))
LIBOBJFILES	=	$(addprefix $(LIBOBJDIR)/,$(LIBOBJ))
BINSRC		=	$(addprefix $(SRCDIR)/,$(BINOBJ:.o=.cpp))
EXAMPLESRC	=	$(addprefix $(SRCDIR)/,$(EXAMPLEOBJ:.o=.cpp))
LIBSRC		=	$(addprefix $(SRCDIR)/,$(LIBOBJ:.o=.cpp))
LIBSRCHEADER	=	$(addprefix $(SRCDIR)/,$(LIBHEADER))
ALLSRC		=	$(BINSRC) $(EXAMPLESRC) $(LIBSRC) $(LIBSRCHEADER)

#-----------------------------------------------------------------------------
# External Libraries
#-----------------------------------------------------------------------------

GMPDEP	:=	$(SRCDIR)/depend.gmp
GMPSRC	:=	$(shell cat $(GMPDEP))
ifeq ($(GMP),true)
ifeq ($(LEGACY),false)
CPPFLAGS	+=	-DSOPLEX_WITH_GMP
LDFLAGS		+=	-lgmp # todo: move this as GMP_LDFLAGS to submakefiles
endif
endif

ZLIBDEP		:=	$(SRCDIR)/depend.zlib
ZLIBSRC		:=	$(shell cat $(ZLIBDEP))
ifeq ($(ZLIB_LDFLAGS),)
ZLIB		=	false
endif
ifeq ($(ZLIB),true)
CPPFLAGS	+=	-DSOPLEX_WITH_ZLIB $(ZLIB_FLAGS)
LDFLAGS		+=	$(ZLIB_LDFLAGS)
endif

EGLIBDEP	:=	$(SRCDIR)/depend.eglib
EGLIBSRC	:=	$(shell cat $(EGLIBDEP))
ifeq ($(EGLIB),true)
ifeq ($(LEGACY),false)
CPPFLAGS	+=	-DSOPLEX_WITH_EGLIB -I$(LIBDIR)/eglib.$(OSTYPE).$(ARCH).$(COMP)/include
LDFLAGS		+=	$(LIBDIR)/eglib.$(OSTYPE).$(ARCH).$(COMP)/lib/EGlib.a
SOFTLINKS	+=	$(LIBDIR)/eglib.$(OSTYPE).$(ARCH).$(COMP)
LINKSINFO	+=	"\n  -> \"eglib.$(OSTYPE).$(ARCH).$(COMP)\" is a directory containing the EGlib installation, i.e., \"eglib.$(OSTYPE).$(ARCH).$(COMP)/include/EGlib.h\" and \"eglib.$(OSTYPE).$(ARCH).$(COMP)/lib/EGlib.a\" should exist.\n"
endif
endif

#-----------------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------------

ifeq ($(VERBOSE),false)
.SILENT:	$(LIBLINK) $(LIBSHORTLINK) $(BINLINK) $(BINSHORTLINK) $(BINFILE) example $(EXAMPLEOBJFILES) $(LIBFILE) $(BINOBJFILES) $(LIBOBJFILES)
MAKE		+= -s
endif

.PHONY: all
all:		makelibfile
		@-$(MAKE) $(BINFILE) $(LIBLINK) $(LIBSHORTLINK) $(BINLINK) $(BINSHORTLINK)

.PHONY: preprocess
preprocess:	checkdefines
ifneq ($(SOFTLINKS),)
		@$(SHELL) -ec 'if test ! -e $(LINKSMARKERFILE) ; \
			then \
				echo "-> generating necessary links" ; \
				$(MAKE) -j1 $(LINKSMARKERFILE) ; \
			fi'
endif
		@-$(MAKE) touchexternal

$(LIBLINK) $(LIBSHORTLINK):	$(LIBFILE)
		@rm -f $@
		cd $(dir $@) && $(LN_s) $(notdir $(LIBFILE)) $(notdir $@)

$(BINLINK) $(BINSHORTLINK):	$(BINFILE)
		@rm -f $@
		cd $(dir $@) && $(LN_s) $(notdir $(BINFILE)) $(notdir $@)

$(BINFILE):	$(LIBOBJFILES) $(BINOBJFILES) | $(BINDIR) $(BINOBJDIR)
		@echo "-> linking $@"
		-$(LINKCXX) $(BINOBJFILES) $(LIBOBJFILES) \
		$(LDFLAGS) $(LINKCXX_o)$@ \
		|| ($(MAKE) errorhints && false)

.PHONY: example
example:	$(LIBOBJFILES) $(EXAMPLEOBJFILES) | $(BINDIR) $(EXAMPLEOBJDIR)
		@echo "-> linking $(EXAMPLEFILE)"
		-$(LINKCXX) $(EXAMPLEOBJFILES) $(LIBOBJFILES) \
		$(LDFLAGS) $(LINKCXX_o)$(EXAMPLEFILE) \
		|| ($(MAKE) errorhints && false)

.PHONY: makelibfile
makelibfile:	preprocess
		@-$(MAKE) $(LIBFILE)

$(LIBFILE):	$(LIBOBJFILES) | $(LIBDIR) $(LIBOBJDIR)
		@echo "-> generating library $@"
		-rm -f $(LIBFILE)
		$(LIBBUILD) $(LIBBUILDFLAGS) $(LIBBUILD_o)$@ $(LIBOBJFILES) $(REPOSIT)
ifneq ($(RANLIB),)
		$(RANLIB) $@
endif

# include target to detect the current git hash
-include make/local/make.detectgithash

# this empty target is needed for the SoPlex release versions
githash::	# do not remove the double-colon

# include local targets 
-include make/local/make.targets

# include install targets
-include make/make.install

.PHONY: lint
lint:		$(BINSRC) $(LIBSRC)
		-rm -f lint.out
ifeq ($(FILES),)
		$(LINT) lint/$(NAME).lnt +os\(lint.out\) -u -zero -Isrc -I/usr/include -e322 -UNDEBUG $^
else
		$(LINT) lint/$(NAME).lnt +os\(lint.out\) -u -zero -Isrc -I/usr/include -e322 -UNDEBUG $(FILES)
endif

.PHONY: doc
doc:		
		cd doc; $(DOXY) $(NAME).dxy

.PHONY: test
test:		#$(BINFILE)
ifeq ($(LEGACY),false)
		cd check; ./test.sh $(TEST) ../$(BINFILE) $(SETTINGS) $(TIME) $(RESDIR)
endif
ifeq ($(LEGACY),true)
		cd check; ./check_legacy.sh $(TEST).test ../$(BINFILE) '$(ALGO)' $(LIMIT)
endif
.PHONY: check
check:	#$(BINFILE)
ifeq ($(LEGACY),false)
		cd check; ./check.sh ../$(BINFILE) $(RESDIR)
endif
ifeq ($(LEGACY),true)
		cd check; ./check_legacy.sh $(TEST).test ../$(BINFILE) '$(ALGO)' $(LIMIT)
endif

valgrind-check:	$(BINFILE)
		cd check; \
		./valgrind.sh $(TEST).test ../$(BINFILE) $(LIMIT) \
		"$(VALGRIND) $(VFLAGS)" $(VSUPPNAME)

memory_exception_test: $(BINFILE)
		cd check; \
		./exception.sh $(TEST).test ../$(BINFILE) $(LIMIT) \
		"$(VALGRIND) $(VFLAGS)" $(VSUPPNAME)

.PHONY: cleanbin
cleanbin:	| $(BINDIR)
		@echo "remove binary $(BINFILE)"
		@-rm -f $(BINFILE) $(BINLINK) $(BINSHORTLINK)

.PHONY: cleanlib
cleanlib:	| $(LIBDIR)
		@echo "remove library $(LIBFILE)" 
		@-rm -f $(LIBFILE) $(LIBLINK) $(LIBSHORTLINK)

.PHONY: clean
clean:          cleanlib cleanbin | $(LIBOBJDIR) $(BINOBJDIR) $(OBJDIR)
		@echo "remove objective files" 
ifneq ($(LIBOBJDIR),)
		@-rm -f $(LIBOBJDIR)/*.o && rmdir $(LIBOBJDIR)
endif
ifneq ($(BINOBJDIR),)
		@-rm -f $(BINOBJDIR)/*.o && rmdir $(BINOBJDIR)
endif
ifneq ($(OBJDIR),)
		@-rm -f $(LASTSETTINGS)
		@-rmdir $(OBJDIR)
endif
		@-rm -f $(EXAMPLEFILE)

vimtags:
		-ctags -o TAGS src/*.cpp src/*.h

etags:
		-ctags -e -o TAGS src/*.cpp src/*.h

$(OBJDIR):	
		@-mkdir -p $(OBJDIR)

$(BINOBJDIR):	| $(OBJDIR)
		@-mkdir -p $(BINOBJDIR)

$(LIBOBJDIR):	| $(OBJDIR)
		@-mkdir -p $(LIBOBJDIR)

$(BINDIR):
		@-mkdir -p $(BINDIR)

$(LIBDIR):
		@-mkdir -p $(LIBDIR)

.PHONY: depend
depend:
		$(SHELL) -ec '$(DCXX) $(DFLAGS) $(CPPFLAGS) \
		$(BINSRC:.o=.cpp) \
		| sed '\''s|^\([0-9A-Za-z_]\{1,\}\)\.o|$$\(BINOBJDIR\)/\1.o|g'\'' \
		>$(DEPEND)'
		$(SHELL) -ec '$(DCXX) $(DFLAGS) $(CPPFLAGS) \
		$(EXAMPLESRC:.o=.cpp) \
		| sed '\''s|^\([0-9A-Za-z_]\{1,\}\)\.o|$$\(BINOBJDIR\)/\1.o|g'\'' \
		>>$(DEPEND)'
		$(SHELL) -ec '$(DCXX) $(DFLAGS) $(CPPFLAGS) \
		$(LIBSRC:.o=.cpp) \
		| sed '\''s|^\([0-9A-Za-z_]\{1,\}\)\.o|$$\(LIBOBJDIR\)/\1.o|g'\'' \
		>>$(DEPEND)'
		@echo `grep -l "SOPLEX_WITH_GMP" $(ALLSRC)` >$(GMPDEP)
		@echo `grep -l "SOPLEX_WITH_ZLIB" $(ALLSRC)` >$(ZLIBDEP)
		@echo `grep -l "SOPLEX_WITH_EGLIB" $(ALLSRC)` >$(EGLIBDEP)
		@echo `grep -l "SOPLEX_LEGACY" $(ALLSRC)` >$(LEGACYDEP)
		@echo `grep -l "DISABLE_VERBOSITY" $(ALLSRC)` >$(PARASCIPDEP)

-include	$(DEPEND)

$(BINOBJDIR)/%.o:	$(SRCDIR)/%.cpp
		@-mkdir -p $(BINOBJDIR)
		@echo "-> compiling $@"
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(BINOFLAGS) $(CXX_c)$< $(CXX_o)$@

$(LIBOBJDIR)/%.o:	$(SRCDIR)/%.cpp
		@-mkdir -p $(LIBOBJDIR)
		@echo "-> compiling $@"
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LIBOFLAGS) $(CXX_c)$< $(CXX_o)$@


-include $(LASTSETTINGS)

.PHONY: touchexternal
touchexternal:	$(GMPDEP) $(ZLIBDEP) $(EGLIBDEP) $(PARASCIPDEP) $(LEGACYDEP) | $(OBJDIR)
ifneq ($(SPXGITHASH),$(LAST_SPXGITHASH))
		@-$(MAKE) githash
endif
		@$(SHELL) -ec 'if test ! -e $(SPXGITHASHFILE) ; \
			then \
				echo "-> generating $(SPXGITHASHFILE)" ; \
				$(MAKE) githash ; \
			fi'
ifneq ($(GMP),$(LAST_GMP))
		@-touch $(GMPSRC)
endif
ifneq ($(ZLIB),$(LAST_ZLIB))
		@-touch $(ZLIBSRC)
endif
ifneq ($(EGLIB),$(LAST_EGLIB))
		@-touch $(EGLIBSRC)
endif
ifneq ($(PARASCIP),$(LAST_PARASCIP))
		@-touch $(PARASCIPSRC)
endif
ifneq ($(LEGACY),$(LAST_LEGACY))
		@-touch $(LEGACYSRC)
endif
ifneq ($(SHARED),$(LAST_SHARED))
		@-touch $(LIBSRC)
		@-touch $(BINSRC)
endif
ifneq ($(USRCXXFLAGS),$(LAST_USRCXXFLAGS))
		@-touch $(LIBSRC)
		@-touch $(BINSRC)
endif
ifneq ($(USRCPPFLAGS),$(LAST_USRCPPFLAGS))
		@-touch $(LIBSRC)
		@-touch $(BINSRC)
endif
ifneq ($(USRLDFLAGS),$(LAST_USRLDFLAGS))
		@-touch -c $(EXAMPLEOBJFILES) $(BINOBJFILES) $(LIBOBJFILES)
endif
ifneq ($(USRARFLAGS),$(LAST_USRARFLAGS))
		@-touch -c $(EXAMPLEOBJFILES) $(BINOBJFILES) $(LIBOBJFILES)
endif
		@-rm -f $(LASTSETTINGS)
		@echo "LAST_SPXGITHASH=$(SPXGITHASH)" >> $(LASTSETTINGS)
		@echo "LAST_GMP=$(GMP)" >> $(LASTSETTINGS)
		@echo "LAST_ZLIB=$(ZLIB)" >> $(LASTSETTINGS)
		@echo "LAST_EGLIB=$(EGLIB)" >> $(LASTSETTINGS)
		@echo "LAST_PARASCIP=$(PARASCIP)" >> $(LASTSETTINGS)
		@echo "LAST_LEGACY=$(LEGACY)" >> $(LASTSETTINGS)
		@echo "LAST_SHARED=$(SHARED)" >> $(LASTSETTINGS)
		@echo "LAST_USRCXXFLAGS=$(USRCXXFLAGS)" >> $(LASTSETTINGS)
		@echo "LAST_USRCPPFLAGS=$(USRCPPFLAGS)" >> $(LASTSETTINGS)
		@echo "LAST_USRLDFLAGS=$(USRLDFLAGS)" >> $(LASTSETTINGS)
		@echo "LAST_USRARFLAGS=$(USRARFLAGS)" >> $(LASTSETTINGS)
		@echo "LAST_USRDFLAGS=$(USRDFLAGS)" >> $(LASTSETTINGS)

$(LINKSMARKERFILE):
		@$(MAKE) links

.PHONY: links
links:		| $(LIBDIR) echosoftlinks $(SOFTLINKS)
		@rm -f $(LINKSMARKERFILE)
		@echo "this is only a marker" > $(LINKSMARKERFILE)

.PHONY: echosoftlinks
echosoftlinks:
		@echo
		@echo "- Current settings: OSTYPE=$(OSTYPE) ARCH=$(ARCH) COMP=$(COMP) SUFFIX=$(LINKLIBSUFFIX) EGLIB=$(EGLIB)"
		@echo
		@echo "* SoPlex needs some softlinks to external programs."
		@echo "* Please insert the paths to the corresponding directories/libraries below."
		@echo "* The links will be installed in the 'lib' directory."
		@echo "* For more information and if you experience problems see the INSTALL file."
		@echo
		@echo -e $(LINKSINFO)

.PHONY: $(SOFTLINKS)
$(SOFTLINKS):
ifeq ($(MAKESOFTLINKS), true)
		@$(SHELL) -ec 'if test ! -e $@ ; \
			then \
				DIRNAME=`dirname $@` ; \
				BASENAMEA=`basename $@ .$(STATICLIBEXT)` ; \
				BASENAMESO=`basename $@ .$(SHAREDLIBEXT)` ; \
				echo ; \
				echo "- preparing missing soft-link \"$@\":" ; \
				if test -e $$DIRNAME/$$BASENAMEA.$(SHAREDLIBEXT) ; \
				then \
					echo "* this soft-link is not necessarily needed since \"$$DIRNAME/$$BASENAMEA.$(SHAREDLIBEXT)\" already exists - press return to skip" ; \
				fi ; \
				if test -e $$DIRNAME/$$BASENAMESO.$(STATICLIBEXT) ; \
				then \
					echo "* this soft-link is not necessarily needed since \"$$DIRNAME/$$BASENAMESO.$(STATICLIBEXT)\" already exists - press return to skip" ; \
				fi ; \
				echo "> Enter soft-link target file or directory for \"$@\" (return if not needed): " ; \
				echo -n "> " ; \
				cd $$DIRNAME ; \
				eval $(READ) TARGET ; \
				cd $(SPXDIR) ; \
				if test "$$TARGET" != "" ; \
				then \
					echo "-> creating softlink \"$@\" -> \"$$TARGET\"" ; \
					rm -f $@ ; \
					$(LN_s) $$TARGET $@ ; \
				else \
					echo "* skipped creation of softlink \"$@\". Call \"make links\" if needed later." ; \
				fi ; \
				echo ; \
			fi'
endif

.PHONY: checkdefines
checkdefines:
ifneq ($(GMP),true)
ifneq ($(GMP),false)
		$(error invalid GMP flag selected: GMP=$(GMP). Possible options are: true false)
endif
endif
ifneq ($(ZLIB),true)
ifneq ($(ZLIB),false)
		$(error invalid ZLIB flag selected: ZLIB=$(ZLIB). Possible options are: true false)
endif
endif
ifneq ($(EGLIB),true)
ifneq ($(EGLIB),false)
		$(error invalid EGLIB flag selected: EGLIB=$(EGLIB). Possible options are: true false)
endif
endif
ifneq ($(PARASCIP),true)
ifneq ($(PARASCIP),false)
		$(error invalid PARASCIP flag selected: PARASCIP=$(PARASCIP). Possible options are: true false)
endif
endif
ifneq ($(LEGACY),true)
ifneq ($(LEGACY),false)
		$(error invalid LEGACY flag selected: LEGACY=$(LEGACY). Possible options are: true false)
endif
endif

.PHONY: errorhints
errorhints:
ifeq ($(ZLIB),true)
		@echo "build failed with ZLIB=true: if ZLIB is not available, try building with ZLIB=false"
endif
ifeq ($(GMP),true)
		@echo "build failed with GMP=true: if GMP is not available, try building with GMP=false"
endif
ifeq ($(EGLIB),true)
		@echo "build failed with EGLIB=true: if EGlib is not available, try building with EGLIB=false"
endif

# --- EOF ---------------------------------------------------------------------
# DO NOT DELETE
