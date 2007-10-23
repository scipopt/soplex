# $Id: Makefile,v 1.85 2007/10/23 09:21:42 bzfberth Exp $
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*   File....: Makefile                                                      *
#*   Name....: SoPlex Makefile                                               *
#*   Author..: Thorsten Koch                                                 *
#*   Copyright by Author, All rights reserved                                *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

.PHONY:		all depend clean distclean lint doc check

VERSION		:=	1.3.3

ARCH            :=      $(shell uname -m | \
                        sed \
			-e s/sun../sparc/ \
			-e s/i.86/x86/ \
	                -e s/i86pc/x86/ \
			-e s/IP../mips/ \
			-e s/9000..../hppa/ \
			-e s/Power\ Macintosh/ppc/ \
			-e s/00........../pwr4/)
OSTYPE		:=	$(shell uname -s | tr '[:upper:]' '[:lower:]' | sed -e s/cygwin.*/cygwin/ -e s/irix../irix/ )
HOSTNAME	:=	$(shell uname -n | tr '[:upper:]' '[:lower:]')

VERBOSE		=	false
OPT		=	opt
LIBEXT		=	a
TEST		=	quick
ALGO		=	1 2 3 4 5 6
LIMIT		=	#

ZLIB		=	true

COMP		=	gnu
CXX		=	g++
DCXX		=	g++
LINT		=	flexelint
AR		=	ar
RANLIB		=	ranlib
DOXY		=	doxygen
VALGRIND	=	valgrind

CPPFLAGS	=	-Isrc
CXXFLAGS	=	-O
BINOFLAGS	=	
LIBOFLAGS	=	
LDFLAGS		=	-lm
ARFLAGS		=	cr
DFLAGS		=	-MM
VFLAGS		=	--tool=memcheck --leak-check=yes --show-reachable=yes #--gen-suppressions=yes

SRCDIR		=	src
BINDIR		=	bin
LIBDIR		=	lib
NAME		=	soplex
FLAGS		=       #
LIBOBJ		= 	changesoplex.o didxset.o \
			dsvector.o dvector.o enter.o factor.o \
			forest.o idxset.o leave.o lpcolset.o lprowset.o \
			lprow.o message.o mpsinput.o nameset.o \
			slufactor.o solve.o soplex.o \
			spxbasis.o spxbounds.o spxchangebasis.o \
			spxequilisc.o spxdantzigpr.o spxdefaultrt.o \
			spxdefines.o spxdesc.o spxdevexpr.o \
			spxfastrt.o spxfileio.o spxgeometsc.o \
			spxharrisrt.o spxhybridpr.o spxid.o spxio.o \
			spxlp.o spxlpfread.o spxmainsm.o spxmpsread.o \
			spxmpswrite.o spxlpfwrite.o \
			spxout.o spxparmultpr.o spxquality.o \
			spxscaler.o spxshift.o spxsolver.o spxsolve.o \
			spxstarter.o spxsteeppr.o spxsumst.o spxvecs.o \
			spxvectorst.o spxweightpr.o spxweightst.o \
			ssvector.o svector.o svset.o timer.o \
			tracemethod.o unitvector.o update.o updatevector.o \
			vector.o vsolve.o \
			gzstream.o
BINOBJ		=	example.o
CHANGEBINOBJ	=	exercise_LP_changes.o
EXCEPTIONBINOBJ	=	status_exception_test.o
REPOSIT		=	# template repository, explicitly empty  #spxproof.o 

BASE		=	$(OSTYPE).$(ARCH).$(COMP).$(OPT)

LASTSETTINGS	=	$(OBJDIR)/make.lastsettings

#------------------------------------------------------------------------------
#--- NOTHING TO CHANGE FROM HERE ON -------------------------------------------
#------------------------------------------------------------------------------

GCCWARN		=	-Wall -W -Wpointer-arith -Wno-unknown-pragmas \
			-Wcast-align -Wwrite-strings -Wconversion \
			-Wctor-dtor-privacy -Wnon-virtual-dtor -Wreorder \
			-Woverloaded-virtual -Wsign-promo -Wsynth -Wundef \
			-Wcast-qual -Wold-style-cast -Wshadow 
#			-Weffc++ -Wredundant-decls    
# gcc 2.xx -Wmissing-declarations -Wbad-function-cast 

#GCCWARN =
#-----------------------------------------------------------------------------
include make/make.$(BASE)
-include make/local/make.$(HOSTNAME)
-include make/local/make.$(HOSTNAME).$(COMP)
-include make/local/make.$(HOSTNAME).$(COMP).$(OPT)
#-----------------------------------------------------------------------------

BINNAME		=	$(NAME)-$(VERSION).$(BASE)
LIBNAME		=	$(NAME)-$(VERSION).$(BASE)
BINFILE		=	$(BINDIR)/$(BINNAME)
CHANGEBINFILE   =	$(BINDIR)/exercise_LP_changes.$(BASE)
EXCEPTIONBINFILE=	$(BINDIR)/status_exception_test.$(BASE)
TESTEXBINFILE	=	$(BINDIR)/mem_exception_test.$(BASE)
LIBFILE		=	$(LIBDIR)/lib$(LIBNAME).$(LIBEXT)
LIBLINK		=	$(LIBDIR)/lib$(NAME).$(BASE).$(LIBEXT)
BINLINK		=	$(BINDIR)/$(NAME).$(BASE)
BINSHORTLINK	=	$(BINDIR)/$(NAME)
DEPEND		=	src/depend

# potential valgrind suppression file name
VSUPPNAME	= 	$(OSTYPE).$(ARCH).$(COMP).supp

OBJDIR		=	obj/O.$(BASE)
BINOBJDIR	=	$(OBJDIR)/bin
LIBOBJDIR	=	$(OBJDIR)/lib
BINOBJFILES	=	$(addprefix $(BINOBJDIR)/,$(BINOBJ))
TESTEXBINOBJFILES	=	$(addprefix $(BINOBJDIR)/,$(BINOBJ))
CHANGEBINOBJFILES =	$(addprefix $(BINOBJDIR)/,$(CHANGEBINOBJ))
EXCEPTIONBINOBJFILES =	$(addprefix $(BINOBJDIR)/,$(EXCEPTIONBINOBJ))
LIBOBJFILES	=	$(addprefix $(LIBOBJDIR)/,$(LIBOBJ))
BINSRC		=	$(addprefix $(SRCDIR)/,$(BINOBJ:.o=.cpp))
LIBSRC		=	$(addprefix $(SRCDIR)/,$(LIBOBJ:.o=.cpp))

ZLIBDEP		:=	$(SRCDIR)/depend.zlib
ZLIBSRC		:=	$(shell cat $(ZLIBDEP))
ifeq ($(ZLIB_LDFLAGS),)
ZLIB		=	false
endif
ifeq ($(ZLIB),true)
CPPFLAGS	+=	-DWITH_ZLIB $(ZLIB_FLAGS)
LDFLAGS		+=	$(ZLIB_LDFLAGS)
endif


ifeq ($(VERBOSE),false)
.SILENT:	$(LIBLINK) $(BINLINK) $(BINSHORTLINK) $(BINFILE) $(LIBFILE) $(CHANGEBINFILE) $(BINOBJFILES) $(LIBOBJFILES)
endif

all:		$(LIBFILE) $(BINFILE) $(LIBLINK) $(BINLINK) $(BINSHORTLINK)

$(LIBLINK):	$(LIBFILE)
		@rm -f $@
		cd $(dir $@) && ln -s $(notdir $(LIBFILE)) $(notdir $@)

$(BINLINK) $(BINSHORTLINK):	$(BINFILE)
		@rm -f $@
		cd $(dir $@) && ln -s $(notdir $(BINFILE)) $(notdir $@)

$(BINFILE):	$(BINDIR) $(BINOBJDIR) $(LIBFILE) $(BINOBJFILES)
		@echo "-> linking $@"
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(BINOBJFILES) \
		-L$(LIBDIR) -l$(LIBNAME) $(LDFLAGS) -o $@

$(LIBFILE):	$(LIBDIR) $(LIBOBJDIR) touchexternal $(LIBOBJFILES) 
		@echo "-> generating library $@"
		-rm -f $(LIBFILE)
		$(AR) $(ARFLAGS) $@ $(LIBOBJFILES) $(REPOSIT)
		$(RANLIB) $@

# build test binaries
$(TESTEXBINFILE):	_$(BINDIR) _$(BINOBJDIR) $(LIBFILE) $(TESTEXBINOBJFILES)
			$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(BINOBJFILES) \
			-L$(LIBDIR) -l$(LIBNAME) $(LDFLAGS) -o $@


$(CHANGEBINFILE): $(BINDIR) $(BINOBJDIR) $(LIBFILE) $(CHANGEBINOBJFILES)
		@echo "-> linking $@"
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(CHANGEBINOBJFILES) \
		-L$(LIBDIR) -l$(LIBNAME) $(LDFLAGS) -o $@

$(EXCEPTIONBINFILE): $(BINDIR) $(BINOBJDIR) $(LIBFILE) $(EXCEPTIONBINOBJFILES)
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(EXCEPTIONBINOBJFILES) \
		-L$(LIBDIR) -l$(LIBNAME) $(LDFLAGS) -o $@

lint:		$(BINSRC) $(LIBSRC)
		$(LINT) lint/$(NAME).lnt -os\(lint.out\) \
		$(CPPFLAGS) -UNDEBUG $^

doc:		
		cd doc; $(DOXY) $(NAME).dxy

change_exerciser: $(CHANGEBINFILE)

status_exception_test: $(EXCEPTIONBINFILE)
		cd bin; \
		../$(EXCEPTIONBINFILE)

all:		$(BINFILE) $(CHANGEBINFILE)

check:		#$(BINFILE)
		cd check; ./check.sh $(TEST).test ../$(BINFILE) '$(ALGO)' $(LIMIT)

valgrind-check:	$(BINFILE)
		cd check; \
		./valgrind.sh $(TEST).test ../$(BINFILE) '$(ALGO)' '$(LIMIT)' \
		"$(VALGRIND) $(VFLAGS)" $(VSUPPNAME)

memory_exception_test: $(BINFILE)
		cd check; \
		./exception.sh $(TEST).test ../$(BINFILE) '$(ALGO)' '$(LIMIT)' \
		"$(VALGRIND) $(VFLAGS)" $(VSUPPNAME)

clean:
		-rm -rf $(OBJDIR)/* $(BINFILE) $(LIBFILE) $(LIBLINK) $(BINLINK) $(BINSHORTLINK)

distclean:	clean
		-rm -rf obj/* $(LIBDIR)/lib$(NAME).* $(BINDIR)/$(NAME).* 

vimtags:
		-ctags -o TAGS src/*.cpp src/*.h

etags:
		-ctags -e -o TAGS src/*.cpp src/*.h

$(OBJDIR):	
		@-mkdir -p $(OBJDIR)

$(BINOBJDIR):	$(OBJDIR)
		@-mkdir -p $(BINOBJDIR)

$(LIBOBJDIR):	$(OBJDIR)
		@-mkdir -p $(LIBOBJDIR)

$(BINDIR):
		@-mkdir -p $(BINDIR)

$(LIBDIR):
		@-mkdir -p $(LIBDIR)

depend:
		$(SHELL) -ec '$(DCXX) $(DFLAGS) $(CPPFLAGS) \
		$(BINSRC:.o=.cpp) \
		| sed '\''s|^\([0-9A-Za-z]\{1,\}\)\.o|$$\(BINOBJDIR\)/\1.o|g'\'' \
		>$(DEPEND)'
		$(SHELL) -ec '$(DCXX) $(DFLAGS) $(CPPFLAGS) \
		$(LIBSRC:.o=.cpp) \
		| sed '\''s|^\([0-9A-Za-z]\{1,\}\)\.o|$$\(LIBOBJDIR\)/\1.o|g'\'' \
		>>$(DEPEND)'
		@echo `grep -l "WITH_ZLIB" $(SRCDIR)/*` >$(ZLIBDEP)

-include	$(DEPEND)

$(BINOBJDIR)/%.o:	$(SRCDIR)/%.cpp
		@echo "-> compiling $@"
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(BINOFLAGS) -c $< -o $@

$(LIBOBJDIR)/%.o:	$(SRCDIR)/%.cpp
		@echo "-> compiling $@"
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LIBOFLAGS) -c $< -o $@


-include $(LASTSETTINGS)

.PHONY: touchexternal
touchexternal:	$(ZLIBDEP)
ifneq ($(ZLIB),$(LAST_ZLIB))
		@-touch $(ZLIBSRC)
endif
		@-rm -f $(LASTSETTINGS)
		@echo "LAST_ZLIB=$(ZLIB)" >> $(LASTSETTINGS)


# --- EOF ---------------------------------------------------------------------
