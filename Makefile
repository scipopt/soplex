# $Id: Makefile,v 1.67 2005/12/12 19:16:58 bzforlow Exp $
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*   File....: Makefile                                                      *
#*   Name....: SoPlex Makefile                                               *
#*   Author..: Thorsten Koch                                                 *
#*   Copyright by Author, All rights reserved                                *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

.PHONY:		depend clean distclean lint doc lib check all

ARCH            :=      $(shell uname -m | \
                        sed \
			-e s/sun../sparc/ \
			-e s/i.86/x86/ \
			-e s/IP../mips/ \
			-e s/9000..../hppa/ \
			-e s/Power\ Macintosh/ppc/ \
			-e s/00........../pwr4/)
OSTYPE		:=	$(shell uname -s | tr '[:upper:]' '[:lower:]' | sed -e s/cygwin.*/cygwin/ -e s/irix../irix/ )
HOSTNAME	:=	$(shell uname -n | tr '[:upper:]' '[:lower:]')

VERBOSE		=	false
OPT		=	opt
LINK		=	static
LIBEXT		=	a
TEST		=	quick
ALGO		=	1 2 3 4 5 6
LIMIT		=	#

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
VFLAGS		=	--tool=memcheck --leak-check=yes --show-reachable=yes

SRCDIR		=	src
BINDIR		=	bin
LIBDIR		=	lib
NAME		=	soplex
FLAGS		=       #
LIBOBJ		= 	changesoplex.o didxset.o \
			dsvector.o dvector.o enter.o factor.o \
			forest.o idxset.o leave.o lpcolset.o lprowset.o \
			lprow.o message.o mpsinput.o nameset.o \
			slufactor.o solve.o soplex.o spxaggregatesm.o \
			spxbasis.o spxbounds.o spxchangebasis.o \
			spxequilisc.o spxdantzigpr.o spxdefaultrt.o \
			spxdefines.o spxdesc.o spxdevexpr.o \
			spxfastrt.o spxfileio.o spxgeneralsm.o spxgeometsc.o \
			spxharrisrt.o spxhybridpr.o spxid.o spxintervalsm.o spxio.o \
			spxlp.o spxlpfread.o spxmpsread.o spxmpswrite.o \
			spxout.o spxparmultpr.o spxquality.o spxredundantsm.o \
			spxscaler.o spxshift.o spxsolver.o spxsolve.o \
			spxstarter.o spxsteeppr.o spxsumst.o spxvecs.o \
			spxvectorst.o spxweightpr.o spxweightst.o \
			ssvector.o svector.o svset.o timer.o \
			tracemethod.o unitvector.o update.o updatevector.o \
			vector.o vsolve.o \
			gzstream.o
BINOBJ		=	example.o
CHANGEBINOBJ	=	exercise_LP_changes.o
REPOSIT		=	# template repository, explicitly empty  #spxproof.o 

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
include make/make.$(OSTYPE).$(ARCH).$(COMP).$(OPT).$(LINK)
-include make/local/make.$(HOSTNAME)
#-----------------------------------------------------------------------------

BINNAME		=	$(NAME).$(OSTYPE).$(ARCH).$(COMP).$(OPT).$(LINK)
LIBNAME		=	$(NAME).$(OSTYPE).$(ARCH).$(COMP).$(OPT).$(LINK)
BINFILE		=	$(BINDIR)/$(BINNAME)
CHANGEBINFILE   =	$(BINDIR)/exercise_LP_changes.$(OSTYPE).$(ARCH).$(COMP).$(OPT).$(LINK)
LIBFILE		=	$(LIBDIR)/lib$(LIBNAME).$(LIBEXT)
DEPEND		=	src/depend

# potential valgrind suppression file name
VSUPPNAME	= 	$(OSTYPE).$(ARCH).$(COMP).supp

OBJDIR		=	obj/O.$(OSTYPE).$(ARCH).$(COMP).$(OPT).$(LINK)
BINOBJDIR	=	$(OBJDIR)/bin
LIBOBJDIR	=	$(OBJDIR)/lib
BINOBJFILES	=	$(addprefix $(BINOBJDIR)/,$(BINOBJ))
CHANGEBINOBJFILES =	$(addprefix $(BINOBJDIR)/,$(CHANGEBINOBJ))
LIBOBJFILES	=	$(addprefix $(LIBOBJDIR)/,$(LIBOBJ))
BINSRC		=	$(addprefix $(SRCDIR)/,$(BINOBJ:.o=.cpp))
LIBSRC		=	$(addprefix $(SRCDIR)/,$(LIBOBJ:.o=.cpp))

$(BINFILE):	_$(BINDIR) _$(BINOBJDIR) $(LIBFILE) $(BINOBJFILES)
		@echo "-> linking $@"
ifeq ($(VERBOSE), true)
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(BINOBJFILES) \
		-L$(LIBDIR) -l$(LIBNAME) $(LDFLAGS) -o $@
else
		@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(BINOBJFILES) \
		-L$(LIBDIR) -l$(LIBNAME) $(LDFLAGS) -o $@
endif

$(LIBFILE):	_$(LIBDIR) _$(LIBOBJDIR) $(LIBOBJFILES) 
		@echo "-> generating library $@"
ifeq ($(VERBOSE), true)
		-rm -f $(LIBFILE)
		$(AR) $(ARFLAGS) $@ $(LIBOBJFILES) $(REPOSIT)
		$(RANLIB) $@
else
		@-rm -f $(LIBFILE)
		@$(AR) $(ARFLAGS) $@ $(LIBOBJFILES) $(REPOSIT)
		@$(RANLIB) $@
endif

$(CHANGEBINFILE): _$(BINDIR) _$(BINOBJDIR) $(LIBFILE) $(CHANGEBINOBJFILES)
		@echo "-> linking $@"
ifeq ($(VERBOSE), true)
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(CHANGEBINOBJFILES) \
		-L$(LIBDIR) -l$(LIBNAME) $(LDFLAGS) -o $@
else
		@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(CHANGEBINOBJFILES) \
		-L$(LIBDIR) -l$(LIBNAME) $(LDFLAGS) -o $@
endif

lint:		$(BINSRC) $(LIBSRC)
		$(LINT) lint/soplex.lnt -os\(lint.out\) \
		$(CPPFLAGS) -UNDEBUG $^

doc:		
		cd doc; $(DOXY) soplex.dxy

lib:		$(LIBFILE)

change_exerciser: $(CHANGEBINFILE)

all:		$(BINFILE) $(CHANGEBINFILE)

check:		$(BINFILE)
		cd check; ./check.sh $(TEST).test ../$(BINFILE) '$(ALGO)' $(LIMIT)

valgrind-check:	$(BINFILE)
		cd check; \
		./valgrind.sh $(TEST).test ../$(BINFILE) '$(ALGO)' '$(LIMIT)' \
		"$(VALGRIND) $(VFLAGS)" $(VSUPPNAME)

clean:
		-rm -rf $(OBJDIR)/* $(BINFILE) $(LIBFILE)

distclean:
		-rm -rf obj/* lib/libsoplex.* bin/soplex.* 

vimtags:
		-ctags -o TAGS src/*.cpp src/*.h

etags:
		-ctags -e -o TAGS src/*.cpp src/*.h

_$(OBJDIR):	
		@-mkdir -p $(OBJDIR)

_$(BINOBJDIR):	_$(OBJDIR)
		@-mkdir -p $(BINOBJDIR)

_$(LIBOBJDIR):	_$(OBJDIR)
		@-mkdir -p $(LIBOBJDIR)

_$(BINDIR):
		@-mkdir -p $(BINDIR)

_$(LIBDIR):
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

-include	$(DEPEND)

$(BINOBJDIR)/%.o:	$(SRCDIR)/%.cpp
		@echo "-> compiling $@"
ifeq ($(VERBOSE), true)
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(BINOFLAGS) -c $< -o $@
else
		@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(BINOFLAGS) -c $< -o $@
endif

$(LIBOBJDIR)/%.o:	$(SRCDIR)/%.cpp
		@echo "-> compiling $@"
ifeq ($(VERBOSE), true)
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LIBOFLAGS) -c $< -o $@
else
		@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LIBOFLAGS) -c $< -o $@
endif

# --- EOF ---------------------------------------------------------------------
