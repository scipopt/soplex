# $Id: Makefile,v 1.34 2002/03/01 13:15:29 bzfpfend Exp $
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*   File....: Makefile                                                      *
#*   Name....: SoPlex Makefile                                               *
#*   Author..: Thorsten Koch                                                 *
#*   Copyright by Author, All rights reserved                                *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

.PHONY:		depend clean distclean lint doc lib
.PHONY:		tests quick check mittel infeas cover

ARCH            :=      $(shell uname -m | \
                        sed \
			-e s/sun../sparc/ \
			-e s/i.86/x86/ \
			-e s/IP../mips/ \
			-e s/9000..../hppa/)
OSTYPE		:=	$(shell uname -s | tr A-Z a-z)

OPT		=	opt

COMP		=	gnu
CXX		=	g++
DCXX		=	g++
LINT		=	flexelint
AR		=	ar
RANLIB		=	ranlib
DOXY		=	doxygen

CPPFLAGS	=	-Isrc
CXXFLAGS	=	-O
LDFLAGS		=	-lm
ARFLAGS		=	cr
DFLAGS		=	-MM

SRCDIR		=	src
BINDIR		=	bin
LIBDIR		=	lib
NAME		=	soplex
FLAGS		=       #
LIBOBJ		= 	changesoplex.o didxset.o \
			dsvector.o dvector.o enter.o factor.o \
			forest.o idxset.o leave.o lpcolset.o lprowset.o \
			lprow.o message.o mpsinput.o nameset.o real.o \
			slufactor.o solve.o soplex.o spxaggregatesm.o \
			spxbasis.o spxbounds.o spxchangebasis.o \
			spxdefaultpr.o spxdefaultrt.o spxdesc.o spxdevexpr.o \
			spxfastrt.o spxgeneralsm.o \
			spxharrisrt.o spxhybridpr.o spxio.o spxlp.o \
			spxlpfread.o spxmpsread.o \
			spxparmultpr.o spxredundantsm.o spxrem1sm.o \
			spxscale.o spxshift.o spxsolver.o spxsolve.o \
			spxstarter.o spxsteeppr.o spxsumst.o spxvecs.o \
			spxvectorst.o spxweightpr.o spxweightst.o \
			ssvector.o subsvector.o svector.o \
			svset.o timer.o unitvector.o update.o updatevector.o \
			vector.o vsolve.o

OBJECT		=	example.o

#------------------------------------------------------------------------------
#--- NOTHING TO CHANGE FROM HERE ON -------------------------------------------
#------------------------------------------------------------------------------

GCCWARN		=	-Wall -W -Wpointer-arith -Wbad-function-cast \
			-Wcast-align -Wwrite-strings -Wconversion \
			-Wstrict-prototypes -Wmissing-prototypes \
			-Wmissing-declarations -Wno-unknown-pragmas \
			-Wctor-dtor-privacy -Wnon-virtual-dtor -Wreorder \
			-Woverloaded-virtual -Wsign-promo -Wsynth -Wundef \
			-Wcast-qual -Wold-style-cast -Wshadow 
#			-Weffc++ -Wredundant-decls    

#GCCWARN =
#-----------------------------------------------------------------------------
include make/make.$(OSTYPE).$(ARCH).$(COMP).$(OPT)
#-----------------------------------------------------------------------------

TARGET		=	$(NAME).$(OSTYPE).$(ARCH).$(COMP).$(OPT)
LIBRARY		=	$(LIBDIR)/lib$(NAME).$(OSTYPE).$(ARCH).$(COMP).$(OPT).a
BINARY		=	$(BINDIR)/$(TARGET)
DEPEND		=	src/depend

OBJDIR		=	obj/O.$(OSTYPE).$(ARCH).$(COMP).$(OPT)
OBJXXX		=	$(addprefix $(OBJDIR)/,$(OBJECT))
LIBXXX		=	$(addprefix $(OBJDIR)/,$(LIBOBJ))
OBJSRC		=	$(addprefix $(SRCDIR)/,$(OBJECT:.o=.cpp))
LIBSRC		=	$(addprefix $(SRCDIR)/,$(LIBOBJ:.o=.cpp))

$(BINARY):	_$(OBJDIR) _$(BINDIR) $(OBJXXX) $(LIBRARY) 
		$(CXX) $(CXXFLAGS) $(OBJXXX) \
		-L$(LIBDIR) -l$(TARGET) $(LDFLAGS) -o $@

$(LIBRARY):	_$(LIBDIR) $(LIBXXX) 
		-rm -f $(LIBRARY)
		$(AR) $(ARFLAGS) $@ $(LIBXXX) 
		$(RANLIB) $@

lint:		$(OBJSRC) $(LIBSRC)
		$(LINT) lint/soplex.lnt -os\(lint.out\) \
		$(CPPFLAGS) -UNDEBUG $^

doc:		
		cd doc; $(DOXY) soplex.dxy

lib:		$(LIBRARY)

tests:		check infeas mittel

check:		
		cd check; ./check.sh netlib.test $(BINARY)

quick:		
		cd check; ./check.sh quick.test $(BINARY)

mittel:		
		cd check; ./check.sh mittelmann.test $(BINARY)

infeas:		
		cd check; ./check.sh infeas.test $(BINARY)

cover:
		cd check; ./cover.sh cover.test $(BINARY)

clean:
		-rm -rf $(OBJDIR)/* $(LIBRARY) $(BINARY)

distclean:
		-rm -rf obj/* lib/libsoplex.* bin/soplex.* 

_$(OBJDIR):	
		-mkdir -p $(OBJDIR)

_$(LIBDIR):
		-mkdir -p $(LIBDIR)

_$(BINDIR):
		-mkdir -p $(BINDIR)

depend:
		$(SHELL) -ec '$(DCXX) $(DFLAGS) $(CPPFLAGS) \
		$(OBJSRC:.o=.cpp) $(LIBSRC:.o=.cpp) \
		| sed '\''s|^\([0-9A-z]\{1,\}\)\.o|$$\(OBJDIR\)/\1.o|g'\'' \
		>$(DEPEND)'

include		$(DEPEND)

$(OBJDIR)/%.o:	$(SRCDIR)/%.cpp
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# --- EOF ---------------------------------------------------------------------


