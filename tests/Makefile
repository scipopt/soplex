#
# General setup.
#

NAME		=	soplex
VERSION		:=	1.4.1

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

OPT		=	opt
LINK		=	static

BASE		=	$(OSTYPE).$(ARCH).$(COMP).$(OPT)

COMP		=	gnu
CXX		=	g++
DCXX		=	g++
LIBEXT		=	a

SRCDIR		=	../src
BINDIR		=	../bin
LIBDIR		=	../lib

CPPFLAGS	=	-I$(SRCDIR)
CXXFLAGS	=	-g
BINOFLAGS	=	
LIBOFLAGS	=	
LDFLAGS		=	-lm -lz -static

LIBNAME		=	$(NAME)-$(VERSION).$(BASE)
LIBFILE		=	$(LIBDIR)/lib$(LIBNAME).$(LIBEXT)

#
# Setup for test binaries.
#
CHANGEBINOBJ	=	exercise_LP_changes.o
CHANGEBINFILE   =	exercise_LP_changes.$(BASE)

EXCEPTIONBINOBJ	=	status_exception_test.o
EXCEPTIONBINFILE=	status_exception_test.$(BASE)

TESTCHANGEELEMENTOBJ	=	testChangeElement.o
TESTCHANGEELEMENTBIN	=	testChangeElement.$(BASE)

TESTCHANGEOBJECTIVEOBJ	=	testChangeObjective.o
TESTCHANGEOBJECTIVEBIN	=	testChangeObjective.$(BASE)

TESTMEMORYOBJ	=	testMemory.o
TESTMEMORYBIN	=	testMemory.$(BASE)


TESTCOPYCONSASSIGNMENTFUNCOBJ    =  testCopyconsAssignmentfunc.o
TESTCOPYCONSASSIGNMENTFUNCBIN    =  testCopyconsAssignmentfunc.$(BASE)


#------------------------------------------------------------------------------
#--- NOTHING TO CHANGE FROM HERE ON -------------------------------------------
#------------------------------------------------------------------------------

GCCWARN		=	-Wall -W -Wpointer-arith -Wno-unknown-pragmas \
			-Wcast-align -Wwrite-strings -Wconversion \
			-Wctor-dtor-privacy -Wnon-virtual-dtor -Wreorder \
			-Woverloaded-virtual -Wsign-promo -Wsynth -Wundef \
			-Wcast-qual -Wold-style-cast -Wshadow 

#-----------------------------------------------------------------------------
#include make/make.$(OSTYPE).$(ARCH).$(COMP).$(OPT).$(LINK)
#-----------------------------------------------------------------------------

.PHONY: change_exerciser status_exception_test testChangeElement testChangeObjective testCopyconsAssignmentfunc testMemory  clean all

all: change_exerciser status_exception_test testChangeElement testCopyconsAssignmentfunc testMemory

change_exerciser: $(CHANGEBINFILE)

status_exception_test: $(EXCEPTIONBINFILE)

testChangeElement: $(TESTCHANGEELEMENTBIN)

testChangeObjective: $(TESTCHANGEOBJECTIVEBIN)

testCopyconsAssignmentfunc: $(TESTCOPYCONSASSIGNMENTFUNCBIN)

testMemory: $(TESTMEMORYBIN)


$(CHANGEBINFILE): $(LIBFILE) $(CHANGEBINOBJ)
		@echo "-> linking $@"
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(CHANGEBINOBJ) \
		-L$(LIBDIR) -l$(LIBNAME) $(LDFLAGS) -o $@

$(EXCEPTIONBINFILE): $(LIBFILE) $(EXCEPTIONBINOBJ)
		@echo "-> linking $@"
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(EXCEPTIONBINOBJ) \
		-L$(LIBDIR) -l$(LIBNAME) $(LDFLAGS) -o $@

$(TESTCHANGEELEMENTBIN): $(LIBFILE) $(TESTCHANGEELEMENTOBJ)
		@echo "-> linking $@"
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(TESTCHANGEELEMENTOBJ) \
		-L$(LIBDIR) -l$(LIBNAME) $(LDFLAGS) -o $@

$(TESTCHANGEOBJECTIVEBIN): $(LIBFILE) $(TESTCHANGEOBJECTIVEOBJ)
		@echo "-> linking $@"
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(TESTCHANGEOBJECTIVEOBJ) \
		-L$(LIBDIR) -l$(LIBNAME) $(LDFLAGS) -o $@

$(TESTCOPYCONSASSIGNMENTFUNCBIN): $(LIBFILE) $(TESTCOPYCONSASSIGNMENTFUNCOBJ)
		@echo "-> linking $@"
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(TESTCOPYCONSASSIGNMENTFUNCOBJ) \
		-L$(LIBDIR) -l$(LIBNAME) $(LDFLAGS) -o $@

$(TESTMEMORYBIN): $(LIBFILE) $(TESTMEMORYOBJ)
		@echo "-> linking $@"
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(TESTMEMORYOBJ) \
		-L$(LIBDIR) -l$(LIBNAME) $(LDFLAGS) -o $@

%.o:	%.cpp
		@echo "-> compiling $@"
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(BINOFLAGS) -c $< -o $@

clean:
		-rm -rf *.o $(EXCEPTIONBINFILE) $(CHANGEBINFILE) $(TESTCHANGEELEMENTBIN) $(TESTCHANGEOBJECTIVEBIN) \
		$(TESTCOPYCONSASSIGNMENTFUNCBIN) $(TESTMEMORYBIN)

