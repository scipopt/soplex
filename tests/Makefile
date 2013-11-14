#-----------------------------------------------------------------------------
# detect host architecture
#-----------------------------------------------------------------------------
include ../make/make.detecthost

#
# General setup.
#

NAME		=	soplex
VERSION		:=	1.8.0.1

OPT		=	opt
LINK		=	static
COMP		=	gnu
GMP		=	false

BASE		=	$(OSTYPE).$(ARCH).$(COMP).$(OPT)

CXX		=	g++
DCXX		=	g++
LIBEXT		=	a

SRCDIR		=	../src
BINDIR		=	../bin
LIBDIR		=	../lib

DEPEND		=	dependencies

CPPFLAGS	=	-I$(SRCDIR)
CXXFLAGS	=	-g
DFLAGS		=	-MM
BINOFLAGS	=	
LIBOFLAGS	=	
LDFLAGS		=	-lm -lz #-static
GMP_FLAGS	=
GMP_LDFLAGS	=	-lgmpxx -lgmp -lz

LIBNAME		=	$(NAME)-$(VERSION).$(BASE)
LIBFILE		=	$(LIBDIR)/lib$(LIBNAME).$(LIBEXT)

ifeq ($(GMP),true)
CPPFLAGS	+=	-DSOPLEX_WITH_GMP $(GMP_FLAGS)
LDFLAGS		+=	$(GMP_LDFLAGS)
endif

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

TESTRATIONALOBJ	=	testRational.o
TESTRATIONALBIN	=	testRational.$(BASE)

TESTCOPYCONSASSIGNMENTFUNCOBJ    =  testCopyconsAssignmentfunc.o
TESTCOPYCONSASSIGNMENTFUNCBIN    =  testCopyconsAssignmentfunc.$(BASE)

ALLOBJ		= $(CHANGEBINOBJ) $(EXCEPTIONBINOBJ) $(TESTCHANGEOBJECTIVEOBJ) $(TESTCHANGEOBJECTIVEOBJ) $(TESTMEMORYOBJ) $(TESTRATIONALOBJ) $(TESTCOPYCONSASSIGNMENTFUNCOBJ)

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

.PHONY: change_exerciser status_exception_test testChangeElement testChangeObjective testCopyconsAssignmentfunc testMemory testRational  clean all

all: change_exerciser status_exception_test testChangeElement testCopyconsAssignmentfunc testMemory testRational

change_exerciser: $(CHANGEBINFILE)

status_exception_test: $(EXCEPTIONBINFILE)

testChangeElement: $(TESTCHANGEELEMENTBIN)

testChangeObjective: $(TESTCHANGEOBJECTIVEBIN)

testCopyconsAssignmentfunc: $(TESTCOPYCONSASSIGNMENTFUNCBIN)

testMemory: $(TESTMEMORYBIN)

testRational: $(TESTRATIONALBIN)


$(CHANGEBINFILE): $(LIBFILE) $(CHANGEBINOBJ)
		@echo "-> linking $@"
ifeq ($(VERBOSE), true)
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(CHANGEBINOBJ) \
		-L$(LIBDIR) -l$(LIBNAME) $(LDFLAGS) -o $@
else
		@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(CHANGEBINOBJ) \
		-L$(LIBDIR) -l$(LIBNAME) $(LDFLAGS) -o $@
endif

$(EXCEPTIONBINFILE): $(LIBFILE) $(EXCEPTIONBINOBJ)
		@echo "-> linking $@"
ifeq ($(VERBOSE), true)
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(EXCEPTIONBINOBJ) \
		-L$(LIBDIR) -l$(LIBNAME) $(LDFLAGS) -o $@
else
		@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(EXCEPTIONBINOBJ) \
		-L$(LIBDIR) -l$(LIBNAME) $(LDFLAGS) -o $@
endif

$(TESTCHANGEELEMENTBIN): $(LIBFILE) $(TESTCHANGEELEMENTOBJ)
		@echo "-> linking $@"
ifeq ($(VERBOSE), true)
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(TESTCHANGEELEMENTOBJ) \
		-L$(LIBDIR) -l$(LIBNAME) $(LDFLAGS) -o $@
else
		@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(TESTCHANGEELEMENTOBJ) \
		-L$(LIBDIR) -l$(LIBNAME) $(LDFLAGS) -o $@
endif

$(TESTCHANGEOBJECTIVEBIN): $(LIBFILE) $(TESTCHANGEOBJECTIVEOBJ)
		@echo "-> linking $@"
ifeq ($(VERBOSE), true)
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(TESTCHANGEOBJECTIVEOBJ) \
		-L$(LIBDIR) -l$(LIBNAME) $(LDFLAGS) -o $@
else
		@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(TESTCHANGEOBJECTIVEOBJ) \
		-L$(LIBDIR) -l$(LIBNAME) $(LDFLAGS) -o $@
endif

$(TESTCOPYCONSASSIGNMENTFUNCBIN): $(LIBFILE) $(TESTCOPYCONSASSIGNMENTFUNCOBJ)
		@echo "-> linking $@"
ifeq ($(VERBOSE), true)
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(TESTCOPYCONSASSIGNMENTFUNCOBJ) \
		-L$(LIBDIR) -l$(LIBNAME) $(LDFLAGS) -o $@
else
		@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(TESTCOPYCONSASSIGNMENTFUNCOBJ) \
		-L$(LIBDIR) -l$(LIBNAME) $(LDFLAGS) -o $@
endif

$(TESTMEMORYBIN): $(LIBFILE) $(TESTMEMORYOBJ)
		@echo "-> linking $@"
ifeq ($(VERBOSE), true)
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(TESTMEMORYOBJ) \
		-L$(LIBDIR) -l$(LIBNAME) $(LDFLAGS) -o $@
else
		@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(TESTMEMORYOBJ) \
		-L$(LIBDIR) -l$(LIBNAME) $(LDFLAGS) -o $@
endif

$(TESTRATIONALBIN): $(LIBFILE) $(TESTRATIONALOBJ)
		@echo "-> linking $@"
ifeq ($(VERBOSE), true)
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(TESTRATIONALOBJ) \
		-L$(LIBDIR) -l$(LIBNAME) $(LDFLAGS) -o $@
else
		@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(TESTRATIONALOBJ) \
		-L$(LIBDIR) -l$(LIBNAME) $(LDFLAGS) -o $@
endif

.PHONY: $(DEPEND)
depend:
		$(SHELL) -ec '$(DCXX) $(DFLAGS) $(CPPFLAGS) \
		$(ALLOBJ:.o=.cpp) \
		| sed '\''s|^\([0-9A-Za-z]\{1,\}\)\.o|\1.o|g'\'' \
		>$(DEPEND)'

-include	$(DEPEND)


%.o:	%.cpp
		@echo "-> compiling $@"
ifeq ($(VERBOSE), true)
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(BINOFLAGS) -c $< -o $@
else
		@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(BINOFLAGS) -c $< -o $@
endif

clean:
		-rm -rf *.o $(EXCEPTIONBINFILE) $(CHANGEBINFILE) $(TESTCHANGEELEMENTBIN) $(TESTCHANGEOBJECTIVEBIN) \
		$(TESTCOPYCONSASSIGNMENTFUNCBIN) $(TESTMEMORYBIN) $(TESTRATIONALBIN)

