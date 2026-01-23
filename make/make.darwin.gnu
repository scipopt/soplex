LDFLAGS		+=	-lm
ARFLAGS		=	crs
ZLIB_FLAGS	=
ZLIB_LDFLAGS 	=	-lz
CXXFLAGS        +=      -ffp-contract=off -std=c++14

ifeq ($(LTO),true)
  # for GCC < 10, use just -flto, otherwise use -flto=auto, which should give faster link times
  # -fno-fat-lto-objects (since GCC 5) should improve compilation time a bit
  GCCVERSION := $(shell $(CXX) -dumpversion | cut -f1 -d.)
  LTOFLAG := $(word $(shell expr \( $(GCCVERSION) \>= 10 \) + 1), -flto -flto=auto)

  CXXFLAGS	+=	$(LTOFLAG) -fno-fat-lto-objects
  LDFLAGS	+=	$(LTOFLAG)
ifeq ($(SHARED),true)
  LIBBUILDFLAGS +=	$(LTOFLAG)
endif
endif
