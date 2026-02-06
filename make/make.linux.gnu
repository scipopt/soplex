LDFLAGS		+=	-lm
ARFLAGS		=	crs
ZLIB_FLAGS	=
ZLIB_LDFLAGS 	=	-lz
CXXFLAGS        +=      $(GCCWARN) -ffp-contract=off -std=c++14

# ThreadSanitizer (https://github.com/google/sanitizers/wiki/ThreadSanitizerCppManual)
ifeq ($(SANITIZE),thread)
  SANITIZERFLAGS = -g -fsanitize=thread
endif

# AddressSanitizer (https://github.com/google/sanitizers/wiki/AddressSanitizer)
ifeq ($(SANITIZE),address)
  SANITIZERFLAGS = -g -fsanitize=address
endif

ifeq ($(SANITIZE),memory)
  $(warning Memory Sanitizer not available with GCC)
endif

# UndefinedBehaviorSanitizer if SANITIZE is true, thread, address, or memory
ifneq ($(filter $(SANITIZE),true thread address memory),)
  SANITIZERFLAGS += -g -fsanitize=undefined -fsanitize=float-cast-overflow -fsanitize=float-divide-by-zero
endif

CXXFLAGS += $(SANITIZERFLAGS)
ifeq ($(SHARED),true)
  LIBBUILDFLAGS += $(SANITIZERFLAGS)
endif
LDFLAGS += $(SANITIZERFLAGS)

ifeq ($(LTO),true)
  # for GCC < 10, use just -flto, otherwise use -flto=auto, which should give faster link times
  # -fno-fat-lto-objects (since GCC 5) should improve compilation time a bit
  GCCVERSION := $(shell $(CXX) -dumpversion | cut -f1 -d.)
  LTOFLAG := $(word $(shell expr \( $(GCCVERSION) \>= 10 \) + 1), -flto -flto=auto)

  CXXFLAGS	+=	$(LTOFLAG) -fno-fat-lto-objects
  LDFLAGS	+=	$(LTOFLAG) -Wno-stringop-overflow -Wno-alloc-size-larger-than
ifeq ($(SHARED),true)
  LIBBUILDFLAGS +=	$(LTOFLAG) -Wno-stringop-overflow -Wno-alloc-size-larger-than
endif
endif
