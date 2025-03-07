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
