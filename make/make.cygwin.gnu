LDFLAGS		+=	-lm
ARFLAGS		=	crs
ZLIB_FLAGS	=
ZLIB_LDFLAGS 	=	-lz
CXXFLAGS        +=      $(GCCWARN) -ffp-contract=off -std=c++14
