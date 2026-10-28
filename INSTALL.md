Installation                     {#INSTALL}
============

There are two ways to compile the code: CMake and plain Makefiles. Using CMake
is highly recommended.


CMake Build System
------------------

[CMake](https://cmake.org/) is a build system generator that can create, e.g.,
Makefiles for UNIX and macOS or Visual Studio project files for Windows.

CMake provides an
[extensive documentation](https://cmake.org/cmake/help/latest/manual/cmake.1.html)
explaining available features and use cases as well as an
[FAQ section](https://cmake.org/Wiki/CMake_FAQ). These are the usual steps on a
Linux or macOS system:

    mkdir build
    cd build
    cmake <path/to/SoPlex>
    make

    # optional: run a quick test
    make test

    #optional: run a (slow) memory check
    ctest -T MemCheck

    # optional: install SoPlex executable, library, and headers
    make install

CMake uses an out-of-source build, i.e., compiled binaries and object files are
separated from the source tree and located in another directory, e.g, `build`.
From within this directory, run `cmake <path/to/SoPlex>` to configure your build,
followed by `make` to compile the code according to the current configuration.

Afterwards, successive calls to `make` are going to recompile modified source code,
without requiring another call to `cmake`. The generated executable and libraries are
put in directories `bin` and `lib` respectively.


### Modifying a CMake configuration

There are several options that can be passed to the `cmake <path/to/SoPlex>` call to
modify how the code is built. For all of these options and parameters you have to
use `-D<Parameter_name>=<value>`. Following a list of available options, for the full
list run `cmake <path/to/SoPlex> -LH`:

| CMake option         | Available values; default 1st| Makefile equivalent      | Enabled feature |
|----------------------|------------------------------|--------------------------|---------|
| CMAKE_BUILD_TYPE     | Release, Debug, ...          | OPT=[opt, dbg]           | compile with optimization flags or with debug information |
| BUILD_SHARED_LIBS    | on, off                      | --                       | build shared libraries |
| SHARED               | auto, on, off                | SHARED=[on, off]         | build shared libraries (auto lets BUILD_SHARED_LIBS decide) |
| GMP                  | auto, on, off                | GMP=[true, false]        | exact arithmetic via GMP |
| MPFR                 | auto, on, off                | MPFR=[true, false]       | variable precision via MPFR |
| QUADMATH             | auto, on, off                | QUADMATH=[true, false]   | quadruple precision |
| BOOST                | auto, on, off                | BOOST=[true, false]      | necessary for higher precision arithmetic |
| PAPILO               | auto, on, off                | --                       | PaPILO for presolving |
| SYSTEM_FMT           | off, on, auto                | SYSTEM_FMT=[true, false] | use system's fmt library instead of the bundled one |
| LTO                  | auto, on, off                | LTO=[true, false]        | compile with link-time optimization |
| ZLIB                 | auto, on, off                | ZLIB=[true, false]       | zlib compression for files |
| COVERAGE             | on, off                      | --                       |         |
| SANITIZE             | on, off, thread, memory, address | SANITIZE=...         | which sanitizer to enable if available |
| CMAKE_INSTALL_PREFIX | \<path\>                     | INSTALLDIR=\<path\>      |         |
| GMP_DIR              | \<path/to/GMP\>              | --                       |         |
| ..._DIR              | \<custom/path/to/package\>   | --                       |         |

Parameters can be set all at once or in subsequent calls to `cmake` - extending
or modifying the existing configuration.

### Installation

CMake uses a default directory for installation, e.g., `/usr/local` on Linux.
This can be modified using `-DCMAKE_INSTALL_PREFIX=<custom/path>`.


Makefile
--------

The plain Makefile system only reliably works on UNIX systems:

Description                | Command
---------------------------|---------
On systems with GNU g++    | `make COMP=gnu OPT=opt` (default)
Linux/x86 with Intel C++   | `make COMP=intel OPT=opt`
Tru64 with Compaq C++      | `make COMP=compaq OPT=opt`
Solaris with SUN Forte C++ | `make COMP=sun OPT=opt`
IRIX with SGI Mips Pro C++ | `make COMP=sgi OPT=opt`
HP-UX with HP aCC          | `make COMP=hp OPT=opt`
AIX with VisualAge C++     | `make COMP=ibm OPT=opt`

Then type `make COMP=<as before> OPT=<as before> test`. This should report no fails.

If ZLIB is not available, building may fail. In this case try

    make COMP=<as before> OPT=<as before> ZLIB=false

which will deactivate the possibility to read gzipped LP and MPS files.

### Boost support

Boost is required for higher precision and rational solving methods.
By default, building with boost is enabled. If you wish to only build the SoPlex library
- using cmake, set `cmake -DBOOST=off`.
- using make, use `make BOOST=false makelibfile`.

### GMP support

For using SoPlex as an exact rational LP solver, SoPlex must be compiled
with support for the [GNU Multiple Precision library](http://www.gmplib.org/)
for this.  If GMP is not available, you can deactivate it by building with

    make COMP=<as before> OPT=<as before> GMP=false.

If you use a different build system than the provided Makefile and want to
build with GMP support, you need to define `SOPLEX_HAS_GMP` for the preprocessor
and link with the GMP callable library.

Note for building SCIP with SoPlex:  If SoPlex was built with GMP, then SCIP
also needs to be built with GMP (default).

### Webdemo

SoPlex can be compiled into a single HTML file and used in a browser. This
requires `emscripten` to be installed from their website (unfortunately, e.g.
`sudo apt install emscripten` in Ubuntu Linux is broken):

    https://emscripten.org/docs/getting_started/downloads.html

Then, run

    sh build_webdemo.sh

This will create the file `build_webdemo/bin/soplex.html`. For fast edit
iterations run

    find src | entr -rs 'make -C build_webdemo soplex; echo'

This will rebuild `soplex.html` every time a source file is modified (e.g.
from Visual Studio Code).

Using the Library
=================

Examples on how to use the SoPlex library are provided in the files
`src/soplexmain.cpp` and `src/example.cpp`.
