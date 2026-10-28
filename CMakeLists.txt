# Overall structure
#
# POLICIES
# PROJECT DEFINITIONS
#    Version
#    Project declaration
#    Minimum dependency versions
#    Directories and search paths
#    Packaging
# OPTIONS
#    Declarations with preprocessing/printing
# BASIC SETTINGS
#    Shared/static libraries
#    Build type
#    LTO
#    Githash
# DEPENDENCIES
#    Search and process
#    Check for position-independent code
# TOOLS
#    Code coverage
#    Global compilation settings
# SUBDIRECTORIES
# INSTALLATION
#    Licenses
#    Config(-version) files
#    Build-tree targets
# [LOGGING]
#    [Non-default parameters]

cmake_minimum_required(VERSION 3.11)

##### POLICIES #####

# INTERPROCEDURAL_OPTIMIZATION is enforced when enabled.
if(POLICY CMP0069)
   cmake_policy(SET CMP0069 NEW)
endif()

# FindBoost is removed in version CMake 3.30.
if(POLICY CMP0167)
   cmake_policy(SET CMP0167 NEW)
endif()

##### PROJECT DEFINITIONS #####
  
set(SOPLEX_VERSION_MAJOR 9)
set(SOPLEX_VERSION_MINOR 0)
set(SOPLEX_VERSION_PATCH 0)
set(SOPLEX_VERSION_API 22)

# Project
project(SOPLEX
  VERSION ${SOPLEX_VERSION_MAJOR}.${SOPLEX_VERSION_MINOR}.${SOPLEX_VERSION_PATCH}
  LANGUAGES CXX
  )

# Minimum versions of dependencies.
set(GMP_MIN_VERSION "5.0.0")
set(MPFR_MIN_VERSION "3.1.1")
set(BOOST_MIN_VERSION "1.68.0")
set(BOOST_RECOMMENDED_VERSION "1.70.0")

# Set installation directories according to GNU coding standards.
include(GNUInstallDirs)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

# Path to CMake modules such as FindGMP.cmake.
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/Modules")

# Packaging
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE")
set(CPACK_PACKAGE_VERSION_MAJOR "${SOPLEX_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${SOPLEX_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${SOPLEX_VERSION_PATCH}")
set(CPACK_PACKAGE_VENDOR "Zuse Institute Berlin")
include(CPack)

##### OPTIONS #####

option(BUILD_PACKAGE "are we building a package" OFF)

option(BUILD_SHARED_LIBS "Build shared libraries" ON)

set(COVERAGE_CTEST_ARGS "" CACHE STRING "additional ctest arguments for coverage")

option(GITHASH_BRANCH "append the current branch to the githash" OFF)

set(GMP AUTO CACHE STRING "should SoPlex link to GMP")
set_property(CACHE GMP PROPERTY STRINGS AUTO ON OFF)
string(TOUPPER "${GMP}" GMP)

set(QUADMATH AUTO CACHE STRING "should SoPlex link to Quadmath")
set_property(CACHE QUADMATH PROPERTY STRINGS AUTO ON OFF)
string(TOUPPER "${QUADMATH}" QUADMATH)

set(MPFR AUTO CACHE STRING "should SoPlex link to MPFR")
set_property(CACHE MPFR PROPERTY STRINGS AUTO ON OFF)
string(TOUPPER "${MPFR}" MPFR)

option(EMSCRIPTEN_HTML "Emscripten HTML output" OFF)

option(SANITIZE "should sanitizers be enabled if available" OFF)
option(COVERAGE "enable coverage support" OFF)

set(PAPILO AUTO CACHE STRING "should PaPILO library be linked")
set_property(CACHE PAPILO PROPERTY STRINGS AUTO ON OFF)
string(TOUPPER "${PAPILO}" PAPILO)

set(BOOST AUTO CACHE STRING "should Boost libraries be used")
set_property(CACHE BOOST PROPERTY STRINGS AUTO ON OFF)
string(TOUPPER "${BOOST}" BOOST)

set(SYSTEM_FMT OFF CACHE STRING "Use system FMT library instead of bundled headers")
set_property(CACHE SYSTEM_FMT PROPERTY STRINGS OFF ON AUTO)
string(TOUPPER "${SYSTEM_FMT}" SYSTEM_FMT)

set(LTO OFF CACHE STRING "enable Link Time Optimization = Interprocedural Optimization")
set_property(CACHE LTO PROPERTY STRINGS AUTO ON OFF)
string(TOUPPER "${LTO}" LTO)

set(ZLIB AUTO CACHE STRING "should SoPlex link to ZLIB")
set_property(CACHE ZLIB PROPERTY STRINGS AUTO ON OFF)
string(TOUPPER "${ZLIB}" ZLIB)


##### BASIC SETTINGS #####

if(BUILD_SHARED_LIBS)
   message(STATUS "Building shared libraries.")
   set(SOPLEX_HAS_PIC 1)
else()
   message(STATUS "Building static libraries.")
   set(SOPLEX_HAS_PIC 0)
endif()

# Default build type is 'Release'.
if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
elseif(NOT CMAKE_BUILD_TYPE MATCHES "Release|Debug|RelWithDebInfo|MinSizeRel")
   message(FATAL_ERROR "Option CMAKE_BUILD_TYPE has invalid value \"${CMAKE_BUILD_TYPE}\"! Possible are `Release`, `Debug`, `RelWithDebInfo` and `MinSizeRel`.")
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Enable link-time optimization = Interprocedural Optimization
if(LTO)
   include(CheckIPOSupported)
   check_ipo_supported(RESULT LTO_AVAILABLE
      LANGUAGES CXX
      OUTPUT lto_error)
   if(LTO_AVAILABLE)
      message(STATUS "LTO is enabled.")

      # With link-time optimization enabled, we get additional warnings at link time on
      # - some presumed malloc() calls with size 18446744073709551608 ((uint)-1)
      # - some presumed memmove() calls in boost::container::short_vector via papilo to a 0-length memory
      # we disable these warnings (need to pass them to the compiler when linking)
      if((CMAKE_CXX_COMPILER_ID MATCHES "GNU") OR (CMAKE_CXX_COMPILER_ID MATCHES "CLANG"))
         add_link_options(-Wno-alloc-size-larger-than -Wno-stringop-overflow)
      endif()
   elseif(LTO STREQUAL "AUTO")
      message(STATUS "LTO is not available: ${lto_error}")
   else()
      message(FATAL_ERROR "LTO is not available: ${lto_error}")
   endif()
else()
   message(STATUS "LTO is disabled.")
   set(LTO_AVAILABLE 0)
endif()

# Query git hash and create a target that creates a preprocessor macro for it.
include(GitHash)
githash_query(SOPLEX_GITHASH "${PROJECT_SOURCE_DIR}/src/soplex/git_hash.cpp" "SPX_GITHASH" ${GITHASH_BRANCH})
githash_add_target(soplex_update_githash "${PROJECT_SOURCE_DIR}/src/soplex/git_hash.cpp" "SPX_GITHASH" ${GITHASH_BRANCH})
if(SOPLEX_GITHASH)
   message(STATUS "Githash: ${SOPLEX_GITHASH}")
endif()

# Remove config.h that might have been generated by make.
add_custom_target(soplex_cleanup_config ALL
   COMMAND ${CMAKE_COMMAND} -E remove -f "${PROJECT_SOURCE_DIR}/src/soplex/config.h")


##### DEPENDENCIES #####

# Search for ZLIB.
if(ZLIB)
   message(STATUS "Searching for ZLIB.")
   if(ZLIB STREQUAL "AUTO")
      find_package(ZLIB)
   else()
      find_package(ZLIB REQUIRED)
   endif()
   if(ZLIB_FOUND)
      set(SOPLEX_HAS_ZLIB 1)
   else()
      set(SOPLEX_HAS_ZLIB 0)
   endif()
else()
   message(STATUS "ZLIB is disabled.")
   set(SOPLEX_HAS_ZLIB 0)
endif()

# Search for GMP.
if(GMP)
   if(GMP STREQUAL "AUTO")
      find_package(GMP ${GMP_MIN_VERSION})
   else()
      find_package(GMP ${GMP_MIN_VERSION} REQUIRED)
   endif()
   if(GMP_FOUND)
      set(SOPLEX_HAS_GMP 1)
   else()
      set(SOPLEX_HAS_GMP 0)
   endif()
else()
   message(STATUS "GMP is disabled.")
   set(SOPLEX_HAS_GMP 0)
endif()

# Support legacy FindGMP.cmake versions.
if(SOPLEX_HAS_GMP AND NOT TARGET GMP::GMP)
   if(TARGET GMP::gmp)
      message(WARNING "Found GMP version ${GMP_VERSION} has no GMP::GMP target; interfacing GMP::gmp target instead.")
      add_library(GMP::GMP INTERFACE IMPORTED)
      target_link_libraries(GMP::GMP INTERFACE GMP::gmp)
   elseif(TARGET GMP)
      message(WARNING "Found GMP version ${GMP_VERSION} has no GMP::GMP target; interfacing GMP target instead.")
      add_library(GMP::GMP INTERFACE IMPORTED)
      target_link_libraries(GMP::GMP INTERFACE GMP)
   elseif(TARGET gmp)
      message(WARNING "Found GMP version ${GMP_VERSION} has no GMP::GMP target; interfacing gmp target instead.")
      add_library(GMP::GMP INTERFACE IMPORTED)
      target_link_libraries(GMP::GMP INTERFACE gmp)
   elseif(GMP_INCLUDE_DIRS AND GMP_LIBRARIES)
      message(WARNING "Found GMP version ${GMP_VERSION} has no GMP::GMP target; using legacy variables GMP_INCLUDE_DIRS=\"${GMP_INCLUDE_DIRS}\" and GMP_LIBRARIES=\"${GMP_LIBRARIES}\" instead.")
      # GMP_LIBRARIES may be a list, so we need to create each library seperately and glue them via an INTERFACE library.
      list(LENGTH GMP_LIBRARIES num_gmp_lib)
      math(EXPR last_gmp_lib ${num_gmp_lib}-1)
      add_library(GMP_libgmp INTERFACE IMPORTED)
      foreach(id_gmp_lib RANGE ${last_gmp_lib})
         add_library(GMP::libgmp_${id_gmp_lib} UNKNOWN IMPORTED)
         list(GET GMP_LIBRARIES ${id_gmp_lib} the_gmp_lib)
         set_target_properties(GMP::libgmp_${id_gmp_lib} PROPERTIES
            IMPORTED_LOCATION "${the_gmp_lib}"
            INTERFACE_INCLUDE_DIRECTORIES "${GMP_INCLUDE_DIRS}")
         target_link_libraries(GMP_libgmp INTERFACE GMP::libgmp_${id_gmp_lib})
      endforeach()
      add_library(GMP::GMP INTERFACE IMPORTED)
      target_link_libraries(GMP::GMP INTERFACE GMP_libgmp)
   else()
      message(FATAL_ERROR "Found GMP version ${GMP_VERSION} has no GMP::GMP target; unable to recover!")
   endif()
endif()

# Search for Quadmath.
if(QUADMATH)
   if(QUADMATH STREQUAL "AUTO")
      find_package(Quadmath)
   else()
      find_package(Quadmath REQUIRED)
   endif()
   if(Quadmath_FOUND)
      set(SOPLEX_HAS_QUADMATH 1)
   else()
      set(SOPLEX_HAS_QUADMATH 0)
   endif()
else()
   message(STATUS "Quadmath is disabled.")
   set(SOPLEX_HAS_QUADMATH 0)
endif()

# Support for legacy FindQuadmath.cmake version.
if(SOPLEX_HAS_QUADMATH AND NOT TARGET Quadmath::Quadmath)
   if(TARGET Quadmath::quadmath)
      message(WARNING "Found Quadmath version has no Quadmath::Quadmath target; interfacing Quadmath::quadmath target instead.")
      add_library(Quadmath::Quadmath INTERFACE IMPORTED)
      target_link_libraries(Quadmath::Quadmath INTERFACE Quadmath::quadmath)
   else()
      message(FATAL_ERROR "Found Quadmath version has no Quadmath::Quadmath target; unable to recover!")
   endif()
endif()

# Search for MPFR.  
set(SOPLEX_HAS_MPFR 0)
if(MPFR)
   if(MPFR STREQUAL "AUTO")
      find_package(MPFR ${MPFR_MIN_VERSION})
   else()
      find_package(MPFR ${MPFR_MIN_VERSION} REQUIRED)
   endif()
   if(MPFR_FOUND)
      set(SOPLEX_HAS_MPFR 1)
   endif()
else()
   message(STATUS "MPFR is disabled.")
endif()

# Support legacy FindMPFR.cmake versions.
if(SOPLEX_HAS_MPFR AND NOT TARGET MPFR::MPFR)
   if(TARGET MPFR::mpfr)
      message(WARNING "Found MPFR version ${MPFR_VERSION} has no MPFR::MPFR target; interfacing MPFR::mpfr instead.")
      add_library(MPFR::MPFR INTERFACE IMPORTED)
      target_link_libraries(MPFR::MPFR INTERFACE MPFR::mpfr)
   elseif(TARGET MPFR)
      message(WARNING "Found MPFR version ${MPFR_VERSION} has no MPFR::MPFR target; interfacing MPFR instead.")
      add_library(MPFR::MPFR INTERFACE IMPORTED)
      target_link_libraries(MPFR::MPFR INTERFACE MPFR)
   elseif(TARGET mpfr)
      message(WARNING "Found MPFR version ${MPFR_VERSION} has no MPFR::MPFR target; interfacing mpfr instead.")
      add_library(MPFR::MPFR INTERFACE IMPORTED)
      target_link_libraries(MPFR::MPFR INTERFACE mpfr)
   elseif(MPFR_INCLUDE_DIRS AND MPFR_LIBRARIES)
      message(WARNING "Found MPFR version ${MPFR_VERSION} has no MPFR::MPFR target; using legacy variables MPFR_INCLUDE_DIRS=\"${MPFR_INCLUDE_DIRS}\" and MPFR_LIBRARIES=\"${MPFR_LIBRARIES}\" instead.")
      # MPFR_LIBRARIES may be a list, so we need to create each library seperately and glue them via an INTERFACE library.
      list(LENGTH MPFR_LIBRARIES num_mpfr_lib)
      math(EXPR last_mpfr_lib ${num_mpfr_lib}-1)
      add_library(MPFR_libmpfr INTERFACE IMPORTED)
      foreach(id_mpfr_lib RANGE ${last_mpfr_lib})
         add_library(MPFR::libmpfr_${id_mpfr_lib} UNKNOWN IMPORTED)
         list(GET MPFR_LIBRARIES ${id_mpfr_lib} the_mpfr_lib)
         set_target_properties(MPFR::libmpfr_${id_mpfr_lib} PROPERTIES
            IMPORTED_LOCATION "${the_mpfr_lib}"
            INTERFACE_INCLUDE_DIRECTORIES "${MPFR_INCLUDE_DIRS}")
         target_link_libraries(MPFR_libmpfr INTERFACE MPFR::libmpfr_${id_mpfr_lib})
      endforeach()
      add_library(MPFR::MPFR INTERFACE IMPORTED)
      target_link_libraries(MPFR::MPFR INTERFACE MPFR_libmpfr)
   else()
      message(FATAL_ERROR "Found MPFR version ${MPFR_VERSION} has no MPFR::MPFR target; unable to recover!")
   endif()
endif()

# Search for PaPILO.
if(PAPILO)
   if(PAPILO STREQUAL "AUTO")
      find_package(PAPILO ${PAPILO_MIN_VERSION} CONFIG HINTS "${PAPILO_DIR}")
   else()
      find_package(PAPILO ${PAPILO_MIN_VERSION} REQUIRED CONFIG HINTS "${PAPILO_DIR}")
   endif()
   if(PAPILO_FOUND)
      set(SOPLEX_HAS_PAPILO 1)
      if(PAPILO_VERSION VERSION_LESS 4)
         message(STATUS "Found PaPILO version ${PAPILO_VERSION}: ${PAPILO_CONFIG}")
      endif()
   else()
      set(SOPLEX_HAS_PAPILO 0)
   endif()
else()
   message(STATUS "PaPILO is disabled.")
   set(SOPLEX_HAS_PAPILO 0)
endif()

# Support for legacy PaPILO version.
if(SOPLEX_HAS_PAPILO AND NOT TARGET PAPILO::PAPILO)
   if(TARGET papilo)
      message(WARNING "Found PaPILO version ${PAPILO_VERSION} has no PAPILO::PAPILO target; interfacing papilo target instead.")
      add_library(PAPILO::PAPILO INTERFACE IMPORTED)
      target_link_libraries(PAPILO::PAPILO INTERFACE papilo)
   else()
      message(FATAL_ERROR "Found PaPILO version ${PAPILO_VERSION} has no PAPILO::PAPILO target; unable to recover!")
   endif()

   # Legacy PaPILO may depend on Quadmath::quadmath, but we may have found Quadmath::Quadmath.
   if(NOT TARGET Quadmath::quadmath AND TARGET Quadmath::Quadmath)
      message(WARNING "Creating target alias Quadmath::quadmath -> Quadmath::Quadmath to satisfy legacy PaPILO.")
      add_library(Quadmath::quadmath INTERFACE IMPORTED)
      target_link_libraries(Quadmath::quadmath INTERFACE Quadmath::Quadmath)
   endif()
endif()

# Search for Boost.
set(SOPLEX_HAS_BOOST 0)
set(SOPLEX_HAS_CPPRAT 0)
set(SOPLEX_HAS_CPPMPF 0)
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.15)
   cmake_policy(SET CMP0093 NEW) # Boost_VERSION shall be in x.y.z format.
endif()
if(BOOST)
   message(STATUS "Searching for Boost.")
   if(BOOST STREQUAL "AUTO")
      find_package(Boost ${BOOST_MIN_VERSION})
   else()
      find_package(Boost ${BOOST_MIN_VERSION} REQUIRED)
   endif()
   if(Boost_FOUND)
      if(NOT DEFINED Boost_VERSION_STRING)
         math(EXPR _boost_major "${Boost_VERSION} / 100000")
         math(EXPR _boost_minor "${Boost_VERSION} / 100 % 1000")
         math(EXPR _boost_patch "${Boost_VERSION} % 100")
         set(Boost_VERSION_STRING "${_boost_major}.${_boost_minor}.${_boost_patch}")
      endif()
      message(STATUS "Found Boost version ${Boost_VERSION_STRING}.")
      if(Boost_VERSION_STRING VERSION_LESS BOOST_RECOMMENDED_VERSION)
         message(WARNING "Compiling with Boost versions below 1.70 may result in a significant performance degradation for high precision / exact arithmetic!")
      endif()

      set(SOPLEX_HAS_BOOST 1)
      if(NOT SOPLEX_HAS_GMP)
         set(SOPLEX_HAS_CPPRAT 1)
      endif()
      if(NOT SOPLEX_HAS_MPFR)
         set(SOPLEX_HAS_CPPMPF 1)
      endif()
   else()
      if(GMP AND NOT GMP STREQUAL "AUTO")
         message(FATAL_ERROR "GMP=on and BOOST=auto marks Boost as required!")
      endif()
      if(MPFR AND NOT MPFR STREQUAL "AUTO")
         message(FATAL_ERROR "MPFR=on and BOOST=auto marks Boost as required!")
      endif()
      if(SOPLEX_HAS_QUADMATH)
         message(STATUS "Disabling Quadmath since we don't need it.")
         set(SOPLEX_HAS_QUADMATH 0) # Disable Quadmath again if Boost is not found.
      endif()
      if(SOPLEX_HAS_GMP)
         message(STATUS "Disabling GMP since we don't need it.")
         set(SOPLEX_HAS_GMP 0) # Disable GMP again if Boost is not found.
      endif()
      if(SOPLEX_HAS_MPFR)
         message(STATUS "Disabling MPFR since we don't need it.")
         set(SOPLEX_HAS_MPFR 0) # Disable MPFR again if Boost is not found.
      endif()
   endif()
else()
   message(STATUS "Boost is disabled.")
   if(GMP AND NOT GMP STREQUAL "AUTO")
      message(FATAL_ERROR "Conflicting options: BOOST=off and GMP=on.")
   elseif(QUADMATH AND NOT QUADMATH STREQUAL "AUTO")
      message(FATAL_ERROR "Conflicting options: BOOST=off and QUADMATH=on.")
   elseif(MPFR AND NOT MPFR STREQUAL "AUTO")
      message(FATAL_ERROR "Conflicting options: BOOST=off and MPFR=on.")
   endif()
   if(SOPLEX_HAS_QUADMATH)
      message(STATUS "Disabling Quadmath since we don't need it.")
      set(SOPLEX_HAS_QUADMATH 0) # Disable Quadmath again if Boost is not found.
   endif()
   if(SOPLEX_HAS_GMP)
      message(STATUS "Disabling GMP since we don't need it.")
      set(SOPLEX_HAS_GMP 0) # Disable GMP again if Boost is not found.
   endif()
   if(SOPLEX_HAS_MPFR)
      message(STATUS "Disabling MPFR since we don't need it.")
      set(SOPLEX_HAS_MPFR 0) # Disable MPFR again if Boost is not found.
   endif()
endif()

# Output supported precisions.
set(soplex_precision "double")
if(SOPLEX_HAS_QUADMATH)
   set(soplex_precision "${soplex_precision}, quad")
endif()
if(SOPLEX_HAS_BOOST AND SOPLEX_HAS_MPFR)
   set(soplex_precision "${soplex_precision}, variable (MPFR)")
elseif(SOPLEX_HAS_BOOST AND SOPLEX_HAS_CPPMPF)
   set(soplex_precision "${soplex_precision}, high (Boost)")
endif()
if(SOPLEX_HAS_BOOST AND SOPLEX_HAS_GMP)
   set(soplex_precision "${soplex_precision}, rational (GMP)")
elseif(SOLPEX_HAS_BOOST AND SOPLEX_HAS_CPPRAT)
   set(soplex_precision "${soplex_precision}, rational (Boost)")
endif()
message(STATUS "Number precisions supported by executable: ${soplex_precision}")

# FMT
set(SOPLEX_HAS_SYSTEM_FMT 0)
if(SYSTEM_FMT)
   if(SYSTEM_FMT STREQUAL "AUTO")
      find_package(fmt CONFIG)
   else()
      find_package(fmt CONFIG REQUIRED)
   endif()
   if(fmt_FOUND)
      set(SOPLEX_HAS_SYSTEM_FMT 1)
   endif()
endif()
if(SOPLEX_HAS_SYSTEM_FMT)
   message(STATUS "Using system FMT library.")
else()
   message(STATUS "Using bundled FMT library.")
endif()

# Check if dependencies are compiled with position-independent code.
if(BUILD_SHARED_LIBS AND SOPLEX_HAS_PAPILO AND (DEFINED PAPILO_HAS_PIC) AND (NOT PAPILO_HAS_PIC))
   message(FATAL_ERROR "Cannot link shared SoPlex library to a PaPILO library without position-independent code!")
endif()


##### TOOLS #####

# Enable coverage support.
if(COVERAGE)
    include(CodeCoverage)
    APPEND_COVERAGE_COMPILER_FLAGS()
    set(COVERAGE_EXCLUDES '/usr*')
    SETUP_TARGET_FOR_COVERAGE(NAME coverage
       EXECUTABLE bash -c "ctest ${COVERAGE_CTEST_ARGS}"
       DEPENDENCIES soplex)
endif()

# Flags for release with debug information config.
set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_DEBUG} ${CMAKE_C_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS_RELEASE}")

# Set the correct rpath for OS X
set(CMAKE_MACOSX_RPATH ON)

# Add search path such that installed libsoplex is found.
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")

# Each directory used during linking is appended to the rpath of installed libraries.
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Use C++14 standard.
set(CMAKE_CXX_STANDARD 14)

# Set function visibility default to hidden.
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN 1)

if(NOT "${CMAKE_CXX_COMPILE_FEATURES}" MATCHES "cxx_thread_local")
   add_definitions(-DTHREADLOCAL=)
endif()

# Disable fused floating point contraction to enhance reproducibility across compilers and architectures.
if(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
   # TODO: manually export symbols
   set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
   add_compile_options(/fp:precise /bigobj)
   add_definitions(/wd4100)
   add_definitions(/wd4661)
   add_definitions(-D_CRT_SECURE_NO_WARNINGS)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
   add_compile_options(-fp-model=precise)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
   # require at least gcc 5
   if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5)
      message(WARNING "GCC version not supported, should be at least 5!")
   endif()
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-psabi")
   add_compile_options(-ffp-contract=off)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
   add_compile_options(-ffp-contract=off)
endif()

# Enable all warnings and explicitly switch on colored output for ninja.
if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
   set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -Wno-strict-overflow")
   if(CMAKE_GENERATOR STREQUAL "Ninja")
      set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -fdiagnostics-color=always")
   endif()
endif()


##### SUBDIRECTORIES #####

# Main code.
add_subdirectory(src)

# Create tests.
add_subdirectory(tests/c_interface)
add_subdirectory(check)
enable_testing()

##### INSTALLATION #####

# Install license files.
install(FILES "${PROJECT_SOURCE_DIR}/LICENSE"
   DESTINATION "${CMAKE_INSTALL_DATADIR}/licenses/soplex")

if(NOT SOPLEX_HAS_SYSTEM_FMT)
   install(FILES "${PROJECT_SOURCE_DIR}/src/soplex/external/fmt/LICENSE.rst"
      DESTINATION "${CMAKE_INSTALL_DATADIR}/licenses/soplex/fmt")
endif()

# Create soplex-config-version.cmake.
include(CMakePackageConfigHelpers)
write_basic_package_version_file("${CMAKE_BINARY_DIR}/soplex-config-version.cmake"
   VERSION ${SOPLEX_VERSION_MAJOR}.${SOPLEX_VERSION_MINOR}.${SOPLEX_VERSION_PATCH}
   COMPATIBILITY SameMajorVersion
)

# Unless we build a package, we populate search directories for dependencies.
if(NOT BUILD_PACKAGE AND ZLIB_DIR)
   set(SOPLEX_CONFIG_ZLIB_DIR "if(NOT ZLIB_FOUND AND NOT ZLIB_DIR)\n   set(ZLIB_DIR \"${ZLIB_DIR}\")\nendif()")
endif()
if(NOT BUILD_PACKAGE AND GMP_DIR)
   set(SOPLEX_CONFIG_GMP_DIR "if(NOT GMP_FOUND AND NOT GMP_DIR)\n   set(GMP_DIR \"${GMP_DIR}\")\nendif()")
endif()
if(NOT BUILD_PACKAGE AND QUADMATH_DIR)
   set(SOPLEX_CONFIG_QUADMATH_DIR "if(NOT Quadmath_FOUND AND NOT QUADMATH_DIR)\n   set(QUADMATH_DIR \"${QUADMATH_DIR}\")\nendif()")
endif()
if(NOT BUILD_PACKAGE AND MPFR_DIR)
   set(SOPLEX_CONFIG_MPFR_DIR "if(NOT MPFR_FOUND AND NOT MPFR_DIR)\n   set(MPFR_DIR \"${MPFR_DIR}\")\nendif()")
endif()
if(NOT BUILD_PACKAGE AND PAPILO_DIR)
   set(SOPLEX_CONFIG_PAPILO_DIR "if(NOT PAPILO_FOUND AND NOT PAPILO_DIR)\n   set(PAPILO_DIR \"${PAPILO_DIR}\")\nendif()")
endif()
if(NOT BUILD_PACKAGE AND Boost_DIR)
   set(SOPLEX_CONFIG_BOOST_DIR "if(NOT Boost_FOUND AND NOT Boost_DIR)\n   set(Boost_DIR \"${Boost_DIR}\")\nendif()")
endif()

# Configure soplex-config.cmake for the build-tree:
# We add the CMake module path from the sources.
set(EXTRA_CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/Modules")
configure_file("${PROJECT_SOURCE_DIR}/soplex-config.cmake.in"
   "${CMAKE_BINARY_DIR}/soplex-config.cmake"
   @ONLY)

# Configure soplex-config.cmake for the installation.
# We add the current directory of the installed soplex-config.cmake.
set(EXTRA_CMAKE_MODULE_PATH "\${CMAKE_CURRENT_LIST_DIR}")
configure_file(${PROJECT_SOURCE_DIR}/soplex-config.cmake.in
   "${PROJECT_BINARY_DIR}/${CMAKE_FILES_DIRECTORY}/soplex-config.cmake"
   @ONLY)

# Install CMake config and modules.
install(FILES
   "${PROJECT_BINARY_DIR}/${CMAKE_FILES_DIRECTORY}/soplex-config.cmake"
   "${CMAKE_BINARY_DIR}/soplex-config-version.cmake"
   $<$<BOOL:${SOPLEX_HAS_GMP}>:${PROJECT_SOURCE_DIR}/cmake/Modules/FindGMP.cmake>
   $<$<BOOL:${SOPLEX_HAS_QUADMATH}>:${PROJECT_SOURCE_DIR}/cmake/Modules/FindQuadmath.cmake>
   $<$<BOOL:${SOPLEX_HAS_MPFR}>:${PROJECT_SOURCE_DIR}/cmake/Modules/FindMPFR.cmake>
   DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/soplex")

# Install the soplex-targets export set.
install(EXPORT soplex-targets
   NAMESPACE SOPLEX::
   DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/soplex")

# Export targets for the build-tree.
export(EXPORT soplex-targets
   NAMESPACE SOPLEX::
   FILE "${CMAKE_BINARY_DIR}/soplex-targets.cmake")

# If we include this directory then we mark SoPlex as found.
set(SOPLEX_FOUND 1)
