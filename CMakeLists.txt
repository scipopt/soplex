cmake_minimum_required(VERSION 3.3.0)

set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_DEBUG} ${CMAKE_C_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS_RELEASE}")

project(SOPLEX
    LANGUAGES CXX)

set(SOPLEX_VERSION_MAJOR 4)
set(SOPLEX_VERSION_MINOR 1)
set(SOPLEX_VERSION_PATCH 0)
set(SOPLEX_VERSION_SUB 0)
set(SOPLEX_VERSION_API 10)

set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/COPYING")
set(CPACK_PACKAGE_VERSION_MAJOR "${SOPLEX_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${SOPLEX_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${SOPLEX_VERSION_PATCH}")
set(CPACK_PACKAGE_VENDOR "Zuse Institute Berlin")
include(CPack)

option(ZLIB "Use ZLIB" ON)
option(GMP "Use GMP" ON)
option(STATIC_GMP "Prefer static GMP lib" OFF)
option(Boost "Use Boost (required to build the binary). Disable if you only want to build libsoplex." OFF)
option(MPFR "Use MPFR" OFF)
option(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS "Export all symbols into the DLL" ON)
option(MT "use static runtime libraries for Visual Studio compiler" OFF)
option(SANITIZE_ADDRESS "should the address sanitizer be enabled in debug mode if available" OFF)
option(SANITIZE_MEMORY "should the memory sanitizer be enabled in debug mode if available" OFF)
option(SANITIZE_UNDEFINED "should the undefined behavior sanitizer be enabled in debug mode if available" ON)
option(SANITIZE_THREAD "should the thread sanitizer be enabled in debug mode if available" OFF)
option(COVERAGE "enable coverage support" OFF)
SET(COVERAGE_CTEST_ARGS "" CACHE STRING "additional ctest arguments for coverage")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# path to e.g. findGMP module
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules)

# make 'Release' the default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# set the correct rpath for OS X
set(CMAKE_MACOSX_RPATH ON)

# use C++11 standard
set(CMAKE_CXX_STANDARD 11)

# set function visibility default to hidden
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN 1)

# create a target for updating the current git hash
file(WRITE ${CMAKE_BINARY_DIR}/soplex_update_githash.cmake "
find_program(GIT git)
if(EXISTS \${DST})
   file(STRINGS \${DST} GITHASH_OLD)
   string(REGEX REPLACE \"#define SPX_GITHASH \\\"(.*)\\\"\" \"\\\\1\" GITHASH_OLD \${GITHASH_OLD})
endif()
if((GIT) AND (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git))
   execute_process(
      COMMAND \${GIT} describe --always --dirty
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      OUTPUT_VARIABLE GITHASH OUTPUT_STRIP_TRAILING_WHITESPACE)
   string(REGEX REPLACE \"^.*-g\" \"\" GITHASH \${GITHASH})
   if(NOT \${GITHASH} STREQUAL \"\${GITHASH_OLD}\")
      file(WRITE \${DST} \"#define SPX_GITHASH \\\"\${GITHASH}\\\"\n\")
   endif()
else()
   set(GITHASH \${GITHASH_OLD})
endif()
message(STATUS \"Git hash: \" \${GITHASH})
")
add_custom_target(soplex_update_githash
                  COMMAND ${CMAKE_COMMAND} -DDST=${PROJECT_SOURCE_DIR}/src/soplex/git_hash.cpp
                                           -P ${CMAKE_BINARY_DIR}/soplex_update_githash.cmake)

if(NOT "${CMAKE_CXX_COMPILE_FEATURES}" MATCHES "cxx_thread_local")
    add_definitions(-DTHREADLOCAL=)
endif()

# enable coverage support
if(COVERAGE)
    include(CodeCoverage)
    APPEND_COVERAGE_COMPILER_FLAGS()
    set(COVERAGE_EXCLUDES '/usr*')
    SETUP_TARGET_FOR_COVERAGE(NAME coverage
                          EXECUTABLE bash -c "ctest ${COVERAGE_CTEST_ARGS}"
                          DEPENDENCIES soplex)
endif()

if(ZLIB)
    find_package(ZLIB)
endif()
if(ZLIB_FOUND)
    add_definitions(-DSOPLEX_WITH_ZLIB)
    set(libs ${libs} ${ZLIB_LIBRARIES})
    include_directories(${ZLIB_INCLUDE_DIRS})
endif()

if(GMP)
    find_package(GMP)
endif()
if(GMP_FOUND)
    add_definitions(-DSOPLEX_WITH_GMP)
    include_directories(${GMP_INCLUDE_DIRS})
    set(libs ${libs} ${GMP_LIBRARIES})
endif()

if(MPFR)
  find_package(MPFR)
endif()
if(MPFR_FOUND)
    message("SoPlex with Boost MPFR libraries.")
    add_definitions(-DSOPLEX_WITH_MPFR)
    include_directories(${MPFR_INCLUDE_DIRS})
    set(libs ${libs} ${MPFR_LIBRARIES})
else()
    message("SoPlex with Boost CPP multiprecision libraries.")
    add_definitions(-DSOPLEX_WITH_CPPMPF)
endif()



# Made with: https://stackoverflow.com/a/3917033/4223038
if(Boost)
    # Technically the version 1.70.0 is not needed for program_option, but for the
    # multiprecision, but it's header only, this is one way to guarantee that the
    # user has the right version.   
    find_package(Boost 1.58.0 COMPONENTS program_options REQUIRED) # Refer to comments in cmake/Modules/FindBoost.cmake
    add_definitions(-DSOPLEX_WITH_BOOST)
    if(Boost_FOUND)
        include_directories(${Boost_INCLUDE_DIRS})
    else()
        message(WARNING "Boost could not be found. If you do not inted to use \
        the multiprecision feature and you only require the soplex library, \
        consider building with -DBoost=off.")
    endif()
else()
    message(WARNING "When building without Boost, only the SoPlex library \
    will be built. Boost with program_options is required to build the binary.")
endif()


# enable all warnings and explicitly switch on colored output for ninja
if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -lboost_program_options -Wall -Wno-strict-overflow")
    if(CMAKE_GENERATOR STREQUAL "Ninja")
        set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -fdiagnostics-color=always")
    endif()
endif()
if(MSVC)
#    add_definitions(/W4)
    add_definitions(/wd4100)
    add_definitions(/wd4661)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
    add_compile_options(/bigobj)

    # set MT or MTd flag when using the static VS runtime libraries
    if(MT)
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")
    endif()
endif()

add_subdirectory(src)
add_subdirectory(check)

enable_testing()
