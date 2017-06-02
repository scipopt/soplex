cmake_minimum_required(VERSION 2.6)

set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_DEBUG} ${CMAKE_C_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS_RELEASE}")

project(SoPlex)

set(SoPlex_VERSION_MAJOR 3)
set(SoPlex_VERSION_MINOR 0)
set(SoPlex_VERSION_PATCH 1)
set(SoPlex_VERSION_SUB 0)
set(SoPlex_VERSION_API 0)

set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/COPYING")
set(CPACK_PACKAGE_VERSION_MAJOR "${SoPlex_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${SoPlex_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${SoPlex_VERSION_PATCH}")
set(CPACK_PACKAGE_VENDOR "Zuse Institute Berlin")
include(CPack)

option(ZLIB "Use ZLIB" ON)
option(GMP "Use GMP" ON)
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
option(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS "Export all symbols into the DLL" ON)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${SoPlex_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${SoPlex_BINARY_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${SoPlex_BINARY_DIR})

# update the cmake submodule
if(EXISTS ${PROJECT_SOURCE_DIR}/.git AND NOT EXISTS ${PROJECT_SOURCE_DIR}/cmake/.git)
   find_program(GIT git)
   if(GIT)
      execute_process(
         COMMAND ${GIT} submodule update --init ${PROJECT_SOURCE_DIR}/cmake)
   endif()
endif()

# path to e.g. findGMP module
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules)

# make 'Release' the default build type
if(CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE Release)
endif()

# use C++11 standard
set(CMAKE_CXX_STANDARD 11)

#create target to update the githash
file(WRITE ${CMAKE_BINARY_DIR}/update_githash.cmake "
find_program(GIT git)
if(GIT)
   execute_process(
      COMMAND \${GIT} describe --always --dirty
      OUTPUT_VARIABLE GITHASH OUTPUT_STRIP_TRAILING_WHITESPACE)
   if(GITHASH STREQUAL \"\")
      set(GITHASH \"000000\")
   endif()
   string(REGEX REPLACE \"^.*-g\" \"\" GITHASH \${GITHASH})
else()
   set(GITHASH \"000000\")
endif()
file(WRITE \${DST} \"#define SPX_GITHASH \\\"\${GITHASH}\\\"\")
message(STATUS \"Git hash: \" \${GITHASH})
")
# encode the current git hash into the binary
add_custom_target(update_githash
                  COMMAND ${CMAKE_COMMAND} -DDST=${PROJECT_SOURCE_DIR}/src/git_hash.cpp
                                           -P ${CMAKE_BINARY_DIR}/update_githash.cmake)

option(LEGACY "Build SoPlex in legacy mode" OFF)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.7")
        set(LEGACY ON)
    endif()
endif()

if(LEGACY)
    add_definitions(-DSOPLEX_LEGACY)
endif()

if(ZLIB)
    find_package(ZLIB QUIET)
endif()
if(ZLIB_FOUND)
    add_definitions(-DSOPLEX_WITH_ZLIB)
    set(libs ${libs} ${ZLIB_LIBRARIES})
    include_directories(${ZLIB_INCLUDE_DIRS})
endif()

if(GMP AND NOT LEGACY)
    find_package(GMP QUIET)
endif()
if(GMP_FOUND)
    add_definitions(-DSOPLEX_WITH_GMP)
    include_directories(${GMP_INCLUDE_DIRS})
    set(libs ${libs} ${GMP_LIBRARIES})
endif()

# enable all warnings and explicitly switch on colored output for ninja
if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -Wno-strict-overflow")
    if(CMAKE_GENERATOR STREQUAL "Ninja")
        set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -fdiagnostics-color=always")
    endif()
endif()
if(MSVC)
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} /W4")
endif()

add_subdirectory(src)
add_subdirectory(check)

enable_testing()