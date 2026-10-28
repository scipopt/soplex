# Overall structure:
#  
# Configure build-dependent header file.
# [Define common compilation settings.]
# Define targets along with properties, dependencies, includes, linked libraries, definitions.
# Define target aliases.
# Add sanitizers.
# Install headers.
# Install targets.
# [Generate export header file.]

# Create build/soplex/config.h, which contains all definitions based on the configuration.
if(BUILD_SHARED_LIBS)
   set(SOPLEX_STATIC_LIB 0)
else()
   set(SOPLEX_STATIC_LIB 1)
endif()
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/soplex/config.h.in"
   "${PROJECT_BINARY_DIR}/soplex/config.h"
   @ONLY)

# Create a SoPlex library.
add_library(libsoplex
   soplex_interface.cpp
   soplex/didxset.cpp
   soplex/idxset.cpp
   soplex/mpsinput.cpp
   soplex/nameset.cpp
   soplex/spxdefines.cpp
   soplex/spxgithash.cpp
   soplex/spxid.cpp
   soplex/spxout.cpp
   soplex/usertimer.cpp
   soplex/wallclocktimer.cpp
   )
set_target_properties(libsoplex PROPERTIES
   PREFIX ""
   CXX_VISIBILITY_PRESET default
   VERSION ${SOPLEX_VERSION_MAJOR}.${SOPLEX_VERSION_MINOR}.${SOPLEX_VERSION_PATCH}
   SOVERSION ${SOPLEX_VERSION_MAJOR}.${SOPLEX_VERSION_MINOR}
   INTERPROCEDURAL_OPTIMIZATION ${LTO_AVAILABLE})
add_dependencies(libsoplex soplex_update_githash)
target_include_directories(libsoplex BEFORE PUBLIC
   $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
   $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
   $<INSTALL_INTERFACE:include>
   )
target_link_libraries(libsoplex
   $<$<BOOL:${SOPLEX_HAS_ZLIB}>:ZLIB::ZLIB>
   $<$<BOOL:${SOPLEX_HAS_BOOST}>:Boost::boost>
   $<$<BOOL:${SOPLEX_HAS_GMP}>:GMP::GMP>
   $<$<BOOL:${SOPLEX_HAS_QUADMATH}>:Quadmath::Quadmath>
   $<$<BOOL:${SOPLEX_HAS_PAPILO}>:PAPILO::PAPILO>
   $<$<BOOL:${SOPLEX_HAS_MPFR}>:MPFR::MPFR>
   )

# The SOPLEX_BUILD_SHARED_LIB and SOPLEX_BUILD_STATIC_LIB macros indicate that the shared / static
# library is built.
if(BUILD_SHARED_LIBS)
   target_compile_definitions(libsoplex PRIVATE SOPLEX_BUILD_SHARED_LIB)
endif()

if(MSVC)
   target_compile_options(libsoplex PUBLIC "/Zc:__cplusplus")
endif()

# Create the SoPlex binary.
add_executable(soplex soplexmain.cpp)
target_link_libraries(soplex LINK_PUBLIC libsoplex)
if(EMSCRIPTEN AND EMSCRIPTEN_HTML)
   set_target_properties(soplex PROPERTIES
      LINK_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/soplex_webdemo_shell.html")
   set(CMAKE_EXECUTABLE_SUFFIX ".html")
endif()
set_target_properties(soplex PROPERTIES
   INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}"
   INTERPROCEDURAL_OPTIMIZATION ${LTO_AVAILABLE})

# Create example binary.
add_executable(example EXCLUDE_FROM_ALL example.cpp)
target_link_libraries(example libsoplex)
set_target_properties(example PROPERTIES
   INTERPROCEDURAL_OPTIMIZATION ${LTO_AVAILABLE})

# Create aliases internally in case SoPlex is included via `add_subdirectory`.
add_library(SOPLEX::SOPLEX ALIAS libsoplex)
add_library(SOPLEX::libsoplex ALIAS libsoplex)
add_executable(SOPLEX::soplex ALIAS soplex)

# Add sanitizers.
if(SANITIZE)
   find_package(Sanitizers)
   add_sanitizers(libsoplex)
   add_sanitizers(soplex)
endif(SANITIZE)

# Install public header files.
foreach(header_file
      soplex.h
      soplex.hpp
      soplex_interface.h
      soplex/array.h
      soplex/basevectors.h
      soplex/changesoplex.hpp
      soplex/classarray.h
      soplex/classset.h
      soplex/clufactor.h
      soplex/clufactor.hpp
      soplex/clufactor_rational.h
      soplex/clufactor_rational.hpp
      soplex/cring.h
      soplex/dataarray.h
      soplex/datahashtable.h
      soplex/datakey.h
      soplex/dataset.h
      soplex/didxset.h
      soplex/dsvectorbase.h
      soplex/dsvector.h
      soplex/enter.hpp
      soplex/exceptions.h
      soplex/fmt.hpp
      soplex/idlist.h
      soplex/idxset.h
      soplex/islist.h
      soplex/leave.hpp
      soplex/lpcolbase.h
      soplex/lpcol.h
      soplex/lpcolsetbase.h
      soplex/lpcolset.h
      soplex/lprowbase.h
      soplex/lprow.h
      soplex/lprowsetbase.h
      soplex/lprowset.h
      soplex/mpsinput.h
      soplex/nameset.h
      soplex/notimer.h
      soplex/random.h
      soplex/rational.h
      soplex/ratrecon.h
      soplex/ratrecon.hpp
      soplex/slinsolver.h
      soplex/slinsolver_rational.h
      soplex/slufactor.h
      soplex/slufactor.hpp
      soplex/slufactor_rational.h
      soplex/slufactor_rational.hpp
      soplex/solbase.h
      soplex/sol.h
      soplex/solverational.hpp
      soplex/solvereal.hpp
      soplex/sorter.h
      soplex/spxalloc.h
      soplex/spxautopr.h
      soplex/spxautopr.hpp
      soplex/spxbasis.h
      soplex/spxbasis.hpp
      soplex/spxboundflippingrt.h
      soplex/spxboundflippingrt.hpp
      soplex/spxbounds.hpp
      soplex/spxchangebasis.hpp
      soplex/spxdantzigpr.h
      soplex/spxdantzigpr.hpp
      soplex/spxdefaultrt.h
      soplex/spxdefaultrt.hpp
      soplex/spxdefines.h
      soplex/spxdefines.hpp
      soplex/spxdesc.hpp
      soplex/spxdevexpr.h
      soplex/spxdevexpr.hpp
      soplex/spxequilisc.h
      soplex/spxequilisc.hpp
      soplex/spxfastrt.h
      soplex/spxfastrt.hpp
      soplex/spxfileio.h
      soplex/spxfileio.hpp
      soplex/spxgeometsc.h
      soplex/spxgeometsc.hpp
      soplex/spxgithash.h
      soplex/spxharrisrt.h
      soplex/spxharrisrt.hpp
      soplex/spxhybridpr.h
      soplex/spxhybridpr.hpp
      soplex/spxid.h
      soplex/spxleastsqsc.h
      soplex/spxleastsqsc.hpp
      soplex/spxlpbase.h
      soplex/spxlpbase_rational.hpp
      soplex/spxlpbase_real.hpp
      soplex/spxlp.h
      soplex/spxmainsm.h
      soplex/spxmainsm.hpp
      soplex/spxout.h
      soplex/spxparmultpr.h
      soplex/spxparmultpr.hpp
      soplex/spxpapilo.h
      soplex/spxpricer.h
      soplex/spxquality.hpp
      soplex/spxratiotester.h
      soplex/spxscaler.h
      soplex/spxscaler.hpp
      soplex/spxshift.hpp
      soplex/spxsimplifier.h
      soplex/spxsolve.hpp
      soplex/spxsolver.h
      soplex/spxsolver.hpp
      soplex/spxstarter.h
      soplex/spxstarter.hpp
      soplex/spxsteepexpr.h
      soplex/spxsteeppr.h
      soplex/spxsteeppr.hpp
      soplex/spxsumst.h
      soplex/spxsumst.hpp
      soplex/spxvecs.hpp
      soplex/spxvectorst.h
      soplex/spxvectorst.hpp
      soplex/spxweightpr.h
      soplex/spxweightpr.hpp
      soplex/spxweightst.h
      soplex/spxweightst.hpp
      soplex/spxwritestate.hpp
      soplex/ssvectorbase.h
      soplex/ssvector.h
      soplex/stablesum.h
      soplex/statistics.h
      soplex/statistics.hpp
      soplex/svectorbase.h
      soplex/svector.h
      soplex/svsetbase.h
      soplex/svset.h
      soplex/timerfactory.h
      soplex/timer.h
      soplex/testsoplex.hpp
      soplex/unitvectorbase.h
      soplex/unitvector.h
      soplex/updatevector.h
      soplex/updatevector.hpp
      soplex/usertimer.h
      soplex/validation.h
      soplex/validation.hpp
      soplex/vectorbase.h
      soplex/vector.h
      soplex/dvector.h
      soplex/wallclocktimer.h
   )
   get_filename_component(directory "${header_file}" DIRECTORY)
   install(
      FILES "${PROJECT_SOURCE_DIR}/src/${header_file}"
      DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${directory}")
endforeach()

# Install fmt.
if(NOT SOPLEX_HAS_SYSTEM_FMT)
   foreach(header_file
         soplex/external/fmt/chrono.h
         soplex/external/fmt/color.h
         soplex/external/fmt/compile.h
         soplex/external/fmt/core.h
         soplex/external/fmt/format-inl.h
         soplex/external/fmt/format.h
         soplex/external/fmt/locale.h
         soplex/external/fmt/os.h
         soplex/external/fmt/ostream.h
         soplex/external/fmt/posix.h
         soplex/external/fmt/printf.h
         soplex/external/fmt/ranges.h
      )
      get_filename_component(directory "${header_file}" DIRECTORY)
      install(
         FILES "${PROJECT_SOURCE_DIR}/src/${header_file}"
         DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${directory}")
   endforeach()
endif()

# Install config.h
install(
   FILES "${PROJECT_BINARY_DIR}/soplex/config.h"
   DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/soplex/")

# Install zstr headers.
if(SOPLEX_HAS_ZLIB)
   install(FILES
      soplex/external/zstr/zstr.hpp
      soplex/external/zstr/strict_fstream.hpp
      DESTINATION "include/soplex/external/zstr")
   install(FILES ${PROJECT_SOURCE_DIR}/src/soplex/external/zstr/License.txt
      DESTINATION "${CMAKE_INSTALL_DATADIR}/licenses/soplex/zstr")
endif()

# Install targets, and record them in soplex-targets export set.
install(TARGETS soplex libsoplex
   EXPORT soplex-targets
   LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
   ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
   RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
   INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
   )

# Install the debug program database files.
if(MSVC)
   install(FILES $<TARGET_PDB_FILE:soplex> DESTINATION ${CMAKE_INSTALL_BINDIR} OPTIONAL)
endif()

