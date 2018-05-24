include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(sources
    soplex/changesoplex.cpp
    soplex/clufactor.cpp
    soplex/clufactor_rational.cpp
    soplex/didxset.cpp
    soplex/enter.cpp
    soplex/gzstream.cpp
    soplex/idxset.cpp
    soplex/leave.cpp
    soplex/mpsinput.cpp
    soplex/nameset.cpp
    soplex/rational.cpp
    soplex/ratrecon.cpp
    soplex/slufactor.cpp
    soplex/slufactor_rational.cpp
    soplex/solvedbds.cpp
    soplex/solverational.cpp
    soplex/solvereal.cpp
    soplex/spxautopr.cpp
    soplex/spxbasis.cpp
    soplex/spxboundflippingrt.cpp
    soplex/spxbounds.cpp
    soplex/spxchangebasis.cpp
    soplex/spxdantzigpr.cpp
    soplex/spxdefaultrt.cpp
    soplex/spxdefines.cpp
    soplex/spxdesc.cpp
    soplex/spxdevexpr.cpp
    soplex/spxequilisc.cpp
    soplex/spxfastrt.cpp
    soplex/spxfileio.cpp
    soplex/spxgeometsc.cpp
    soplex/spxgithash.cpp
    soplex/spxharrisrt.cpp
    soplex/spxhybridpr.cpp
    soplex/spxid.cpp
    soplex/spxleastsqsc.cpp
    soplex/spxlpbase_rational.cpp
    soplex/spxlpbase_real.cpp
    soplex/spxmainsm.cpp
    soplex/spxout.cpp
    soplex/spxparmultpr.cpp
    soplex/spxquality.cpp
    soplex/spxscaler.cpp
    soplex/spxshift.cpp
    soplex/spxsolve.cpp
    soplex/spxsolver.cpp
    soplex/spxstarter.cpp
    soplex/spxsteeppr.cpp
    soplex/spxsumst.cpp
    soplex/spxvecs.cpp
    soplex/spxvectorst.cpp
    soplex/spxweightpr.cpp
    soplex/spxweightst.cpp
    soplex/spxwritestate.cpp
    soplex/statistics.cpp
    soplex/testsoplex.cpp
    soplex/updatevector.cpp
    soplex/usertimer.cpp
    soplex/validation.cpp
    soplex/wallclocktimer.cpp
    soplex.cpp)

set(headers
    soplex/array.h
    soplex/basevectors.h
    soplex/classarray.h
    soplex/clufactor.h
    soplex/clufactor_rational.h
    soplex/cring.h
    soplex/dataarray.h
    soplex/datahashtable.h
    soplex/datakey.h
    soplex/dataset.h
    soplex/didxset.h
    soplex/dsvectorbase.h
    soplex/dsvector.h
    soplex/dvectorbase.h
    soplex/dvector.h
    soplex/exceptions.h
    soplex/gzstream.h
    soplex/idlist.h
    soplex/idxset.h
    soplex/islist.h
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
    soplex/slinsolver.h
    soplex/slinsolver_rational.h
    soplex/slufactor.h
    soplex/slufactor_rational.h
    soplex/solbase.h
    soplex/sol.h
    soplex/sorter.h
    soplex/spxalloc.h
    soplex/spxautopr.h
    soplex/spxbasis.h
    soplex/spxboundflippingrt.h
    soplex/spxdantzigpr.h
    soplex/spxdefaultrt.h
    soplex/spxdefines.h
    soplex/spxdevexpr.h
    soplex/spxequilisc.h
    soplex/spxfastrt.h
    soplex/spxfileio.h
    soplex/spxgeometsc.h
    soplex/spxgithash.h
    soplex/spxharrisrt.h
    soplex/spxhybridpr.h
    soplex/spxid.h
    soplex/spxleastsqsc.h
    soplex/spxlpbase.h
    soplex/spxlp.h
    soplex/spxmainsm.h
    soplex/spxout.h
    soplex/spxparmultpr.h
    soplex/spxpricer.h
    soplex/spxratiotester.h
    soplex/spxscaler.h
    soplex/spxsimplifier.h
    soplex/spxsolver.h
    soplex/spxstarter.h
    soplex/spxsteepexpr.h
    soplex/spxsteeppr.h
    soplex/spxsumst.h
    soplex/spxvectorst.h
    soplex/spxweightpr.h
    soplex/spxweightst.h
    soplex/ssvectorbase.h
    soplex/ssvector.h
    soplex/statistics.h
    soplex/svectorbase.h
    soplex/svector.h
    soplex/svsetbase.h
    soplex/svset.h
    soplex/timerfactory.h
    soplex/timer.h
    soplex/unitvectorbase.h
    soplex/unitvector.h
    soplex/updatevector.h
    soplex/usertimer.h
    soplex/validation.h
    soplex/vectorbase.h
    soplex/vector.h
    soplex/wallclocktimer.h)

    # create soplex library with pic
    add_library(libsoplex-pic STATIC ${sources})
    set_target_properties(libsoplex-pic PROPERTIES
    POSITION_INDEPENDENT_CODE on
    VERSION ${SOPLEX_VERSION_MAJOR}.${SOPLEX_VERSION_MINOR}.${SOPLEX_VERSION_PATCH}.${SOPLEX_VERSION_SUB}
    SOVERSION ${SOPLEX_VERSION_MAJOR}.${SOPLEX_VERSION_MINOR})
add_dependencies(libsoplex-pic soplex_update_githash)
target_link_libraries(libsoplex-pic ${libs})

# create soplex library without pic
add_library(libsoplex STATIC ${sources})
set_target_properties(libsoplex PROPERTIES
    VERSION ${SOPLEX_VERSION_MAJOR}.${SOPLEX_VERSION_MINOR}.${SOPLEX_VERSION_PATCH}.${SOPLEX_VERSION_SUB}
    SOVERSION ${SOPLEX_VERSION_MAJOR}.${SOPLEX_VERSION_MINOR})
add_dependencies(libsoplex soplex_update_githash)
target_link_libraries(libsoplex ${libs})

# create soplex binary using library without pic
add_executable(soplex soplexmain.cpp)
target_link_libraries(soplex libsoplex)

if(CMAKE_BUILD_TYPE EQUAL "Debug")
    find_package(Sanitizers)
    add_sanitizers(soplex)
endif()

add_executable(example EXCLUDE_FROM_ALL example.cpp)
target_link_libraries(example libsoplex)

# set the install rpath to the installed destination
set_target_properties(soplex PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
# on UNIX system the 'lib' prefix is automatically added
set_target_properties(libsoplex PROPERTIES
    OUTPUT_NAME "soplex"
    MACOSX_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set_target_properties(libsoplex-pic PROPERTIES
    OUTPUT_NAME "soplex-pic"
    MACOSX_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# install the header files of soplex
install(FILES ${headers} DESTINATION include/soplex)
install(FILES soplex.h DESTINATION include)

# install the binary and the library to appropriate lcoations and add them to an export group
install(TARGETS soplex libsoplex libsoplex-pic EXPORT soplex-targets
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
        INCLUDES DESTINATION include)

# Add library targets to the build-tree export set
export(TARGETS libsoplex libsoplex-pic
  FILE "${CMAKE_BINARY_DIR}/soplex-targets.cmake")

#configure the config file for the build tree
set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/src")
configure_file(${PROJECT_SOURCE_DIR}/soplex-config.cmake.in
  "${CMAKE_BINARY_DIR}/soplex-config.cmake" @ONLY)

#configure the config file for the install
set(CONF_INCLUDE_DIRS "\${CMAKE_CURRENT_LIST_DIR}/../../../include")
configure_file(${PROJECT_SOURCE_DIR}/soplex-config.cmake.in
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/soplex-config.cmake" @ONLY)

# install the targets of the soplex export group and the config file so that other projects
# can link easily against soplex
install(EXPORT soplex-targets FILE soplex-targets.cmake DESTINATION lib/cmake/soplex)
install(FILES "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/soplex-config.cmake" DESTINATION lib/cmake/soplex)
