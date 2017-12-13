set(sources
    changesoplex.cpp
    clufactor.cpp
    clufactor_rational.cpp
    didxset.cpp
    enter.cpp
    gzstream.cpp
    idxset.cpp
    leave.cpp
    mpsinput.cpp
    nameset.cpp
    rational.cpp
    ratrecon.cpp
    slufactor.cpp
    slufactor_rational.cpp
    solvedbds.cpp
    solverational.cpp
    solvereal.cpp
    soplex.cpp
    soplexlegacy.cpp
    spxautopr.cpp
    spxbasis.cpp
    spxboundflippingrt.cpp
    spxbounds.cpp
    spxchangebasis.cpp
    spxdantzigpr.cpp
    spxdefaultrt.cpp
    spxdefines.cpp
    spxdesc.cpp
    spxdevexpr.cpp
    spxequilisc.cpp
    spxfastrt.cpp
    spxfileio.cpp
    spxgeometsc.cpp
    spxgithash.cpp
    spxharrisrt.cpp
    spxhybridpr.cpp
    spxid.cpp
    spxleastsqsc.cpp
    spxlpbase_rational.cpp
    spxlpbase_real.cpp
    spxmainsm.cpp
    spxout.cpp
    spxparmultpr.cpp
    spxquality.cpp
    spxscaler.cpp
    spxshift.cpp
    spxsolve.cpp
    spxsolver.cpp
    spxstarter.cpp
    spxsteeppr.cpp
    spxsumst.cpp
    spxvecs.cpp
    spxvectorst.cpp
    spxweightpr.cpp
    spxweightst.cpp
    spxwritestate.cpp
    statistics.cpp
    testsoplex.cpp
    updatevector.cpp
    usertimer.cpp
    validation.cpp
    wallclocktimer.cpp)

set(headers
    array.h
    basevectors.h
    classarray.h
    clufactor.h
    clufactor_rational.h
    cring.h
    dataarray.h
    datahashtable.h
    datakey.h
    dataset.h
    didxset.h
    dsvectorbase.h
    dsvector.h
    dvectorbase.h
    dvector.h
    exceptions.h
    gzstream.h
    idlist.h
    idxset.h
    islist.h
    lpcolbase.h
    lpcol.h
    lpcolsetbase.h
    lpcolset.h
    lprowbase.h
    lprow.h
    lprowsetbase.h
    lprowset.h
    mpsinput.h
    nameset.h
    notimer.h
    random.h
    rational.h
    ratrecon.h
    slinsolver.h
    slinsolver_rational.h
    slufactor.h
    slufactor_rational.h
    solbase.h
    sol.h
    soplex.h
    soplexlegacy.h
    sorter.h
    spxalloc.h
    spxautopr.h
    spxbasis.h
    spxboundflippingrt.h
    spxdantzigpr.h
    spxdefaultrt.h
    spxdefines.h
    spxdevexpr.h
    spxequilisc.h
    spxfastrt.h
    spxfileio.h
    spxgeometsc.h
    spxgithash.h
    spxharrisrt.h
    spxhybridpr.h
    spxid.h
    spxleastsqsc.h
    spxlpbase.h
    spxlp.h
    spxmainsm.h
    spxout.h
    spxparmultpr.h
    spxpricer.h
    spxratiotester.h
    spxscaler.h
    spxsimplifier.h
    spxsolver.h
    spxstarter.h
    spxsteepexpr.h
    spxsteeppr.h
    spxsumst.h
    spxvectorst.h
    spxweightpr.h
    spxweightst.h
    ssvectorbase.h
    ssvector.h
    statistics.h
    svectorbase.h
    svector.h
    svsetbase.h
    svset.h
    timerfactory.h
    timer.h
    unitvectorbase.h
    unitvector.h
    updatevector.h
    usertimer.h
    validation.h
    vectorbase.h
    vector.h
    wallclocktimer.h)

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

if(NOT LEGACY EQUAL ON)
    add_executable(example EXCLUDE_FROM_ALL example.cpp)
    target_link_libraries(example libsoplex)
endif()

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
install(FILES ${headers} DESTINATION include)

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
