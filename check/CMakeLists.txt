include(CTest)

#
# define the set of feasible instances
#
set(feasibleInstances
    adlittle
    afiro
    agg
    beaconfd
    blend
    bore3d
    brandy
    capri
    etamacro
    finnis
    grow7
    israel
    kb2
    lotfi
    recipe
    sc105
    sc205
    sc50a
    sc50b
    scagr25
    scagr7
    scfxm1
    scorpion
    scrs8
    scsd1
    seba
    share1b
    share2b
    shell
    vtp-base
    )

#
# define the set of infeasible instances
#
set(infeasibleInstances
    bgetam
    box1
    ex72a
    forest6
    galenet
    gams10am
    klein1
    refinery
    woodinfe)

#
# define the set of unbounded instances
#
set(unboundedInstances
    gas11)

#
# define the set of settings
#
set(settings
    devex
    steep
    exact
    default
    )

#
# define a macro to add tests
#
# add_instancetests takes an instance group and a status
# that SoPlex should report as arguments
#
macro(add_instancetests instances status)
#
# loop over the instances
#
foreach(instance ${${instances}})
    #
    # loop over all settings
    #
    foreach(setting ${settings})

    #
    # specify the instance and the settings load command for soplex
    #
    set(inst "${PROJECT_SOURCE_DIR}/check/instances/${instance}.mps")
    if($setting STREQUAL "")
        set(loadCommand "")
    else()
        set(loadCommand "--loadset=${PROJECT_SOURCE_DIR}/settings/${setting}.set")
    endif()

    #
    # add a test and test properties
    #
    add_test(NAME ${instance}-${setting} COMMAND soplex ${inst} $loadCommand)
    set_tests_properties (${instance}-${setting} PROPERTIES PASS_REGULAR_EXPRESSION "problem is solved \\[${status}\\]")
    endforeach(setting)
endforeach(instance)
endmacro(add_instancetests)

#
# add tests for infeasible, feasible and unbounded instances
#
add_instancetests( feasibleInstances "optimal" )
add_instancetests( infeasibleInstances "infeasible" )
add_instancetests( unboundedInstances "unbounded" )