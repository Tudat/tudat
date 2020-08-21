#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

set(ESTIMATIONSETUPSDIR "${SIMULATIONSETUPDIR}/EstimationSetup")
set(ENVIRONMENTSETUPDIR "${SIMULATIONSETUPDIR}/EnvironmentSetup")
set(PROPAGATIONSETUPSDIR "${SIMULATIONSETUPDIR}/PropagationSetup")

# Add subdirectories.
if( BUILD_WITH_ESTIMATION_TOOLS )
    add_subdirectory("${SRCROOT}${ESTIMATIONSETUPSDIR}")
endif( )
add_subdirectory("${SRCROOT}${ENVIRONMENTSETUPDIR}")
add_subdirectory("${SRCROOT}${PROPAGATIONSETUPSDIR}")

# Get target properties for static libraries.
get_target_property(SIMULATION_ENVIRONMENT_SETUP_SOURCES tudat_environment_setup SOURCES)
if( BUILD_WITH_ESTIMATION_TOOLS )
    get_target_property(SIMULATION_ESTIMATION_SETUP_SOURCES tudat_estimation_setup SOURCES)
endif( )
get_target_property(SIMULATION_PROPAGATION_SETUP_SOURCES tudat_propagation_setup SOURCES)

# Add unit tests.
if(USE_CSPICE)

    tudat_add_test(
      TARGET test_EnvironmentCreation
      SOURCES "${SRCROOT}${SIMULATIONSETUPDIR}/UnitTests/unitTestEnvironmentModelSetup.cpp"
      LINK_LIBRARIES ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

    tudat_add_test(
      TARGET test_AccelerationModelCreation
      SOURCES "${SRCROOT}${SIMULATIONSETUPDIR}/UnitTests/unitTestAccelerationModelSetup.cpp"
      LINK_LIBRARIES ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

endif()
