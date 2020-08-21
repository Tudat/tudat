 #    Copyright (c) 2010-2019, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.
 #

# Set the source files.
set(JSONINTERFACE_SOURCES
    "${SRCROOT}${JSONINTERFACEDIR}/Support/keys.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Support/deserialization.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Support/path.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Support/valueAccess.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/spice.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/body.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/atmosphere.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/ephemeris.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/aerodynamics.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/gravityField.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/rotationModel.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/shapeModel.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/radiationPressure.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/gravityFieldVariation.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/groundStations.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Propagation/variable.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Propagation/termination.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Propagation/acceleration.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Propagation/massRateModel.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Propagation/thrust.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Propagation/torque.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Mathematics/interpolation.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Propagation/export.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Support/options.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/jsonInterface.cpp"
)

# Set the header files.
set(JSONINTERFACE_HEADERS
    "${SRCROOT}${JSONINTERFACEDIR}/Support/errorHandling.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Support/keys.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Support/deserialization.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Support/path.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Support/utilities.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Support/valueAccess.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Support/valueConversions.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/spice.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/body.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/atmosphere.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/ephemeris.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/aerodynamics.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/gravityField.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/rotationModel.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/shapeModel.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/radiationPressure.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/gravityFieldVariation.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/groundStations.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Propagation/variable.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Propagation/state.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Propagation/propagator.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Propagation/termination.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Propagation/referenceFrames.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Propagation/acceleration.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Propagation/thrust.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Propagation/massRateModel.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Propagation/torque.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Mathematics/integrator.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Mathematics/interpolation.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Propagation/export.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Support/options.h"
    "${SRCROOT}${JSONINTERFACEDIR}/UnitTests/unitTestSupport.h"
    "${SRCROOT}${JSONINTERFACEDIR}/jsonInterface.h"
)

if( BUILD_WITH_ESTIMATION_TOOLS )
    set(JSONINTERFACE_SOURCES
        ${JSONINTERFACE_SOURCES}
        "${SRCROOT}${JSONINTERFACEDIR}/Estimation/parameter.cpp"
        "${SRCROOT}${JSONINTERFACEDIR}/Estimation/observation.cpp"
#        "${SRCROOT}${JSONINTERFACEDIR}/Estimation/orbitDetermination.cpp"
        "${SRCROOT}${JSONINTERFACEDIR}/jsonInterfaceVariational.cpp"
#        "${SRCROOT}${JSONINTERFACEDIR}/jsonEstimationInterface.cpp"
    )

    # Set the header files.
    set(JSONINTERFACE_HEADERS
        ${JSONINTERFACE_HEADERS}
        "${SRCROOT}${JSONINTERFACEDIR}/Estimation/parameter.h"
        "${SRCROOT}${JSONINTERFACEDIR}/Estimation/observation.h"
#        "${SRCROOT}${JSONINTERFACEDIR}/Estimation/orbitDetermination.h"
        "${SRCROOT}${JSONINTERFACEDIR}/jsonInterfaceVariational.h"
#        "${SRCROOT}${JSONINTERFACEDIR}/jsonEstimationInterface.h"
    )
endif( )

# Add static libraries.
add_library(json_interface_library STATIC ${JSONINTERFACE_SOURCES} ${JSONINTERFACE_HEADERS})
setup_tudat_library_target(json_interface_library "")

if( BUILD_WITH_ESTIMATION_TOOLS )
   list(APPEND JSON_PROPAGATION_LIBRARIES json_interface_library ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})
else( )
   list(APPEND JSON_PROPAGATION_LIBRARIES json_interface_library ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})
endif( )

# Set the global macros for setting up targets.
macro(setup_executable_target target_name)
  set_property(TARGET ${target_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY "${BINROOT}")
  install(TARGETS ${target_name} RUNTIME DESTINATION "${BINROOT}")
endmacro()


# Add JSON-based application.

add_executable(json_interface "${SRCROOT}${JSONINTERFACEDIR}/jsonInterfaceExecute.cpp")
setup_executable_target(json_interface)
target_link_libraries(json_interface ${JSON_PROPAGATION_LIBRARIES})

# Add unit tests.
set(JSON_TESTS_DIR "${SRCROOT}${JSONINTERFACEDIR}/UnitTests")

# Acceleration
tudat_add_test(
  TARGET test_JsonInterfaceAcceleration
  SOURCES "${JSON_TESTS_DIR}/unitTestAcceleration.cpp"
  LINK_LIBRARIES ${JSON_PROPAGATION_LIBRARIES})

# Aerodynamics
tudat_add_test(
  TARGET test_JsonInterfaceAerodynamics
  SOURCES "${JSON_TESTS_DIR}/unitTestAerodynamics.cpp"
  LINK_LIBRARIES ${JSON_PROPAGATION_LIBRARIES})

# Atmosphere
tudat_add_test(
  TARGET test_JsonInterfaceAtmosphere
  SOURCES "${JSON_TESTS_DIR}/unitTestAtmosphere.cpp"
  LINK_LIBRARIES ${JSON_PROPAGATION_LIBRARIES})

# Body
tudat_add_test(
  TARGET test_JsonInterfaceBody
  SOURCES "${JSON_TESTS_DIR}/unitTestBody.cpp"
  LINK_LIBRARIES ${JSON_PROPAGATION_LIBRARIES})

# Deserialization
tudat_add_test(
  TARGET test_JsonInterfaceDeserialization
  SOURCES "${JSON_TESTS_DIR}/unitTestDeserialization.cpp"
  LINK_LIBRARIES ${JSON_PROPAGATION_LIBRARIES})

# Ephemeris
tudat_add_test(
  TARGET test_JsonInterfaceEphemeris
  SOURCES "${JSON_TESTS_DIR}/unitTestEphemeris.cpp"
  LINK_LIBRARIES ${JSON_PROPAGATION_LIBRARIES})

# Export
tudat_add_test(
  TARGET test_JsonInterfaceExport
  SOURCES "${JSON_TESTS_DIR}/unitTestExport.cpp"
  LINK_LIBRARIES ${JSON_PROPAGATION_LIBRARIES})

# GravityField
tudat_add_test(
  TARGET test_JsonInterfaceGravityField
  SOURCES "${JSON_TESTS_DIR}/unitTestGravityField.cpp"
  LINK_LIBRARIES ${JSON_PROPAGATION_LIBRARIES})

# GravityFieldVariation
tudat_add_test(
  TARGET test_JsonInterfaceGravityFieldVariation
  SOURCES "${JSON_TESTS_DIR}/unitTestGravityFieldVariation.cpp"
  LINK_LIBRARIES ${JSON_PROPAGATION_LIBRARIES})

# GroundStation
tudat_add_test(
  TARGET test_JsonInterfaceGroundStation
  SOURCES "${JSON_TESTS_DIR}/unitTestGroundStation.cpp"
  LINK_LIBRARIES ${JSON_PROPAGATION_LIBRARIES})

# Integrator
tudat_add_test(
  TARGET test_JsonInterfaceIntegrator
  SOURCES "${JSON_TESTS_DIR}/unitTestIntegrator.cpp"
  LINK_LIBRARIES ${JSON_PROPAGATION_LIBRARIES})

# Interpolation
tudat_add_test(
  TARGET test_JsonInterfaceInterpolation
  SOURCES "${JSON_TESTS_DIR}/unitTestInterpolation.cpp"
  LINK_LIBRARIES ${JSON_PROPAGATION_LIBRARIES})

# MassRateModel
tudat_add_test(
  TARGET test_JsonInterfaceMassRateModel
  SOURCES "${JSON_TESTS_DIR}/unitTestMassRateModel.cpp"
  LINK_LIBRARIES ${JSON_PROPAGATION_LIBRARIES})

# Propagator
tudat_add_test(
  TARGET test_JsonInterfacePropagator
  SOURCES "${JSON_TESTS_DIR}/unitTestPropagator.cpp"
  LINK_LIBRARIES ${JSON_PROPAGATION_LIBRARIES})

# RadiationPressure
tudat_add_test(
  TARGET test_JsonInterfaceRadiationPressure
  SOURCES "${JSON_TESTS_DIR}/unitTestRadiationPressure.cpp"
  LINK_LIBRARIES ${JSON_PROPAGATION_LIBRARIES})

# ReferenceFrames
tudat_add_test(
  TARGET test_JsonInterfaceReferenceFrames
  SOURCES "${JSON_TESTS_DIR}/unitTestReferenceFrames.cpp"
  LINK_LIBRARIES ${JSON_PROPAGATION_LIBRARIES})

# RotationModel
tudat_add_test(
  TARGET test_JsonInterfaceRotationModel
  SOURCES "${JSON_TESTS_DIR}/unitTestRotationModel.cpp"
  LINK_LIBRARIES ${JSON_PROPAGATION_LIBRARIES})

# ShapeModel
tudat_add_test(
  TARGET test_JsonInterfaceShapeModel
  SOURCES "${JSON_TESTS_DIR}/unitTestShapeModel.cpp"
  LINK_LIBRARIES ${JSON_PROPAGATION_LIBRARIES})

# SimulationSingleSatellite
tudat_add_test(
  TARGET test_JsonInterfaceSimulationSingleSatellite
  SOURCES "${JSON_TESTS_DIR}/unitTestSimulationSingleSatellite.cpp"
  LINK_LIBRARIES ${JSON_PROPAGATION_LIBRARIES})

# SimulationSinglePerturbedSatellite
tudat_add_test(
  TARGET test_JsonInterfaceSimulationSinglePerturbedSatellite
  SOURCES "${JSON_TESTS_DIR}/unitTestSimulationSinglePerturbedSatellite.cpp"
  LINK_LIBRARIES ${JSON_PROPAGATION_LIBRARIES})

# SimulationInnerSolarSystem
tudat_add_test(
  TARGET test_JsonInterfaceSimulationInnerSolarSystem
  SOURCES "${JSON_TESTS_DIR}/unitTestSimulationInnerSolarSystem.cpp"
  LINK_LIBRARIES ${JSON_PROPAGATION_LIBRARIES})

# SimulationGalileoConstellation
tudat_add_test(
  TARGET test_JsonInterfaceSimulationGalileoConstellation
  SOURCES "${JSON_TESTS_DIR}/unitTestSimulationGalileoConstellation.cpp"
  LINK_LIBRARIES ${JSON_PROPAGATION_LIBRARIES})

# SimulationThrustAlongVelocityVector
tudat_add_test(
  TARGET test_JsonInterfaceSimulationThrustAlongVelocityVector
  SOURCES "${JSON_TESTS_DIR}/unitTestSimulationThrustAlongVelocityVector.cpp"
  LINK_LIBRARIES ${JSON_PROPAGATION_LIBRARIES})

# SimulationThrustAccelerationFromFile
tudat_add_test(
  TARGET test_JsonInterfaceSimulationThrustAccelerationFromFile
  SOURCES "${JSON_TESTS_DIR}/unitTestSimulationThrustAccelerationFromFile.cpp"
  LINK_LIBRARIES ${JSON_PROPAGATION_LIBRARIES})

# Spice
tudat_add_test(
  TARGET test_JsonInterfaceSpice
  SOURCES "${JSON_TESTS_DIR}/unitTestSpice.cpp"
  LINK_LIBRARIES ${JSON_PROPAGATION_LIBRARIES})

# State
tudat_add_test(
  TARGET test_JsonInterfaceState
  SOURCES "${JSON_TESTS_DIR}/unitTestState.cpp"
  LINK_LIBRARIES ${JSON_PROPAGATION_LIBRARIES})

# Termination
tudat_add_test(
  TARGET test_JsonInterfaceTermination
  SOURCES "${JSON_TESTS_DIR}/unitTestTermination.cpp"
  LINK_LIBRARIES ${JSON_PROPAGATION_LIBRARIES})

# Thrust
tudat_add_test(
  TARGET test_JsonInterfaceThrust
  SOURCES "${JSON_TESTS_DIR}/unitTestThrust.cpp"
  LINK_LIBRARIES ${JSON_PROPAGATION_LIBRARIES})

# Torque
tudat_add_test(
  TARGET test_JsonInterfaceTorque
  SOURCES "${JSON_TESTS_DIR}/unitTestTorque.cpp"
  LINK_LIBRARIES ${JSON_PROPAGATION_LIBRARIES})

# Variable
tudat_add_test(
  TARGET test_JsonInterfaceVariable
  SOURCES "${JSON_TESTS_DIR}/unitTestVariable.cpp"
  LINK_LIBRARIES ${JSON_PROPAGATION_LIBRARIES})

if( BUILD_WITH_ESTIMATION_TOOLS )

  tudat_add_test(
    TARGET test_JsonInterfaceObservation
    SOURCES "${JSON_TESTS_DIR}/unitTestObservation.cpp"
    LINK_LIBRARIES ${JSON_PROPAGATION_LIBRARIES})

  tudat_add_test(
    TARGET test_JsonInterfaceParameter
    SOURCES "${JSON_TESTS_DIR}/unitTestParameter.cpp"
    LINK_LIBRARIES ${JSON_PROPAGATION_LIBRARIES})

  tudat_add_test(
    TARGET test_JsonInterfaceSimulationSingleSatelliteVariational
    SOURCES "${JSON_TESTS_DIR}/unitTestSimulationSingleSatelliteVariational.cpp"
    LINK_LIBRARIES ${JSON_PROPAGATION_LIBRARIES})

  #tudat_add_test(
  #  TARGET test_JsonInterfaceSimulationSingleSatelliteEstimation
  #  SOURCES "${JSON_TESTS_DIR}/unitTestSimulationSingleSatelliteEstimation.cpp"
  #  LINK_LIBRARIES ${JSON_PROPAGATION_LIBRARIES})

endif( )


