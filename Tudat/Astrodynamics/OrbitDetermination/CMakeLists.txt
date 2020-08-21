 #    Copyright (c) 2010-2019, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.

set(ACCELERATIONPARTIALSDIR "${ORBITDETERMINATIONDIR}/AccelerationPartials")
set(TORQUEPARTIALSDIR "${ORBITDETERMINATIONDIR}/RotationalDynamicsPartials")
set(ESTIMATABLEPARAMETERSDIR "${ORBITDETERMINATIONDIR}/EstimatableParameters")
set(OBSERVATIONPARTIALSDIR "${ORBITDETERMINATIONDIR}/ObservationPartials")
set(LIGHTTIMECORRECTIONPARTIALSDIR "${ORBITDETERMINATIONDIR}/LightTimeCorrectionPartials")

# Add subdirectories.
add_subdirectory("${SRCROOT}${ACCELERATIONPARTIALSDIR}")
add_subdirectory("${SRCROOT}${TORQUEPARTIALSDIR}")
add_subdirectory("${SRCROOT}${ESTIMATABLEPARAMETERSDIR}")
add_subdirectory("${SRCROOT}${OBSERVATIONPARTIALSDIR}")

# Get target properties for static libraries.
get_target_property(ACCELERATIONPARTIALSSOURCES tudat_acceleration_partials SOURCES)
get_target_property(TORQUEPARTIALSSOURCES tudat_torque_partials SOURCES)
get_target_property(ESTIMATABLEPARAMETERSOURCES tudat_estimatable_parameters SOURCES)
get_target_property(OBSERVATIONPARTIALSSOURCES tudat_observation_partials SOURCES)

# Set the source files.
set(ORBIT_DETERMINATION_SOURCES
  "${SRCROOT}${ORBITDETERMINATIONDIR}/stateDerivativePartial.cpp"
  "${SRCROOT}${ORBITDETERMINATIONDIR}/podInputOutputTypes.cpp"
)

# Set the header files.
set(ORBIT_DETERMINATION_HEADERS
  "${SRCROOT}${ORBITDETERMINATIONDIR}/stateDerivativePartial.h"
  "${SRCROOT}${ORBITDETERMINATIONDIR}/podInputOutputTypes.h"
)


# Add static libraries.
add_library(tudat_orbit_determination STATIC ${ORBIT_DETERMINATION_SOURCES} ${ORBIT_DETERMINATION_HEADERS})
setup_tudat_library_target(tudat_orbit_determination "${SRCROOT}{ORBITDETERMINATIONDIR}")

# Add unit tests.
set(ORBDET_TEST_DIR "${SRCROOT}${ORBITDETERMINATIONDIR}/UnitTests")
if( BUILD_WITH_PROPAGATION_TESTS )

  tudat_add_test(
    TARGET test_EstimationFromPositionDoubleDouble
    SOURCES "${ORBDET_TEST_DIR}/unitTestEstimationFromIdealDataDoubleDouble.cpp"
    LINK_LIBRARIES ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

  tudat_add_test(
    TARGET test_EstimationInput
    SOURCES "${ORBDET_TEST_DIR}/unitTestEstimationInput.cpp"
    LINK_LIBRARIES ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

  tudat_add_test(
    TARGET test_TidalPropertyEstimation
    SOURCES "${ORBDET_TEST_DIR}/unitTestTidalPropertyEstimation.cpp"
    LINK_LIBRARIES ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

  tudat_add_test(
    TARGET test_MultiArcStateEstimation
    SOURCES "${ORBDET_TEST_DIR}/unitTestMultiArcStateEstimation.cpp"
    LINK_LIBRARIES ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

  tudat_add_test(
    TARGET test_HybridArcStateEstimation
    SOURCES "${ORBDET_TEST_DIR}/unitTestHybridArcStateEstimation.cpp"
    LINK_LIBRARIES ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

  tudat_add_test(
    TARGET test_ParameterInfluenceDetermination
    SOURCES "${ORBDET_TEST_DIR}/unitTestParameterInfluenceDetermination.cpp"
    LINK_LIBRARIES ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

  tudat_add_test(
    TARGET test_BiasEstimation
    SOURCES "${ORBDET_TEST_DIR}/unitTestBiasEstimation.cpp"
    LINK_LIBRARIES ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

  tudat_add_test(
    TARGET test_RotationalStateEstimation
    SOURCES "${ORBDET_TEST_DIR}/unitTestRotationalStateEstimation.cpp"
    LINK_LIBRARIES ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

  tudat_add_test(
    TARGET test_ArcwiseEnvironmentEstimation
    SOURCES "${ORBDET_TEST_DIR}/unitTestArcwiseEnvironmentParameters.cpp"
    LINK_LIBRARIES ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

  tudat_add_test(
    TARGET test_DesaturationDeltaVsEstimation
    SOURCES "${ORBDET_TEST_DIR}/unitTestDesaturationDeltaVsEstimation.cpp"
    LINK_LIBRARIES ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

  tudat_add_test(
    TARGET test_FullPlanetaryRotationModelEstimation
    SOURCES "${ORBDET_TEST_DIR}/unitTestFullPlanetaryRotationModelEstimation.cpp"
    LINK_LIBRARIES ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

  if( BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS )

    tudat_add_test(
      TARGET test_EstimationFromPositionDoubleLongDouble
      SOURCES "${ORBDET_TEST_DIR}/unitTestEstimationFromIdealDataDoubleLongDouble.cpp"
      LINK_LIBRARIES ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

    tudat_add_test(
      TARGET test_EstimationFromPositionTimeDouble
      SOURCES "${ORBDET_TEST_DIR}/unitTestEstimationFromIdealDataTimeDouble.cpp"
      LINK_LIBRARIES ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

    tudat_add_test(
      TARGET test_EstimationFromPositionTimeLongDouble
      SOURCES "${ORBDET_TEST_DIR}/unitTestEstimationFromIdealDataTimeLongDouble.cpp"
      LINK_LIBRARIES ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

  endif( )

endif( )

