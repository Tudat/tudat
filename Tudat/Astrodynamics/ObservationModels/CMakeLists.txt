 #    Copyright (c) 2010-2019, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.

set(OBSERVABLECORRECTIONSDIR "${OBSERVATIONMODELSDIR}/ObservableCorrections")


# Set the source files.
set(OBSERVATION_MODELS_SOURCES
  "${SRCROOT}${OBSERVATIONMODELSDIR}/lightTimeSolution.cpp"
  "${SRCROOT}${OBSERVATIONMODELSDIR}/observableTypes.cpp"
  "${SRCROOT}${OBSERVATIONMODELSDIR}/linkTypeDefs.cpp"
  "${SRCROOT}${OBSERVATIONMODELSDIR}/observationViabilityCalculator.cpp"
  "${SRCROOT}${OBSERVATIONMODELSDIR}/observationManager.cpp"
  "${SRCROOT}${OBSERVATIONMODELSDIR}/observationSimulator.cpp"
  "${SRCROOT}${OBSERVATIONMODELSDIR}/observationModel.cpp"
  "${SRCROOT}${OBSERVATIONMODELSDIR}/simulateObservations.cpp"
  "${SRCROOT}${OBSERVABLECORRECTIONSDIR}/firstOrderRelativisticCorrection.cpp"
)

# Set the header files.
set(OBSERVATION_MODELS_HEADERS
  "${SRCROOT}${OBSERVATIONMODELSDIR}/angularPositionObservationModel.h"
  "${SRCROOT}${OBSERVATIONMODELSDIR}/lightTimeSolution.h"
  "${SRCROOT}${OBSERVATIONMODELSDIR}/linkTypeDefs.h"
  "${SRCROOT}${OBSERVATIONMODELSDIR}/observableTypes.h"
  "${SRCROOT}${OBSERVATIONMODELSDIR}/observationModel.h"
  "${SRCROOT}${OBSERVATIONMODELSDIR}/observationManager.h"
  "${SRCROOT}${OBSERVATIONMODELSDIR}/observationSimulator.h"
  "${SRCROOT}${OBSERVATIONMODELSDIR}/observationBias.h"
  "${SRCROOT}${OBSERVATIONMODELSDIR}/nWayRangeObservationModel.h"
  "${SRCROOT}${OBSERVATIONMODELSDIR}/oneWayRangeObservationModel.h"
  "${SRCROOT}${OBSERVATIONMODELSDIR}/oneWayDopplerObservationModel.h"
  "${SRCROOT}${OBSERVATIONMODELSDIR}/oneWayDifferencedRangeRateObservationModel.h"
  "${SRCROOT}${OBSERVATIONMODELSDIR}/simulateObservations.h"
  "${SRCROOT}${OBSERVATIONMODELSDIR}/positionObservationModel.h"
  "${SRCROOT}${OBSERVATIONMODELSDIR}/velocityObservationModel.h"
  "${SRCROOT}${OBSERVATIONMODELSDIR}/twoWayDopplerObservationModel.h"
  "${SRCROOT}${OBSERVATIONMODELSDIR}/eulerAngleObservationModel.h"
  "${SRCROOT}${OBSERVATIONMODELSDIR}/UnitTests/testLightTimeCorrections.h"
  "${SRCROOT}${OBSERVATIONMODELSDIR}/observationViabilityCalculator.h"
  "${SRCROOT}${OBSERVABLECORRECTIONSDIR}/lightTimeCorrection.h"
  "${SRCROOT}${OBSERVABLECORRECTIONSDIR}/firstOrderRelativisticCorrection.h"
)


# Add static libraries.
add_library(tudat_observation_models STATIC ${OBSERVATION_MODELS_SOURCES} ${OBSERVATION_MODELS_HEADERS})
setup_tudat_library_target(tudat_observation_models "${SRCROOT}${OBSERVATIONMODELSDIR}")

# Add unit tests.
set(OBSMDL_TEST_DIR "${SRCROOT}${OBSERVATIONMODELSDIR}/UnitTests")
if(USE_CSPICE)

  tudat_add_test(
    TARGET test_LightTime
    SOURCES "${OBSMDL_TEST_DIR}/unitTestLightTimeSolution.cpp"
    LINK_LIBRARIES tudat_observation_models
                   tudat_ephemerides
                   tudat_basic_astrodynamics
                   tudat_spice_interface
                   ${SPICE_LIBRARIES}
                   ${Boost_LIBRARIES})

  tudat_add_test(
    TARGET test_AngularPositionModel
    SOURCES "${OBSMDL_TEST_DIR}/unitTestAngularPositionModel.cpp"
    LINK_LIBRARIES ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

  tudat_add_test(
    TARGET test_OneWayRangeModel
    SOURCES "${OBSMDL_TEST_DIR}/unitTestOneWayRangeObservationModel.cpp"
    LINK_LIBRARIES ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

  tudat_add_test(
    TARGET test_OneWayDifferencedRangeModel
    SOURCES "${OBSMDL_TEST_DIR}/unitTestOneWayDifferencedRangeModel.cpp"
    LINK_LIBRARIES ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

  tudat_add_test(
    TARGET test_DopplerModels
    SOURCES "${OBSMDL_TEST_DIR}/unitTestDopplerModels.cpp"
    LINK_LIBRARIES ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

  tudat_add_test(
    TARGET test_PositionObservationModel
    SOURCES "${OBSMDL_TEST_DIR}/unitTestPositionObservationModel.cpp"
    LINK_LIBRARIES ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

  tudat_add_test(
    TARGET test_EulerAngleObservationModel
    SOURCES "${OBSMDL_TEST_DIR}/unitTestEulerAngleObservationModel.cpp"
    LINK_LIBRARIES ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

  tudat_add_test(
    TARGET test_VelocityObservationModel
    SOURCES "${OBSMDL_TEST_DIR}/unitTestVelocityObservationModel.cpp"
    LINK_LIBRARIES ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

  tudat_add_test(
    TARGET test_ObservationNoiseSimulation
    SOURCES "${OBSMDL_TEST_DIR}/unitTestSimulatedObservationNoise.cpp"
    LINK_LIBRARIES ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

  tudat_add_test(
    TARGET test_NWayRangeObservationModel
    SOURCES "${OBSMDL_TEST_DIR}/unitTestNWayRangeObservationModel.cpp"
    LINK_LIBRARIES ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

  tudat_add_test(
    TARGET test_ObservationViabilityCalculators
    SOURCES "${OBSMDL_TEST_DIR}/unitTestObservationViabilityCalculators.cpp"
    LINK_LIBRARIES ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

endif( )
