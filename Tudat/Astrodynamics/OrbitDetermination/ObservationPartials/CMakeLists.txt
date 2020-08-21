 #    Copyright (c) 2010-2019, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.

# Set the source files.
set(OBSERVATION_PARTIALS_SOURCES
  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/angularPositionPartial.cpp"
  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/rotationMatrixPartial.cpp"
  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/positionPartials.cpp"
  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/oneWayRangePartial.cpp"
  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/twoWayDopplerPartial.cpp"
  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/nWayRangePartial.cpp"
  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/differencedOneWayRangeRatePartial.cpp"
  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/oneWayDopplerPartial.cpp"
  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/observationPartial.cpp"
  "${SRCROOT}${LIGHTTIMECORRECTIONPARTIALSDIR}/lightTimeCorrectionPartial.cpp"
  "${SRCROOT}${LIGHTTIMECORRECTIONPARTIALSDIR}/firstOrderRelativisticPartial.cpp"
)

# Set the header files.
set(OBSERVATION_PARTIALS_HEADERS
  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/angularPositionPartial.h"
  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/rotationMatrixPartial.h"
  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/positionPartials.h"
  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/oneWayRangePartial.h"
  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/nWayRangePartial.h"
  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/differencedOneWayRangeRatePartial.h"
  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/oneWayDopplerPartial.h"
  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/twoWayDopplerPartial.h"
  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/observationPartial.h"
  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/eulerAngleObservablePartials.h"
  "${SRCROOT}${LIGHTTIMECORRECTIONPARTIALSDIR}/lightTimeCorrectionPartial.h"
  "${SRCROOT}${LIGHTTIMECORRECTIONPARTIALSDIR}/firstOrderRelativisticPartial.h"
)

# Set the header files.
set(OBSERVATION_PARTIALS_TEST_HEADERS
  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/UnitTests/observationPartialTestFunctions.h"
  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/UnitTests/numericalObservationPartial.h"
)

set(OBSERVATION_PARTIALS_TEST_SOURCES
  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/UnitTests/observationPartialTestFunctions.cpp"
  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/UnitTests/numericalObservationPartial.cpp"
)
# Add static libraries.
add_library(tudat_observation_partials STATIC ${OBSERVATION_PARTIALS_SOURCES} ${OBSERVATION_PARTIALS_HEADERS})
setup_tudat_library_target(tudat_observation_partials "${SRCROOT}{OBSERVATIONPARTIALSDIR}")

add_library(tudat_observation_partials_test_functions STATIC ${OBSERVATION_PARTIALS_TEST_SOURCES} ${OBSERVATION_PARTIALS_TEST_HEADERS})
setup_tudat_library_target(tudat_observation_partials "${SRCROOT}{OBSERVATIONPARTIALSDIR}/UnitTests")

# Add unit tests.
set(OBSPARTIALS_TEST_DIR "${SRCROOT}${OBSERVATIONPARTIALSDIR}/UnitTests")

tudat_add_test(
  TARGET test_OneWayRangePartials
  SOURCES "${OBSPARTIALS_TEST_DIR}/unitTestOneWayRangePartials.cpp"
  LINK_LIBRARIES tudat_observation_partials_test_functions
                 ${TUDAT_ESTIMATION_LIBRARIES}
                 ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_NWayRangePartials
  SOURCES "${OBSPARTIALS_TEST_DIR}/unitTestNWayRangePartials.cpp"
  LINK_LIBRARIES tudat_observation_partials_test_functions
                 ${TUDAT_ESTIMATION_LIBRARIES}
                 ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_OneWayDifferencedRangePartials
  SOURCES "${OBSPARTIALS_TEST_DIR}/unitTestDifferencedOneWayRangePartials.cpp"
  LINK_LIBRARIES tudat_observation_partials_test_functions
                 ${TUDAT_ESTIMATION_LIBRARIES}
                 ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_OneWayDopplerPartials
  SOURCES "${OBSPARTIALS_TEST_DIR}/unitTestOneWayDopplerPartials.cpp"
  LINK_LIBRARIES tudat_observation_partials_test_functions
                 ${TUDAT_ESTIMATION_LIBRARIES}
                 ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_TwoWayDopplerPartials
  SOURCES "${OBSPARTIALS_TEST_DIR}/unitTestTwoWayDopplerPartials.cpp"
  LINK_LIBRARIES tudat_observation_partials_test_functions
                 ${TUDAT_ESTIMATION_LIBRARIES}
                 ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_AngularPositionPartials
  SOURCES "${OBSPARTIALS_TEST_DIR}/unitTestAngularPositionPartials.cpp"
  LINK_LIBRARIES tudat_observation_partials_test_functions
                 ${TUDAT_ESTIMATION_LIBRARIES}
                 ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_PositionPartials
  SOURCES "${OBSPARTIALS_TEST_DIR}/unitTestPositionPartials.cpp"
  LINK_LIBRARIES tudat_observation_partials_test_functions
                 ${TUDAT_ESTIMATION_LIBRARIES}
                 ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_EulerAngleObservationPartials
  SOURCES "${OBSPARTIALS_TEST_DIR}/unitTestEulerAngleObservationPartials.cpp"
  LINK_LIBRARIES tudat_observation_partials_test_functions
                 ${TUDAT_ESTIMATION_LIBRARIES}
                 ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_RotationMatrixPartials
  SOURCES "${OBSPARTIALS_TEST_DIR}/unitTestRotationMatrixPartials.cpp"
  LINK_LIBRARIES ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_FullRotationPartials
  SOURCES "${OBSPARTIALS_TEST_DIR}/unitTestFullRotationPartials.cpp"
  LINK_LIBRARIES ${TUDAT_PROPAGATION_LIBRARIES}
                 ${TUDAT_ESTIMATION_LIBRARIES}
                 ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_LightTimePartials
  SOURCES "${SRCROOT}${LIGHTTIMECORRECTIONPARTIALSDIR}/UnitTests/unitTestLightTimePartials.cpp"
  LINK_LIBRARIES tudat_observation_partials_test_functions 
                 ${TUDAT_ESTIMATION_LIBRARIES}
                 ${Boost_LIBRARIES})
