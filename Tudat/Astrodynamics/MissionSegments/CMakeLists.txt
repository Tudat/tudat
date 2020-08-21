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
set(MISSIONSEGMENTS_SOURCES
  "${SRCROOT}${MISSIONSEGMENTSDIR}/escapeAndCapture.cpp"
  "${SRCROOT}${MISSIONSEGMENTSDIR}/gravityAssist.cpp"
  "${SRCROOT}${MISSIONSEGMENTSDIR}/improvedInversePolynomialWall.cpp"
  "${SRCROOT}${MISSIONSEGMENTSDIR}/lambertTargeterIzzo.cpp"
  "${SRCROOT}${MISSIONSEGMENTSDIR}/lambertTargeterGooding.cpp"
  "${SRCROOT}${MISSIONSEGMENTSDIR}/lambertRoutines.cpp"
  "${SRCROOT}${MISSIONSEGMENTSDIR}/multiRevolutionLambertTargeterIzzo.cpp"
  "${SRCROOT}${MISSIONSEGMENTSDIR}/oscillatingFunctionNovak.cpp"
  "${SRCROOT}${MISSIONSEGMENTSDIR}/zeroRevolutionLambertTargeterIzzo.cpp"
)

# Set the header files.
set(MISSIONSEGMENTS_HEADERS
  "${SRCROOT}${MISSIONSEGMENTSDIR}/escapeAndCapture.h"
  "${SRCROOT}${MISSIONSEGMENTSDIR}/gravityAssist.h"
  "${SRCROOT}${MISSIONSEGMENTSDIR}/improvedInversePolynomialWall.h"
  "${SRCROOT}${MISSIONSEGMENTSDIR}/lambertTargeter.h"
  "${SRCROOT}${MISSIONSEGMENTSDIR}/lambertTargeterIzzo.h"
  "${SRCROOT}${MISSIONSEGMENTSDIR}/lambertTargeterGooding.h"
  "${SRCROOT}${MISSIONSEGMENTSDIR}/lambertRoutines.h"
  "${SRCROOT}${MISSIONSEGMENTSDIR}/multiRevolutionLambertTargeterIzzo.h"
  "${SRCROOT}${MISSIONSEGMENTSDIR}/oscillatingFunctionNovak.h"
  "${SRCROOT}${MISSIONSEGMENTSDIR}/zeroRevolutionLambertTargeterIzzo.h"
)

# Add static libraries.
add_library(tudat_mission_segments STATIC ${MISSIONSEGMENTS_SOURCES} ${MISSIONSEGMENTS_HEADERS})
setup_tudat_library_target(tudat_mission_segments "${SRCROOT}${MISSIONSEGMENTSDIR}")

# Add unit tests.
set(MISS_SEG_TEST_DIR "${SRCROOT}${MISSIONSEGMENTSDIR}/UnitTests")

tudat_add_test(
  TARGET test_EscapeAndCapture
  SOURCES "${MISS_SEG_TEST_DIR}/unitTestEscapeAndCapture.cpp"
  LINK_LIBRARIES tudat_mission_segments
                 tudat_gravitation
                 tudat_basic_astrodynamics
                 tudat_basic_mathematics
                 tudat_root_finders
                 tudat_input_output
                 ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_GravityAssist
  SOURCES "${MISS_SEG_TEST_DIR}/unitTestGravityAssist.cpp"
  LINK_LIBRARIES tudat_gravitation
                 tudat_mission_segments
                 tudat_basic_astrodynamics
                 tudat_basic_mathematics
                 tudat_root_finders
                 tudat_input_output
                 ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_LambertTargeterIzzo
  SOURCES "${MISS_SEG_TEST_DIR}/unitTestLambertTargeterIzzo.cpp"
  LINK_LIBRARIES tudat_mission_segments
                 tudat_root_finders
                 tudat_basic_astrodynamics
                 tudat_basic_mathematics
                 ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_LambertTargeterGooding
  SOURCES "${MISS_SEG_TEST_DIR}/unitTestLambertTargeterGooding.cpp"
  LINK_LIBRARIES tudat_mission_segments
                 tudat_root_finders
                 tudat_basic_astrodynamics
                 tudat_basic_mathematics
                 ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_LambertTargeter
  SOURCES "${MISS_SEG_TEST_DIR}/unitTestLambertTargeter.cpp"
  LINK_LIBRARIES tudat_mission_segments
                 tudat_root_finders
                 tudat_basic_astrodynamics
                 tudat_basic_mathematics
                 ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_LambertRoutines
  SOURCES "${MISS_SEG_TEST_DIR}/unitTestLambertRoutines.cpp"
  LINK_LIBRARIES tudat_mission_segments
                 tudat_root_finders
                 tudat_basic_astrodynamics
                 tudat_basic_mathematics
                 ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_ZeroRevolutionLambertTargeterIzzo
  SOURCES "${MISS_SEG_TEST_DIR}/unitTestZeroRevolutionLambertTargeterIzzo.cpp"
  LINK_LIBRARIES tudat_mission_segments
                 tudat_root_finders
                 tudat_basic_astrodynamics
                 tudat_basic_mathematics
                 ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_MultiRevolutionLambertTargeterIzzo
  SOURCES "${MISS_SEG_TEST_DIR}/unitTestMultiRevolutionLambertTargeterIzzo.cpp"
  LINK_LIBRARIES tudat_mission_segments
                 tudat_root_finders
                 tudat_basic_astrodynamics
                 tudat_basic_mathematics
                 ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_MathematicalShapeFunctions
  SOURCES "${MISS_SEG_TEST_DIR}/unitTestMathematicalShapeFunctions.cpp"
  LINK_LIBRARIES tudat_mission_segments tudat_basic_mathematics ${Boost_LIBRARIES})
