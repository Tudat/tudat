#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.

# Set the source files.
set(EARTH_ORIENTATION_SOURCES
  "${SRCROOT}${EARTHORIENTATIONDIR}/earthOrientationCalculator.cpp"
  "${SRCROOT}${EARTHORIENTATIONDIR}/terrestrialTimeScaleConverter.cpp"
  "${SRCROOT}${EARTHORIENTATIONDIR}/eopReader.cpp"
  "${SRCROOT}${EARTHORIENTATIONDIR}/polarMotionCalculator.cpp"
  "${SRCROOT}${EARTHORIENTATIONDIR}/precessionNutationCalculator.cpp"
  "${SRCROOT}${EARTHORIENTATIONDIR}/readAmplitudeAndArgumentMultipliers.cpp"
  "${SRCROOT}${EARTHORIENTATIONDIR}/UnitTests/sofaEarthOrientationCookbookExamples.cpp"
  "${SRCROOT}${EARTHORIENTATIONDIR}/shortPeriodEarthOrientationCorrectionCalculator.cpp"
)

# Set the header files.
set(EARTH_ORIENTATION_HEADERS
  "${SRCROOT}${EARTHORIENTATIONDIR}/earthOrientationCalculator.h"
  "${SRCROOT}${EARTHORIENTATIONDIR}/terrestrialTimeScaleConverter.h"
  "${SRCROOT}${EARTHORIENTATIONDIR}/eopReader.h"
  "${SRCROOT}${EARTHORIENTATIONDIR}/polarMotionCalculator.h"
  "${SRCROOT}${EARTHORIENTATIONDIR}/precessionNutationCalculator.h"
  "${SRCROOT}${EARTHORIENTATIONDIR}/readAmplitudeAndArgumentMultipliers.h"
  "${SRCROOT}${EARTHORIENTATIONDIR}/shortPeriodEarthOrientationCorrectionCalculator.h"
  "${SRCROOT}${EARTHORIENTATIONDIR}/UnitTests/sofaEarthOrientationCookbookExamples.h"
)

# Add static libraries.
add_library(tudat_earth_orientation STATIC ${EARTH_ORIENTATION_SOURCES} ${EARTH_ORIENTATION_HEADERS})
setup_tudat_library_target(tudat_earth_orientation "${SRCROOT}${EARTHORIENTATIONDIR}")

# Add unit tests.
set(EO_TEST_DIR "${SRCROOT}${EARTHORIENTATIONDIR}/UnitTests")

tudat_add_test(
  TARGET test_EarthOrientationCalculator
  SOURCES "${EO_TEST_DIR}/unitTestEarthOrientationCalculator.cpp"
  LINK_LIBRARIES tudat_earth_orientation
                 tudat_spice_interface
                 tudat_sofa_interface
                 tudat_interpolators
                 tudat_basic_astrodynamics
                 tudat_basic_mathematics
                 tudat_input_output
                 tudat_input_output
                 sofa cspice
                 ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_EopReader
  SOURCES "${EO_TEST_DIR}/unitTestEopReader.cpp"
  LINK_LIBRARIES tudat_earth_orientation
                 tudat_sofa_interface
                 tudat_interpolators
                 tudat_basic_astrodynamics
                 tudat_basic_mathematics
                 tudat_input_output sofa
                 ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_PolarMotion
  SOURCES "${EO_TEST_DIR}/unitTestPolarMotionCalculator.cpp"
  LINK_LIBRARIES tudat_earth_orientation
                 tudat_sofa_interface
                 tudat_interpolators
                 tudat_basic_astrodynamics
                 tudat_basic_mathematics
                 tudat_input_output sofa
                 ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_TimeScaleConverter
  SOURCES "${EO_TEST_DIR}/unitTestTimeScaleConverter.cpp"
  LINK_LIBRARIES tudat_earth_orientation
                 tudat_sofa_interface
                 tudat_interpolators
                 tudat_basic_astrodynamics
                 tudat_basic_mathematics
                 tudat_input_output
                 sofa
                 ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_ShortPeriodEopCorrections
  SOURCES "${EO_TEST_DIR}/unitTestShortPeriodEopCorrections.cpp"
  LINK_LIBRARIES tudat_earth_orientation
                 tudat_sofa_interface
                 tudat_basic_astrodynamics
                 tudat_basic_astrodynamics
                 tudat_basic_mathematics
                 tudat_input_output
                 sofa
                 ${Boost_LIBRARIES})


