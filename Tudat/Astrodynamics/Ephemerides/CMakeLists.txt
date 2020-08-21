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
set(EPHEMERIDES_SOURCES
  "${SRCROOT}${EPHEMERIDESDIR}/approximatePlanetPositionsBase.cpp"
  "${SRCROOT}${EPHEMERIDESDIR}/approximatePlanetPositions.cpp"
  "${SRCROOT}${EPHEMERIDESDIR}/approximatePlanetPositionsCircularCoplanar.cpp"
  "${SRCROOT}${EPHEMERIDESDIR}/ephemeris.cpp"
  "${SRCROOT}${EPHEMERIDESDIR}/rotationalEphemeris.cpp"
  "${SRCROOT}${EPHEMERIDESDIR}/cartesianStateExtractor.cpp"
  "${SRCROOT}${EPHEMERIDESDIR}/keplerStateExtractor.cpp"
  "${SRCROOT}${EPHEMERIDESDIR}/keplerEphemeris.cpp"
  "${SRCROOT}${EPHEMERIDESDIR}/rotationalEphemeris.cpp"
  "${SRCROOT}${EPHEMERIDESDIR}/simpleRotationalEphemeris.cpp"
  "${SRCROOT}${EPHEMERIDESDIR}/tabulatedEphemeris.cpp"
  "${SRCROOT}${EPHEMERIDESDIR}/frameManager.cpp"
  "${SRCROOT}${EPHEMERIDESDIR}/compositeEphemeris.cpp"
  "${SRCROOT}${EPHEMERIDESDIR}/tabulatedRotationalEphemeris.cpp"
  "${SRCROOT}${EPHEMERIDESDIR}/synchronousRotationalEphemeris.cpp"
  "${SRCROOT}${EPHEMERIDESDIR}/fullPlanetaryRotationModel.cpp"
    "${SRCROOT}${EPHEMERIDESDIR}/tleEphemeris.cpp"
)

# Set the header files.
set(EPHEMERIDES_HEADERS
  "${SRCROOT}${EPHEMERIDESDIR}/approximatePlanetPositionsBase.h"
  "${SRCROOT}${EPHEMERIDESDIR}/approximatePlanetPositions.h"
  "${SRCROOT}${EPHEMERIDESDIR}/approximatePlanetPositionsCircularCoplanar.h"
  "${SRCROOT}${EPHEMERIDESDIR}/approximatePlanetPositionsDataContainer.h"
  "${SRCROOT}${EPHEMERIDESDIR}/ephemeris.h"
  "${SRCROOT}${EPHEMERIDESDIR}/constantEphemeris.h"
  "${SRCROOT}${EPHEMERIDESDIR}/cartesianStateExtractor.h"
  "${SRCROOT}${EPHEMERIDESDIR}/keplerStateExtractor.h"
  "${SRCROOT}${EPHEMERIDESDIR}/keplerEphemeris.h"
  "${SRCROOT}${EPHEMERIDESDIR}/rotationalEphemeris.h"
  "${SRCROOT}${EPHEMERIDESDIR}/constantRotationalEphemeris.h"
  "${SRCROOT}${EPHEMERIDESDIR}/simpleRotationalEphemeris.h"
  "${SRCROOT}${EPHEMERIDESDIR}/tabulatedEphemeris.h"
  "${SRCROOT}${EPHEMERIDESDIR}/frameManager.h"
  "${SRCROOT}${EPHEMERIDESDIR}/itrsToGcrsRotationModel.h"
  "${SRCROOT}${EPHEMERIDESDIR}/compositeEphemeris.h"
  "${SRCROOT}${EPHEMERIDESDIR}/constantEphemeris.h"
  "${SRCROOT}${EPHEMERIDESDIR}/constantRotationalEphemeris.h"
  "${SRCROOT}${EPHEMERIDESDIR}/multiArcEphemeris.h"
  "${SRCROOT}${EPHEMERIDESDIR}/tabulatedRotationalEphemeris.h"
  "${SRCROOT}${EPHEMERIDESDIR}/fullPlanetaryRotationModel.h"
  "${SRCROOT}${EPHEMERIDESDIR}/synchronousRotationalEphemeris.h"
    "${SRCROOT}${EPHEMERIDESDIR}/tleEphemeris.h")

# Add static libraries.
add_library(tudat_ephemerides STATIC ${EPHEMERIDES_SOURCES} ${EPHEMERIDES_HEADERS})
setup_tudat_library_target(tudat_ephemerides "${SRCROOT}${EPHEMERIDESDIR}")

# Add unit tests.
set(EPHEM_TEST_DIR "${SRCROOT}${EPHEMERIDESDIR}/UnitTests")

tudat_add_test(
  TARGET test_ApproximatePlanetPositions
  SOURCES "${EPHEM_TEST_DIR}/unitTestApproximatePlanetPositions.cpp"
  LINK_LIBRARIES tudat_ephemerides
                 tudat_gravitation
                 tudat_basic_astrodynamics
                 tudat_basic_mathematics
                 tudat_input_output
                 tudat_root_finders
                 ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_TabulatedEphemeris
  SOURCES "${EPHEM_TEST_DIR}/unitTestTabulatedEphemeris.cpp"
  LINK_LIBRARIES tudat_ephemerides
                 tudat_interpolators
                 tudat_basic_astrodynamics
                 tudat_basic_mathematics
                 ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_CartesianStateExtractor
  SOURCES "${EPHEM_TEST_DIR}/unitTestCartesianStateExtractor.cpp"
  LINK_LIBRARIES tudat_input_output tudat_ephemerides ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_KeplerStateExtractor
  SOURCES "${EPHEM_TEST_DIR}/unitTestKeplerStateExtractor.cpp"
  LINK_LIBRARIES tudat_ephemerides
                 tudat_input_output
                 tudat_basic_astrodynamics
                 tudat_basic_mathematics
                 tudat_root_finders
                 ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_RotationalEphemeris
  SOURCES "${EPHEM_TEST_DIR}/unitTestRotationalEphemeris.cpp"
  LINK_LIBRARIES tudat_ephemerides
                 tudat_reference_frames
                 tudat_input_output
                 tudat_basic_astrodynamics
                 tudat_basic_mathematics
                 tudat_spice_interface
                 ${TUDAT_EXTERNAL_LIBRARIES}
                 ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_SimpleRotationalEphemeris
  SOURCES "${EPHEM_TEST_DIR}/unitTestSimpleRotationalEphemeris.cpp"
  LINK_LIBRARIES tudat_ephemerides
                 tudat_reference_frames
                 tudat_input_output
                 tudat_basic_astrodynamics
                 tudat_basic_mathematics
                 ${Boost_LIBRARIES})

if(USE_CSPICE)

  tudat_add_test(
    TARGET test_FrameManager
    SOURCES "${EPHEM_TEST_DIR}/unitTestFrameManager.cpp"
    LINK_LIBRARIES tudat_ephemerides
                   tudat_reference_frames
                   ${TUDAT_EXTERNAL_INTERFACE_LIBRARIES}
                   tudat_input_output
                   tudat_basic_astrodynamics
                   tudat_basic_mathematics
                   ${TUDAT_EXTERNAL_LIBRARIES}
                   ${Boost_LIBRARIES})

  tudat_add_test(
    TARGET test_CompositeEphemeris
    SOURCES "${EPHEM_TEST_DIR}/unitTestCompositeEphemeris.cpp"
    LINK_LIBRARIES ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

  tudat_add_test(
    TARGET test_TabulatedRotationalEphemeris
    SOURCES "${EPHEM_TEST_DIR}/unitTestTabulatedRotationalEphemeris.cpp"
    LINK_LIBRARIES tudat_ephemerides
                   tudat_reference_frames
                   tudat_input_output
                   tudat_basic_astrodynamics
                   tudat_basic_mathematics
                   tudat_spice_interface
                   ${TUDAT_EXTERNAL_LIBRARIES}
                   ${Boost_LIBRARIES})

  tudat_add_test(
    TARGET test_SynchronousRotationalEphemeris
    SOURCES "${EPHEM_TEST_DIR}/unitTestSynchronousRotationalEphemeris.cpp"
    LINK_LIBRARIES ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

  tudat_add_test(
    TARGET test_TwoLineElementsEphemeris
    SOURCES "${EPHEM_TEST_DIR}/unitTestTwoLineElementsEphemeris.cpp"
    LINK_LIBRARIES ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

endif()

tudat_add_test(
  TARGET test_KeplerEphemeris
  SOURCES "${EPHEM_TEST_DIR}/unitTestKeplerEphemeris.cpp"
  LINK_LIBRARIES tudat_ephemerides
                 tudat_reference_frames
                 tudat_input_output
                 tudat_basic_astrodynamics
                 tudat_basic_mathematics
                 ${Boost_LIBRARIES})

if(USE_SOFA)
  tudat_add_test(
    TARGET test_GcrsToItrsRotation
    SOURCES "${EPHEM_TEST_DIR}/unitTestItrsToGcrsRotationModel.cpp"
    LINK_LIBRARIES tudat_ephemerides
                   tudat_spice_interface
                   tudat_reference_frames
                   tudat_earth_orientation
                   tudat_sofa_interface
                   tudat_basic_astrodynamics
                   tudat_basic_mathematics
                   tudat_input_output
                   cspice sofa
                   ${Boost_LIBRARIES})
endif( )

tudat_add_test(
  TARGET test_PlanetaryRotationModel
  SOURCES "${EPHEM_TEST_DIR}/unitTestPlanetaryRotationModel.cpp"
  LINK_LIBRARIES ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})
