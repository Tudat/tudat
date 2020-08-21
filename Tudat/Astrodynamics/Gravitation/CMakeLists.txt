 #    Copyright (c) 2010-2019, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.
 #
 #    Notes
 #
 #

# Set the source files.
set(GRAVITATION_SOURCES
  "${SRCROOT}${GRAVITATIONDIR}/basicSolidBodyTideGravityFieldVariations.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/gravityFieldVariations.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/centralGravityModel.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/centralJ2GravityModel.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/centralJ2J3GravityModel.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/centralJ2J3J4GravityModel.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/gravityFieldModel.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/jacobiEnergy.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/librationPoint.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/sphericalHarmonicsGravityModel.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/sphericalHarmonicsGravityField.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/thirdBodyPerturbation.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/timeDependentSphericalHarmonicsGravityField.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/unitConversionsCircularRestrictedThreeBodyProblem.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/UnitTests/planetTestData.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/unitConversionsCircularRestrictedThreeBodyProblem.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/triAxialEllipsoidGravity.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/tabulatedGravityFieldVariations.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/mutualSphericalHarmonicGravityModel.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/secondDegreeGravitationalTorque.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/directTidalDissipationAcceleration.cpp"
)

# Set the header files.
set(GRAVITATION_HEADERS
  "${SRCROOT}${GRAVITATIONDIR}/basicSolidBodyTideGravityFieldVariations.h"
  "${SRCROOT}${GRAVITATIONDIR}/gravityFieldVariations.h"
  "${SRCROOT}${GRAVITATIONDIR}/centralGravityModel.h"
  "${SRCROOT}${GRAVITATIONDIR}/centralJ2GravityModel.h"
  "${SRCROOT}${GRAVITATIONDIR}/centralJ2J3GravityModel.h"
  "${SRCROOT}${GRAVITATIONDIR}/centralJ2J3J4GravityModel.h"
  "${SRCROOT}${GRAVITATIONDIR}/gravityFieldModel.h"
  "${SRCROOT}${GRAVITATIONDIR}/jacobiEnergy.h"
  "${SRCROOT}${GRAVITATIONDIR}/librationPoint.h"
  "${SRCROOT}${GRAVITATIONDIR}/sphericalHarmonicsGravityModel.h"
  "${SRCROOT}${GRAVITATIONDIR}/sphericalHarmonicsGravityModelBase.h"
  "${SRCROOT}${GRAVITATIONDIR}/sphericalHarmonicsGravityField.h"
  "${SRCROOT}${GRAVITATIONDIR}/thirdBodyPerturbation.h"
  "${SRCROOT}${GRAVITATIONDIR}/timeDependentSphericalHarmonicsGravityField.h"
  "${SRCROOT}${GRAVITATIONDIR}/unitConversionsCircularRestrictedThreeBodyProblem.h"
  "${SRCROOT}${GRAVITATIONDIR}/UnitTests/planetTestData.h"
  "${SRCROOT}${GRAVITATIONDIR}/triAxialEllipsoidGravity.h"
  "${SRCROOT}${GRAVITATIONDIR}/tabulatedGravityFieldVariations.h"
  "${SRCROOT}${GRAVITATIONDIR}/mutualSphericalHarmonicGravityModel.h"
  "${SRCROOT}${GRAVITATIONDIR}/secondDegreeGravitationalTorque.h"
  "${SRCROOT}${GRAVITATIONDIR}/directTidalDissipationAcceleration.h"
  "${SRCROOT}${GRAVITATIONDIR}/sphericalHarmonicGravitationalTorque.h"
)

# Add static libraries.
add_library(tudat_gravitation STATIC ${GRAVITATION_SOURCES} ${GRAVITATION_HEADERS})
setup_tudat_library_target(tudat_gravitation "${SRCROOT}${GRAVITATIONDIR}")

# Add unit tests.
set(GRAV_TEST_DIR "${SRCROOT}${GRAVITATIONDIR}/UnitTests")

tudat_add_test(
  TARGET test_SphericalHarmonicsGravityField
  SOURCES "${GRAV_TEST_DIR}/unitTestSphericalHarmonicsGravityField.cpp"
  LINK_LIBRARIES tudat_gravitation
                 tudat_basic_astrodynamics
                 tudat_basic_mathematics
                 ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_GravitationalForce
  SOURCES "${GRAV_TEST_DIR}/unitTestGravitationalForce.cpp"
  LINK_LIBRARIES tudat_gravitation tudat_basic_astrodynamics ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_CentralAndZonalGravityModel
  SOURCES "${GRAV_TEST_DIR}/unitTestCentralAndZonalGravityModel.cpp"
  LINK_LIBRARIES tudat_gravitation
                   tudat_basic_astrodynamics
                   tudat_basic_astrodynamics
                   tudat_numerical_integrators
                   ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_UnitConversionsCircularRestrictedThreeBodyProblem
  SOURCES "${GRAV_TEST_DIR}/unitTestUnitConversionsCircularRestrictedThreeBodyProblem.cpp"
  LINK_LIBRARIES tudat_gravitation ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_JacobiEnergy
  SOURCES "${GRAV_TEST_DIR}/unitTestJacobiEnergy.cpp"
  LINK_LIBRARIES tudat_gravitation
                   tudat_mission_segments
                   tudat_root_finders
                   tudat_basic_astrodynamics
                   tudat_input_output
                   ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_LibrationPoints
  SOURCES "${GRAV_TEST_DIR}/unitTestLibrationPoints.cpp"
  LINK_LIBRARIES tudat_ephemerides
                   tudat_root_finders
                   tudat_gravitation
                   tudat_input_output
                   ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_SphericalHarmonicsGravityModel
  SOURCES "${GRAV_TEST_DIR}/unitTestSphericalHarmonicsGravityModel.cpp"
  LINK_LIBRARIES tudat_gravitation
                   tudat_basic_astrodynamics
                   tudat_basic_mathematics
                   ${Boost_LIBRARIES} )

tudat_add_test(
  TARGET test_ThirdBodyPerturbation
  SOURCES "${GRAV_TEST_DIR}/unitTestThirdBodyPerturbation.cpp"
  LINK_LIBRARIES tudat_gravitation
                 tudat_basic_astrodynamics
                 tudat_basic_mathematics
                 ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_TriAxialEllipsoidGravity
  SOURCES "${GRAV_TEST_DIR}/unitTestTriAxialEllipsoidGravity.cpp"
  LINK_LIBRARIES tudat_gravitation tudat_basic_mathematics ${Boost_LIBRARIES} )

tudat_add_test(
  TARGET test_DirectTidalDissipationAcceleration
  SOURCES "${GRAV_TEST_DIR}/unitTestDirectTidalDissipation.cpp"
  LINK_LIBRARIES ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

if(USE_CSPICE)

  tudat_add_test(
    TARGET test_GravityFieldVariations
    SOURCES "${GRAV_TEST_DIR}/unitTestGravityFieldVariations.cpp"
    DEPENDS ${SPICE_LIBRARIES}
    LINK_LIBRARIES tudat_gravitation
                   tudat_interpolators
                   tudat_basic_mathematics
                   tudat_spice_interface
                   ${SPICE_LIBRARIES}
                   ${Boost_LIBRARIES} )

  tudat_add_test(
    TARGET test_MutualSphericalHarmonicsGravityModel
    SOURCES "${GRAV_TEST_DIR}/unitTestMutualSphericalHarmonicAcceleration.cpp"
    LINK_LIBRARIES ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

  tudat_add_test(
    TARGET test_GravitationalTorques
    SOURCES "${GRAV_TEST_DIR}/unitTestGravitationalTorques.cpp"
    LINK_LIBRARIES ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

endif()
