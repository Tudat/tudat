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
set(BASICASTRODYNAMICS_SOURCES
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/accelerationModelTypes.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/accelerationModel.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/attitudeElementConversions.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/clohessyWiltshirePropagator.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/geodeticCoordinateConversions.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/missionGeometry.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/modifiedEquinoctialElementConversions.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/timeConversions.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/astrodynamicsFunctions.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/physicalConstants.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/bodyShapeModel.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/sphericalStateConversions.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/unifiedStateModelQuaternionElementConversions.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/unifiedStateModelModifiedRodriguesParameterElementConversions.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/unifiedStateModelExponentialMapElementConversions.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/stateRepresentationConversions.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/torqueModelTypes.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/torqueModel.cpp"
)

# Set the header files.
set(BASICASTRODYNAMICS_HEADERS
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/accelerationModelTypes.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/accelerationModel.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/attitudeElementConversions.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/celestialBodyConstants.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/convertMeanToEccentricAnomalies.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/clohessyWiltshirePropagator.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/customTorque.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/geodeticCoordinateConversions.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/keplerPropagator.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/missionGeometry.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/modifiedEquinoctialElementConversions.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/stateVectorIndices.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/timeConversions.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/testAccelerationModels.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/testBody.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/keplerPropagatorTestData.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/astrodynamicsFunctions.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/orbitalElementConversions.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/physicalConstants.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/unitConversions.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/bodyShapeModel.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/oblateSpheroidBodyShapeModel.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/sphericalBodyShapeModel.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/massRateModel.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/sphericalStateConversions.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/unifiedStateModelQuaternionElementConversions.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/unifiedStateModelModifiedRodriguesParameterElementConversions.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/unifiedStateModelExponentialMapElementConversions.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/stateRepresentationConversions.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/empiricalAcceleration.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/torqueModel.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/torqueModelTypes.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/dissipativeTorqueModel.h"
)

# Add static libraries.
add_library(tudat_basic_astrodynamics STATIC ${BASICASTRODYNAMICS_SOURCES} ${BASICASTRODYNAMICS_HEADERS})
setup_tudat_library_target(tudat_basic_astrodynamics "${SRCROOT}${BASICASTRODYNAMICSDIR}")

# Add unit tests.
set(BASICASTRODYNAMICS_UNITTESTS
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/unitTestBasicAstrodynamics.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/unitTestAstrodynamicsFunctions.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/unitTestOrbitalElementConversions.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/unitTestPhysicalConstants.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/unitTestUnitConversions.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/unitTestUnifiedStateModelQuaternionsElementConversions.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/unitTestUnifiedStateModelModifiedRodriguesParametersElementConversions.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/unitTestUnifiedStateModelExponentialMapElementConversions.cpp"
)

set(BASIC_ASTRO_TEST_DIR "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests")

tudat_add_test(
  TARGET test_AstrodynamicsFunctions
  SOURCES "${BASIC_ASTRO_TEST_DIR}/unitTestAstrodynamicsFunctions.cpp"
  LINK_LIBRARIES tudat_basic_astrodynamics tudat_basic_mathematics ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_OrbitalElementConversions
  SOURCES "${BASIC_ASTRO_TEST_DIR}/unitTestOrbitalElementConversions.cpp"
  LINK_LIBRARIES tudat_basic_astrodynamics tudat_basic_mathematics ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_PhysicalConstants
  SOURCES "${BASIC_ASTRO_TEST_DIR}/unitTestPhysicalConstants.cpp"
  LINK_LIBRARIES tudat_basic_astrodynamics tudat_basic_mathematics ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_UnitConversions
  SOURCES "${BASIC_ASTRO_TEST_DIR}/unitTestUnitConversions.cpp"
  LINK_LIBRARIES tudat_basic_astrodynamics tudat_basic_mathematics ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_MeanToEccentricAnomalyConversion
  SOURCES "${BASIC_ASTRO_TEST_DIR}/unitTestConvertMeanAnomalyToEccentricAnomaly.cpp"
  LINK_LIBRARIES tudat_basic_astrodynamics
                 tudat_basic_mathematics
                 tudat_root_finders
                 ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_MeanToHyperbolicEccentricAnomalyConversion
  SOURCES "${BASIC_ASTRO_TEST_DIR}/unitTestConvertMeanAnomalyToHyperbolicEccentricAnomaly.cpp"
  LINK_LIBRARIES tudat_basic_astrodynamics
                 tudat_basic_mathematics
                 tudat_root_finders
                 ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_KeplerPropagator
  SOURCES "${BASIC_ASTRO_TEST_DIR}/unitTestKeplerPropagator.cpp"
  LINK_LIBRARIES tudat_input_output
                 tudat_gravitation
                 tudat_basic_astrodynamics
                 tudat_basic_mathematics
                 tudat_root_finders
                 ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_AccelerationModel
  SOURCES "${BASIC_ASTRO_TEST_DIR}/unitTestAccelerationModel.cpp"
  LINK_LIBRARIES tudat_basic_astrodynamics ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_ClohessyWiltshirePropagator
  SOURCES "${BASIC_ASTRO_TEST_DIR}/unitTestClohessyWiltshirePropagator.cpp"
  LINK_LIBRARIES tudat_basic_astrodynamics ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_MissionGeometry
  SOURCES "${BASIC_ASTRO_TEST_DIR}/unitTestMissionGeometry.cpp"
  LINK_LIBRARIES tudat_basic_astrodynamics ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_TimeConversions
  SOURCES "${BASIC_ASTRO_TEST_DIR}/unitTestTimeConversions.cpp"
  LINK_LIBRARIES tudat_basic_astrodynamics ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_CelestialBodyConstants
  SOURCES "${BASIC_ASTRO_TEST_DIR}/unitTestCelestialBodyConstants.cpp"
  LINK_LIBRARIES tudat_basic_astrodynamics ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_ModifiedEquinoctialElementConversions
  SOURCES "${BASIC_ASTRO_TEST_DIR}/unitTestModifiedEquinoctialElementConversions.cpp"
  LINK_LIBRARIES tudat_basic_astrodynamics tudat_basic_mathematics ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_GeodeticCoordinateConversions
  SOURCES "${BASIC_ASTRO_TEST_DIR}/unitTestGeodeticCoordinateConversions.cpp"
  LINK_LIBRARIES tudat_basic_astrodynamics tudat_basic_mathematics ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_StateElementConversions
  SOURCES "${BASIC_ASTRO_TEST_DIR}/unitTestStateConversions.cpp"
  LINK_LIBRARIES tudat_basic_astrodynamics tudat_basic_mathematics ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_SphericalOrbitalStateConversions
  SOURCES "${BASIC_ASTRO_TEST_DIR}/unitTestSphericalOrbitStateConversions.cpp"
  LINK_LIBRARIES tudat_basic_astrodynamics
                 tudat_reference_frames
                 tudat_basic_mathematics
                 ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_BodyShapeModels
  SOURCES "${BASIC_ASTRO_TEST_DIR}/unitTestBodyShapeModels.cpp"
  LINK_LIBRARIES tudat_basic_astrodynamics tudat_basic_mathematics ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_UnifiedStateModelQuaternionElementConversions
  SOURCES "${BASIC_ASTRO_TEST_DIR}/unitTestUnifiedStateModelQuaternionElementConversions.cpp"
  LINK_LIBRARIES tudat_basic_astrodynamics tudat_basic_mathematics ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_UnifiedStateModelMRPElementConversions
  SOURCES "${BASIC_ASTRO_TEST_DIR}/unitTestUnifiedStateModelMRPElementConversions.cpp"
  LINK_LIBRARIES tudat_basic_astrodynamics tudat_basic_mathematics ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_UnifiedStateModelEMElementConversions
  SOURCES "${BASIC_ASTRO_TEST_DIR}/unitTestUnifiedStateModelEMElementConversions.cpp"
  LINK_LIBRARIES tudat_basic_astrodynamics tudat_basic_mathematics ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_EmpiricalAccelerations
  SOURCES "${BASIC_ASTRO_TEST_DIR}/unitTestEmpiricalAcceleration.cpp"
  LINK_LIBRARIES ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})
