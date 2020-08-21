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
set(AERODYNAMICS_SOURCES
  "${SRCROOT}${AERODYNAMICSDIR}/aerodynamicAcceleration.cpp"
  "${SRCROOT}${AERODYNAMICSDIR}/aerodynamicCoefficientGenerator.cpp"
  "${SRCROOT}${AERODYNAMICSDIR}/aerodynamicTorque.cpp"
  "${SRCROOT}${AERODYNAMICSDIR}/aerodynamicForce.cpp"
  "${SRCROOT}${AERODYNAMICSDIR}/aerodynamics.cpp"
  "${SRCROOT}${AERODYNAMICSDIR}/customConstantTemperatureAtmosphere.cpp"
  "${SRCROOT}${AERODYNAMICSDIR}/exponentialAtmosphere.cpp"
  "${SRCROOT}${AERODYNAMICSDIR}/hypersonicLocalInclinationAnalysis.cpp"
  "${SRCROOT}${AERODYNAMICSDIR}/tabulatedAtmosphere.cpp"
  "${SRCROOT}${AERODYNAMICSDIR}/flightConditions.cpp"
  "${SRCROOT}${AERODYNAMICSDIR}/trimOrientation.cpp"
  "${SRCROOT}${AERODYNAMICSDIR}/equilibriumWallTemperature.cpp"
)

# Set the header files.
set(AERODYNAMICS_HEADERS
  "${SRCROOT}${AERODYNAMICSDIR}/aerodynamicAcceleration.h"
  "${SRCROOT}${AERODYNAMICSDIR}/aerodynamicTorque.h"
  "${SRCROOT}${AERODYNAMICSDIR}/aerodynamicCoefficientGenerator.h"
  "${SRCROOT}${AERODYNAMICSDIR}/aerodynamicCoefficientInterface.h"
  "${SRCROOT}${AERODYNAMICSDIR}/aerodynamicAcceleration.h"
  "${SRCROOT}${AERODYNAMICSDIR}/aerodynamicForce.h"
  "${SRCROOT}${AERODYNAMICSDIR}/aerodynamics.h"
  "${SRCROOT}${AERODYNAMICSDIR}/atmosphereModel.h"
  "${SRCROOT}${AERODYNAMICSDIR}/customConstantTemperatureAtmosphere.h"
  "${SRCROOT}${AERODYNAMICSDIR}/exponentialAtmosphere.h"
  "${SRCROOT}${AERODYNAMICSDIR}/hypersonicLocalInclinationAnalysis.h"
  "${SRCROOT}${AERODYNAMICSDIR}/tabulatedAtmosphere.h"
  "${SRCROOT}${AERODYNAMICSDIR}/standardAtmosphere.h"
  "${SRCROOT}${AERODYNAMICSDIR}/customAerodynamicCoefficientInterface.h"
  "${SRCROOT}${AERODYNAMICSDIR}/controlSurfaceAerodynamicCoefficientInterface.h"
  "${SRCROOT}${AERODYNAMICSDIR}/flightConditions.h"
  "${SRCROOT}${AERODYNAMICSDIR}/UnitTests/testApolloCapsuleCoefficients.h"
  "${SRCROOT}${AERODYNAMICSDIR}/flightConditions.cpp"
  "${SRCROOT}${AERODYNAMICSDIR}/aerodynamicGuidance.h"
  "${SRCROOT}${AERODYNAMICSDIR}/equilibriumWallTemperature.h"
  "${SRCROOT}${AERODYNAMICSDIR}/windModel.h"
)

if(USE_NRLMSISE00)
  set(AERODYNAMICS_SOURCES "${AERODYNAMICS_SOURCES}"
    "${SRCROOT}${AERODYNAMICSDIR}/nrlmsise00Atmosphere.cpp"
    "${SRCROOT}${AERODYNAMICSDIR}/nrlmsise00InputFunctions.cpp")
  set(AERODYNAMICS_HEADERS "${AERODYNAMICS_HEADERS}"
    "${SRCROOT}${AERODYNAMICSDIR}/nrlmsise00Atmosphere.h"
    "${SRCROOT}${AERODYNAMICSDIR}/nrlmsise00InputFunctions.h")
endif( )

# Add static libraries
add_library(tudat_aerodynamics STATIC ${AERODYNAMICS_SOURCES} ${AERODYNAMICS_HEADERS})
setup_tudat_library_target(tudat_aerodynamics "${SRCROOT}{AERODYNAMICSDIR}")

# Add unit tests
set(AERO_TEST_DIR "${SRCROOT}${AERODYNAMICSDIR}/UnitTests")

tudat_add_test(
  TARGET test_AerodynamicMomentAndAerodynamicForce
  SOURCES "${AERO_TEST_DIR}/unitTestAerodynamicMomentAndAerodynamicForce.cpp"
  LINK_LIBRARIES ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_AerodynamicsNamespace
  SOURCES "${AERO_TEST_DIR}/unitTestAerodynamicsNamespace.cpp"
  LINK_LIBRARIES tudat_aerodynamics tudat_basic_mathematics ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_AerodynamicCoefficientGenerator
  SOURCES "${AERO_TEST_DIR}/unitTestCoefficientGenerator.cpp"
  LINK_LIBRARIES tudat_aerodynamics tudat_geometric_shapes tudat_interpolators tudat_basic_mathematics ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_ExponentialAtmosphere
  SOURCES "${AERO_TEST_DIR}/unitTestExponentialAtmosphere.cpp"
  LINK_LIBRARIES tudat_aerodynamics tudat_basic_mathematics ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_CustomConstantTemperatureAtmosphere
  SOURCES "${AERO_TEST_DIR}/unitTestCustomConstantTemperatureAtmosphere.cpp"
  LINK_LIBRARIES tudat_aerodynamics tudat_basic_mathematics ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_TabulatedAtmosphere
  SOURCES "${AERO_TEST_DIR}/unitTestTabulatedAtmosphere.cpp"
  LINK_LIBRARIES tudat_aerodynamics
                 tudat_input_output
                 tudat_interpolators
                 tudat_basic_mathematics
                 tudat_basic_astrodynamics
                 ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_TabulatedAerodynamicCoefficients
  SOURCES "${AERO_TEST_DIR}/unitTestTabulatedAerodynamicCoefficients.cpp"
  LINK_LIBRARIES ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_HeatTransfer
  SOURCES "${AERO_TEST_DIR}/unitTestHeatTransfer.cpp"
  LINK_LIBRARIES tudat_aerodynamics
                 tudat_root_finders
                 tudat_basic_mathematics
                 tudat_basic_astrodynamics
                 ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_ControlSurfaceIncrements
  SOURCES "${AERO_TEST_DIR}/unitTestControlSurfaceIncrements.cpp"
  LINK_LIBRARIES ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_AerodynamicCoefficientsFromFile
  SOURCES "${AERO_TEST_DIR}/unitTestAerodynamicCoefficientsFromFile.cpp"
  LINK_LIBRARIES ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_WindModel
  SOURCES "${AERO_TEST_DIR}/unitTestWindModel.cpp"
  LINK_LIBRARIES ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

if(USE_NRLMSISE00)
    tudat_add_test(
      TARGET test_NRLMSISE00Atmosphere
      SOURCES "${AERO_TEST_DIR}/unitTestNRLMSISE00Atmosphere.cpp"
      LINK_LIBRARIES tudat_aerodynamics
                     tudat_interpolators
                     tudat_basic_mathematics
                     nrlmsise00
                     tudat_input_output
                     tudat_basic_astrodynamics
                     ${Boost_LIBRARIES})
endif( )
