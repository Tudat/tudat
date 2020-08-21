 #    Copyright (c) 2010-2019, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.
 #

# Add source files.
set(NUMERICALINTEGRATORS_SOURCES
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/rungeKuttaCoefficients.cpp"
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/createNumericalIntegrator.cpp"
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/bulirschStoerVariableStepsizeIntegrator.cpp"
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/UnitTests/numericalIntegratorTests.cpp"
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/numericalIntegrator.cpp"
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/reinitializableNumericalIntegrator.cpp"
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/rungeKutta4Integrator.cpp"
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/rungeKuttaVariableStepSizeIntegrator.cpp"
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/adamsBashforthMoultonIntegrator.cpp"
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/bulirschStoerVariableStepsizeIntegrator.cpp"
)

# Add header files.
set(NUMERICALINTEGRATORS_HEADERS
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/adamsBashforthMoultonIntegrator.h"
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/createNumericalIntegrator.h"
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/bulirschStoerVariableStepsizeIntegrator.h"
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/euler.h"
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/numericalIntegrator.h"
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/reinitializableNumericalIntegrator.h"
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/rungeKutta4Integrator.h"
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/rungeKuttaCoefficients.h"
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/rungeKuttaVariableStepSizeIntegrator.h"
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/UnitTests/burdenAndFairesNumericalIntegratorTest.h"
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/UnitTests/numericalIntegratorTests.h"
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/UnitTests/numericalIntegratorTestFunctions.h"
)

# Add static libraries.
add_library(tudat_numerical_integrators STATIC ${NUMERICALINTEGRATORS_SOURCES} ${NUMERICALINTEGRATORS_HEADERS})
setup_tudat_library_target(tudat_numerical_integrators "${SRCROOT}${NUMERICALINTEGRATORSDIR}")

# Add unit tests.
set(MATHS_DEINT_TEST_DIR "${SRCROOT}${NUMERICALINTEGRATORSDIR}/UnitTests")

tudat_add_test(
  TARGET test_AdamsBashforthMoultonIntegrator
  SOURCES "${MATHS_DEINT_TEST_DIR}/unitTestAdamsBashforthMoultonIntegrator.cpp"
  LINK_LIBRARIES tudat_numerical_integrators tudat_input_output ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_EulerIntegrator
  SOURCES "${MATHS_DEINT_TEST_DIR}/unitTestEulerIntegrator.cpp"
  LINK_LIBRARIES tudat_numerical_integrators tudat_input_output ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_NumericalIntegrator
  SOURCES "${MATHS_DEINT_TEST_DIR}/unitTestNumericalIntegrator.cpp"
  LINK_LIBRARIES tudat_numerical_integrators ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_RungeKutta4Integrator
  SOURCES "${MATHS_DEINT_TEST_DIR}/unitTestRungeKutta4Integrator.cpp"
  LINK_LIBRARIES tudat_numerical_integrators tudat_input_output ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_RungeKuttaVariableStepSizeIntegrator
  SOURCES "${MATHS_DEINT_TEST_DIR}/unitTestRungeKuttaVariableStepSizeIntegrator.cpp"
  LINK_LIBRARIES tudat_numerical_integrators ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_RungeKuttaCoefficients
  SOURCES "${MATHS_DEINT_TEST_DIR}/unitTestRungeKuttaCoefficients.cpp"
  LINK_LIBRARIES tudat_numerical_integrators ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_RungeKuttaFehlberg45Integrator
  SOURCES "${MATHS_DEINT_TEST_DIR}/unitTestRungeKuttaFehlberg45Integrator.cpp"
  LINK_LIBRARIES tudat_numerical_integrators tudat_input_output ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_RungeKuttaFehlberg56Integrator
  SOURCES "${MATHS_DEINT_TEST_DIR}/unitTestRungeKuttaFehlberg56Integrator.cpp"
  LINK_LIBRARIES tudat_numerical_integrators tudat_input_output ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_RungeKuttaFehlberg78Integrator
  SOURCES "${MATHS_DEINT_TEST_DIR}/unitTestRungeKuttaFehlberg78Integrator.cpp"
  LINK_LIBRARIES tudat_numerical_integrators tudat_input_output ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_RungeKutta87DormandPrinceIntegrator
  SOURCES "${MATHS_DEINT_TEST_DIR}/unitTestRungeKutta87DormandPrinceIntegrator.cpp"
  LINK_LIBRARIES tudat_numerical_integrators tudat_input_output ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_BulirschStoerVariableStepSizeIntegrator
  SOURCES "${MATHS_DEINT_TEST_DIR}/unitTestBulirschStoerVariableStepSizeIntegrator.cpp"
  LINK_LIBRARIES tudat_numerical_integrators tudat_input_output ${Boost_LIBRARIES})
