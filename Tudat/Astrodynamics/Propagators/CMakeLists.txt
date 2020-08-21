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
set(PROPAGATORS_SOURCES
  "${SRCROOT}${PROPAGATORSDIR}/nBodyStateDerivative.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/nBodyCowellStateDerivative.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/nBodyEnckeStateDerivative.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/nBodyGaussKeplerStateDerivative.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/nBodyGaussModifiedEquinoctialStateDerivative.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/dynamicsStateDerivativeModel.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/nBodyUnifiedStateModelQuaternionsStateDerivative.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/nBodyUnifiedStateModelModifiedRodriguesParametersStateDerivative.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/nBodyUnifiedStateModelExponentialMapStateDerivative.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/variationalEquations.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/stateTransitionMatrixInterface.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/environmentUpdateTypes.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/singleStateTypeDerivative.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/rotationalMotionStateDerivative.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/rotationalMotionQuaternionsStateDerivative.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/rotationalMotionModifiedRodriguesParametersStateDerivative.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/rotationalMotionExponentialMapStateDerivative.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/stateDerivativeCircularRestrictedThreeBodyProblem.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/integrateEquations.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/dynamicsStateDerivativeModel.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/propagateCovariance.cpp"
)

# Add header files.
set(PROPAGATORS_HEADERS
  "${SRCROOT}${PROPAGATORSDIR}/centralBodyData.h"
  "${SRCROOT}${PROPAGATORSDIR}/nBodyStateDerivative.h"
  "${SRCROOT}${PROPAGATORSDIR}/nBodyCowellStateDerivative.h"
  "${SRCROOT}${PROPAGATORSDIR}/nBodyEnckeStateDerivative.h"
  "${SRCROOT}${PROPAGATORSDIR}/nBodyGaussKeplerStateDerivative.h"
  "${SRCROOT}${PROPAGATORSDIR}/nBodyGaussModifiedEquinoctialStateDerivative.h"
  "${SRCROOT}${PROPAGATORSDIR}/nBodyUnifiedStateModelQuaternionsStateDerivative.h"
  "${SRCROOT}${PROPAGATORSDIR}/nBodyUnifiedStateModelModifiedRodriguesParametersStateDerivative.h"
  "${SRCROOT}${PROPAGATORSDIR}/nBodyUnifiedStateModelExponentialMapStateDerivative.h"
  "${SRCROOT}${PROPAGATORSDIR}/dynamicsStateDerivativeModel.h"
  "${SRCROOT}${PROPAGATORSDIR}/singleStateTypeDerivative.h"
  "${SRCROOT}${PROPAGATORSDIR}/integrateEquations.h"
  "${SRCROOT}${PROPAGATORSDIR}/bodyMassStateDerivative.h"
  "${SRCROOT}${PROPAGATORSDIR}/variationalEquations.h"
  "${SRCROOT}${PROPAGATORSDIR}/stateTransitionMatrixInterface.h"
  "${SRCROOT}${PROPAGATORSDIR}/environmentUpdateTypes.h"
  "${SRCROOT}${PROPAGATORSDIR}/customStateDerivative.h"
  "${SRCROOT}${PROPAGATORSDIR}/rotationalMotionStateDerivative.h"
  "${SRCROOT}${PROPAGATORSDIR}/rotationalMotionQuaternionsStateDerivative.h"
  "${SRCROOT}${PROPAGATORSDIR}/rotationalMotionModifiedRodriguesParametersStateDerivative.h"
  "${SRCROOT}${PROPAGATORSDIR}/rotationalMotionExponentialMapStateDerivative.h"
  "${SRCROOT}${PROPAGATORSDIR}/stateDerivativeCircularRestrictedThreeBodyProblem.h"
  "${SRCROOT}${PROPAGATORSDIR}/getZeroProperModeRotationalInitialState.h"
  "${SRCROOT}${PROPAGATORSDIR}/propagateCovariance.h"
)

# Add static libraries.
add_library(tudat_propagators STATIC ${PROPAGATORS_SOURCES} ${PROPAGATORS_HEADERS})
setup_tudat_library_target(tudat_propagators "${SRCROOT}${PROPAGATORSDIR}")


include(tudatLinkLibraries)

# Add unit tests.
set(PROPAGATORS_TEST_DIR "${SRCROOT}${PROPAGATORSDIR}/UnitTests")

tudat_add_test(
  TARGET test_CentralBodyData
  SOURCES "${PROPAGATORS_TEST_DIR}/unitTestCentralBodyData.cpp"
  LINK_LIBRARIES tudat_propagators tudat_basic_mathematics ${Boost_LIBRARIES})

if(USE_CSPICE)

  if( BUILD_WITH_PROPAGATION_TESTS )

    tudat_add_test(
      TARGET test_CowellStateDerivative
      SOURCES "${PROPAGATORS_TEST_DIR}/unitTestCowellStateDerivative.cpp"
      LINK_LIBRARIES ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

    tudat_add_test(
      TARGET test_EnckeStateDerivative
      SOURCES "${PROPAGATORS_TEST_DIR}/unitTestEnckeStateDerivative.cpp"
      LINK_LIBRARIES ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

    tudat_add_test(
      TARGET test_GaussStateDerivative
      SOURCES "${PROPAGATORS_TEST_DIR}/unitTestGaussStateDerivative.cpp"
      LINK_LIBRARIES ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

    tudat_add_test(
      TARGET test_UnifiedStateModelStateDerivative
      SOURCES "${PROPAGATORS_TEST_DIR}/unitTestUnifiedStateModelStateDerivative.cpp"
      LINK_LIBRARIES ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

    if( BUILD_WITH_ESTIMATION_TOOLS )

      tudat_add_test(
        TARGET test_SequentialVariationEquationIntegration
        SOURCES "${PROPAGATORS_TEST_DIR}/unitTestSequentialVariationalEquationIntegration.cpp"
        LINK_LIBRARIES ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

      tudat_add_test(
        TARGET test_VariationalEquations
        SOURCES "${PROPAGATORS_TEST_DIR}/unitTestVariationalEquationPropagation.cpp"
        LINK_LIBRARIES ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

      tudat_add_test(
        TARGET test_MultiArcVariationalEquations
        SOURCES "${PROPAGATORS_TEST_DIR}/unitTestMultiArcVariationalEquationPropagation.cpp"
        LINK_LIBRARIES ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

      tudat_add_test(
        TARGET test_HybridArcVariationalEquations
        SOURCES "${PROPAGATORS_TEST_DIR}/unitTestHybridArcVariationalEquationPropagation.cpp"
        LINK_LIBRARIES ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

      tudat_add_test(
        TARGET test_DependentVariableOutput
        SOURCES "${PROPAGATORS_TEST_DIR}/unitTestDependentVariableOutput.cpp"
        LINK_LIBRARIES ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

    endif( )

    tudat_add_test(
      TARGET test_EnvironmentModelComputations
      SOURCES "${PROPAGATORS_TEST_DIR}/unitTestEnvironmentUpdater.cpp"
      LINK_LIBRARIES ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

    tudat_add_test(
      TARGET test_BodyMassPropagation
      SOURCES "${PROPAGATORS_TEST_DIR}/unitTestBodyMassPropagation.cpp"
      LINK_LIBRARIES ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

    tudat_add_test(
      TARGET test_MultiTypeStatePropagation
      SOURCES "${PROPAGATORS_TEST_DIR}/unitTestMultiTypeStatePropagation.cpp"
      LINK_LIBRARIES ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

    tudat_add_test(
      TARGET test_StoppingConditions
      SOURCES "${PROPAGATORS_TEST_DIR}/unitTestStoppingConditions.cpp"
      LINK_LIBRARIES ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

    tudat_add_test(
      TARGET test_CustomStatePropagation
      SOURCES "${PROPAGATORS_TEST_DIR}/unitTestCustomStatePropagation.cpp"
      LINK_LIBRARIES ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

    tudat_add_test(
      TARGET test_MultiArcDynamics
      SOURCES "${PROPAGATORS_TEST_DIR}/unitTestMultiArcDynamics.cpp"
      LINK_LIBRARIES ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

    tudat_add_test(
      TARGET test_HybridArcDynamics
      SOURCES "${PROPAGATORS_TEST_DIR}/unitTestHybridArcDynamics.cpp"
      LINK_LIBRARIES ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

    tudat_add_test(
      TARGET test_PropagationTerminationReason
      SOURCES "${PROPAGATORS_TEST_DIR}/unitTestPropagationTerminationReason.cpp"
      LINK_LIBRARIES ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

    tudat_add_test(
      TARGET test_RotationalDynamicsPropagation
      SOURCES "${PROPAGATORS_TEST_DIR}/unitTestRotationalDynamicsPropagator.cpp"
      LINK_LIBRARIES ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

    tudat_add_test(
      TARGET test_CR3BPPropagation
      SOURCES "${PROPAGATORS_TEST_DIR}/unitTestCR3BPPropagation.cpp"
      LINK_LIBRARIES ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

    tudat_add_test(
      TARGET test_PropagationTerminationCheckOnFinalStep
      SOURCES "${PROPAGATORS_TEST_DIR}/uniTestPropagationTerminationCheckOnFinalStep.cpp"
      LINK_LIBRARIES ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

    tudat_add_test(
      TARGET test_GlobalFrameOrigin
      SOURCES "${PROPAGATORS_TEST_DIR}/unitTestGlobalFrameOrigin.cpp"
      LINK_LIBRARIES ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

    tudat_add_test(
      TARGET test_ForwardBackwardIntegration
      SOURCES "${PROPAGATORS_TEST_DIR}/unitTestForwardsBackwardsIntegration.cpp"
      LINK_LIBRARIES ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

    tudat_add_test(
      TARGET test_ExactTermination
      SOURCES "${PROPAGATORS_TEST_DIR}/unitTestExactTermination.cpp"
      LINK_LIBRARIES ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

  endif( )

endif( )

tudat_add_test(
  TARGET test_StateDerivativeRestrictedThreeBodyProblem
  SOURCES "${PROPAGATORS_TEST_DIR}/unitTestStateDerivativeCircularRestrictedThreeBodyProblem.cpp"
  LINK_LIBRARIES tudat_mission_segments
                 tudat_root_finders
                 tudat_propagators
                 tudat_numerical_integrators
                 tudat_basic_astrodynamics
                 tudat_input_output
                 ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_FullPropagationRestrictedThreeBodyProblem
  SOURCES "${PROPAGATORS_TEST_DIR}/unitTestFullPropagationCircularRestrictedThreeBodyProblem.cpp"
  LINK_LIBRARIES ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_FullPropagationLambertTargeter
  SOURCES "${PROPAGATORS_TEST_DIR}/unitTestFullPropagationLambertTargeter.cpp"
  LINK_LIBRARIES tudat_trajectory_design
                 tudat_mission_segments
                 tudat_ephemerides
                 tudat_basic_astrodynamics
                 tudat_basic_mathematics
                 ${TUDAT_PROPAGATION_LIBRARIES}
                 ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_FullPropagationPatchedConicsTrajectory
  SOURCES "${PROPAGATORS_TEST_DIR}/unitTestFullPropagationPatchedConicsTrajectory.cpp"
  LINK_LIBRARIES tudat_trajectory_design
                 tudat_mission_segments
                 tudat_ephemerides
                 tudat_basic_astrodynamics
                 tudat_basic_mathematics
                 ${TUDAT_PROPAGATION_LIBRARIES}
                 ${Boost_LIBRARIES})


