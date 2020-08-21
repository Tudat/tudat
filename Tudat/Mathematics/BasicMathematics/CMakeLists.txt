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
set(BASICMATHEMATICS_SOURCES
  "${SRCROOT}${BASICMATHEMATICSDIR}/coordinateConversions.cpp"
  "${SRCROOT}${BASICMATHEMATICSDIR}/legendrePolynomials.cpp"
  "${SRCROOT}${BASICMATHEMATICSDIR}/nearestNeighbourSearch.cpp"
  "${SRCROOT}${BASICMATHEMATICSDIR}/numericalDerivative.cpp"
  "${SRCROOT}${BASICMATHEMATICSDIR}/sphericalHarmonics.cpp"
  "${SRCROOT}${BASICMATHEMATICSDIR}/rotationAboutArbitraryAxis.cpp"
  "${SRCROOT}${BASICMATHEMATICSDIR}/basicMathematicsFunctions.cpp"
  "${SRCROOT}${BASICMATHEMATICSDIR}/coordinateConversions.cpp"
  "${SRCROOT}${BASICMATHEMATICSDIR}/linearAlgebra.cpp"
  "${SRCROOT}${BASICMATHEMATICSDIR}/leastSquaresEstimation.cpp"
  "${SRCROOT}${BASICMATHEMATICSDIR}/rotationRepresentations.cpp"
)

# Add header files.
set(BASICMATHEMATICS_HEADERS
  "${SRCROOT}${BASICMATHEMATICSDIR}/basicFunction.h"
  "${SRCROOT}${BASICMATHEMATICSDIR}/convergenceException.h"
  "${SRCROOT}${BASICMATHEMATICSDIR}/coordinateConversions.h"
  "${SRCROOT}${BASICMATHEMATICSDIR}/function.h"
  "${SRCROOT}${BASICMATHEMATICSDIR}/functionProxy.h"
  "${SRCROOT}${BASICMATHEMATICSDIR}/legendrePolynomials.h"
  "${SRCROOT}${BASICMATHEMATICSDIR}/linearAlgebra.h"
  "${SRCROOT}${BASICMATHEMATICSDIR}/nearestNeighbourSearch.h"
  "${SRCROOT}${BASICMATHEMATICSDIR}/numericalDerivative.h"
  "${SRCROOT}${BASICMATHEMATICSDIR}/sphericalHarmonics.h"
  "${SRCROOT}${BASICMATHEMATICSDIR}/rotationAboutArbitraryAxis.h"
  "${SRCROOT}${BASICMATHEMATICSDIR}/basicMathematicsFunctions.h"
  "${SRCROOT}${BASICMATHEMATICSDIR}/coordinateConversions.h"
  "${SRCROOT}${BASICMATHEMATICSDIR}/linearAlgebra.h"
  "${SRCROOT}${BASICMATHEMATICSDIR}/mathematicalConstants.h"
  "${SRCROOT}${BASICMATHEMATICSDIR}/leastSquaresEstimation.h"
  "${SRCROOT}${BASICMATHEMATICSDIR}/rotationRepresentations.h"
)

# Add static libraries.
add_library(tudat_basic_mathematics STATIC ${BASICMATHEMATICS_SOURCES} ${BASICMATHEMATICS_HEADERS})
setup_tudat_library_target(tudat_basic_mathematics "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics")

# Add unit tests.
set(BASIC_MATHS_TEST_DIR "${SRCROOT}${BASICMATHEMATICSDIR}/UnitTests")

tudat_add_test(
  TARGET test_MathematicalConstants
  SOURCES "${BASIC_MATHS_TEST_DIR}/unitTestMathematicalConstants.cpp"
  LINK_LIBRARIES tudat_basic_mathematics ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_LinearAlgebra
  SOURCES "${BASIC_MATHS_TEST_DIR}/unitTestLinearAlgebra.cpp"
  LINK_LIBRARIES tudat_basic_mathematics ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_CoordinateConversions
  SOURCES "${BASIC_MATHS_TEST_DIR}/unitTestCoordinateConversions.cpp"
  LINK_LIBRARIES tudat_basic_mathematics ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_NearestNeighbourSearch
  SOURCES "${BASIC_MATHS_TEST_DIR}/unitTestNearestNeighbourSearch.cpp"
  LINK_LIBRARIES tudat_input_output tudat_basic_mathematics ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_NumericalDerivative
  SOURCES "${BASIC_MATHS_TEST_DIR}/unitTestNumericalDerivative.cpp"
  LINK_LIBRARIES tudat_basic_mathematics ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_LegendrePolynomials
  SOURCES "${BASIC_MATHS_TEST_DIR}/unitTestLegendrePolynomials.cpp"
  LINK_LIBRARIES tudat_basic_mathematics ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_SphericalHarmonics
  SOURCES "${BASIC_MATHS_TEST_DIR}/unitTestSphericalHarmonics.cpp"
  LINK_LIBRARIES tudat_basic_mathematics ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_RotationAboutArbitraryAxis
  SOURCES "${BASIC_MATHS_TEST_DIR}/unitTestRotationAboutArbitraryAxis.cpp"
  LINK_LIBRARIES tudat_basic_mathematics ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_RotationPartials
  SOURCES "${BASIC_MATHS_TEST_DIR}/unitTestRotationPartials.cpp"
  LINK_LIBRARIES tudat_reference_frames tudat_basic_mathematics ${Boost_LIBRARIES})
