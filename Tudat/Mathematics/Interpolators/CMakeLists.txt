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
set(INTERPOLATORS_SOURCES
  "${SRCROOT}${MATHEMATICSDIR}/Interpolators/cubicSplineInterpolator.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/Interpolators/linearInterpolator.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/Interpolators/lagrangeInterpolator.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/Interpolators/interpolator.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/Interpolators/multiLinearInterpolator.cpp"
)

# Add header files.
set(INTERPOLATORS_HEADERS
  "${SRCROOT}${MATHEMATICSDIR}/Interpolators/cubicSplineInterpolator.h"
  "${SRCROOT}${MATHEMATICSDIR}/Interpolators/hermiteCubicSplineInterpolator.h"
  "${SRCROOT}${MATHEMATICSDIR}/Interpolators/linearInterpolator.h"
  "${SRCROOT}${MATHEMATICSDIR}/Interpolators/lagrangeInterpolator.h"
  "${SRCROOT}${MATHEMATICSDIR}/Interpolators/interpolator.h"
  "${SRCROOT}${MATHEMATICSDIR}/Interpolators/lookupScheme.h"
  "${SRCROOT}${MATHEMATICSDIR}/Interpolators/multiDimensionalInterpolator.h"
  "${SRCROOT}${MATHEMATICSDIR}/Interpolators/oneDimensionalInterpolator.h"
  "${SRCROOT}${MATHEMATICSDIR}/Interpolators/multiLinearInterpolator.h"
  "${SRCROOT}${MATHEMATICSDIR}/Interpolators/piecewiseConstantInterpolator.h"
  "${SRCROOT}${MATHEMATICSDIR}/Interpolators/jumpDataLinearInterpolator.h"
  "${SRCROOT}${MATHEMATICSDIR}/Interpolators/createInterpolator.h"
)

# Add static libraries.
add_library(tudat_interpolators STATIC ${INTERPOLATORS_SOURCES} ${INTERPOLATORS_HEADERS})
setup_tudat_library_target(tudat_interpolators "${SRCROOT}${MATHEMATICSDIR}/Interpolators")

# Add unit tests.
set(MATHS_INTERP_TEST_DIR "${SRCROOT}${MATHEMATICSDIR}/Interpolators/UnitTests")

tudat_add_test(
  TARGET test_CubicSplineInterpolator
  SOURCES "${MATHS_INTERP_TEST_DIR}/unitTestCubicSplineInterpolator.cpp"
  LINK_LIBRARIES tudat_input_output
                 tudat_interpolators
                 tudat_basic_mathematics
                 ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_HermiteCubicSplineInterpolator
  SOURCES "${MATHS_INTERP_TEST_DIR}/unitTestHermiteCubicSplineInterpolator.cpp"
  LINK_LIBRARIES tudat_input_output
                 tudat_interpolators
                 tudat_basic_mathematics
                 ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_LinearInterpolator
  SOURCES "${MATHS_INTERP_TEST_DIR}/unitTestLinearInterpolator.cpp"
  LINK_LIBRARIES tudat_input_output
                 tudat_interpolators
                 tudat_basic_mathematics
                 ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_MultiLinearInterpolator
  SOURCES "${MATHS_INTERP_TEST_DIR}/unitTestMultiLinearInterpolator.cpp"
  LINK_LIBRARIES tudat_input_output
                 tudat_interpolators
                 tudat_basic_mathematics
                 ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_LagrangeInterpolator
  SOURCES "${MATHS_INTERP_TEST_DIR}/unitTestLagrangeInterpolators.cpp"
  LINK_LIBRARIES tudat_input_output
                 tudat_interpolators
                 tudat_basic_mathematics
                 ${Boost_LIBRARIES})


