#   Copyright (c) 2010-2019, Delft University of Technology
#   All rigths reserved
#
#   This file is part of the Tudat. Redistribution and use in source and
#   binary forms, with or without modification, are permitted exclusively
#   under the terms of the Modified BSD license. You should have received
#   a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

# Add header files.
set(NUMERICAL_QUADRATURE_HEADERS
  "${SRCROOT}${MATHEMATICSDIR}/NumericalQuadrature/numericalQuadrature.h"
  "${SRCROOT}${MATHEMATICSDIR}/NumericalQuadrature/trapezoidQuadrature.h"
  "${SRCROOT}${MATHEMATICSDIR}/NumericalQuadrature/gaussianQuadrature.h"
  "${SRCROOT}${MATHEMATICSDIR}/NumericalQuadrature/createNumericalQuadrature.h"
)

set(NUMERICAL_QUADRATURE_SOURCES
  "${SRCROOT}${MATHEMATICSDIR}/NumericalQuadrature/dummyNumericalQuadrature.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/NumericalQuadrature/gaussianQuadrature.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/NumericalQuadrature/createNumericalQuadrature.cpp"
)


# Add static libraries.
add_library(tudat_numerical_quadrature STATIC ${NUMERICAL_QUADRATURE_SOURCES} ${NUMERICAL_QUADRATURE_HEADERS} )
setup_tudat_library_target(tudat_numerical_quadrature "${SRCROOT}${MATHEMATICSDIR}/NumericalQuadrature")
set_target_properties(tudat_numerical_quadrature PROPERTIES LINKER_LANGUAGE CXX)

# Add unit tests.
set(MATHS_QUAD_TEST_DIR "${SRCROOT}${MATHEMATICSDIR}/NumericalQuadrature/UnitTests")

tudat_add_test(
  TARGET test_TrapezoidalIntegrator
  SOURCES "${MATHS_QUAD_TEST_DIR}/unitTestTrapezoidQuadrature.cpp"
  LINK_LIBRARIES tudat_numerical_quadrature ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_GaussianQuadrature
  SOURCES "${MATHS_QUAD_TEST_DIR}/unitTestGaussianQuadrature.cpp"
  LINK_LIBRARIES tudat_input_output tudat_numerical_quadrature ${Boost_LIBRARIES})
