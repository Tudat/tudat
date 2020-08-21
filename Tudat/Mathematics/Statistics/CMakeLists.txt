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
set(STATISTICS_SOURCES
  "${SRCROOT}${MATHEMATICSDIR}/Statistics/basicStatistics.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/Statistics/simpleLinearRegression.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/Statistics/multiVariateGaussianProbabilityDistributions.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/Statistics/continuousProbabilityDistributions.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/Statistics/boostProbabilityDistributions.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/Statistics/kernelDensityDistribution.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/Statistics/randomSampling.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/Statistics/randomVariableGenerator.cpp"
)

# Add header files.
set(STATISTICS_HEADERS
  "${SRCROOT}${MATHEMATICSDIR}/Statistics/basicStatistics.h"
  "${SRCROOT}${MATHEMATICSDIR}/Statistics/simpleLinearRegression.h"
  "${SRCROOT}${MATHEMATICSDIR}/Statistics/multiVariateGaussianProbabilityDistributions.h"
  "${SRCROOT}${MATHEMATICSDIR}/Statistics/continuousProbabilityDistributions.h"
  "${SRCROOT}${MATHEMATICSDIR}/Statistics/boostProbabilityDistributions.h"
  "${SRCROOT}${MATHEMATICSDIR}/Statistics/kernelDensityDistribution.h"
  "${SRCROOT}${MATHEMATICSDIR}/Statistics/randomSampling.h"
  "${SRCROOT}${MATHEMATICSDIR}/Statistics/randomVariableGenerator.h"
)

# Add static libraries.
add_library(tudat_statistics STATIC ${STATISTICS_SOURCES} ${STATISTICS_HEADERS})
setup_tudat_library_target(tudat_statistics "${SRCROOT}${MATHEMATICSDIR}/Statistics")

# Add unit tests.
set(MATHS_STAT_TEST_DIR "${SRCROOT}${MATHEMATICSDIR}/Statistics/UnitTests")

tudat_add_test(
  TARGET test_SimpleLinearRegression
  SOURCES "${MATHS_STAT_TEST_DIR}/unitTestSimpleLinearRegression.cpp"
  LINK_LIBRARIES tudat_statistics ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_BasicStatistics
  SOURCES "${MATHS_STAT_TEST_DIR}/unitTestBasicStatistics.cpp"
  LINK_LIBRARIES tudat_statistics tudat_basics ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_MultiVariateProbabilityDistributions
  SOURCES "${MATHS_STAT_TEST_DIR}/unitTestMultiVariateGaussianProbabilityDistributions.cpp"
  LINK_LIBRARIES tudat_statistics ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_BoostDistributions
  SOURCES "${MATHS_STAT_TEST_DIR}/unitTestBoostRandomVariables.cpp"
  LINK_LIBRARIES tudat_statistics ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_KernelDensityDistribution
  SOURCES "${MATHS_STAT_TEST_DIR}/unitTestKernelDensityDistribution.cpp"
  LINK_LIBRARIES tudat_statistics ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_RandomSampling
  SOURCES "${MATHS_STAT_TEST_DIR}/unitTestRandomSampling.cpp"
  LINK_LIBRARIES tudat_statistics ${TUDAT_EXTERNAL_LIBRARIES} ${Boost_LIBRARIES})
