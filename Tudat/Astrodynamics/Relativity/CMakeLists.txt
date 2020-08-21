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
set(RELATIVITY_SOURCES
  "${SRCROOT}${RELATIVITYDIR}/relativisticLightTimeCorrection.cpp"
  "${SRCROOT}${RELATIVITYDIR}/relativisticAccelerationCorrection.cpp"
  "${SRCROOT}${RELATIVITYDIR}/relativisticTimeConversion.cpp"
  "${SRCROOT}${RELATIVITYDIR}/metric.cpp"
)

# Set the header files.
set(RELATIVITY_HEADERS
  "${SRCROOT}${RELATIVITYDIR}/relativisticLightTimeCorrection.h"
  "${SRCROOT}${RELATIVITYDIR}/relativisticAccelerationCorrection.h"
  "${SRCROOT}${RELATIVITYDIR}/relativisticTimeConversion.h"
  "${SRCROOT}${RELATIVITYDIR}/metric.h"
)

# Add static libraries.
add_library(tudat_relativity STATIC ${RELATIVITY_SOURCES} ${RELATIVITY_HEADERS})
setup_tudat_library_target(tudat_relativity "${SRCROOT}${RELATIVITYDIR}")

# Add unit tests.
set(EINSTEIN_TEST_DIR "${SRCROOT}${RELATIVITYDIR}/UnitTests")
if( BUILD_WITH_ESTIMATION_TOOLS )
  tudat_add_test(
    TARGET test_ShapiroTimeDelay
    SOURCES "${EINSTEIN_TEST_DIR}/unitTestShapiroTimeDelay.cpp"
    LINK_LIBRARIES tudat_relativity
                   tudat_gravitation
                   tudat_observation_models
                   tudat_ephemerides
                   ${Boost_LIBRARIES})
endif( )

if( BUILD_WITH_PROPAGATION_TESTS )
  tudat_add_test(
    TARGET test_RelativisticAccelerationCorrection
    SOURCES "${EINSTEIN_TEST_DIR}/unitTestRelativisticAccelerationCorrection.cpp"
    LINK_LIBRARIES ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})
endif( )
