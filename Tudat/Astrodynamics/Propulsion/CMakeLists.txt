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
set(PROPULSION_SOURCES
  "${SRCROOT}${PROPULSIONDIR}/thrustGuidance.cpp"
  "${SRCROOT}${PROPULSIONDIR}/thrustFunctions.cpp"
  "${SRCROOT}${PROPULSIONDIR}/costateBasedThrustGuidance.cpp"
)

# Set the header files.
set(PROPULSION_HEADERS
  "${SRCROOT}${PROPULSIONDIR}/thrustAccelerationModel.h"
  "${SRCROOT}${PROPULSIONDIR}/thrustGuidance.h"
  "${SRCROOT}${PROPULSIONDIR}/thrustFunctions.h"
  "${SRCROOT}${PROPULSIONDIR}/thrustMagnitudeWrapper.h"
  "${SRCROOT}${PROPULSIONDIR}/massRateFromThrust.h"
  "${SRCROOT}${PROPULSIONDIR}/costateBasedThrustGuidance.h"
)

# Add static libraries.
add_library(tudat_propulsion STATIC ${PROPULSION_SOURCES} ${PROPULSION_HEADERS})
setup_tudat_library_target(tudat_propulsion "${SRCROOT}${PROPULSIONDIR}")

# Add unit tests.
tudat_add_test(
  TARGET test_ThrustAcceleration
  SOURCES "${SRCROOT}${PROPULSIONDIR}/UnitTests/unitTestThrustAcceleration.cpp"
  LINK_LIBRARIES ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})
