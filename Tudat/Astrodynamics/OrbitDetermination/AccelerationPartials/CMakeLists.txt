 #    Copyright (c) 2010-2019, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.

# Set the source files.
set(ACCELERATION_PARTIALS_SOURCES
  "${SRCROOT}${ACCELERATIONPARTIALSDIR}/aerodynamicAccelerationPartial.cpp"
  "${SRCROOT}${ACCELERATIONPARTIALSDIR}/centralGravityAccelerationPartial.cpp"
  "${SRCROOT}${ACCELERATIONPARTIALSDIR}/numericalAccelerationPartial.cpp"
  "${SRCROOT}${ACCELERATIONPARTIALSDIR}/relativisticAccelerationPartial.cpp"
  "${SRCROOT}${ACCELERATIONPARTIALSDIR}/radiationPressureAccelerationPartial.cpp"
  "${SRCROOT}${ACCELERATIONPARTIALSDIR}/sphericalHarmonicPartialFunctions.cpp"
  "${SRCROOT}${ACCELERATIONPARTIALSDIR}/sphericalHarmonicAccelerationPartial.cpp"
  "${SRCROOT}${ACCELERATIONPARTIALSDIR}/mutualSphericalHarmonicGravityPartial.cpp"
  "${SRCROOT}${ACCELERATIONPARTIALSDIR}/empiricalAccelerationPartial.cpp"
  "${SRCROOT}${ACCELERATIONPARTIALSDIR}/tidalLoveNumberPartialInterface.cpp"
  "${SRCROOT}${ACCELERATIONPARTIALSDIR}/directTidalDissipationAccelerationPartial.cpp"
  "${SRCROOT}${ACCELERATIONPARTIALSDIR}/panelledRadiationPressureAccelerationPartial.cpp"
  "${SRCROOT}${ACCELERATIONPARTIALSDIR}/thrustAccelerationPartial.cpp"
)

# Set the header files.
set(ACCELERATION_PARTIALS_HEADERS
  "${SRCROOT}${ACCELERATIONPARTIALSDIR}/accelerationPartial.h"
  "${SRCROOT}${ACCELERATIONPARTIALSDIR}/aerodynamicAccelerationPartial.h"
  "${SRCROOT}${ACCELERATIONPARTIALSDIR}/thirdBodyGravityPartial.h"
  "${SRCROOT}${ACCELERATIONPARTIALSDIR}/centralGravityAccelerationPartial.h"
  "${SRCROOT}${ACCELERATIONPARTIALSDIR}/numericalAccelerationPartial.h"
  "${SRCROOT}${ACCELERATIONPARTIALSDIR}/relativisticAccelerationPartial.h"
  "${SRCROOT}${ACCELERATIONPARTIALSDIR}/radiationPressureAccelerationPartial.h"
  "${SRCROOT}${ACCELERATIONPARTIALSDIR}/sphericalHarmonicPartialFunctions.h"
  "${SRCROOT}${ACCELERATIONPARTIALSDIR}/sphericalHarmonicAccelerationPartial.h"
  "${SRCROOT}${ACCELERATIONPARTIALSDIR}/mutualSphericalHarmonicGravityPartial.h"
  "${SRCROOT}${ACCELERATIONPARTIALSDIR}/empiricalAccelerationPartial.h"
  "${SRCROOT}${ACCELERATIONPARTIALSDIR}/tidalLoveNumberPartialInterface.h"
  "${SRCROOT}${ACCELERATIONPARTIALSDIR}/directTidalDissipationAccelerationPartial.h"
  "${SRCROOT}${ACCELERATIONPARTIALSDIR}/panelledRadiationPressureAccelerationPartial.h"
  "${SRCROOT}${ACCELERATIONPARTIALSDIR}/thrustAccelerationPartial.h"
)

# Add static libraries.
add_library(tudat_acceleration_partials STATIC ${ACCELERATION_PARTIALS_SOURCES} ${ACCELERATION_PARTIALS_HEADERS})
setup_tudat_library_target(tudat_acceleration_partials "${SRCROOT}{ACCELERATIONPARTIALSDIR}")

# Add unit tests
set(ACCPARTIALS_TEST_DIR "${SRCROOT}${ACCELERATIONPARTIALSDIR}/UnitTests")
if(USE_CSPICE)

  tudat_add_test(
    TARGET test_AccelerationPartials
    SOURCES "${ACCPARTIALS_TEST_DIR}/unitTestAccelerationPartials.cpp"
    LINK_LIBRARIES ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

  tudat_add_test(
    TARGET test_SphericalHarmonicPartials
    SOURCES "${ACCPARTIALS_TEST_DIR}/unitTestSphericalHarmonicPartials.cpp"
    LINK_LIBRARIES ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

  tudat_add_test(
    TARGET test_MutualSphericalHarmonicPartials
    SOURCES "${ACCPARTIALS_TEST_DIR}/unitTestMutualSphericalHarmonicPartials.cpp"
    LINK_LIBRARIES ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

endif( )


