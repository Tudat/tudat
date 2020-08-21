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
set(ELECTROMAGNETISM_SOURCES
  "${SRCROOT}${ELECTROMAGNETISMDIR}/cannonBallRadiationPressureAcceleration.h"
  "${SRCROOT}${ELECTROMAGNETISMDIR}/cannonBallRadiationPressureForce.h"
  "${SRCROOT}${ELECTROMAGNETISMDIR}/lorentzStaticMagneticForce.h"
  "${SRCROOT}${ELECTROMAGNETISMDIR}/lorentzStaticMagneticAcceleration.h"
  "${SRCROOT}${ELECTROMAGNETISMDIR}/radiationPressureInterface.h"
  "${SRCROOT}${ELECTROMAGNETISMDIR}/basicElectroMagnetism.h"
  "${SRCROOT}${ELECTROMAGNETISMDIR}/panelledRadiationPressure.h"
  "${SRCROOT}${ELECTROMAGNETISMDIR}/solarSailAcceleration.h"
  "${SRCROOT}${ELECTROMAGNETISMDIR}/solarSailForce.h"
)

# Set the header files.
set(ELECTROMAGNETISM_HEADERS
  "${SRCROOT}${ELECTROMAGNETISMDIR}/cannonBallRadiationPressureAcceleration.cpp"
  "${SRCROOT}${ELECTROMAGNETISMDIR}/cannonBallRadiationPressureForce.cpp"
  "${SRCROOT}${ELECTROMAGNETISMDIR}/lorentzStaticMagneticForce.cpp"
  "${SRCROOT}${ELECTROMAGNETISMDIR}/lorentzStaticMagneticAcceleration.cpp"
  "${SRCROOT}${ELECTROMAGNETISMDIR}/radiationPressureInterface.cpp"
  "${SRCROOT}${ELECTROMAGNETISMDIR}/panelledRadiationPressure.cpp"
  "${SRCROOT}${ELECTROMAGNETISMDIR}/solarSailAcceleration.cpp"
  "${SRCROOT}${ELECTROMAGNETISMDIR}/solarSailForce.cpp"
)

# Add static libraries.
add_library(tudat_electro_magnetism STATIC ${ELECTROMAGNETISM_SOURCES} ${ELECTROMAGNETISM_HEADERS})
setup_tudat_library_target(tudat_electro_magnetism "${SRCROOT}${ELECTROMAGNETISMDIR}")

# Add unit tests.
set(EM_TEST_DIR "${SRCROOT}${ELECTROMAGNETISMDIR}/UnitTests")

tudat_add_test(
  TARGET test_CannonBallRadiationPressureAccelerationAndForce
  SOURCES "${EM_TEST_DIR}/unitTestCannonBallRadiationPressureAccelerationAndForce.cpp"
  LINK_LIBRARIES tudat_electro_magnetism tudat_basic_astrodynamics ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_LorentzStaticMagneticAccelerationAndForce
  SOURCES "${EM_TEST_DIR}/unitTestLorentzStaticMagneticAccelerationAndForce.cpp"
  LINK_LIBRARIES tudat_electro_magnetism tudat_basic_astrodynamics ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_RadiationPressureInterface
  SOURCES "${EM_TEST_DIR}/unitTestRadiationPressureInterface.cpp"
  LINK_LIBRARIES tudat_electro_magnetism tudat_basic_astrodynamics ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_PanelledRadiationPressure
  SOURCES "${EM_TEST_DIR}/unitTestPanelledRadiationPressure.cpp"
  LINK_LIBRARIES ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_SolarSailAccelerationAndForce
  SOURCES "${EM_TEST_DIR}/unitTestSolarSailAccelerationAndForce.cpp"
  LINK_LIBRARIES tudat_electro_magnetism tudat_basic_astrodynamics ${Boost_LIBRARIES})
