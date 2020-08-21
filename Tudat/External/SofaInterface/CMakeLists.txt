#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.

#Set the source files.
set(SOFAINTERFACE_SOURCES
  "${SRCROOT}${EXTERNALDIR}/SofaInterface/sofaTimeConversions.h"
  "${SRCROOT}${EXTERNALDIR}/SofaInterface/fundamentalArguments.h"
  "${SRCROOT}${EXTERNALDIR}/SofaInterface/earthOrientation.h"

)

# Set the header files.
set(SOFAINTERFACE_HEADERS
"${SRCROOT}${EXTERNALDIR}/SofaInterface/sofaTimeConversions.cpp"
"${SRCROOT}${EXTERNALDIR}/SofaInterface/fundamentalArguments.cpp"
  "${SRCROOT}${EXTERNALDIR}/SofaInterface/earthOrientation.cpp"
)

# Add static libraries.
add_library(tudat_sofa_interface STATIC ${SOFAINTERFACE_SOURCES} ${SOFAINTERFACE_HEADERS})
setup_tudat_library_target(tudat_sofa_interface "${SRCROOT}${SOFAINTERFACEDIR}")

# Add unit tests.
set(SOFAIF_TEST_DIR "${SRCROOT}${EXTERNALDIR}/SofaInterface/UnitTests")

tudat_add_test(
  TARGET test_SofaEarthOrientation
  SOURCES "${SOFAIF_TEST_DIR}/unitTestSofaEarthOrientationInterface.cpp"
  LINK_LIBRARIES tudat_sofa_interface tudat_basic_mathematics sofa ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_SofaFundamentalArguments
  SOURCES "${SOFAIF_TEST_DIR}/unitTestFundamentalArguments.cpp"
  LINK_LIBRARIES tudat_sofa_interface tudat_basic_astrodynamics sofa ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_SofaTimeConversions
  SOURCES "${SOFAIF_TEST_DIR}/unitTestSofaTimeConversions.cpp"
  LINK_LIBRARIES tudat_basic_astrodynamics
                 tudat_basic_mathematics
                 tudat_sofa_interface sofa
                 ${Boost_LIBRARIES})

# End of file.
