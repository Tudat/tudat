 #    Copyright (c) 2010-2018, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.
 #

# Set the source files.
if( USE_PAGMO )
	set(SHAPEBASEDMETHODS_SOURCES
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/baseFunctionsHodographicShaping.cpp"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/createBaseFunctionHodographicShaping.cpp"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/compositeFunctionHodographicShaping.cpp"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/hodographicShaping.cpp"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/baseFunctionsSphericalShaping.cpp"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/compositeFunctionSphericalShaping.cpp"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/sphericalShaping.cpp"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/baseFunctionsExposinsShaping.cpp"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/compositeFunctionExposinsShaping.cpp"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/exposinsShaping.cpp"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/baseFunctionsInvPolyShaping.cpp"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/compositeFunctionInvPolyShaping.cpp"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/invPolyShaping.cpp"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/shapeBasedMethodLeg.cpp"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/invpolyShapingOptimisationSetup.cpp"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/shapeBasedCoastingOptimisationSetup.cpp"
	)
else( )
	set(SHAPEBASEDMETHODS_SOURCES
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/baseFunctionsHodographicShaping.cpp"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/createBaseFunctionHodographicShaping.cpp"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/compositeFunctionHodographicShaping.cpp"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/hodographicShaping.cpp"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/baseFunctionsSphericalShaping.cpp"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/compositeFunctionSphericalShaping.cpp"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/sphericalShaping.cpp"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/baseFunctionsExposinsShaping.cpp"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/compositeFunctionExposinsShaping.cpp"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/exposinsShaping.cpp"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/baseFunctionsInvPolyShaping.cpp"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/compositeFunctionInvPolyShaping.cpp"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/invPolyShaping.cpp"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/shapeBasedMethodLeg.cpp"
	)
endif( )

if( USE_PAGMO )
# Set the header files.
	set(SHAPEBASEDMETHODS_HEADERS
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/baseFunctionsHodographicShaping.h"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/createBaseFunctionHodographicShaping.h"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/compositeFunctionHodographicShaping.h"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/hodographicShaping.h"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/baseFunctionsSphericalShaping.h"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/compositeFunctionSphericalShaping.h"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/sphericalShaping.h"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/baseFunctionsExposinsShaping.h"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/compositeFunctionExposinsShaping.h"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/exposinsShaping.h"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/baseFunctionsInvPolyShaping.h"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/compositeFunctionInvPolyShaping.h"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/invPolyShaping.h"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/shapeBasedMethodLeg.h"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/invpolyShapingOptimisationSetup.h"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/shapeBasedCoastingOptimisationSetup.h"
	)
else( )
	set(SHAPEBASEDMETHODS_HEADERS
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/baseFunctionsHodographicShaping.h"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/createBaseFunctionHodographicShaping.h"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/compositeFunctionHodographicShaping.h"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/hodographicShaping.h"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/baseFunctionsSphericalShaping.h"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/compositeFunctionSphericalShaping.h"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/sphericalShaping.h"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/baseFunctionsExposinsShaping.h"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/compositeFunctionExposinsShaping.h"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/exposinsShaping.h"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/baseFunctionsInvPolyShaping.h"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/compositeFunctionInvPolyShaping.h"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/invPolyShaping.h"
	  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/shapeBasedMethodLeg.h"
	)
endif( )


# Find PaGMO library on local system.
if( USE_PAGMO )
  find_package(PaGMO2)

  # Include PaGMO directories.
  if(NOT APPLE OR APPLE_INCLUDE_FORCE)
    include_directories(SYSTEM AFTER "${PAGMO_INCLUDE_DIR}")
  else( )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${PAGMO_INCLUDE_DIR}\"")
  endif( )
endif( )

# Add static libraries.
add_library(tudat_shape_based_methods STATIC ${SHAPEBASEDMETHODS_SOURCES} ${SHAPEBASEDMETHODS_HEADERS})
setup_tudat_library_target(tudat_shape_based_methods "${SRCROOT}${SHAPEBASEDMETHODSDIR}")

include(tudatLinkLibraries)

## Add unit tests.
add_executable(test_HodographicShaping "${SRCROOT}${SHAPEBASEDMETHODSDIR}/UnitTests/unitTestHodographicShaping.cpp")
setup_custom_test_program(test_HodographicShaping "${SRCROOT}${SHAPEBASEDMETHODSDIR}")
target_link_libraries(test_HodographicShaping tudat_shape_based_methods tudat_low_thrust_direct_methods tudat_numerical_quadrature ${TUDAT_PROPAGATION_LIBRARIES}  ${Boost_LIBRARIES})

add_executable(test_SphericalShaping "${SRCROOT}${SHAPEBASEDMETHODSDIR}/UnitTests/unitTestSphericalShaping.cpp")
setup_custom_test_program(test_SphericalShaping "${SRCROOT}${SHAPEBASEDMETHODSDIR}")
target_link_libraries(test_SphericalShaping tudat_shape_based_methods tudat_low_thrust_direct_methods tudat_numerical_quadrature ${TUDAT_PROPAGATION_LIBRARIES}  ${Boost_LIBRARIES})

add_executable(test_ExposinsShaping "${SRCROOT}${SHAPEBASEDMETHODSDIR}/UnitTests/unitTestExposinsShaping.cpp")
setup_custom_test_program(test_ExposinsShaping "${SRCROOT}${SHAPEBASEDMETHODSDIR}")
target_link_libraries(test_ExposinsShaping tudat_shape_based_methods tudat_low_thrust_direct_methods tudat_numerical_quadrature ${TUDAT_PROPAGATION_LIBRARIES}  ${Boost_LIBRARIES})

add_executable(test_InvPolyShaping "${SRCROOT}${SHAPEBASEDMETHODSDIR}/UnitTests/unitTestInvPolyShaping.cpp")
setup_custom_test_program(test_InvPolyShaping "${SRCROOT}${SHAPEBASEDMETHODSDIR}")
target_link_libraries(test_InvPolyShaping tudat_shape_based_methods tudat_low_thrust_direct_methods tudat_numerical_quadrature ${TUDAT_PROPAGATION_LIBRARIES}  ${Boost_LIBRARIES})

add_executable(test_unitTestInterplanetaryCases "${SRCROOT}${SHAPEBASEDMETHODSDIR}/UnitTests/unitTestInterplanetaryCases.cpp")
setup_custom_test_program(test_unitTestInterplanetaryCases "${SRCROOT}${SHAPEBASEDMETHODSDIR}")
target_link_libraries(test_unitTestInterplanetaryCases tudat_shape_based_methods tudat_low_thrust_direct_methods tudat_numerical_quadrature ${TUDAT_PROPAGATION_LIBRARIES}  ${Boost_LIBRARIES})

add_executable(test_unitTestShapeBasedCoasting "${SRCROOT}${SHAPEBASEDMETHODSDIR}/UnitTests/unitTestShapeBasedCoasting.cpp")
setup_custom_test_program(test_unitTestShapeBasedCoasting "${SRCROOT}${SHAPEBASEDMETHODSDIR}")
target_link_libraries(test_unitTestShapeBasedCoasting tudat_shape_based_methods tudat_low_thrust_direct_methods tudat_numerical_quadrature ${TUDAT_PROPAGATION_LIBRARIES}  ${Boost_LIBRARIES})



