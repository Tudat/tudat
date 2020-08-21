#=============================================================================
#
#    tudat_add_test
#
#    This function facilitates defining a standardized tudat unit test
#
#    Usage:
#        tudat_add_test(
#            TARGET <target_name>
#            [ SOURCES <list> ]
#            [ LINK_LIBRARIES <list> ]
#            [ COMPILE_FLAGS <list> ]
#            [ INCLUDES <list> ]
#            [ DEPENDS <string> ]
#            )
#
#    Input:
#        TARGET        : test target name
#        COMPILE_FLAGS : custom compile flags
#        LINK_FLAGS    : custom link flags
#        LINK_LIBRARIES: libraries to link test with
#        SOURCES       : source files
#        INCLUDES      : include directories
#        DEPENDS       : targets which this test depends on
#
#    Output:
#        Tudat-standardized test with name matching TARGET.
#
#    Example:
#        tudat_add_test(
#            TARGET
#                my_test_application
#            SOURCES
#                file.cpp
#            LINK_LIBRARIES
#                ${Boost_LIBRARIES}
#            )
#
function(tudat_add_test)

  if (NOT TUDAT_DISABLE_TESTS)

    tudat_parse_function_args(
        NAME        tudat_add_test
        ONE_VALUE   TARGET
        MULTI_VALUE SOURCES DEPENDS LINK_LIBRARIES INCLUDES COMPILE_FLAGS
        REQUIRED    TARGET SOURCES
        ARGN        ${ARGN})

    add_executable(${TARGET} ${SOURCES})

    if(DEPENDS)
      add_dependencies(${TARGET} ${DEPENDS})
    endif()

    if(LINK_LIBRARIES)
      target_link_libraries(${TARGET} ${LINK_LIBRARIES})
    endif()

    if(COMPILE_FLAGS)
      target_compile_options(${TARGET} PRIVATE ${COMPILE_FLAGS})
    endif()

    if(INCLUDES)
      target_include_directories(${TARGET} PRIVATE ${INCLUDES})
    endif()

    # Tudat paths
    set_property(TARGET ${TARGET} PROPERTY RUNTIME_OUTPUT_DIRECTORY "${BINROOT}/unit_tests")
    get_property(CUSTOM_TEST_PROGRAM_NAME TARGET ${TARGET} PROPERTY OUTPUT_NAME)

    # Finally, add to CTest
    add_test("${TARGET}" "${BINROOT}/unit_tests/${TARGET}")

  endif()

endfunction()
