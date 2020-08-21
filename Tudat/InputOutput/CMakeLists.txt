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
set(INPUTOUTPUT_SOURCES
  "${SRCROOT}${INPUTOUTPUTDIR}/basicInputOutput.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/dictionaryComparer.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/dictionaryTools.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/fieldValue.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/fixedWidthParser.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/linearFieldTransform.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/missileDatcomData.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/missileDatcomReader.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/parsedDataVectorUtilities.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/separatedParser.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/textParser.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/twoLineElementsFileReader.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/streamFilters.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/parseSolarActivityData.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/extractSolarActivityData.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/solarActivityData.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/multiDimensionalArrayReader.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/aerodynamicCoefficientReader.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/tabulatedAtmosphereReader.cpp"
)

# Add header files.
set(INPUTOUTPUT_HEADERS
  "${SRCROOT}${INPUTOUTPUTDIR}/basicInputOutput.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/dictionaryComparer.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/dictionaryEntry.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/dictionaryTools.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/extractor.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/fieldTransform.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/fieldType.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/fieldValue.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/fixedWidthParser.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/linearFieldTransform.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/missileDatcomData.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/missileDatcomReader.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/parsedDataVectorUtilities.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/parser.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/separatedParser.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/textParser.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/twoLineElementsFileReader.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/basicInputOutput.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/mapTextFileReader.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/matrixTextFileReader.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/streamFilters.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/parseSolarActivityData.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/extractSolarActivityData.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/solarActivityData.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/multiDimensionalArrayReader.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/multiDimensionalArrayWriter.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/aerodynamicCoefficientReader.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/readHistoryFromFile.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/tabulatedAtmosphereReader.h"
)

# Add unit test files.
set(INPUTOUTPUT_UNITTESTS
  "${SRCROOT}${INPUTOUTPUTDIR}/UnitTests/unitTestBasicInputOutput.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/UnitTests/unitTestMatrixTextFileReader.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/UnitTests/unitTestStreamFilters.cpp"
)

# Add static libraries.
add_library(tudat_input_output STATIC ${INPUTOUTPUT_SOURCES} ${INPUTOUTPUT_HEADERS})
setup_tudat_library_target(tudat_input_output "${SRCROOT}${INPUTOUTPUTDIR}")

# Add unit tests.
set(IO_TEST_DIR "${SRCROOT}${INPUTOUTPUTDIR}/UnitTests")

tudat_add_test(
  TARGET test_MapTextFileReader
  SOURCES "${IO_TEST_DIR}/unitTestMapTextFileReader.cpp"
  LINK_LIBRARIES tudat_input_output tudat_basic_astrodynamics ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_MatrixTextFileReader
  SOURCES "${IO_TEST_DIR}/unitTestMatrixTextFileReader.cpp"
  LINK_LIBRARIES tudat_input_output tudat_basic_astrodynamics ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_StreamFilters
  SOURCES "${IO_TEST_DIR}/unitTestStreamFilters.cpp"
  LINK_LIBRARIES tudat_input_output tudat_basic_astrodynamics ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_BasicInputOutput
  SOURCES "${IO_TEST_DIR}/unitTestBasicInputOutput.cpp"
  LINK_LIBRARIES tudat_input_output ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_ParsedDataVectorUtilities
  SOURCES "${IO_TEST_DIR}/unitTestParsedDataVectorUtilities.cpp"
  LINK_LIBRARIES tudat_input_output ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_FieldValue
  SOURCES "${IO_TEST_DIR}/unitTestFieldValue.cpp"
  LINK_LIBRARIES tudat_input_output ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_TextParser
  SOURCES "${IO_TEST_DIR}/unitTestTextParser.cpp"
  LINK_LIBRARIES tudat_input_output ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_MissileDatcomReader
  SOURCES "${IO_TEST_DIR}/unitTestMissileDatcomReader.cpp"
  LINK_LIBRARIES tudat_input_output ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_MissileDatcomData
  SOURCES "${IO_TEST_DIR}/unitTestMissileDatcomData.cpp"
  LINK_LIBRARIES tudat_input_output ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_DictionaryInputSystem
  SOURCES "${IO_TEST_DIR}/unitTestDictionaryInputSystem.cpp"
  LINK_LIBRARIES tudat_input_output ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_Extractor
  SOURCES "${IO_TEST_DIR}/unitTestExtractor.cpp"
  LINK_LIBRARIES tudat_input_output ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_FixedWidthParser
  SOURCES "${IO_TEST_DIR}/unitTestFixedWidthParser.cpp"
  LINK_LIBRARIES tudat_input_output ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_SeparatedParser
  SOURCES "${IO_TEST_DIR}/unitTestSeparatedParser.cpp"
  LINK_LIBRARIES tudat_input_output ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_LinearFieldTransform
  SOURCES "${IO_TEST_DIR}/unitTestLinearFieldTransform.cpp"
  LINK_LIBRARIES tudat_input_output ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_SolarActivityData
  SOURCES "${IO_TEST_DIR}/unitTestSolarActivityData.cpp"
  LINK_LIBRARIES tudat_input_output tudat_basic_astrodynamics ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_MultiArrayReader
  SOURCES "${IO_TEST_DIR}/unitTestMultiArrayReader.cpp"
  LINK_LIBRARIES tudat_input_output tudat_basic_astrodynamics ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_MultiArrayWriter
  SOURCES "${IO_TEST_DIR}/unitTestMultiArrayWriter.cpp"
  LINK_LIBRARIES tudat_input_output tudat_basic_astrodynamics ${Boost_LIBRARIES})

tudat_add_test(
  TARGET test_AerodynamicCoefficientReader
  SOURCES "${IO_TEST_DIR}/unitTestAerodynamicCoefficientReader.cpp"
  LINK_LIBRARIES tudat_input_output
                 tudat_basic_astrodynamics
                 tudat_basics
                 ${Boost_LIBRARIES})
