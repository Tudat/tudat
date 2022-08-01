/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>
#include "tudat/math/interpolators/linearInterpolator.h"

#include <Eigen/Core>

#include "tudat/basics/testMacros.h"
#include "tudat/io/matrixTextFileReader.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/math/basic/mathematicalConstants.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_linear_interpolation )

// Test implementation of linear interpolation function for data vectors.
BOOST_AUTO_TEST_CASE( test_linearInterpolation_vector )
{
    // Set vectors of data.
    const Eigen::Vector3d sortedIndependentVariables( 0.0, 1.0, 3.0 );
    const Eigen::Vector3d associatedDependentVariables( -20.0, 20.0, 21.0 );

    // Test 1: Test vector data with expected result of 0.0.
    {
        // Set target independent value in vector data.
        const double targetIndependentVariableValue = 0.5;

        // Compute interpolation.
        const double interpolatedValue = interpolators::computeLinearInterpolation(
                    sortedIndependentVariables, associatedDependentVariables,
                    targetIndependentVariableValue );

        // Verify that interpolated value corresponds to expected value (0.0).
        BOOST_CHECK_SMALL( std::fabs( interpolatedValue - 0.0 ),
                           std::numeric_limits< double >::min( ) );
    }

    // Test 2: Test vector data with expected result of 20.5.

    // Set target independent value in vector data.
    const double targetIndependentVariableValue = 2.0;

    // Compute interpolation.
    const double interpolatedValue = interpolators::computeLinearInterpolation(
                sortedIndependentVariables, associatedDependentVariables,
                targetIndependentVariableValue );

    // Verify that interpolated value corresponds to expected value (20.5).
    BOOST_CHECK_SMALL( std::fabs( interpolatedValue - 20.5 ),
                       std::numeric_limits< double >::min( ) );
}

// Test linear interpolation with map of vectors with keys as independent variable.
BOOST_AUTO_TEST_CASE( test_linearInterpolation_map )
{
    // Declare map of data and vectors for map value.
    std::map < double, Eigen::VectorXd > sortedIndepedentAndDependentVariables;
    const Eigen::Vector3d vectorOne( 10.0, -10.0, 70.0 );
    const Eigen::Vector3d vectorTwo( 20.0, -5.0, 80.0 );
    const Eigen::Vector3d vectorThree( 30.0, 60.0, 90.0 );

    // Set map values in map using vector data.
    sortedIndepedentAndDependentVariables[ 0.0 ] = vectorOne;
    sortedIndepedentAndDependentVariables[ 1.0 ] = vectorTwo;
    sortedIndepedentAndDependentVariables[ 2.0 ] = vectorThree;

    // Set target independent variable value for interpolation.
    const double targetIndependentVariableValue = 1.5;

    // Compute interpolation.
    Eigen::Vector3d interpolatedVector = interpolators::computeLinearInterpolation(
                sortedIndepedentAndDependentVariables,
                targetIndependentVariableValue );

    // Check that interpolated values correspond to expected elements [25, 27.5, 85].
    BOOST_CHECK_SMALL( std::fabs( interpolatedVector( 0 ) - 25.0 ),
                       std::numeric_limits< double >::min( ) );

    BOOST_CHECK_SMALL( std::fabs( interpolatedVector( 1 ) - 27.5 ),
                       std::numeric_limits< double >::min( ) );

    BOOST_CHECK_SMALL( std::fabs( interpolatedVector( 2 ) - 85.0 ),
                       std::numeric_limits< double >::min( ) );
}

// Test linear interpolation from benchmark values generetd by Matlab, interpolating the
// error function.
BOOST_AUTO_TEST_CASE( test_linearInterpolation_matlab_compare )
{
    using namespace interpolators;

    // Load input data used for generating matlab interpolation.
    Eigen::MatrixXd inputData = input_output::readMatrixFromFile(
                paths::getTudatTestDataPath( ) +
                "/interpolator_test_input_data.dat","," );

    // Put data in STL vectors.
    std::vector< double > independentVariableValues;
    std::vector< double > dependentVariableValues;

    for ( int i = 0; i < inputData.rows( ); i++ )
    {
        independentVariableValues.push_back( inputData( i, 0 ) );
        dependentVariableValues.push_back( inputData( i, 1 ) );
    }

    // Create linear interpolator using hunting algorithm.
    LinearInterpolatorDouble linearInterpolator(
                independentVariableValues, dependentVariableValues, huntingAlgorithm );

    // Load points at which interpolator is to be evaluated and data generated by Matlab.
    Eigen::MatrixXd benchmarkData = input_output::readMatrixFromFile(
                paths::getTudatTestDataPath( ) +
                "/linear_interpolator_test_output_data.dat",
                "," );

    // Perform interpolation for required data points.
    Eigen::VectorXd outputData = Eigen::VectorXd( benchmarkData.rows( ) );
    for ( int i = 0; i < outputData.rows( ); i++ )
    {
        outputData[ i ] = linearInterpolator.interpolate( benchmarkData( i, 0 ) );
    }

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( benchmarkData.block( 0, 1, benchmarkData.rows( ), 1 ),
                                       outputData, 1.0E-13 );

    // Create linear interpolator, now with nearest neighbvour search.
    linearInterpolator = LinearInterpolatorDouble(
                independentVariableValues, dependentVariableValues, binarySearch );

    // Perform interpolation for required data points.
    for ( int i = 0; i < outputData.rows( ); i++ )
    {
        outputData[ i ] = linearInterpolator.interpolate( benchmarkData( i, 0 ) );
    }

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( benchmarkData.block( 0, 1, benchmarkData.rows( ), 1 ),
                                       outputData, 1.0E-13 );
}

// Test linear interpolation outside of independent variable range
BOOST_AUTO_TEST_CASE( test_linearInterpolation_boundary_case )
{
    using namespace interpolators;

    // Load input data used for generating matlab interpolation.
    Eigen::MatrixXd inputData = input_output::readMatrixFromFile(
                paths::getTudatTestDataPath( ) +
                "/interpolator_test_input_data.dat","," );

    // Put data in STL vectors.
    std::vector< double > independentVariableValues;
    std::vector< double > dependentVariableValues;

    for ( int i = 0; i < inputData.rows( ); i++ )
    {
        independentVariableValues.push_back( inputData( i, 0 ) );
        dependentVariableValues.push_back( inputData( i, 1 ) );
    }

    // Create linear interpolator using hunting algorithm.
    double valueOffset = 2.0;
    double valueBelowMinimumValue = independentVariableValues[ 0 ] - valueOffset;
    double valueAboveMaximumValue = independentVariableValues[ inputData.rows( ) - 1 ] + valueOffset;
    double interpolatedValue = TUDAT_NAN, expectedValue = TUDAT_NAN;
    bool exceptionIsCaught = false;

    for( unsigned int i = 0; i < 7; i++ )
    {
        LinearInterpolatorDouble linearInterpolator(
                    independentVariableValues, dependentVariableValues, huntingAlgorithm,
                    static_cast< BoundaryInterpolationType >( i ) );

        if( static_cast< BoundaryInterpolationType >( i ) == throw_exception_at_boundary )
        {
            try
            {
                linearInterpolator.interpolate( valueBelowMinimumValue );
            }
            catch( std::runtime_error const& )
            {
                exceptionIsCaught = true;
            }
            BOOST_CHECK_EQUAL( exceptionIsCaught, true );

            exceptionIsCaught = false;
            try
            {
                linearInterpolator.interpolate( valueAboveMaximumValue );
            }
            catch( std::runtime_error const& )
            {
                exceptionIsCaught = true;
            }
            BOOST_CHECK_EQUAL( exceptionIsCaught, true );
        }
        else if( ( static_cast< BoundaryInterpolationType >( i ) == use_boundary_value ) ||
                 ( static_cast< BoundaryInterpolationType >( i ) == use_boundary_value_with_warning ) )
        {
            interpolatedValue = linearInterpolator.interpolate( valueBelowMinimumValue );
            BOOST_CHECK_CLOSE_FRACTION( interpolatedValue, dependentVariableValues.at( 0 ), 1.0E-15 );

            interpolatedValue = linearInterpolator.interpolate( valueAboveMaximumValue );
            BOOST_CHECK_CLOSE_FRACTION( interpolatedValue, dependentVariableValues.at( inputData.rows( ) - 1 ), 1.0E-15 );

        }
        else if( ( static_cast< BoundaryInterpolationType >( i ) == extrapolate_at_boundary ) ||
                 ( static_cast< BoundaryInterpolationType >( i ) == extrapolate_at_boundary_with_warning ) )
        {
            interpolatedValue = linearInterpolator.interpolate( valueBelowMinimumValue );
            expectedValue = dependentVariableValues.at( 0 ) - valueOffset * (
                        dependentVariableValues.at( 1 ) - dependentVariableValues.at( 0 ) ) / (
                        independentVariableValues.at( 1 ) - independentVariableValues.at( 0 ) );
            BOOST_CHECK_CLOSE_FRACTION( interpolatedValue, expectedValue, 1.0E-15 );

            interpolatedValue = linearInterpolator.interpolate( valueAboveMaximumValue );
            expectedValue = dependentVariableValues.at( inputData.rows( ) - 1 ) + valueOffset * (
                        dependentVariableValues.at( inputData.rows( ) - 1 ) - dependentVariableValues.at( inputData.rows( ) - 2 ) ) / (
                        independentVariableValues.at( inputData.rows( ) - 1 ) - independentVariableValues.at( inputData.rows( ) - 2 ) );
            BOOST_CHECK_CLOSE_FRACTION( interpolatedValue, expectedValue, 1.0E-15 );
        }
        else if( ( static_cast< BoundaryInterpolationType >( i ) == use_default_value ) ||
                 ( static_cast< BoundaryInterpolationType >( i ) == use_default_value_with_warning ) )
        {
            interpolatedValue = linearInterpolator.interpolate( valueBelowMinimumValue );
            BOOST_CHECK_CLOSE_FRACTION( interpolatedValue, 0.0, 1.0E-15 );

            interpolatedValue = linearInterpolator.interpolate( valueAboveMaximumValue );
            BOOST_CHECK_CLOSE_FRACTION( interpolatedValue, 0.0, 1.0E-15 );
        }
    }
}

// Test linear interpolation outside of independent variable range with default extrapolation value
BOOST_AUTO_TEST_CASE( test_linearInterpolation_boundary_case_extrapolation_default_value )
{
    using namespace interpolators;

    // Load input data used for generating matlab interpolation.
    Eigen::MatrixXd inputData = input_output::readMatrixFromFile(
                paths::getTudatTestDataPath( ) +
                "/interpolator_test_input_data.dat","," );

    // Put data in STL vectors.
    std::vector< double > independentVariableValues;

    for ( int i = 0; i < inputData.rows( ); i++ )
    {
        independentVariableValues.push_back( inputData( i, 0 ) );
    }

    // Create linear interpolator using hunting algorithm.
    double valueOffset = 2.0;
    double valueBelowMinimumValue = independentVariableValues[ 0 ] - valueOffset;
    double valueAboveMaximumValue = independentVariableValues[ inputData.rows( ) - 1 ] + valueOffset;

    // Test with long double
    {
        // Put data in STL vectors.
        std::vector< long double > dependentVariableValues;
        for ( int i = 0; i < inputData.rows( ); i++ )
        {
            dependentVariableValues.push_back( inputData( i, 1 ) );
        }
        long double interpolatedValue;

        for( unsigned int i = 5; i < 7; i++ )
        {
            LinearInterpolator< double, long double > linearInterpolator(
                        independentVariableValues, dependentVariableValues, huntingAlgorithm,
                        static_cast< BoundaryInterpolationType >( i ) );

            interpolatedValue = linearInterpolator.interpolate( valueBelowMinimumValue );
            BOOST_CHECK_CLOSE_FRACTION( interpolatedValue, 0.0L, 1.0E-15 );

            interpolatedValue = linearInterpolator.interpolate( valueAboveMaximumValue );
            BOOST_CHECK_CLOSE_FRACTION( interpolatedValue, 0.0L, 1.0E-15 );
        }
    }

    // Test with Eigen::Vector3d
    {
        // Put data in STL vectors.
        std::vector< Eigen::Vector3d > dependentVariableValues;
        Eigen::Vector3d tempData = Eigen::Vector3d::Zero( );
        for ( int i = 0; i < inputData.rows( ); i++ )
        {
            tempData[ 0 ] = inputData( i, 1 );
            dependentVariableValues.push_back( tempData );
        }
        Eigen::Vector3d interpolatedValue;

        for( unsigned int i = 5; i < 7; i++ )
        {
            LinearInterpolator< double, Eigen::Vector3d > linearInterpolator(
                        independentVariableValues, dependentVariableValues, huntingAlgorithm,
                        static_cast< BoundaryInterpolationType >( i ) );

            interpolatedValue = linearInterpolator.interpolate( valueBelowMinimumValue );
            BOOST_CHECK_SMALL( ( interpolatedValue - Eigen::Vector3d::Zero( ) ).norm( ), 1.0E-15 );

            interpolatedValue = linearInterpolator.interpolate( valueAboveMaximumValue );
            BOOST_CHECK_SMALL( ( interpolatedValue - Eigen::Vector3d::Zero( ) ).norm( ), 1.0E-15 );
        }
    }

    // Test with Eigen::Matrix3d
    {
        // Put data in STL vectors.
        std::vector< Eigen::Matrix3d > dependentVariableValues;
        Eigen::Matrix3d tempData = Eigen::Matrix3d::Zero( );
        for ( int i = 0; i < inputData.rows( ); i++ )
        {
            tempData( 0, 2 ) = inputData( i, 1 );
            dependentVariableValues.push_back( tempData );
        }
        Eigen::Matrix3d interpolatedValue;

        for( unsigned int i = 5; i < 7; i++ )
        {
            LinearInterpolator< double, Eigen::Matrix3d > linearInterpolator(
                        independentVariableValues, dependentVariableValues, huntingAlgorithm,
                        static_cast< BoundaryInterpolationType >( i ) );

            interpolatedValue = linearInterpolator.interpolate( valueBelowMinimumValue );
            BOOST_CHECK_SMALL( ( interpolatedValue - Eigen::Matrix3d::Zero( ) ).norm( ), 1.0E-15 );

            interpolatedValue = linearInterpolator.interpolate( valueAboveMaximumValue );
            BOOST_CHECK_SMALL( ( interpolatedValue - Eigen::Matrix3d::Zero( ) ).norm( ), 1.0E-15 );
        }
    }
}

// Test linear interpolation outside of independent variable range with extrapolation value given by user
BOOST_AUTO_TEST_CASE( test_linearInterpolation_boundary_case_extrapolation_user_value )
{
    using namespace interpolators;

    // Load input data used for generating matlab interpolation.
    Eigen::MatrixXd inputData = input_output::readMatrixFromFile(
                paths::getTudatTestDataPath( ) +
                "/interpolator_test_input_data.dat","," );

    // Put data in STL vectors.
    std::vector< double > independentVariableValues;

    for ( int i = 0; i < inputData.rows( ); i++ )
    {
        independentVariableValues.push_back( inputData( i, 0 ) );
    }

    // Create linear interpolator using hunting algorithm.
    double valueOffset = 2.0;
    double valueBelowMinimumValue = independentVariableValues[ 0 ] - valueOffset;
    double valueAboveMaximumValue = independentVariableValues[ inputData.rows( ) - 1 ] + valueOffset;

    // Test with long double
    {
        // Put data in STL vectors.
        std::vector< long double > dependentVariableValues;

        for ( int i = 0; i < inputData.rows( ); i++ )
        {
            dependentVariableValues.push_back( inputData( i, 1 ) );
        }
        long double interpolatedValue;
        long double extrapolationValue = 1.0L;

        for( unsigned int i = 5; i < 7; i++ )
        {
            LinearInterpolator< double, long double > linearInterpolator(
                        independentVariableValues, dependentVariableValues, huntingAlgorithm,
                        static_cast< BoundaryInterpolationType >( i ), extrapolationValue );

            interpolatedValue = linearInterpolator.interpolate( valueBelowMinimumValue );
            BOOST_CHECK_CLOSE_FRACTION( interpolatedValue, extrapolationValue, 1.0E-15 );

            interpolatedValue = linearInterpolator.interpolate( valueAboveMaximumValue );
            BOOST_CHECK_CLOSE_FRACTION( interpolatedValue, extrapolationValue, 1.0E-15 );
        }
    }

    // Test with Eigen::Vector3d
    {
        // Put data in STL vectors.
        std::vector< Eigen::Vector3d > dependentVariableValues;
        Eigen::Vector3d tempData = Eigen::Vector3d::Zero( );
        for ( int i = 0; i < inputData.rows( ); i++ )
        {
            tempData[ 0 ] = inputData( i, 1 );
            dependentVariableValues.push_back( tempData );
        }
        Eigen::Vector3d interpolatedValue;
        Eigen::Vector3d extrapolationValue;
        extrapolationValue[ 0 ] = 1.5;
        extrapolationValue[ 1 ] = -3;
        extrapolationValue[ 2 ] = 0;

        for( unsigned int i = 5; i < 7; i++ )
        {
            LinearInterpolator< double, Eigen::Vector3d > linearInterpolator(
                        independentVariableValues, dependentVariableValues, huntingAlgorithm,
                        static_cast< BoundaryInterpolationType >( i ), extrapolationValue );

            interpolatedValue = linearInterpolator.interpolate( valueBelowMinimumValue );
            BOOST_CHECK_SMALL( ( interpolatedValue - extrapolationValue ).norm( ), 1.0E-15 );

            interpolatedValue = linearInterpolator.interpolate( valueAboveMaximumValue );
            BOOST_CHECK_SMALL( ( interpolatedValue - extrapolationValue ).norm( ), 1.0E-15 );
        }
    }

    // Test with Eigen::Matrix3d
    {
        // Put data in STL vectors.
        std::vector< Eigen::Matrix3d > dependentVariableValues;
        Eigen::Matrix3d tempData = Eigen::Matrix3d::Zero( );
        for ( int i = 0; i < inputData.rows( ); i++ )
        {
            tempData( 0, 2 ) = inputData( i, 1 );
            dependentVariableValues.push_back( tempData );
        }
        Eigen::Matrix3d interpolatedValue;
        Eigen::Matrix3d extrapolationValue = Eigen::Matrix3d::Random( );

        for( unsigned int i = 5; i < 7; i++ )
        {
            LinearInterpolator< double, Eigen::Matrix3d > linearInterpolator(
                        independentVariableValues, dependentVariableValues, huntingAlgorithm,
                        static_cast< BoundaryInterpolationType >( i ), extrapolationValue );

            interpolatedValue = linearInterpolator.interpolate( valueBelowMinimumValue );
            BOOST_CHECK_SMALL( ( interpolatedValue - extrapolationValue ).norm( ), 1.0E-15 );

            interpolatedValue = linearInterpolator.interpolate( valueAboveMaximumValue );
            BOOST_CHECK_SMALL( ( interpolatedValue - extrapolationValue ).norm( ), 1.0E-15 );
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
