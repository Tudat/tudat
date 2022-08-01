/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *
 *    Notes:
 *      Test tolerance was set at 5.0e-15 (or 5.0e-7 for floats) instead of epsilon due to
 *      rounding errors in Eigen types with entries over a number of orders of magnitude,
 *      presumably causing the observed larger than epsilon relative differences between
 *      expected and computed values.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <vector>

#include <boost/bind/bind.hpp>
using namespace boost::placeholders;

#include <boost/make_shared.hpp>
#include <memory>
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

#include <Eigen/Core>

#include "tudat/basics/testMacros.h"

#include "tudat/astro/basic_astro/accelerationModel.h"
#include "tudat/astro/basic_astro/testAccelerationModels.h"
#include "tudat/astro/basic_astro/testBody.h"
#include "tudat/basics/basicTypedefs.h"

namespace tudat
{
namespace unit_tests
{

using basic_astrodynamics::AccelerationModel;
using basic_astrodynamics::updateAndGetAcceleration;

BOOST_AUTO_TEST_SUITE( test_accelerationModel )

//! Test whether DerivedAccelerationModel (acceleration data type=Eigen::Vector3d) functions
//! correctly.
BOOST_AUTO_TEST_CASE( test_derived3dAccelerationModel )
{
    using basic_astrodynamics::AccelerationModel3dPointer;

    // Shortcuts.
    typedef TestBody< 3, double > TestBody3d;
    typedef std::shared_ptr< TestBody3d > TestBody3dPointer;
    typedef DerivedAccelerationModel< > DerivedAccelerationModel3d;

    // Create body with initial state and time.
    TestBody3dPointer body = std::make_shared< TestBody3d >(
                ( Eigen::Vector6d( ) << 1.1, 2.2, 3.3, -0.1, 0.2, 0.3 ).finished( ), 2.0 );

    // Create acceleration model using DerivedAccelerationModel class, and pass pointers to
    // functions in body.
    AccelerationModel3dPointer accelerationModel3d
            = std::make_shared< DerivedAccelerationModel3d >(
                std::bind( &TestBody3d::getCurrentPosition, body ),
                std::bind( &TestBody3d::getCurrentTime, body ) );

    // Declare container of computed accelerations.
    std::vector< Eigen::Vector3d > computedAccelerations( 3 );

    // Get acceleration vector before members are updated.
    computedAccelerations.at( 0 ) = accelerationModel3d->getAcceleration( );

    // Update time and state.
    body->setCurrentTimeAndState(
                -1.1, ( Eigen::Vector6d( )
                        << -0.45, 10.63, -9.81, 0.11, 0.22, 0.33 ).finished( ) );

    // Update acceleration model members.
    accelerationModel3d->updateMembers( );

    // Get acceleration vector, now after members have been updated.
    computedAccelerations.at( 1 ) = accelerationModel3d->getAcceleration( );

    // Update time and state.
    body->setCurrentTimeAndState(
                4.6, ( Eigen::Vector6d( )
                       << -87.685, 101.44, -1.38, -0.12, 0.23, -0.34 ).finished( ) );

    // Update and get acceleration with single function.
    computedAccelerations.at( 2 ) = updateAndGetAcceleration( accelerationModel3d );

    // Set expected accelerations.
    std::vector< Eigen::Vector3d > expectedAccelerations;
    expectedAccelerations.push_back( Eigen::Vector3d(   1.1,     2.2,   3.3  ) / (  2.0 *  2.0 ) );
    expectedAccelerations.push_back( Eigen::Vector3d(  -0.45,   10.63, -9.81 ) / ( -1.1 * -1.1 ) );
    expectedAccelerations.push_back( Eigen::Vector3d( -87.685, 101.44, -1.38 ) / (  4.6 *  4.6 ) );

    // Check that the acceleration vectors before and after the update match expected values.
    for ( unsigned int i = 0; i < computedAccelerations.size( ); i++ )
    {
        TUDAT_CHECK_MATRIX_BASE( computedAccelerations.at( i ), expectedAccelerations.at( i ) )
                BOOST_CHECK_CLOSE_FRACTION( computedAccelerations.at( i ).coeff( row, col ),
                                   expectedAccelerations.at( i ).coeff( row, col ), 5.0e-15 );
    }
}



////! Test whether DerivedAccelerationModel (acceleration data type=Eigen::Vector1i) functions
////! correctly.
//BOOST_AUTO_TEST_CASE( test_derived1iAccelerationModel )
//{
//    // Shortcuts.
//    typedef TestBody< 1, int > TestBody1i;
//    typedef std::shared_ptr< TestBody1i > TestBody1iPointer;
//    typedef AccelerationModel< int > AccelerationModel1i;
//    typedef std::shared_ptr< AccelerationModel1i > AccelerationModel1iPointer;
//    typedef DerivedAccelerationModel< int, int, int > DerivedAccelerationModel1i;

//    // Create body with initial state and time.
//    TestBody1iPointer body = std::make_shared< TestBody1i >( Eigen::Vector2i( 20, 1 ), 2 );

//    // Create acceleration model using DerivedAccelerationModel class, and pass pointers to
//    // functions in body.
//    AccelerationModel1iPointer accelerationModeli1
//            = std::make_shared< DerivedAccelerationModel1i >(
//                std::bind( ( &TestBody1i::getCurrentPosition ), body ),
//                std::bind( ( &TestBody1i::getCurrentTime ), body ) );

//    // Declare container of computed accelerations.
//    std::vector< int > computedAccelerations( 3 );

//    // Get scalar acceleration before members are updated.
//    computedAccelerations.at( 0 ) = accelerationModeli1->getAcceleration( );

//    // Update time and state.
//    body->setCurrentTimeAndState( -3, Eigen::Vector2i( 90, 3 ) );

//    // Update acceleration model members.
//    accelerationModeli1->updateMembers( );

//    // Get scalar acceleration, now after members have been updated.
//    computedAccelerations.at( 1 ) = accelerationModeli1->getAcceleration( );

//    // Update time and state.
//    body->setCurrentTimeAndState( 4, Eigen::Vector2i( 16, -4 ) );

//    // Update and get scalar acceleration with single function.
//    computedAccelerations.at( 2 ) = updateAndGetAcceleration( accelerationModeli1 );

//    // Set expected accelerations.
//    const std::vector< int > expectedAccelerations = { 20 / ( 2 * 2 ), 90 / ( -3 * -3 ), 16 / ( 4 * 4 ) };

//    // Check that the acceleration vectors before and after the update match expected values.
//    for ( unsigned int i = 0; i < computedAccelerations.size( ); i++ )
//    {
//        BOOST_CHECK_EQUAL( computedAccelerations.at( i ), expectedAccelerations.at( i ) );
//    }
//}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
