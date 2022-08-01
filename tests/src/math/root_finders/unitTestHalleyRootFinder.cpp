/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Notes
 *      This unit test does not include testing with the function
 *      testFunctionWithLargeRootDifferences, because this test function does not provide a second
 *      derivative.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/bind/bind.hpp>
using namespace boost::placeholders;

#include <boost/make_shared.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "tudat/math/root_finders/halleyRootFinder.h"
#include "tudat/math/root_finders/testFunction1.h"
#include "tudat/math/root_finders/testFunction2.h"
#include "tudat/math/root_finders/testFunction3.h"
#include "tudat/math/root_finders/testFunctionWithZeroRoot.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( testsuite_rootfinders )

using namespace tudat;
using namespace root_finders;


//! Check if Halley method converges on test function #1 (TestFunction1).
BOOST_AUTO_TEST_CASE( test_halleyRootFinder_testFunction1 )
{
    // Create object containing the test functions.
    std::shared_ptr< TestFunction1 > testFunction = std::make_shared< TestFunction1 >( 2 );

    // The termination condition.
    HalleyRootFinder< double >::TerminationFunction terminationConditionFunction =
            std::bind( &RootAbsoluteToleranceTerminationCondition< double >::checkTerminationCondition,
                         std::make_shared< RootAbsoluteToleranceTerminationCondition< double > >(
                             testFunction->getTrueRootAccuracy( ) ), std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5 );

    // Test Halley object.
    HalleyRootFinder< double > halleyRootFinder( terminationConditionFunction );

    // Let the Halley method search for the root.
    const double root = halleyRootFinder.execute( testFunction, testFunction->getInitialGuess( ) );

    // Check if the result is within the requested accuracy.
    BOOST_CHECK_CLOSE_FRACTION( root, testFunction->getTrueRootLocation( ), 1.0e-15 );
    BOOST_CHECK_LT( testFunction->evaluate( root ), testFunction->getTrueRootAccuracy( ) );
}

//! Check if Halley method converges on test function #2 (TestFunction2).
BOOST_AUTO_TEST_CASE( test_halleyRootFinder_testFunction2 )
{
    // Create object containing the test functions.
    std::shared_ptr< TestFunction2 > testFunction = std::make_shared< TestFunction2 >( 2 );

    // The termination condition.
    HalleyRootFinder< >::TerminationFunction terminationConditionFunction =
            std::bind( &RootAbsoluteToleranceTerminationCondition< double >::checkTerminationCondition,
                         std::make_shared< RootAbsoluteToleranceTerminationCondition< double > >(
                             testFunction->getTrueRootAccuracy( ) ), std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5 );

    // Test Halley object.
    HalleyRootFinder< > halleyRootFinder( terminationConditionFunction );

    // Let the Halley method search for the root.
    const double root = halleyRootFinder.execute( testFunction, testFunction->getInitialGuess( ) );

    // Check if the result is within the requested accuracy.
    BOOST_CHECK_CLOSE_FRACTION( root, testFunction->getTrueRootLocation( ), 1.0e-15 );
    BOOST_CHECK_LT( testFunction->evaluate( root ), testFunction->getTrueRootAccuracy( ) );
}

//! Check if Halley method converges on test function #3 (TestFunction3).
BOOST_AUTO_TEST_CASE( test_halleyRootFinder_testFunction3 )
{
    // Create object containing the test functions.
    std::shared_ptr< TestFunction3 > testFunction = std::make_shared< TestFunction3 >( 2 );

    // The termination condition.
    HalleyRootFinder< >::TerminationFunction terminationConditionFunction =
            std::bind( &RootAbsoluteToleranceTerminationCondition< double >::checkTerminationCondition,
                         std::make_shared< RootAbsoluteToleranceTerminationCondition< double > >(
                             testFunction->getTrueRootAccuracy( ) ), std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5 );

    // Test Halley object.
    HalleyRootFinder< > halleyRootFinder( terminationConditionFunction );

    // Let the Halley method search for the root.
    const double root = halleyRootFinder.execute( testFunction, testFunction->getInitialGuess( ) );

    // Check if the result is within the requested accuracy.
    BOOST_CHECK_CLOSE_FRACTION( root, testFunction->getTrueRootLocation( ), 1.0e-15 );
    BOOST_CHECK_LT( testFunction->evaluate( root ), testFunction->getTrueRootAccuracy( ) );
}

//! Check if Halley method converges on function with zero root (testFunctionWithZeroRoot).
// Not the best test case. Inheritance from old code. Not really relevant anymore. The basic
// idea is that Halley's method should work for both a function that becomes zero, as well as for a
// function that does not become zero. A better case should be written.
BOOST_AUTO_TEST_CASE( test_halleyRootFinder_testFunctionWithZeroRoot )
{
    // Create object containing the test functions.
    std::shared_ptr< TestFunctionWithZeroRoot > testFunction =
            std::make_shared< TestFunctionWithZeroRoot >( 2 );

    // The termination condition.
    HalleyRootFinder< >::TerminationFunction terminationConditionFunction
            = std::bind( &RootAbsoluteToleranceTerminationCondition< double >::
                           checkTerminationCondition,
                           std::make_shared< RootAbsoluteToleranceTerminationCondition< double > >(
                               1.0e-150 ), std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5 );

    // Test Halley object.
    HalleyRootFinder< > halleyRootFinder( terminationConditionFunction );

    // Let the Halley method search for the root.
    const double root = halleyRootFinder.execute( testFunction, testFunction->getInitialGuess( ) );

    // Check if the result is within the requested accuracy.
    BOOST_CHECK_SMALL( root, 1.0e-100 );
    BOOST_CHECK_SMALL( testFunction->evaluate( root ), 1.0e-200 );
}

BOOST_AUTO_TEST_SUITE_END( ) // testsuite_rootfinders

} // namespace unit_tests
} // namespace tudat
