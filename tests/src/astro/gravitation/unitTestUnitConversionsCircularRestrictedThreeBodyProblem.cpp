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
 *      Koon, W.S., et al. Dynamical Systems, the Three-Body Problem and Space Mission Design,
 *          2006.
 *      Jet propulsion Laboratory, NASA. Planets and Pluto: Physical Characteristics.
 *          http://ssd.jpl.nasa.gov/?planet_phys_par, last updated: 5th November, 2008, last
 *          accessed: 22nd November, 2011.
 *      Jet propulsion Laboratory, NASA. astro Constants.
 *          http://ssd.jpl.nasa.gov/?constants, last updated: 6th September, 2011, last accessed:
 *          22nd November, 2011.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <cmath>
#include <limits>

#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/math/basic/mathematicalConstants.h"

#include "tudat/astro/gravitation/unitConversionsCircularRestrictedThreeBodyProblem.h"

namespace tudat
{
namespace unit_tests
{

//! Test if Cartesian state conversion for CRTBP is computed correctly.
BOOST_AUTO_TEST_CASE( testCartesianStateConversionCircularRestrictedThreeBodyProblem )
{
    // Test 1: conversion of dimensionless to dimensional Cartesian state vector for test particle
    // in Sun-Jupiter system [m, m/s] (Koon, 2006).
    {
        // Set distance between Sun and Jupiter [m] [NASA, 2010].
        double distanceSunJupiter = 7.784e11;

        // Set magnitude of velocity of Sun-Jupiter system [m/s].
        double velocitySunJupiter = 13.102e3;

        // Set gravitational parameter of Jupiter [m^3 s^-2].
        double gravitationalParameterJupiter = 6.67259e-11 * 1898.13e24;

        // Set gravitational parameter of Sun [m^3 s^-2].
        double gravitationalParameterSun = 1.32712440018e20;

        // Set mass parameter of Sun-Jupiter system.
        double massParameter = gravitationalParameterJupiter
                / ( gravitationalParameterJupiter + gravitationalParameterSun );

        // Set dimensionless state of the Sun.
        using orbital_element_conversions::xCartesianPositionIndex;
        Eigen::VectorXd dimensionlessStateOfSun = Eigen::VectorXd::Zero( 6 );
        dimensionlessStateOfSun( xCartesianPositionIndex ) = -massParameter;

        // Declare and set dimensionless state of Jupiter.
        Eigen::VectorXd dimensionlessStateOfJupiter = Eigen::VectorXd::Zero( 6 );
        dimensionlessStateOfJupiter( xCartesianPositionIndex ) = 1.0 - massParameter;

        // Set dimensionless state of test particle in Sun-Jupiter system.
        using orbital_element_conversions::yCartesianVelocityIndex;
        Eigen::VectorXd dimensionlessStateOfTestParticle = Eigen::VectorXd::Zero( 6 );
        dimensionlessStateOfTestParticle( xCartesianPositionIndex ) = 1.0;
        dimensionlessStateOfTestParticle( yCartesianVelocityIndex ) = 1.0;

        // Convert dimensionless to dimensional Cartesian state vector for test particle [m, m/s].
        Eigen::VectorXd computedDimensionalStateOfTestParticle( 6 );
        computedDimensionalStateOfTestParticle =
                circular_restricted_three_body_problem::convertDimensionlessCartesianStateToDimensionalUnits(
                    dimensionlessStateOfTestParticle, gravitationalParameterSun,
                    gravitationalParameterJupiter, distanceSunJupiter );

        // Check if computed position and velocity matches expected value.
        BOOST_CHECK_CLOSE_FRACTION( distanceSunJupiter,
                                    computedDimensionalStateOfTestParticle( xCartesianPositionIndex ),
                                    std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_CLOSE_FRACTION( velocitySunJupiter,
                                    computedDimensionalStateOfTestParticle( yCartesianVelocityIndex ),
                                    1.0e-2 );
    }
}

//! Test if time conversion for CRTBP is computed correctly.
BOOST_AUTO_TEST_CASE( testTimeConversionCircularRestrictedThreeBodyProblem )
{
    using mathematical_constants::PI;

    // Test 1: conversion of dimensionless to dimensional time for test particle in Sun-Jupiter
    // system [s] (Koon, 2006).
    {
        // Set distance between Sun and Jupiter [m] (NASA, 2010).
        double distanceSunJupiter = 7.784e11;

        // Set gravitational parameter of Jupiter [m^3 s^-2].
        double gravitationalParameterJupiter = 6.67259e-11 * 1898.13e24;

        // Set gravitational parameter of Sun [m^3 s^-2].
        double gravitationalParameterSun = 1.32712440018e20;

        // Compute dimensional time for one complete orbit.
        double computedDimensionalTime =
                circular_restricted_three_body_problem::convertDimensionlessTimeToDimensionalTime(
                    2.0 * PI, gravitationalParameterSun, gravitationalParameterJupiter,
                    distanceSunJupiter );

        // Set expected orbital period of system [s].
        double expectedOrbitalPeriodOfSunJupiterSystem = 3.733e08;

        // Check if computed orbital period matches expected value.
        BOOST_CHECK_CLOSE_FRACTION( expectedOrbitalPeriodOfSunJupiterSystem,
                                    computedDimensionalTime, 1.0e-2 );
    }
}

} // namespace unit_tests
} // namespace tudat
