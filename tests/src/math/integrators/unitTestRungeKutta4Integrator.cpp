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
 *      Burden, R.L., Faires, J.D. Numerical Analysis, 7th Edition, Books/Cole, 2001.
 *      The MathWorks, Inc. Symbolic Math Toolbox, 2012.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <string>

#include <boost/make_shared.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "tudat/io/basicInputOutput.h"
#include "tudat/io/matrixTextFileReader.h"

#include "tudat/math/integrators/numericalIntegrator.h"
#include "tudat/math/integrators/reinitializableNumericalIntegrator.h"
#include "tudat/math/integrators/rungeKutta4Integrator.h"
#include "tudat/support/numericalIntegratorTests.h"
#include "tudat/math/integrators/numericalIntegratorTestFunctions.h"

namespace tudat
{
namespace unit_tests
{

using numerical_integrators::RungeKutta4IntegratorXd;
using numerical_integrators::NumericalIntegratorXdPointer;
using numerical_integrators::ReinitializableNumericalIntegratorXdPointer;

using numerical_integrator_test_functions::computeNonAutonomousModelStateDerivative;

BOOST_AUTO_TEST_SUITE( test_runge_kutta_4_integrator )

//! Test Runge-Kutta 4 integrator using benchmark data from (The MathWorks, 2012).
BOOST_AUTO_TEST_CASE( testRungeKutta4IntegratorUsingMatlabData )
{
    using namespace numerical_integrator_tests;

    // Read in benchmark data (generated using Symbolic Math Toolbox in Matlab
    // (The MathWorks, 2012)). This data is generated using the RK4 numerical integrator.
    const std::string pathToForwardIntegrationOutputFile = paths::getTudatTestDataPath( )
            + "/matlabOutputRungeKutta4Forwards.txt";
    const std::string pathToBackwardIntegrationOutputFile = paths::getTudatTestDataPath( )
            + "/matlabOutputRungeKutta4Backwards.txt";
    const std::string pathToDiscreteEventIntegrationOutputFile = paths::getTudatTestDataPath( )
            + "/matlabOutputRungeKutta4DiscreteEvent.txt";

    // Store benchmark data in matrix.
    const Eigen::MatrixXd matlabForwardIntegrationData =
            input_output::readMatrixFromFile( pathToForwardIntegrationOutputFile, "," );
    const Eigen::MatrixXd matlabBackwardIntegrationData =
            input_output::readMatrixFromFile( pathToBackwardIntegrationOutputFile, "," );
    const Eigen::MatrixXd matlabDiscreteEventIntegrationData =
            input_output::readMatrixFromFile( pathToDiscreteEventIntegrationOutputFile, "," );

    // Case 1: Execute integrateTo() to integrate one step forward in time.
    {
        // Declare integrator with all necessary settings.
        NumericalIntegratorXdPointer integrator
                = std::make_shared< RungeKutta4IntegratorXd >(
                    &computeNonAutonomousModelStateDerivative,
                    matlabForwardIntegrationData( FIRST_ROW, TIME_COLUMN_INDEX ),
                    ( Eigen::VectorXd( 1 )
                      << matlabForwardIntegrationData( FIRST_ROW,
                                                       STATE_COLUMN_INDEX ) ).finished( ) );

        executeOneIntegrateToStep( matlabForwardIntegrationData, 1.0e-15, integrator );
    }

    // Case 2: Execute performIntegrationStep() to perform multiple integration steps until final
    //         time.
    {
        // Declare integrator with all necessary settings.
        NumericalIntegratorXdPointer integrator
                = std::make_shared< RungeKutta4IntegratorXd >(
                    &computeNonAutonomousModelStateDerivative,
                    matlabForwardIntegrationData( FIRST_ROW, TIME_COLUMN_INDEX ),
                    ( Eigen::VectorXd( 1 )
                      << matlabForwardIntegrationData( FIRST_ROW,
                                                       STATE_COLUMN_INDEX ) ).finished( ) );

        performIntegrationStepToSpecifiedTime( matlabForwardIntegrationData,
                                               1.0e-15,  1.0e-14, integrator );
    }

    // Case 3: Execute performIntegrationStep() to perform multiple integration steps until initial
    //         time (backwards).
    {
        // Declare integrator with all necessary settings.
        NumericalIntegratorXdPointer integrator
                = std::make_shared< RungeKutta4IntegratorXd >(
                    &computeNonAutonomousModelStateDerivative,
                    matlabBackwardIntegrationData( FIRST_ROW, TIME_COLUMN_INDEX ),
                    ( Eigen::VectorXd( 1 )
                      << matlabBackwardIntegrationData( FIRST_ROW,
                                                        STATE_COLUMN_INDEX ) ).finished( ) );

        performIntegrationStepToSpecifiedTime( matlabBackwardIntegrationData,
                                               1.0e-15, 1.0e-14, integrator );
    }

    // Case 4: Execute integrateTo() to integrate to specified time in one step.
    {
        // Declare integrator with all necessary settings.
        NumericalIntegratorXdPointer integrator
                = std::make_shared< RungeKutta4IntegratorXd >(
                    &computeNonAutonomousModelStateDerivative,
                    matlabForwardIntegrationData( FIRST_ROW, TIME_COLUMN_INDEX ),
                    ( Eigen::VectorXd( 1 )
                      << matlabForwardIntegrationData( FIRST_ROW,
                                                       STATE_COLUMN_INDEX ) ).finished( ) );

        executeIntegrateToToSpecifiedTime( matlabForwardIntegrationData, 1.0e-14, integrator,
                                           matlabForwardIntegrationData(
                                               matlabForwardIntegrationData.rows( ) - 1,
                                               TIME_COLUMN_INDEX ) );
    }

    // Case 5: Execute performIntegrationstep() to integrate to specified time in multiple steps,
    //         including discrete events.
    {
        // Declare integrator with all necessary settings.
        ReinitializableNumericalIntegratorXdPointer integrator
                = std::make_shared< RungeKutta4IntegratorXd >(
                    &computeNonAutonomousModelStateDerivative,
                    matlabDiscreteEventIntegrationData( FIRST_ROW, TIME_COLUMN_INDEX ),
                    ( Eigen::VectorXd( 1 )
                      << matlabDiscreteEventIntegrationData( FIRST_ROW,
                                                             STATE_COLUMN_INDEX ) ).finished( ) );

        performIntegrationStepToSpecifiedTimeWithEvents( matlabDiscreteEventIntegrationData,
                                                         1.0e-15, 1.0e-13, integrator );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
