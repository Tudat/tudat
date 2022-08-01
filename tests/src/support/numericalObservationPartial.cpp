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

#include "tudat/simulation/propagation_setup/dynamicsSimulator.h"
#include "tudat/support/numericalObservationPartial.h"

namespace tudat
{

namespace observation_partials
{


void emptyVoidFunction( ){ }


//! Function to compute numerical partial derivative of double observable w.r.t. double parameter.
Eigen::Matrix< double, 1, 1 > calculateNumericalObservationParameterPartial(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter,
        const double parameterPerturbation,
        std::function< double( const double ) > observationFunction,
        const double evaluationTime,
        std::function< void( ) > updateFunction )
{

    double unperturbedParameterValue = parameter->getParameterValue( );

    parameter->setParameterValue( unperturbedParameterValue + parameterPerturbation );
    updateFunction( );
    Eigen::Matrix< double, 1, 1 > upPerturbedValue;
    upPerturbedValue( 0, 0 ) = observationFunction( evaluationTime );

    parameter->setParameterValue( unperturbedParameterValue - parameterPerturbation );
    updateFunction( );
    Eigen::Matrix< double, 1, 1 > downPerturbedValue;
    downPerturbedValue( 0, 0 ) = observationFunction( evaluationTime );

    parameter->setParameterValue( unperturbedParameterValue );
    updateFunction( );

    return ( upPerturbedValue - downPerturbedValue ) / ( 2.0 * parameterPerturbation );
}

//! Function to compute numerical partial derivative of vector observable w.r.t. double parameter.
Eigen::Matrix< double, Eigen::Dynamic, 1 > calculateNumericalObservationParameterPartial(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter,
        const double parameterPerturbation,
        std::function< Eigen::VectorXd( const double ) > observationFunction,
        const double evaluationTime,
        std::function< void( ) > updateFunction )
{

    double unperturbedParameterValue = parameter->getParameterValue( );

    parameter->setParameterValue( unperturbedParameterValue + parameterPerturbation );
    updateFunction( );
    Eigen::Matrix< double, Eigen::Dynamic, 1 >  upPerturbedValue = observationFunction( evaluationTime );

    parameter->setParameterValue( unperturbedParameterValue - parameterPerturbation );
    updateFunction( );
    Eigen::Matrix< double, Eigen::Dynamic, 1 >  downPerturbedValue = observationFunction( evaluationTime );

    parameter->setParameterValue( unperturbedParameterValue );
    updateFunction( );

    return ( upPerturbedValue - downPerturbedValue ) / ( 2.0 * parameterPerturbation );
}

//! Function to compute numerical partial derivative of vector observable w.r.t. vector parameter.
Eigen::MatrixXd calculateNumericalObservationParameterPartial(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter,
        const Eigen::VectorXd parameterPerturbation,
        std::function< Eigen::VectorXd( const double ) > observationFunction,
        const double evaluationTime,
        std::function< void( ) > updateFunction )
{
    Eigen::MatrixXd parameterPartial = Eigen::MatrixXd::Zero( observationFunction( evaluationTime ).rows( ),
                                                              parameter->getParameterSize( ) );

    Eigen::VectorXd unperturbedParameterValue = parameter->getParameterValue( );
    Eigen::VectorXd perturbedParameterValue;

    for( int i = 0; i < parameter->getParameterSize( ); i++ )
    {
        perturbedParameterValue = unperturbedParameterValue;
        perturbedParameterValue( i ) += parameterPerturbation( i );
        parameter->setParameterValue( perturbedParameterValue );
        updateFunction( );
        Eigen::VectorXd upPerturbedValue = observationFunction( evaluationTime );

        perturbedParameterValue = unperturbedParameterValue;
        perturbedParameterValue( i ) -= parameterPerturbation( i );
        parameter->setParameterValue( perturbedParameterValue );
        updateFunction( );
        Eigen::VectorXd downPerturbedValue = observationFunction( evaluationTime );

        parameterPartial.block( 0, i, downPerturbedValue.rows( ), 1 ) = ( upPerturbedValue - downPerturbedValue ) /
                ( 2.0 * parameterPerturbation( i ) );
    }

    parameter->setParameterValue( unperturbedParameterValue );
    updateFunction( );

    return parameterPartial;
}


}

}
