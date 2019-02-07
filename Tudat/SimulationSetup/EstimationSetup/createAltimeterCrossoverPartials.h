/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEALTIMETERCROSSOVERPARTIALS_H
#define TUDAT_CREATEALTIMETERCROSSOVERPARTIALS_H

#include <vector>
#include <map>

#include <memory>

#include <Eigen/Core>

#include "Tudat/Mathematics/Interpolators/interpolator.h"

#include "Tudat/SimulationSetup/EstimationSetup/createCartesianStatePartials.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/altimeterCrossoverPartial.h"
#include "Tudat/SimulationSetup/EstimationSetup/createLightTimeCorrectionPartials.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/initialTranslationalState.h"
#include "Tudat/Astrodynamics/ObservationModels/linkTypeDefs.h"


namespace tudat
{

namespace observation_partials
{

//! Function to generate angular position partial wrt a position of a body.
/*!
 *  Function to generate angular position partial wrt a position of a body, for a single link ends (which must contain a
 *  transmitter and receiever  linkEndType).
 *  \param altimeterCrossoverLinkEnds Link ends (transmitter and receiever) for which partials are to be calculated
 *  (i.e. for which angular position observation observations are to be processed).
 *  \param bodyMap List of all bodies, for creating angular position partial.
 *  \param bodyToEstimate Name of body wrt position of which a partial is to be created.
 *  \param altimeterCrossoverScaler Object scale position partials to angular position partials for current link ends.
 *  \param lightTimeCorrectionPartialObjects List of light time correction partials to be used (empty by default)
 *  \return Angular position partial object wrt a current position of a body (is nullptr if no parameter dependency exists).
 */
std::shared_ptr< AltimeterCrossoverPartial > createAltimeterCrossoverPartialWrtBodyPosition(
        const observation_models::LinkEnds altimeterCrossoverLinkEnds,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::string bodyToEstimate,
        const std::shared_ptr< AltimeterCrossoverScaling > altimeterCrossoverScaler,
        const std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > >&
        lightTimeCorrectionPartialObjects =
        std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > >( ) );

//! Function to generate angular position partial wrt a single  parameter.
/*!
 *  Function to generate angular position partial wrt a single  parameter, for a single link ends (which must contain a
 *  transmitter and receiever linkEndType).
 *  \tparam ParameterType Type of parameter (double for size 1, VectorXd for larger size).
 *  \param altimeterCrossoverLinkEnds Link ends (transmitter and receiever) for which angular position partials are to be
 *  calculated (i.e. for which  angular position observations are to be processed).
 *  \param bodyMap List of all bodies, for creating angular position partial.
 *  \param parameterToEstimate Object of current parameter that is to be estimated.
 *  \param altimeterCrossoverScaler Object scale position partials to angular position partials for current link ends.
 *  \param lightTimeCorrectionPartialObjects List of light time correction partials to be used (empty by default)
 *  \return Angular position partial object wrt a single parameter (is nullptr if no parameter dependency exists).
 */
template< typename ParameterType >
std::shared_ptr< AltimeterCrossoverPartial > createAltimeterCrossoverPartialWrtParameter(
        const observation_models::LinkEnds altimeterCrossoverLinkEnds,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< estimatable_parameters::EstimatableParameter< ParameterType > > parameterToEstimate,
        const std::shared_ptr< AltimeterCrossoverScaling > altimeterCrossoverScaler,
        const std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > >&
        lightTimeCorrectionPartialObjects =
        std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > >( ) )
{
    std::shared_ptr< AltimeterCrossoverPartial > altimeterCrossoverPartial;

//    {
//        std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > > positionPartials =
//                createCartesianStatePartialsWrtParameter( altimeterCrossoverLinkEnds, bodyMap, parameterToEstimate );

//        std::shared_ptr< AltimeterCrossoverPartial > testAltimeterCrossoverPartial = std::make_shared< AltimeterCrossoverPartial >(
//                    altimeterCrossoverScaler, positionPartials, parameterToEstimate->getParameterName( ),
//                    lightTimeCorrectionPartialObjects );

//        // Create angular position partials if any position partials are created (i.e. if any dependency exists).
//        if( positionPartials.size( ) > 0 || testAltimeterCrossoverPartial->getNumberOfLighTimeCorrectionPartialsFunctions( ) )
//        {
//            altimeterCrossoverPartial = testAltimeterCrossoverPartial;
//        }
//    }
    return altimeterCrossoverPartial;
}

//! Function to generate angular position partials and associated scaler for single link end.
/*!
 *  Function to generate angular position partials and associated scaler for all parameters that are to be estimated,
 *  for a single link ends.
 *  The set of parameters and bodies that are to be estimated, as well as the set of link ends
 *  (each of which must contain a transmitter and receiever linkEndType) that are to be used.
 *  \param altimeterCrossoverLinkEnds Link ends (transmitter and receiever) for which angular position partials are to be
 *  calculated (i.e. for which angular position observations are to be processed).
 *  \param bodyMap List of all bodies, for creating angular position partials.
 *  \param parametersToEstimate Set of parameters that are to be estimated (in addition to initial states of
 *  requested bodies)
 *  \param lightTimeCorrections List of light time correction partials to be used (empty by default)
 *  \return Set of observation partials with associated indices in complete vector of parameters that are estimated,
 *  representing all  necessary angular position partials of a single link end, and AltimeterCrossoverScaling, object, used for
 *  scaling the position partial members of all AltimeterCrossoverPartials in link end.
 */
template< typename ParameterType >
std::pair< std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< 1 > > >,
std::shared_ptr< PositionPartialScaling > >
createAltimeterCrossoverPartials(
        const observation_models::LinkEnds altimeterCrossoverLinkEnds,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
        const std::vector< std::shared_ptr< observation_models::LightTimeCorrection > >& lightTimeCorrections =
        std::vector< std::shared_ptr< observation_models::LightTimeCorrection > >( ) )
{
    std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > > lightTimeCorrectionPartialObjects;
    if( lightTimeCorrections.size( ) > 0 )
    {
        lightTimeCorrectionPartialObjects = observation_partials::createLightTimeCorrectionPartials( lightTimeCorrections );
    }

    // Create scaling object, to be used for all angular position partials in current link end.
    std::shared_ptr< AltimeterCrossoverScaling > altimeterCrossoverScaling =
            std::make_shared< AltimeterCrossoverScaling >( );

    SingleLinkObservationPartialList altimeterCrossoverPartials;


    // Initialize vector index variables.
    int currentIndex = 0;
    std::pair< int, int > currentPair = std::pair< int, int >( currentIndex, 1 );

    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter<
            Eigen::Matrix< ParameterType, Eigen::Dynamic, 1 > > > > initialDynamicalParameters =
            parametersToEstimate->getEstimatedInitialStateParameters( );

    // Iterate over list of bodies of which the partials of the accelerations acting on them are required.
    for( unsigned int i = 0; i < initialDynamicalParameters.size( ); i++ )
    {
        std::string acceleratedBody;
        if( initialDynamicalParameters.at( i )->getParameterName( ).first == estimatable_parameters::initial_body_state )
        {
            acceleratedBody = initialDynamicalParameters.at( i )->getParameterName( ).second.first;
        }
        else if( initialDynamicalParameters.at( i )->getParameterName( ).first == estimatable_parameters::arc_wise_initial_body_state )
        {
            acceleratedBody = initialDynamicalParameters.at( i )->getParameterName( ).second.first;
        }
        else
        {
            throw std::runtime_error( "Error when making angular position partials, could not identify parameter" );
        }

        // Create position angular position partial for current body
        std::shared_ptr< AltimeterCrossoverPartial > currentAltimeterCrossoverPartial =
                createAltimeterCrossoverPartialWrtBodyPosition(
                    altimeterCrossoverLinkEnds, bodyMap, acceleratedBody, altimeterCrossoverScaling,
                    lightTimeCorrectionPartialObjects );

        // Check if partial is non-nullptr (i.e. whether dependency exists between current angular position and current body)
        if( currentAltimeterCrossoverPartial != nullptr )
        {
            // Add partial to the list.
            currentPair = std::pair< int, int >( currentIndex, 6 );
            altimeterCrossoverPartials[ currentPair ] = currentAltimeterCrossoverPartial;
        }

        // Increment current index by size of body initial state (6).
        currentIndex += 6;
    }

    // Iterate over all double parameters that are to be estimated.
    std::map< int, std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > > doubleParametersToEstimate =
            parametersToEstimate->getDoubleParameters( );
    for( std::map< int, std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > >::iterator
         parameterIterator = doubleParametersToEstimate.begin( );
         parameterIterator != doubleParametersToEstimate.end( ); parameterIterator++ )
    {
        // Create position angular position partial for current parameter
        std::shared_ptr< AltimeterCrossoverPartial > currentAltimeterCrossoverPartial =
                createAltimeterCrossoverPartialWrtParameter(
                    altimeterCrossoverLinkEnds, bodyMap, parameterIterator->second, altimeterCrossoverScaling,
                    lightTimeCorrectionPartialObjects );

        if( currentAltimeterCrossoverPartial != nullptr )
        {
            // Add partial to the list.
            currentPair = std::pair< int, int >( parameterIterator->first, 1 );
            altimeterCrossoverPartials[ currentPair ] = currentAltimeterCrossoverPartial;
        }
    }

    // Iterate over all vector parameters that are to be estimated.
    std::map< int, std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > >
            vectorParametersToEstimate = parametersToEstimate->getVectorParameters( );
    for( std::map< int, std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd  > > >::iterator
         parameterIterator = vectorParametersToEstimate.begin( );
         parameterIterator != vectorParametersToEstimate.end( ); parameterIterator++ )
    {

        // Create position angular position partial for current parameter
        std::shared_ptr< ObservationPartial< 1 > > currentAltimeterCrossoverPartial;

        if( !isParameterObservationLinkProperty( parameterIterator->second->getParameterName( ).first )  )
        {
            currentAltimeterCrossoverPartial = createAltimeterCrossoverPartialWrtParameter(
                        altimeterCrossoverLinkEnds, bodyMap, parameterIterator->second, altimeterCrossoverScaling,
                        lightTimeCorrectionPartialObjects );
        }
        else
        {
            currentAltimeterCrossoverPartial = createObservationPartialWrtLinkProperty< 1 >(
                        altimeterCrossoverLinkEnds, observation_models::angular_position, parameterIterator->second );
        }

        // Check if partial is non-nullptr (i.e. whether dependency exists between current observable and current parameter)
        if( currentAltimeterCrossoverPartial != nullptr )
        {
            // Add partial to the list.
            currentPair = std::pair< int, int >( parameterIterator->first,
                                                 parameterIterator->second->getParameterSize( ) );
            altimeterCrossoverPartials[ currentPair ] = currentAltimeterCrossoverPartial;
        }

    }
    return std::make_pair( altimeterCrossoverPartials, altimeterCrossoverScaling );
}

//! Function to generate angular position partials for all parameters that are to be estimated, for all sets of link ends.
/*!
 *  Function to generate angular position partials for all parameters that are to be estimated, for all sets of link ends.
 *  The angular position partials are generated per set of link ends. The set of parameters and bodies that are to be
 *  estimated, as well as the set of link ends (each of which must contain a transmitter and receiever linkEndType)
 *  that are to be used.
 *  \param linkEnds Vector of all link ends for which angular position partials are to be calculated (i.e. for which one-way
 *  angular position observations are  to be processed).
 *  \param bodyMap List of all bodies, for creating angular position partials.
 *  \param parametersToEstimate Set of parameters that are to be estimated (in addition to initial states
 *  of requested bodies)
 *  \param lightTimeCorrections List of light time correction partials to be used (empty by default)
 *  \return Map of SingleLinkObservationPartialList, representing all necessary angular position partials of a single link
 *  end, and AltimeterCrossoverScaling, object, used for scaling the position partial members of all AltimeterCrossoverPartials in
 *  link end.
 */
template< typename ParameterType >
std::map< observation_models::LinkEnds, std::pair< SingleLinkObservationPartialList,
std::shared_ptr< PositionPartialScaling > > >
createAltimeterCrossoverPartials(
        const std::vector< observation_models::LinkEnds > linkEnds,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
        const std::map< observation_models::LinkEnds,
        std::vector< std::vector< std::shared_ptr< observation_models::LightTimeCorrection > > > >& lightTimeCorrections =
        std::map< observation_models::LinkEnds,
        std::vector< std::vector< std::shared_ptr< observation_models::LightTimeCorrection > > > >( ) )
{
    // Declare return list.
    std::map< observation_models::LinkEnds, std::pair< SingleLinkObservationPartialList,
            std::shared_ptr< PositionPartialScaling > > > altimeterCrossoverPartials;

    // Iterate over all link ends.
    for( unsigned int i = 0; i < linkEnds.size( ); i++ )
    {
        // Check if required link end types are present
        if( ( linkEnds[ i ].count( observation_models::first_arc_body ) == 0 ) ||
                ( linkEnds[ i ].count( observation_models::second_arc_body ) == 0 ) )
        {
            throw std::runtime_error( "Error when making angular position partials, did not find both receiver and transmitter in link ends" );

        }

        std::vector< std::shared_ptr< observation_models::LightTimeCorrection > > singleLinkLightTimeCorrections;
        if( lightTimeCorrections.count( linkEnds.at( i ) ) > 0 )
        {
            if( lightTimeCorrections.at( linkEnds.at( i ) ).size( ) != 1 )
            {
                std::cerr << "Error when making angular position partials, light time corrections for " <<
                           lightTimeCorrections.at( linkEnds.at( i ) ).size( ) << " links found" << std::endl;
            }
            singleLinkLightTimeCorrections = lightTimeCorrections.at( linkEnds.at( i ) ).at( 0 );
        }

        // Create angular position partials for current link ends
        altimeterCrossoverPartials[ linkEnds[ i ] ] = createAltimeterCrossoverPartials(
                    linkEnds[ i ], bodyMap, parametersToEstimate, singleLinkLightTimeCorrections );
    }
    return altimeterCrossoverPartials;
}

}

}

#endif // TUDAT_CREATEALTIMETERCROSSOVERPARTIALS_H
