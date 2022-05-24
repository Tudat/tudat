/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEANGULARPOSITIONPARTIALS_H
#define TUDAT_CREATEANGULARPOSITIONPARTIALS_H

#include <vector>
#include <map>

#include <memory>

#include <Eigen/Core>

#include "tudat/math/interpolators/interpolator.h"

#include "tudat/simulation/estimation_setup/createCartesianStatePartials.h"
#include "tudat/astro/orbit_determination/observation_partials/angularPositionPartial.h"
#include "tudat/simulation/estimation_setup/createLightTimeCorrectionPartials.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/initialTranslationalState.h"
#include "tudat/astro/observation_models/linkTypeDefs.h"


namespace tudat
{

namespace observation_partials
{

//! Function to generate angular position partial wrt a position of a body.
/*!
 *  Function to generate angular position partial wrt a position of a body, for a single link ends (which must contain a
 *  transmitter and receiever  linkEndType).
 *  \param angularPositionLinkEnds Link ends (transmitter and receiever) for which partials are to be calculated
 *  (i.e. for which angular position observation observations are to be processed).
 *  \param bodies List of all bodies, for creating angular position partial.
 *  \param bodyToEstimate Name of body wrt position of which a partial is to be created.
 *  \param angularPositionScaler Object scale position partials to angular position partials for current link ends.
 *  \param lightTimeCorrectionPartialObjects List of light time correction partials to be used (empty by default)
 *  \return Angular position partial object wrt a current position of a body (is nullptr if no parameter dependency exists).
 */
std::shared_ptr< AngularPositionPartial > createAngularPositionPartialWrtBodyPosition(
        const observation_models::LinkEnds angularPositionLinkEnds,
        const simulation_setup::SystemOfBodies& bodies,
        const std::string bodyToEstimate,
        const std::shared_ptr< AngularPositionScaling > angularPositionScaler,
        const std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > >&
        lightTimeCorrectionPartialObjects =
        std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > >( ) );

//! Function to generate angular position partial wrt a single  parameter.
/*!
 *  Function to generate angular position partial wrt a single  parameter, for a single link ends (which must contain a
 *  transmitter and receiever linkEndType).
 *  \tparam ParameterType Type of parameter (double for size 1, VectorXd for larger size).
 *  \param angularPositionLinkEnds Link ends (transmitter and receiever) for which angular position partials are to be
 *  calculated (i.e. for which  angular position observations are to be processed).
 *  \param bodies List of all bodies, for creating angular position partial.
 *  \param parameterToEstimate Object of current parameter that is to be estimated.
 *  \param angularPositionScaler Object scale position partials to angular position partials for current link ends.
 *  \param lightTimeCorrectionPartialObjects List of light time correction partials to be used (empty by default)
 *  \return Angular position partial object wrt a single parameter (is nullptr if no parameter dependency exists).
 */
template< typename ParameterType >
std::shared_ptr< AngularPositionPartial > createAngularPositionPartialWrtParameter(
        const observation_models::LinkEnds angularPositionLinkEnds,
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameter< ParameterType > > parameterToEstimate,
        const std::shared_ptr< AngularPositionScaling > angularPositionScaler,
        const std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > >&
        lightTimeCorrectionPartialObjects =
        std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > >( ) )
{
    std::shared_ptr< AngularPositionPartial > angularPositionPartial;

    {
        std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > > positionPartials =
                createCartesianStatePartialsWrtParameter( angularPositionLinkEnds, bodies, parameterToEstimate );

        std::shared_ptr< AngularPositionPartial > testAngularPositionPartial = std::make_shared< AngularPositionPartial >(
                    angularPositionScaler, positionPartials, parameterToEstimate->getParameterName( ),
                    lightTimeCorrectionPartialObjects );
        // Create angular position partials if any position partials are created (i.e. if any dependency exists).
        if( positionPartials.size( ) > 0 || testAngularPositionPartial->getNumberOfLighTimeCorrectionPartialsFunctions( ) )
        {
            angularPositionPartial = testAngularPositionPartial;
        }
    }
    return angularPositionPartial;
}

//! Function to generate angular position partials and associated scaler for single link end.
/*!
 *  Function to generate angular position partials and associated scaler for all parameters that are to be estimated,
 *  for a single link ends.
 *  The set of parameters and bodies that are to be estimated, as well as the set of link ends
 *  (each of which must contain a transmitter and receiever linkEndType) that are to be used.
 *  \param angularPositionLinkEnds Link ends (transmitter and receiever) for which angular position partials are to be
 *  calculated (i.e. for which angular position observations are to be processed).
 *  \param bodies List of all bodies, for creating angular position partials.
 *  \param parametersToEstimate Set of parameters that are to be estimated (in addition to initial states of
 *  requested bodies)
 *  \param lightTimeCorrections List of light time correction partials to be used (empty by default)
 *  \return Set of observation partials with associated indices in complete vector of parameters that are estimated,
 *  representing all  necessary angular position partials of a single link end, and AngularPositionScaling, object, used for
 *  scaling the position partial members of all AngularPositionPartials in link end.
 */
template< typename ParameterType >
std::pair< std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< 2 > > >,
std::shared_ptr< PositionPartialScaling > >
createAngularPositionPartials(
        const observation_models::LinkEnds angularPositionLinkEnds,
        const simulation_setup::SystemOfBodies& bodies,
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
    std::shared_ptr< AngularPositionScaling > angularPositionScaling =
            std::make_shared< AngularPositionScaling >( );

    SingleLinkObservationTwoPartialList angularPositionPartials;


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
        std::shared_ptr< AngularPositionPartial > currentAngularPositionPartial = createAngularPositionPartialWrtBodyPosition(
                    angularPositionLinkEnds, bodies, acceleratedBody, angularPositionScaling,
                    lightTimeCorrectionPartialObjects );

        // Check if partial is non-nullptr (i.e. whether dependency exists between current angular position and current body)
        if( currentAngularPositionPartial != nullptr )
        {
            // Add partial to the list.
            currentPair = std::pair< int, int >( currentIndex, 6 );
            angularPositionPartials[ currentPair ] = currentAngularPositionPartial;
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
        std::shared_ptr< AngularPositionPartial > currentAngularPositionPartial = createAngularPositionPartialWrtParameter(
                    angularPositionLinkEnds, bodies, parameterIterator->second, angularPositionScaling,
                    lightTimeCorrectionPartialObjects );

        if( currentAngularPositionPartial != nullptr )
        {
            // Add partial to the list.
            currentPair = std::pair< int, int >( parameterIterator->first, 1 );
            angularPositionPartials[ currentPair ] = currentAngularPositionPartial;
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
        std::shared_ptr< ObservationPartial< 2 > > currentAngularPositionPartial;

        if( !isParameterObservationLinkProperty( parameterIterator->second->getParameterName( ).first )  )
        {
            currentAngularPositionPartial = createAngularPositionPartialWrtParameter(
                        angularPositionLinkEnds, bodies, parameterIterator->second, angularPositionScaling,
                        lightTimeCorrectionPartialObjects );
        }
        else
        {
            currentAngularPositionPartial = createObservationPartialWrtLinkProperty< 2 >(
                        angularPositionLinkEnds, observation_models::angular_position, parameterIterator->second );
        }

        // Check if partial is non-nullptr (i.e. whether dependency exists between current observable and current parameter)
        if( currentAngularPositionPartial != nullptr )
        {
            // Add partial to the list.
            currentPair = std::pair< int, int >( parameterIterator->first,
                                                 parameterIterator->second->getParameterSize( ) );
            angularPositionPartials[ currentPair ] = currentAngularPositionPartial;
        }

    }
    return std::make_pair( angularPositionPartials, angularPositionScaling );
}

//! Function to generate angular position partials for all parameters that are to be estimated, for all sets of link ends.
/*!
 *  Function to generate angular position partials for all parameters that are to be estimated, for all sets of link ends.
 *  The angular position partials are generated per set of link ends. The set of parameters and bodies that are to be
 *  estimated, as well as the set of link ends (each of which must contain a transmitter and receiever linkEndType)
 *  that are to be used.
 *  \param linkEnds Vector of all link ends for which angular position partials are to be calculated (i.e. for which one-way
 *  angular position observations are  to be processed).
 *  \param bodies List of all bodies, for creating angular position partials.
 *  \param parametersToEstimate Set of parameters that are to be estimated (in addition to initial states
 *  of requested bodies)
 *  \param lightTimeCorrections List of light time correction partials to be used (empty by default)
 *  \return Map of SingleLinkObservationPartialList, representing all necessary angular position partials of a single link
 *  end, and AngularPositionScaling, object, used for scaling the position partial members of all AngularPositionPartials in
 *  link end.
 */
template< typename ParameterType >
std::map< observation_models::LinkEnds, std::pair< SingleLinkObservationTwoPartialList,
std::shared_ptr< PositionPartialScaling > > >
createAngularPositionPartials(
        const std::vector< observation_models::LinkEnds > linkEnds,
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
        const std::map< observation_models::LinkEnds,
        std::vector< std::vector< std::shared_ptr< observation_models::LightTimeCorrection > > > >& lightTimeCorrections =
        std::map< observation_models::LinkEnds,
        std::vector< std::vector< std::shared_ptr< observation_models::LightTimeCorrection > > > >( ) )
{
    // Declare return list.
    std::map< observation_models::LinkEnds, std::pair< SingleLinkObservationTwoPartialList,
            std::shared_ptr< PositionPartialScaling > > > angularPositionPartials;

    // Iterate over all link ends.
    for( unsigned int i = 0; i < linkEnds.size( ); i++ )
    {
        // Check if required link end types are present
        if( ( linkEnds[ i ].count( observation_models::receiver ) == 0 ) ||
                ( linkEnds[ i ].count( observation_models::transmitter ) == 0 ) )
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
        angularPositionPartials[ linkEnds[ i ] ] = createAngularPositionPartials(
                    linkEnds[ i ], bodies, parametersToEstimate, singleLinkLightTimeCorrections );
    }
    return angularPositionPartials;
}

}

}

#endif // TUDAT_CREATEANGULARPOSITIONPARTIALS_H
