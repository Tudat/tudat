/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/ephemerides/simpleRotationalEphemeris.h"
#include "tudat/simulation/estimation_setup/createCartesianStatePartials.h"


namespace tudat
{

namespace observation_partials
{

//! Function to return partial(s) of position of ground station(s) w.r.t. state of a single body.
std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > > createCartesianStatePartialsWrtBodyState(
        const observation_models::LinkEnds& linkEnds,
        const simulation_setup::SystemOfBodies& bodies,
        const std::string bodyToEstimate )
{
    // Declare data map to return.
    std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > > partialMap;

    // Declare local variable to use in loop
    std::string currentBodyName;

    // Iterate over all like ends.
    for( observation_models::LinkEnds::const_iterator linkEndIterator = linkEnds.begin( );
         linkEndIterator != linkEnds.end( ); linkEndIterator++ )
    {
        // Check if current link end is on body that is requested.
        currentBodyName = linkEndIterator->second.first;
        if( bodyToEstimate == currentBodyName )
        {
            // Create partial
            partialMap[ linkEndIterator->first ] = std::make_shared< CartesianStatePartialWrtCartesianState >( );

        }
    }

    return partialMap;
}

//! Function to return partial(s) of position of ground station(s) w.r.t. rotational state of a single body.
std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > > createCartesianStatePartialsWrtBodyRotationalState(
        const observation_models::LinkEnds& linkEnds,
        const simulation_setup::SystemOfBodies& bodies,
        const std::string& bodyToEstimate )
{
    // Declare data map to return.
    std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > > partialMap;

    // Declare local variable to use in loop
    std::string currentBodyName;

    // Iterate over all like ends.
    for( observation_models::LinkEnds::const_iterator linkEndIterator = linkEnds.begin( );
         linkEndIterator != linkEnds.end( ); linkEndIterator++ )
    {
        // Check if current link end is on body that is requested.
        currentBodyName = linkEndIterator->second.first;
        if( bodyToEstimate == currentBodyName )
        {
            std::shared_ptr< simulation_setup::Body > currentBody = bodies.at( currentBodyName );

            if( currentBody->getGroundStationMap( ).count( linkEndIterator->second.second ) == 0 )
            {
                throw std::runtime_error(
                            "Error when making cartesian state partial w.r.t. rotation parameter, ground station " +
                            linkEndIterator->second.second + " not found on body " + linkEndIterator->second.first );
            }

            // Set ground station position function
            std::function< Eigen::Vector3d( const double ) > groundStationPositionFunction =
                    std::bind( &ground_stations::GroundStationState::getCartesianPositionInTime,
                                 currentBody->getGroundStation( linkEndIterator->second.second )->getNominalStationState( ),
                                 std::placeholders::_1, basic_astrodynamics::JULIAN_DAY_ON_J2000 );

            // Create partial
            partialMap[ linkEndIterator->first ] = std::make_shared< CartesianStatePartialWrtRotationMatrixParameter >(
                        std::make_shared< RotationMatrixPartialWrtRotationalState >(
                            std::bind( &ephemerides::RotationalEphemeris::getRotationToBaseFrame,
                                         currentBody->getRotationalEphemeris( ), std::placeholders::_1 ) ), groundStationPositionFunction );
        }
    }

    return partialMap;
}

//! Function to return partial object(s) of position of reference point w.r.t. a (double) parameter.
std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > > createCartesianStatePartialsWrtParameter(
        const observation_models::LinkEnds linkEnds,
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameterToEstimate )
{
    using namespace ephemerides;

    // Declare data map to return.
    std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > > partialMap;

    // Declare local variable to use in loop
    std::shared_ptr< simulation_setup::Body > currentBody;
    std::string currentBodyName;

    // Iterate over all like ends.
    for( observation_models::LinkEnds::const_iterator linkEndIterator = linkEnds.begin( );
         linkEndIterator != linkEnds.end( ); linkEndIterator++ )
    {
        // Check if current link end body corresponds to body with property to estimate.
        if( ( linkEndIterator->second.first == parameterToEstimate->getParameterName( ).second.first ) &&
                ( linkEndIterator->second.second != "" ) )
        {
            // Set current body name and object.
            currentBodyName = linkEndIterator->second.first;
            currentBody = bodies.at( currentBodyName );

            // Check if parameter is a rotation model property, in which case a CartesianStatePartialWrtRotationMatrixParameter,
            // with the rotation matrix partial created from createRotationMatrixPartialsWrtParameter function.
            if( estimatable_parameters::isParameterRotationMatrixProperty( parameterToEstimate->getParameterName( ).first ) )
            {
                // Set ground station position function
                std::function< Eigen::Vector3d( const double ) > groundStationPositionFunction =
                        std::bind( &ground_stations::GroundStationState::getCartesianPositionInTime,
                                     ( currentBody )->getGroundStation( linkEndIterator->second.second )
                                     ->getNominalStationState( ), std::placeholders::_1, basic_astrodynamics::JULIAN_DAY_ON_J2000 );

                // Create parameter partial object.
                partialMap[ linkEndIterator->first ] = std::make_shared< CartesianStatePartialWrtRotationMatrixParameter >(
                            createRotationMatrixPartialsWrtParameter( bodies, parameterToEstimate ),
                            groundStationPositionFunction );
            }
            else
            {
                // Check which parameter is requested and create position partial if direct dependency between position and
                // parameter exists.
                switch( parameterToEstimate->getParameterName( ).first )
                {
                case estimatable_parameters::gravitational_parameter:
                    break;
                case estimatable_parameters::constant_drag_coefficient:
                    break;
                case estimatable_parameters::radiation_pressure_coefficient:
                    break;
                case estimatable_parameters::ppn_parameter_gamma:
                    break;
                case estimatable_parameters::ppn_parameter_beta:
                    break;
                case estimatable_parameters::equivalence_principle_lpi_violation_parameter:
                    break;
                case estimatable_parameters::mean_moment_of_inertia:
                    break;
                case estimatable_parameters::direct_dissipation_tidal_time_lag:
                    break;
                default:

                    std::string errorMessage =
                            "Parameter " + std::to_string(
                                parameterToEstimate->getParameterName( ).first ) +
                            " not implemented when making position partial";
                    throw std::runtime_error( errorMessage );
                    break;
                }
            }
        }
    }

    return partialMap;
}

//! Function to return partial(s) of position of ground station(s) w.r.t. a (vector) parameter.
std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > > createCartesianStatePartialsWrtParameter(
        const observation_models::LinkEnds linkEnds,
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameterToEstimate )
{
    using namespace ephemerides;

    // Declare data map to return.
    std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > > partialMap;

    // Declare local variable to use in loop
    std::shared_ptr< simulation_setup::Body > currentBody;
    std::string currentBodyName;

    // Iterate over all like ends.
    for( observation_models::LinkEnds::const_iterator linkEndIterator = linkEnds.begin( );
         linkEndIterator != linkEnds.end( ); linkEndIterator++ )
    {
        // Check if current link end body corresponds to body with property to estimate.
        if( linkEndIterator->second.first == parameterToEstimate->getParameterName( ).second.first &&
                linkEndIterator->second.second != "" )
        {
            // Set current body name and object.
            currentBodyName = linkEndIterator->second.first;
            currentBody = bodies.at( currentBodyName );

            // Check if parameter is a rotation model property, in which case a CartesianStatePartialWrtRotationMatrixParameter,
            // with the rotation matrix partial created from createRotationMatrixPartialsWrtParameter function.
            if( estimatable_parameters::isParameterRotationMatrixProperty( parameterToEstimate->getParameterName( ).first ) )
            {
                if( currentBody->getGroundStationMap( ).count( linkEndIterator->second.second ) == 0 )
                {
                    throw std::runtime_error(
                                "Error when making cartesian state partial w.r.t. rotation parameter, ground station " +
                                linkEndIterator->second.second + " not found on body " + linkEndIterator->second.first );
                }

                // Set ground station position function
                std::function< Eigen::Vector3d( const double ) > groundStationPositionFunction =
                        std::bind( &ground_stations::GroundStationState::getCartesianPositionInTime,
                                     currentBody->getGroundStation( linkEndIterator->second.second )
                                     ->getNominalStationState( ), std::placeholders::_1, basic_astrodynamics::JULIAN_DAY_ON_J2000 );

                // Create parameter partial object.
                partialMap[ linkEndIterator->first ] = std::make_shared< CartesianStatePartialWrtRotationMatrixParameter >(
                            createRotationMatrixPartialsWrtParameter( bodies, parameterToEstimate ),
                            groundStationPositionFunction );
            }
            else
            {
                // Check which parameter is requested and create position partial if direct dependency between position and
                // parameter exists.
                switch( parameterToEstimate->getParameterName( ).first )
                {
                case estimatable_parameters::spherical_harmonics_cosine_coefficient_block:
                    break;
                case estimatable_parameters::spherical_harmonics_sine_coefficient_block:
                    break;
                case estimatable_parameters::arc_wise_radiation_pressure_coefficient:
                    break;
                case estimatable_parameters::arc_wise_constant_drag_coefficient:
                    break;
                case estimatable_parameters::constant_additive_observation_bias:
                    break;
                case estimatable_parameters::arcwise_constant_additive_observation_bias:
                    break;
                case estimatable_parameters::constant_relative_observation_bias:
                    break;
                case estimatable_parameters::arcwise_constant_relative_observation_bias:
                    break;
                case estimatable_parameters::empirical_acceleration_coefficients:
                    break;
                case estimatable_parameters::arc_wise_empirical_acceleration_coefficients:
                    break;
                case estimatable_parameters::full_degree_tidal_love_number:
                    break;
                case estimatable_parameters::single_degree_variable_tidal_love_number:
                    break;
                case estimatable_parameters::ground_station_position:

                    // Check if current link end station is same station as that of which position is to be estimated.
                    if( linkEndIterator->second.second == parameterToEstimate->getParameterName( ).second.second )
                    {
                        if( currentBody->getRotationalEphemeris( ) == nullptr )
                        {
                            throw std::runtime_error(
                                        "Warning, body's rotation model is not found when making position w.r.t. ground station position position partial" );
                        }
                        if( currentBody->getGroundStationMap( ).count( linkEndIterator->second.second ) == 0 )
                        {
                                std::runtime_error( "Warning, ground station " + linkEndIterator->second.second +
                                           "not found when making ground station position position partial" );
                        }

                        // Create partial object.
                        partialMap[ linkEndIterator->first ] = std::make_shared< CartesianPartialWrtBodyFixedPosition >(
                                    currentBody->getRotationalEphemeris( ) );
                    }
                    break;
                default:
                    std::string errorMessage =
                            "Parameter " + std::to_string(
                                parameterToEstimate->getParameterName( ).first ) +
                            " not implemented when making position partial";
                    throw std::runtime_error( errorMessage );
                    break;
                }
            }

        }

    }
    return partialMap;
}

//! Function to create partial object(s) of rotation matrix wrt translational state
std::shared_ptr< RotationMatrixPartial > createRotationMatrixPartialsWrtTranslationalState(
        const std::shared_ptr< simulation_setup::Body > currentBody )
{
    std::shared_ptr< RotationMatrixPartial > rotationMatrixPartial;
    if( std::dynamic_pointer_cast< ephemerides::SynchronousRotationalEphemeris >(
                currentBody->getRotationalEphemeris( ) ) != nullptr )
    {
        rotationMatrixPartial = std::make_shared< SynchronousRotationMatrixPartialWrtTranslationalState >(
                    std::dynamic_pointer_cast< ephemerides::SynchronousRotationalEphemeris >(
                        currentBody->getRotationalEphemeris( ) ) );
    }
    return rotationMatrixPartial;
}

//! Function to create partial object(s) of rotation matrix wrt a (double) parameter.
std::shared_ptr< RotationMatrixPartial > createRotationMatrixPartialsWrtParameter(
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameterToEstimate )
{
    using namespace simulation_setup;
    using namespace ephemerides;

    // Declare return object.
    std::shared_ptr< RotationMatrixPartial >  rotationMatrixPartial;

    // Get body for rotation of which partial is to be created.
    std::shared_ptr< Body > currentBody = bodies.at( parameterToEstimate->getParameterName( ).second.first );

    // Check for which rotation model parameter the partial object is to be created.
    switch( parameterToEstimate->getParameterName( ).first )
    {
    case estimatable_parameters::constant_rotation_rate:

        if( std::dynamic_pointer_cast< ephemerides::SimpleRotationalEphemeris >(
                    currentBody->getRotationalEphemeris( ) ) == nullptr )
        {
            throw std::runtime_error( "Warning, body's rotation model is not simple when making position w.r.t. constant rtoation rate partial" ) ;
        }

        // Create rotation matrix partial object
        rotationMatrixPartial = std::make_shared< RotationMatrixPartialWrtConstantRotationRate >(
                    std::dynamic_pointer_cast< SimpleRotationalEphemeris>( currentBody->getRotationalEphemeris( ) ) );
        break;

    case estimatable_parameters::core_factor:

        if( std::dynamic_pointer_cast< ephemerides::PlanetaryRotationModel >(
                    currentBody->getRotationalEphemeris() ) == nullptr )
        {
            std::string errorMessage = "Warning, body's rotation model is not a full planetary rotational model when making"
                                       "position w.r.t. core factor partial";
            throw std::runtime_error( errorMessage );
        }

        // Create rotation matrix partial object
        rotationMatrixPartial = std::make_shared< RotationMatrixPartialWrtCoreFactor >(
                    std::dynamic_pointer_cast< PlanetaryRotationModel >( currentBody->getRotationalEphemeris( ) ) );

        break;

    case estimatable_parameters::free_core_nutation_rate:

        if( std::dynamic_pointer_cast< ephemerides::PlanetaryRotationModel >(
                    currentBody->getRotationalEphemeris() ) == nullptr ){
            std::string errorMessage = "Warning, body's rotation model is not a full planetary rotational model when making"
                                       "position w.r.t. polar motion amplitude partial";
            throw std::runtime_error( errorMessage );
        }

        // Create rotation matrix partial object
        rotationMatrixPartial = std::make_shared< RotationMatrixPartialWrtFreeCoreNutationRate >(
                    std::dynamic_pointer_cast< PlanetaryRotationModel >( currentBody->getRotationalEphemeris() ));

        break;
    case estimatable_parameters::scaled_longitude_libration_amplitude:

        if( std::dynamic_pointer_cast< ephemerides::SynchronousRotationalEphemeris >(
                    currentBody->getRotationalEphemeris() ) == nullptr ){
            std::string errorMessage = "Warning, body's rotation model is not a synchronous rotational model when making"
                                       "position w.r.t. longitude libration amplitude partial";
            throw std::runtime_error( errorMessage );
        }

        // Create rotation matrix partial object
        rotationMatrixPartial = std::make_shared< RotationMatrixPartialWrtScaledLongitudeLibrationAmplitude >(
                    std::dynamic_pointer_cast< SynchronousRotationalEphemeris >( currentBody->getRotationalEphemeris( ) ) );

        break;
    default:
        std::string errorMessage = "Warning, rotation matrix partial not implemented for parameter " +
                std::to_string( parameterToEstimate->getParameterName( ).first );
        throw std::runtime_error( errorMessage );

        break;
    }

    return rotationMatrixPartial;


}

//! Function to create partial object(s) of rotation matrix wrt a (vector) parameter.
std::shared_ptr< RotationMatrixPartial > createRotationMatrixPartialsWrtParameter(
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameterToEstimate )

{
    using namespace simulation_setup;
    using namespace ephemerides;

    // Declare return object.
    std::shared_ptr< RotationMatrixPartial >  rotationMatrixPartial;

    // Get body for rotation of which partial is to be created.
    std::shared_ptr< Body > currentBody = bodies.at( parameterToEstimate->getParameterName( ).second.first );

    // Check for which rotation model parameter the partial object is to be created.
    switch( parameterToEstimate->getParameterName( ).first )
    {
    case estimatable_parameters::rotation_pole_position:


        if( std::dynamic_pointer_cast< ephemerides::SimpleRotationalEphemeris >(
                    currentBody->getRotationalEphemeris( ) ) == nullptr )
        {
            std::string errorMessage = "Warning, body's rotation model is not simple when making position w.r.t. pole position partial";
            throw std::runtime_error( errorMessage );
        }

        // Create rotation matrix partial object
        rotationMatrixPartial = std::make_shared< RotationMatrixPartialWrtPoleOrientation >(
                    std::dynamic_pointer_cast< SimpleRotationalEphemeris>( currentBody->getRotationalEphemeris( ) ) );
        break;

    case estimatable_parameters::periodic_spin_variation:

        if( std::dynamic_pointer_cast< ephemerides::PlanetaryRotationModel >(
                    currentBody->getRotationalEphemeris() ) == nullptr ){
            std::string errorMessage = "Warning, body's rotation model is not a full planetary rotational model when making"
                                       "position w.r.t. periodic spin variation partial";
            throw std::runtime_error( errorMessage );
        }

        // Create rotation matrix partial object
        rotationMatrixPartial = std::make_shared< RotationMatrixPartialWrtPeriodicSpinVariations >(
                    std::dynamic_pointer_cast< PlanetaryRotationModel >( currentBody->getRotationalEphemeris() ));
        break;

    case estimatable_parameters::polar_motion_amplitude:

        if( std::dynamic_pointer_cast< ephemerides::PlanetaryRotationModel >(
                    currentBody->getRotationalEphemeris() ) == nullptr ){
            std::string errorMessage = "Warning, body's rotation model is not a full planetary rotational model when making"
                                       "position w.r.t. polar motion amplitude partial";
            throw std::runtime_error( errorMessage );
        }
        else{

            if ( std::dynamic_pointer_cast< ephemerides::PlanetaryRotationModel >( currentBody->getRotationalEphemeris() )
                 ->getPlanetaryOrientationAngleCalculator()->getXpolarMotionCoefficients().size()
                 != std::dynamic_pointer_cast< ephemerides::PlanetaryRotationModel >( currentBody->getRotationalEphemeris() )
                 ->getPlanetaryOrientationAngleCalculator()->getYpolarMotionCoefficients().size() ){

                throw std::runtime_error( "Error, unconsistent sizes when comparing x and y polar motion"
                                          "amplitude coefficients." );
            }

            // Create rotation matrix partial object
            rotationMatrixPartial = std::make_shared< RotationMatrixPartialWrtPolarMotionAmplitude >(
                        std::dynamic_pointer_cast< PlanetaryRotationModel >( currentBody->getRotationalEphemeris() ));

        }

        break;
    default:
        std::string errorMessage = "Warning, rotation matrix partial not implemented for parameter " +
                std::to_string( parameterToEstimate->getParameterName( ).first );
        throw std::runtime_error( errorMessage );
        break;
    }

    return rotationMatrixPartial;

}

//! Function to create an objects that computes the partial derivatives of a three-dimensional position observable w.r.t.
//! the position of a body.
std::shared_ptr< PositionObervationPartial > createPositionObservablePartialWrtPosition(
        const  observation_models::LinkEnds linkEnds,
        const simulation_setup::SystemOfBodies& bodies,
        const std::string bodyToEstimate,
        const std::shared_ptr< PositionObservationScaling > positionObservableScaler )
{
    std::map<  observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > > positionPartials =
            createCartesianStatePartialsWrtBodyState( linkEnds, bodies, bodyToEstimate );
    std::shared_ptr< PositionObervationPartial > positionObervationPartial;

    if( positionPartials.size( ) > 0 )
    {
        positionObervationPartial = std::make_shared< PositionObervationPartial >(
                    positionObservableScaler, positionPartials, std::make_pair(
                        estimatable_parameters::initial_body_state, std::make_pair( bodyToEstimate, "") ) );
    }

    return positionObervationPartial;
}

//! Function to create an objects that computes the partial derivatives of a three-dimensional position observable w.r.t.
//! the Velocity of a body.
std::shared_ptr< VelocityObervationPartial > createVelocityObservablePartialWrtVelocity(
        const  observation_models::LinkEnds linkEnds,
        const simulation_setup::SystemOfBodies& bodies,
        const std::string bodyToEstimate,
        const std::shared_ptr< VelocityObservationScaling > velocityObservableScaler )
{
    std::map<  observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > > velocityPartials =
            createCartesianStatePartialsWrtBodyState( linkEnds, bodies, bodyToEstimate );
    std::shared_ptr< VelocityObervationPartial > velocityObervationPartial;

    if( velocityPartials.size( ) > 0 )
    {
        velocityObervationPartial = std::make_shared< VelocityObervationPartial >(
                    velocityObservableScaler, velocityPartials, std::make_pair(
                        estimatable_parameters::initial_body_state, std::make_pair( bodyToEstimate, "") ) );
    }

    return velocityObervationPartial;
}

}

}
