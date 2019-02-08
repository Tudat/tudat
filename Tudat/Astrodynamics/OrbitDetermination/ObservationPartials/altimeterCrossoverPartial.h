/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_ALTIMETERCROSSOVERPARTIAL_H
#define TUDAT_ALTIMETERCROSSOVERPARTIAL_H

#include <functional>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/Ephemerides/ephemeris.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/observationPartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/positionPartials.h"
#include "Tudat/Astrodynamics/OrbitDetermination/LightTimeCorrectionPartials/lightTimeCorrectionPartial.h"
#include "Tudat/Astrodynamics/ObservationModels/linkTypeDefs.h"

namespace tudat
{

namespace observation_partials
{


class AltimeterCrossoverScaling: public PositionPartialScaling
{
public:

    //! Destructor
    ~AltimeterCrossoverScaling( ){ }


    void update( const std::vector< Eigen::Vector6d >& linkEndStates,
                 const std::vector< double >& times,
                 const observation_models::LinkEndType fixedLinkEnd,
                 const Eigen::VectorXd currentObservation = Eigen::VectorXd::Constant( 1, TUDAT_NAN ) ){ }
};

//! Class to compute the partial derivatives of a altimeter crossover observation partial.
class AltimeterCrossoverPartial: public ObservationPartial< 1 >
{

public:

    typedef std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > AltimeterCrossoverPartialReturnType;
    typedef std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > SingleAltimeterCrossoverPartialReturnType;

    AltimeterCrossoverPartial(
            const std::shared_ptr< AltimeterCrossoverScaling > altimeterCrossoverScaler,
            const std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > >& positionPartialList,
            const estimatable_parameters::EstimatebleParameterIdentifier parameterIdentifier,
            const std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > >&
            lighTimeCorrectionPartials =
            std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > >( ) ):
        ObservationPartial< 1 >( parameterIdentifier ), altimeterCrossoverScaler_( altimeterCrossoverScaler ),
        positionPartialList_( positionPartialList )
    {
        std::pair< std::function< SingleAltimeterCrossoverPartialReturnType(
                    const std::vector< Eigen::Vector6d >&, const std::vector< double >& ) >,
                bool > lightTimeCorrectionPartial;

//        // Create light time correction partial functions
//        for( unsigned int i = 0; i < lighTimeCorrectionPartials.size( ); i++ )
//        {
//            lightTimeCorrectionPartial = getLightTimeParameterPartialFunction(
//                        parameterIdentifier, lighTimeCorrectionPartials.at( i ) );
//            if( lightTimeCorrectionPartial.second != 0 )
//            {
//                lighTimeCorrectionPartialsFunctions_.push_back( lightTimeCorrectionPartial.first );
//            }
//        }
    }

    //! Destructor.
    ~AltimeterCrossoverPartial( ) { }

    //! Function to calculate the observation partial(s) at required time and state
    /*!
     *  Function to calculate the observation partial(s) at required time and state. State and time
     *  are typically obtained from evaluation of observation model.
     *  \param states Link end states. Index maps to link end for a given ObsevableType through getLinkEndIndex function.
     *  \param times Link end time.
     *  \param linkEndOfFixedTime Link end that is kept fixed when computing the observable.
     *  \param currentObservation Value of the observation for which the partial is to be computed (default NaN for
     *  compatibility purposes)
     *  \return Vector of pairs containing partial values and associated times.
     */
    virtual AltimeterCrossoverPartialReturnType calculatePartial(
            const std::vector< Eigen::Vector6d >& states,
            const std::vector< double >& times,
            const observation_models::LinkEndType linkEndOfFixedTime,
            const Eigen::Vector1d& currentObservation = Eigen::Vector1d::Constant( TUDAT_NAN ) );

//    //! Function to get scaling object used for mapping partials of positions to partials of observable
//    /*!
//     * Function to get scaling object used for mapping partials of positions to partials of observable
//     * \return
//     */
//    std::shared_ptr< AltimeterCrossoverScaling > getaltimeterCrossoverScaler( )
//    {
//        return altimeterCrossoverScaler_;
//    }

//    //! Function to get the number of light-time correction partial functions.
//    /*!
//     * Number of light-time correction partial functions.
//     * \return Number of light-time correction partial functions.
//     */
//    int getNumberOfLighTimeCorrectionPartialsFunctions( )
//    {
//        return lighTimeCorrectionPartialsFunctions_.size( );
//    }

protected:

    //! Scaling object used for mapping partials of positions to partials of observable
    std::shared_ptr< AltimeterCrossoverScaling > altimeterCrossoverScaler_;

    //! List of position partials per link end.
    std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > > positionPartialList_;

    //! Iterator over list of position partials per link end.
    std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > >::iterator positionPartialIterator_;

//    //! List of light-time correction partial functions.
//    std::vector< std::function< SingleAltimeterCrossoverPartialReturnType(
//            const std::vector< Eigen::Vector6d >&, const std::vector< double >& ) > >
//    lighTimeCorrectionPartialsFunctions_;

//    //! List of light-time correction partial objects.
//    std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > > lighTimeCorrectionPartials_;

    //! Pre-declared state variable to be used in calculatePartial function.
    Eigen::Vector6d currentState_;

    //! Pre-declared time variable to be used in calculatePartial function.
    double currentTime_;

};

}

}



#endif // TUDAT_ALTIMETERCROSSOVERPARTIAL_H
