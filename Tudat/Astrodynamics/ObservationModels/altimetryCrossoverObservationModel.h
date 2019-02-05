/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_ALTIMETRYCROSSOVEROBSERVATIONMODEL_H
#define TUDAT_ALTIMETRYCROSSOVEROBSERVATIONMODEL_H

#include <map>

#include <functional>
#include <boost/make_shared.hpp>

#include <Eigen/Geometry>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

#include "Tudat/Astrodynamics/Ephemerides/simpleRotationalEphemeris.h"
#include "Tudat/Astrodynamics/ObservationModels/observationModel.h"
#include "Tudat/Astrodynamics/ObservationModels/lightTimeSolution.h"

namespace tudat
{

namespace observation_models
{

template< typename ObservationScalarType = double,
          typename TimeType = double >
class AltimetryCrossoverObservationModel: public ObservationModel< 1, ObservationScalarType, TimeType >
{
public:    
    typedef Eigen::Matrix< ObservationScalarType, 6, 1 > StateType;
    typedef Eigen::Matrix< ObservationScalarType, 3, 1 > PositionType;


    AltimetryCrossoverObservationModel(
            const std::function< StateType( const TimeType ) > firstArcBodyStateFunction,
            const std::function< StateType( const TimeType ) > secondArcBodyStateFunction,
            const std::string& centralBody,
            const std::map< double, double >& crossoverTimes,
            const std::shared_ptr< ObservationBias< 1 > > observationBiasCalculator = nullptr ):
        ObservationModel< 1, ObservationScalarType, TimeType >( altimetry_crossover, observationBiasCalculator ),
        firstArcBodyStateFunction_( firstArcBodyStateFunction ), secondArcBodyStateFunction_( secondArcBodyStateFunction ),
        centralBody_( centralBody ), crossoverTimes_( crossoverTimes ){ }

    //! Destructor
    ~AltimetryCrossoverObservationModel( ){ }

    Eigen::Matrix< ObservationScalarType, 1, 1 > computeIdealObservations(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime )
    {
        std::vector< double > linkEndTimes;
        std::vector< Eigen::Matrix< double, 6, 1 > > linkEndStates;

        return computeIdealObservationsWithLinkEndData(
                    time, linkEndAssociatedWithTime, linkEndTimes, linkEndStates );
    }

    //! Function to compute one-way range observable without any corrections.
    /*!
     *  Function to compute one-way range  observable without any corrections, i.e. the true physical range as computed
     *  from the defined link ends. Note that this observable does include light-time
     *  corrections, which represent physically true corrections. It does not include e.g. system-dependent measurement
     *  errors, such as biases or clock errors.
     *  The times and states of the link ends are also returned in full precision (determined by class template
     *  arguments). These states and times are returned by reference.
     *  \param time Time at which observable is to be evaluated.
     *  \param linkEndAssociatedWithTime Link end at which given time is valid, i.e. link end for which associated time
     *  is kept constant (to input value)
     *  \param linkEndTimes List of times at each link end during observation.
     *  \param linkEndStates List of states at each link end during observation.
     *  \return Ideal one-way range observable.
     */
    Eigen::Matrix< ObservationScalarType, 1, 1 > computeIdealObservationsWithLinkEndData(
                    const TimeType time,
                    const LinkEndType linkEndAssociatedWithTime,
                    std::vector< double >& linkEndTimes,
                    std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates )
    {
        Eigen::Matrix< ObservationScalarType, 1, 1 > crossoverAltimetryObservation;

        return crossoverAltimetryObservation;
    }

private:

    std::function< StateType( const TimeType ) > firstArcBodyStateFunction_;

    std::function< StateType( const TimeType ) > secondArcBodyStateFunction_;

    std::string centralBody_;

    std::map< double, double > crossoverTimes_;

};

} // namespace observation_models

} // namespace tudat

#endif // TUDAT_ALTIMETRYCROSSOVEROBSERVATIONMODEL_H
