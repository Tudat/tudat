/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


#include <algorithm>

#include <boost/make_shared.hpp>

#include "Tudat/SimulationSetup/EstimationSetup/createAltimeterCrossoverPartials.h"

namespace tudat
{

namespace observation_partials
{

using namespace ephemerides;

//! Function to generate one-way range partial wrt an initial position of a body.
std::shared_ptr< altimeterCrossoverPartial > createaltimeterCrossoverPartialWrtBodyPosition(
        const observation_models::LinkEnds altimeterCrossoverLinkEnds,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::string bodyToEstimate,
        const std::shared_ptr< AltimeterCrossoverScaling > altimeterCrossoverScaler,
        const std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > >&
        lightTimeCorrectionPartialObjects  )
{
    // Create position partials of link ends for current body position
    std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > > positionPartials =
            createCartesianStatePartialsWrtBodyState( altimeterCrossoverLinkEnds, bodyMap, bodyToEstimate );

    // Create one-range partials if any position partials are created (i.e. if any dependency exists).
    std::shared_ptr< altimeterCrossoverPartial > altimeterCrossoverPartial;
    if( positionPartials.size( ) > 0 )
    {
        altimeterCrossoverPartial = std::make_shared< altimeterCrossoverPartial >(
                    altimeterCrossoverScaler, positionPartials, std::make_pair(
                        estimatable_parameters::initial_body_state, std::make_pair( bodyToEstimate, "" ) ),
                    lightTimeCorrectionPartialObjects );
    }

    // Return range partial object (nullptr if no dependency exists).
    return altimeterCrossoverPartial;
}

}

}


