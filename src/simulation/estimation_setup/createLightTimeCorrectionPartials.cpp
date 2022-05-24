/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


#include "tudat/astro/observation_models/corrections/firstOrderRelativisticCorrection.h"
#include "tudat/astro/orbit_determination/observation_partials/firstOrderRelativisticPartial.h"
#include "tudat/simulation/estimation_setup/createLightTimeCorrectionPartials.h"


namespace tudat
{

namespace observation_partials
{

//! Function to create a partial objects from list of light time corrections.
std::vector< std::shared_ptr< LightTimeCorrectionPartial > > createLightTimeCorrectionPartials(
        const std::vector< std::shared_ptr< observation_models::LightTimeCorrection > >& lightTimeCorrectionList )
{
    std::vector< std::shared_ptr< LightTimeCorrectionPartial > > partialList;

    // Iterate over all light time corrections
    for( unsigned int i = 0; i < lightTimeCorrectionList.size( ); i++ )
    {
        // Check type of light time correction
        switch( lightTimeCorrectionList.at( i )->getLightTimeCorrectionType( ) )
        {
        case observation_models::first_order_relativistic:
        {
            std::shared_ptr< observation_models::FirstOrderLightTimeCorrectionCalculator > currentCorrection =
                    std::dynamic_pointer_cast< observation_models::FirstOrderLightTimeCorrectionCalculator >(
                        lightTimeCorrectionList.at( i ) );
            if( currentCorrection == nullptr )
            {
                throw std::runtime_error( "Error when making first order light time correction partial, type id observation_models::first_order_relativistic not consistent with class type." );
            }
            else
            {
                // Create partial of first-order relativistic light-time correction
                partialList.push_back(
                            std::make_shared< FirstOrderRelativisticLightTimeCorrectionPartial >( currentCorrection ) );
            }

            break;
        }
        default:
            break;
        }
    }

    return partialList;
}

}

}

