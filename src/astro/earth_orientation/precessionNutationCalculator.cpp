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

#include <boost/bind/bind.hpp>
using namespace boost::placeholders;


#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/math/interpolators/cubicSplineInterpolator.h"
#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/astro/earth_orientation/precessionNutationCalculator.h"
#include "tudat/interface/sofa/sofaTimeConversions.h"

namespace tudat
{

namespace earth_orientation
{

//! Constructor for (CIO-based) precession-nutation calculation object.
PrecessionNutationCalculator::PrecessionNutationCalculator(
        basic_astrodynamics::IAUConventions precessionNutationTheory,
        std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector2d > >
        dailyCorrectionInterpolator ):
    precessionNutationTheory_( precessionNutationTheory ),
    dailyCorrectionInterpolator_( dailyCorrectionInterpolator )
{
    // Link selected SOFA function wrapper for direct calculation of precession-nutation.
    nominalCipPositionFunction_ =
            std::bind( sofa_interface::getPositionOfCipInGcrs,
                         std::placeholders::_1, basic_astrodynamics::JULIAN_DAY_ON_J2000, precessionNutationTheory );
}

//! Function to calculate the position of CIP in GCRS (CIO-based precession-nutation) and CIO-locator.
std::pair< Eigen::Vector2d, double > PrecessionNutationCalculator::getPositionOfCipInGcrs(
        const double terrestrialTime )
{
    // Calculate current UTC from SOFA.
    double utc = sofa_interface::convertTTtoUTC( terrestrialTime );

    // Call function to compte precession-nutation from UTC and TT.
    return getPositionOfCipInGcrs( terrestrialTime, utc );
}

//! Function to calculate the position of CIP in GCRS (CIO-based precession-nutation) and CIO-locator.
std::pair< Eigen::Vector2d, double > PrecessionNutationCalculator::getPositionOfCipInGcrs(
        const double terrestrialTime,
        const double utc )
{
    // Calculate nominal precession-nutation values.
    std::pair< Eigen::Vector2d, double > nominalCipPosition = nominalCipPositionFunction_( terrestrialTime );

    // Retrieve measured corrections to model.
    Eigen::Vector2d iersCorrections = dailyCorrectionInterpolator_->interpolate( utc );

    // Add nominal values and corrections and return.
    return std::pair< Eigen::Vector2d, double >(
                nominalCipPosition.first + iersCorrections, nominalCipPosition.second );
}

}

}
