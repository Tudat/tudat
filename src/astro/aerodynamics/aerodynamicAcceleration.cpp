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

#include "tudat/astro/aerodynamics/aerodynamicAcceleration.h"

namespace tudat
{
namespace aerodynamics
{

//! Compute the aerodynamic acceleration in same reference frame as input coefficients.
Eigen::Vector3d computeAerodynamicAcceleration( const double dynamicPressure,
                                                const double referenceArea,
                                                const Eigen::Vector3d& aerodynamicCoefficients,
                                                const double vehicleMass )
{
    return computeAerodynamicForce( dynamicPressure, referenceArea, aerodynamicCoefficients )
            / vehicleMass;
}

//! Compute the aerodynamic acceleration in same reference frame as input coefficients.
Eigen::Vector3d computeAerodynamicAcceleration(
        const double dynamicPressure,
        AerodynamicCoefficientInterfacePointer coefficientInterface,
        const double vehicleMass )
{
    return computeAerodynamicForce( dynamicPressure, coefficientInterface ) / vehicleMass;
}


} // namespace aerodynamics
} // namespace tudat
