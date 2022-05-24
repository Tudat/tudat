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

#ifndef TUDAT_SHAPE_BASED_METHOD_H
#define TUDAT_SHAPE_BASED_METHOD_H

#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include <map>
#include "tudat/astro/low_thrust/lowThrustLeg.h"
#include "tudat/math/quadrature/createNumericalQuadrature.h"

using namespace tudat::low_thrust_trajectories;

namespace tudat
{
namespace shape_based_methods
{


// Base class for shape based methods.
class ShapeBasedMethod : public LowThrustLeg
{
public:

    //! Empty constructor.
    ShapeBasedMethod( const Eigen::Vector6d& stateAtDeparture,
                      const Eigen::Vector6d& stateAtArrival,
                      const double timeOfFlight ) :
    LowThrustLeg( stateAtDeparture, stateAtArrival, timeOfFlight, false ){ }

    //! Default destructor.
    virtual ~ShapeBasedMethod( ) { }

    //! Returns initial value of the independent variable.
    virtual double getInitialValueInpendentVariable( ) = 0;

    //! Returns final value of the independent variable.
    virtual double getFinalValueInpendentVariable( ) = 0;

    //! Returns state history.
    void getTrajectory(
            std::vector< double >& epochsVector,
            std::map< double, Eigen::Vector6d >& propagatedTrajectory );



protected:

    //! Numerical quadrature settings, required to compute the time of flight and total deltaV.
    std::shared_ptr< numerical_quadrature::QuadratureSettings< double > > quadratureSettings_;


};


} // namespace shape_based_methods
} // namespace tudat

#endif // TUDAT_SHAPE_BASED_METHOD_H
