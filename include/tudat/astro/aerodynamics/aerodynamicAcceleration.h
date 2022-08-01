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

#ifndef TUDAT_AERODYNAMIC_ACCELERATION_H
#define TUDAT_AERODYNAMIC_ACCELERATION_H

#include <functional>
#include <memory>

#include <Eigen/Core>

#include "tudat/astro/aerodynamics/aerodynamicCoefficientInterface.h"
#include "tudat/astro/aerodynamics/aerodynamicForce.h"
#include "tudat/astro/basic_astro/accelerationModel.h"

namespace tudat
{
namespace aerodynamics
{

//! Compute the aerodynamic acceleration in same reference frame as input coefficients.
/*!
 * This function computes the aerodynamic acceleration. It takes primitive types as arguments to
 * perform the calculations. Therefore, these quantities (dynamic pressure, reference area and
 * aerodynamic coefficients) have to computed before passing them to this function.
 * \param dynamicPressure Dynamic pressure at which the body undergoing the acceleration flies.
 * \param referenceArea Reference area of the aerodynamic coefficients.
 * \param aerodynamicCoefficients Aerodynamic coefficients in right-handed reference frame.
 * \param vehicleMass Mass of vehicle undergoing acceleration.
 * \return Resultant aerodynamic acceleration, given in reference frame in which the
 *          aerodynamic coefficients were given (assuming coefficients in positive direction).
 */
Eigen::Vector3d computeAerodynamicAcceleration( const double dynamicPressure,
                                                const double referenceArea,
                                                const Eigen::Vector3d& aerodynamicCoefficients,
                                                const double vehicleMass );

//! Compute the aerodynamic acceleration in same reference frame as input coefficients.
/*!
 * This function computes the aerodynamic acceleration. It takes the dynamic pressure and an
 * aerodynamic coefficient interface as input. The coefficient interface has to have been
 * updated with current vehicle conditions before being passed to this function. Aerodynamic
 * coefficients and reference area are then retrieved from it.
 * \param dynamicPressure Dynamic pressure at which the body undergoing the acceleration flies.
 * \param coefficientInterface AerodynamicCoefficientInterface class from which reference area
 *          and coefficients are retrieved.
 * \param vehicleMass Mass of vehicle undergoing acceleration.
 * \return Resultant aerodynamic acceleration, given in reference frame in which the
 *          aerodynamic coefficients were given (assuming coefficients in positive direction).
 */
Eigen::Vector3d computeAerodynamicAcceleration(
        const double dynamicPressure,
        AerodynamicCoefficientInterfacePointer coefficientInterface,
        const double vehicleMass );

//! Class for calculation of aerodynamic accelerations.
/*!
 * Class for calculation of aerodynamic accelerations.
 * \sa AccelerationModel.
 */
class AerodynamicAcceleration : public basic_astrodynamics::AccelerationModel< Eigen::Vector3d >
{
private:

    //! Typedef for double-returning function.
    typedef std::function< double ( ) > DoubleReturningFunction;

    //! Typedef for coefficient-returning function.
    typedef std::function< void( Eigen::Vector3d& ) > CoefficientReturningFunction;

public:

    //! Acceleration model constructor, taking constant values of mass and reference area.
    /*!
     * Acceleration model constructor, taking constant values of mass and reference area.
     * \param coefficientFunction Function which retrieves current values of aerodynamic
     *          coefficients.
     * \param densityFunction Function which retrieves current value of the density.
     * \param airSpeedFunction Function which retrieves current value of the airspeed.
     * \param constantMass Value of vehicle mass that is used for all calls of this class.
     * \param constantReferenceArea Value of aerodynamic coefficient reference area that is used
     *          for all calls of this class.
     * \param areCoefficientsInNegativeDirection Boolean that determines whether to invert
     *          direction of aerodynamic coefficients. This is typically done for lift, drag and
     *          side force coefficients that point in negative direction in the local frame
     *          (default true).
     */
    AerodynamicAcceleration( const CoefficientReturningFunction coefficientFunction,
                             const DoubleReturningFunction densityFunction,
                             const DoubleReturningFunction airSpeedFunction,
                             const double constantMass,
                             const double constantReferenceArea,
                             const bool areCoefficientsInNegativeDirection = true ):
        coefficientFunction_( coefficientFunction ),
        densityFunction_( densityFunction ),
        airSpeedFunction_( airSpeedFunction ),
        massFunction_( [ = ]( ){ return constantMass; } ),
        referenceAreaFunction_( [ = ]( ){ return constantReferenceArea; } )
    {
        coefficientMultiplier_ = areCoefficientsInNegativeDirection == true ? -1.0 : 1.0;
    }

    //! Acceleration model constructor.
    /*!
     * Acceleration model constructor, taking function pointers for all member variables.
     * \param coefficientFunction Function which retrieves current values of aerodynamic
     *          coefficients.
     * \param densityFunction Function which retrieves current value of the density.
     * \param airSpeedFunction Function which retrieves current value of the airspeed.
     * \param massFunction Function which retrieves current value of the vehicle mass.
     * \param referenceAreaFunction Function which retrieves current value of the aerodynamic
     *          coefficient reference area.
     * \param areCoefficientsInNegativeDirection Boolean that determines whether to invert
     *          direction of aerodynamic coefficients. This is typically done for lift, drag and
     *          side force coefficients that point in negative direction in the local frame
     *          (default true).
     */
    AerodynamicAcceleration( const CoefficientReturningFunction coefficientFunction,
                             const DoubleReturningFunction densityFunction,
                             const DoubleReturningFunction airSpeedFunction,
                             const DoubleReturningFunction massFunction,
                             const DoubleReturningFunction referenceAreaFunction,
                             const bool areCoefficientsInNegativeDirection = true ):
        coefficientFunction_( coefficientFunction ),
        densityFunction_( densityFunction ),
        airSpeedFunction_( airSpeedFunction ),
        massFunction_( massFunction ),
        referenceAreaFunction_( referenceAreaFunction )
    {
        coefficientMultiplier_ = areCoefficientsInNegativeDirection ? -1.0 : 1.0;
    }

    //! Destructor
    ~AerodynamicAcceleration( ){ }

    //! Update member variables used by the aerodynamic acceleration model.
    /*!
     * Updates member variables used by the aerodynamic acceleration model.
     * Function pointers to retrieve the current values of quantities from which the
     * acceleration is to be calculated are set by constructor. This function calls
     * them to update the associated variables to their current state.
     * \param currentTime Time at which acceleration model is to be updated.
     */
    void updateMembers( const double currentTime = TUDAT_NAN )
    {
        if( !( this->currentTime_ == currentTime ) )
        {
            this->coefficientFunction_( currentForceCoefficients_ );
            currentForceCoefficients_ *= coefficientMultiplier_;
            currentDensity_ = this->densityFunction_( );
            currentMass_ = this->massFunction_( );
            currentAirspeed_ = this->airSpeedFunction_( );
            currentReferenceArea_ = this->referenceAreaFunction_( );

            currentTime_ = currentTime;
            currentAcceleration_ = computeAerodynamicAcceleration(
                                0.5 * currentDensity_ * currentAirspeed_ * currentAirspeed_,
                                currentReferenceArea_, currentForceCoefficients_, currentMass_ );
        }
    }

    //! Function to return current mass of body undergoing acceleration
    /*!
     * Function to return current mass of body undergoing acceleration as set from massFunction_ by updateMembers
     * \return Current mass of body undergoing acceleration
     */
    double getCurrentMass( )
    {
        return currentMass_;
    }


private:

    //! Function to retrieve the current aerodynamic force coefficients.
    const CoefficientReturningFunction coefficientFunction_;

    //! Function to retrieve the current density.
    const DoubleReturningFunction densityFunction_;

    //! Function to retrieve the current airspeed.
    const DoubleReturningFunction airSpeedFunction_;

    //! Function to retrieve the current mass.
    const DoubleReturningFunction massFunction_;

    //! Function to retrieve the current reference area.
    const DoubleReturningFunction referenceAreaFunction_;

    //! Current aerodynamic force coefficients.
    Eigen::Vector3d currentForceCoefficients_;

    //! Current density.
    double currentDensity_;

    //! Current airspeed.
    double currentAirspeed_;

    //! Current mass as set by massFunction_.
    double currentMass_;

    //! Current reference area, as set by referenceAreaFunction_.
    double currentReferenceArea_;

    //! Multiplier to reverse direction of coefficients.
    double coefficientMultiplier_;
};

//! Typedef for shared-pointer to AerodynamicAcceleration object.
typedef std::shared_ptr< AerodynamicAcceleration > AerodynamicAccelerationPointer;

} // namespace aerodynamics
} // namespace tudat

#endif // TUDAT_AERODYNAMIC_ACCELERATION_H
