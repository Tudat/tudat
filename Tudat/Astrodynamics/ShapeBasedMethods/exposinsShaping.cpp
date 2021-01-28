/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */


#include "Tudat/Astrodynamics/BasicAstrodynamics/celestialBodyConstants.h"
#include "Tudat/Basics/timeType.h"
#include "Tudat/Astrodynamics/Propulsion/thrustAccelerationModel.h"
#include "exposinsShaping.h"

namespace tudat
{
namespace shape_based_methods
{


//! Constructur for exposins shaping.
ExposinsShaping::ExposinsShaping( Eigen::Vector6d initialState,
                                    Eigen::Vector6d finalState,
                                    double requiredTimeOfFlight,
                                    int numberOfRevolutions,
                                    simulation_setup::NamedBodyMap& bodyMap,
                                    const std::string bodyToPropagate,
                                    const std::string centralBody,
                                    double windingParameter,
                                    std::shared_ptr< root_finders::RootFinderSettings >& rootFinderSettings,
                                    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings ):
    ShapeBasedMethodLeg( initialState, finalState, requiredTimeOfFlight, bodyMap, bodyToPropagate, centralBody, integratorSettings ),
    numberOfRevolutions_( numberOfRevolutions ),
    bodyMap_( bodyMap ), bodyToPropagate_( bodyToPropagate ), centralBody_( centralBody ),
    windingParameter_( windingParameter ),
    rootFinderSettings_( rootFinderSettings ),
    integratorSettings_( integratorSettings )
{
    // Retrieve gravitational parameter of the central body.
    double centralBodyGravitationalParameter = bodyMap_[ centralBody_ ]->getGravityFieldModel()->getGravitationalParameter( );

    // Normalize the initial state.
    initialState_.segment( 0, 3 ) = initialState.segment( 0, 3 ) / physical_constants::ASTRONOMICAL_UNIT;
    initialState_.segment( 3, 3 ) = initialState.segment( 3, 3 ) * physical_constants::JULIAN_YEAR / physical_constants::ASTRONOMICAL_UNIT;

    // Normalize the final state.
    finalState_.segment( 0, 3 ) = finalState.segment( 0, 3 ) / physical_constants::ASTRONOMICAL_UNIT;
    finalState_.segment( 3, 3 ) = finalState.segment( 3, 3 ) * physical_constants::JULIAN_YEAR / physical_constants::ASTRONOMICAL_UNIT;

    // Normalize the required time of flight.
    requiredTimeOfFlight_ = requiredTimeOfFlight / physical_constants::JULIAN_YEAR;

    // Normalize the gravitational parameter of the central body.
    centralBodyGravitationalParameter_ = centralBodyGravitationalParameter * std::pow( physical_constants::JULIAN_YEAR, 2.0 )
            / std::pow( physical_constants::ASTRONOMICAL_UNIT, 3.0 );


//    std::cout << centralBodyGravitationalParameter << std::endl;
    // Compute initial state in spherical coordinates.
    initialStateSphericalCoordinates_ = coordinate_conversions::convertCartesianToSphericalState( initialState_ );

    // Compute final state in spherical coordinates.
    finalStateSphericalCoordinates_ = coordinate_conversions::convertCartesianToSphericalState( finalState_ );

    if ( initialStateSphericalCoordinates_( 1 ) < 0.0 )
    {
//        std::cout << "initialStateSphericalCoordinates_ HIT" << initialStateSphericalCoordinates_( 1 ) << std::endl;
        initialStateSphericalCoordinates_( 1 ) += 2.0 * mathematical_constants::PI;
    }
    if ( finalStateSphericalCoordinates_( 1 ) < 0.0 )
    {
//        std::cout << "finalStateSphericalCoordinates_ HIT" << finalStateSphericalCoordinates_( 1 ) << std::endl;
        finalStateSphericalCoordinates_( 1 ) += 2.0 * mathematical_constants::PI;
    }

    // Retrieve the initial value of the azimuth angle.
    initialAzimuthAngle_ = initialStateSphericalCoordinates_[ 1 ];

    // Compute final value of the azimuth angle.
    if ( ( finalStateSphericalCoordinates_( 1 ) - initialStateSphericalCoordinates_( 1 ) ) < 0.0 )
    {
//        std::cout << "difference negative HIT" << finalStateSphericalCoordinates_( 1 ) - initialStateSphericalCoordinates_( 1 ) <<  std::endl;
        finalAzimuthAngle_ = finalStateSphericalCoordinates_( 1 ) + 2.0 * mathematical_constants::PI * ( numberOfRevolutions_ + 1.0 );
    }
    else
    {
//        std::cout << "difference positive HIT" << finalStateSphericalCoordinates_( 1 ) - initialStateSphericalCoordinates_( 1 ) <<  std::endl;
        finalAzimuthAngle_ = finalStateSphericalCoordinates_( 1 ) + 2.0 * mathematical_constants::PI * ( numberOfRevolutions_ );
    }


    // Compute initial and final values of the derivative of the azimuth angle w.r.t. time.
    double initialDerivativeAzimuthAngle = initialStateSphericalCoordinates_[ 4 ]
            / ( initialStateSphericalCoordinates_[ 0 ] * std::cos( initialStateSphericalCoordinates_[ 2 ] ) );
    double finalDerivativeAzimuthAngle = finalStateSphericalCoordinates_[ 4 ]
            / ( finalStateSphericalCoordinates_[ 0 ] * std::cos( finalStateSphericalCoordinates_[ 2 ] ) );

    // Compute initial state parametrized by azimuth angle theta.
    initialStateParametrizedByAzimuthAngle_ = ( Eigen::Vector6d() << initialStateSphericalCoordinates_[ 0 ],
            initialStateSphericalCoordinates_[ 1 ],
            initialStateSphericalCoordinates_[ 2 ],
            initialStateSphericalCoordinates_[ 3 ] / initialDerivativeAzimuthAngle,
            initialStateSphericalCoordinates_[ 4 ] / initialDerivativeAzimuthAngle,
            initialStateSphericalCoordinates_[ 5 ] / initialDerivativeAzimuthAngle ).finished();

    // Compute final state parametrized by azimuth angle theta.
    finalStateParametrizedByAzimuthAngle_ = ( Eigen::Vector6d() << finalStateSphericalCoordinates_[ 0 ],
            finalStateSphericalCoordinates_[ 1 ],
            finalStateSphericalCoordinates_[ 2 ],
            finalStateSphericalCoordinates_[ 3 ] / finalDerivativeAzimuthAngle,
            finalStateSphericalCoordinates_[ 4 ] / finalDerivativeAzimuthAngle,
            finalStateSphericalCoordinates_[ 5 ] / finalDerivativeAzimuthAngle ).finished();


    // Initialise coefficients for radial distance and elevation angle functions.
    coefficientsRadialDistanceFunction_.resize( 4 );
    for ( int i = 0 ; i < 4 ; i++ )
    {
        coefficientsRadialDistanceFunction_[ i ] = 1.0;
    }
//    // NOT USED
//    coefficientsElevationAngleFunction_.resize( 4 );
//    for ( int i = 0 ; i < 4 ; i++ )
//    {
//        coefficientsElevationAngleFunction_[ i ] = 1.0;
//    }


    // Define coefficients for radial distance and elevation angle composite functions.
    radialDistanceCompositeFunction_ = std::make_shared< CompositeRadialFunctionExposinsShaping >( coefficientsRadialDistanceFunction_ );


    // /////////////////////////////////////////////////

    // //// ACTUAL ANGLE - account for diff elevation
    // https://en.wikipedia.org/wiki/Great-circle_distance#Computational_formulas
    // output is planar

    //A Travelled azimuth angle
    initialPositionUnit_ = initialState_.segment( 0, 3 ).normalized();
    finalPositionUnit_   = finalState_.segment( 0, 3 ).normalized();
    axisOfRotation_      = initialPositionUnit_.cross(finalPositionUnit_).normalized(); // needs to be unit vector

    travelledAzimuthAngle_ = finalAzimuthAngle_ - initialAzimuthAngle_; // numberOfrevolutions already included in finalAzimuthAngle
    computeTravelledAzimuthAngle();

    initialCylindricalRadius_ = initialStateSphericalCoordinates_[0]; // in-plane cylyndrical radius
    finalCylindricalRadius_   = finalStateSphericalCoordinates_[0];

    // Compute bounds for TOF parameter
    lowerBoundGamma_ = std::atan(windingParameter_/2*(-log(initialCylindricalRadius_/finalCylindricalRadius_)*(1/std::tan((windingParameter_*travelledAzimuthAngle_)/2)) - std::sqrt((2*(1-std::cos(windingParameter_*travelledAzimuthAngle_)))/std::pow(windingParameter_,4) - std::pow(std::log(initialCylindricalRadius_/finalCylindricalRadius_),2))));
    upperBoundGamma_ = std::atan(windingParameter_/2*(-log(initialCylindricalRadius_/finalCylindricalRadius_)*(1/std::tan((windingParameter_*travelledAzimuthAngle_)/2)) + std::sqrt((2*(1-std::cos(windingParameter_*travelledAzimuthAngle_)))/std::pow(windingParameter_,4) - std::pow(std::log(initialCylindricalRadius_/finalCylindricalRadius_),2))));

    // 64 nodes - best accuracy - need at least 32 for smooth results
    quadratureSettings_ = std::make_shared< numerical_quadrature::GaussianQuadratureSettings < double > >( 0, 64 );

    // Compute TOF boundaries from TOF parameter
    requiredGamma_ = lowerBoundGamma_;
    double lowerBoundTOF = computeNormalizedTimeOfFlight();
    requiredGamma_ = upperBoundGamma_;
    double upperBoundTOF = computeNormalizedTimeOfFlight();

    // Initial guess TOF parameter
    requiredGamma_ = (lowerBoundGamma_ + upperBoundGamma_)/2;

    // Check if requested TOF is within TOF bounds
    if ((requiredTimeOfFlight_ - lowerBoundTOF) * (requiredTimeOfFlight_ - upperBoundTOF) <= 0)
    {
        // Find the ciefficients such that the TOF requirement is satisfied
        iterateToMatchRequiredTimeOfFlight( rootFinderSettings_, lowerBoundGamma_, upperBoundGamma_, requiredGamma_ );

        // Store final shape coefficients in correct arrays and variables and initialize them
        updateShapeCoefficients( );

//        std::cout << "requiredGamma_" << requiredGamma_*180/mathematical_constants::PI << std::endl;

//        // Retrieve initial step size.
//        double initialStepSize = integratorSettings->initialTimeStep_;

//        // Vector of azimuth angles at which the time should be computed.
//        Eigen::VectorXd azimuthAnglesToComputeAssociatedEpochs =
//                Eigen::VectorXd::LinSpaced( std::ceil( computeNormalizedTimeOfFlight() * physical_constants::JULIAN_YEAR / initialStepSize ),
//                                            0, travelledAzimuthAngle_ );

//        std::map< double, double > dataToInterpolate;
//        for ( int i = 0 ; i < azimuthAnglesToComputeAssociatedEpochs.size() ; i++ )
//        {
//            dataToInterpolate[ computeCurrentTimeFromAzimuthAngle( azimuthAnglesToComputeAssociatedEpochs[ i ] ) * physical_constants::JULIAN_YEAR ]
//                    = azimuthAnglesToComputeAssociatedEpochs[ i ];
//        }

//        // Create interpolator.
//        std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings =
//                std::make_shared< interpolators::LagrangeInterpolatorSettings >( 10 );

//        interpolator_ = interpolators::createOneDimensionalInterpolator( dataToInterpolate, interpolatorSettings );


        infeasibleTOF_ = false;

    }
    else
    {
        infeasibleTOF_ = true;
    }



}


//! Convert time to independent variable. - TO DO
double ExposinsShaping::convertTimeToIndependentVariable( const double time )
{
//    // Retrieve initial step size.
//    double initialStepSize = integratorSettings_->initialTimeStep_;

//    // Vector of azimuth angles at which the time should be computed.
//    Eigen::VectorXd azimuthAnglesToComputeAssociatedEpochs =
//            Eigen::VectorXd::LinSpaced( std::ceil( computeNormalizedTimeOfFlight( ) * physical_constants::JULIAN_YEAR / initialStepSize ),
//                                        initialAzimuthAngle_, finalAzimuthAngle_ );

//    std::map< double, double > dataToInterpolate;
//    for ( int i = 0 ; i < azimuthAnglesToComputeAssociatedEpochs.size() ; i++ )
//    {
//        dataToInterpolate[ computeCurrentTimeFromAzimuthAngle( azimuthAnglesToComputeAssociatedEpochs[ i ] ) * physical_constants::JULIAN_YEAR ]
//                = azimuthAnglesToComputeAssociatedEpochs[ i ];
//    }

//    // Create interpolator.
//    std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings =
//            std::make_shared< interpolators::LagrangeInterpolatorSettings >( 10 );

//    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > > interpolator =
//            interpolators::createOneDimensionalInterpolator( dataToInterpolate, interpolatorSettings );

//    double independentVariable = interpolator->interpolate( time );

    return  interpolator_->interpolate( time ); //independentVariable;
}

//! Convert independent variable to time. - Done
double ExposinsShaping::convertIndependentVariableToTime( const double independentVariable )
{
    // Define the derivative of the time function w.r.t the azimuth angle theta.
    std::function< double( const double ) > derivativeTimeFunction = [ = ] ( const double currentAzimuthAngle )
    {
        return computeBaseTimeOfFlightFunction( currentAzimuthAngle );
    };

    // Create numerical quadrature from quadrature settings.
    std::shared_ptr< numerical_quadrature::NumericalQuadrature< double, double > > quadrature =
            numerical_quadrature::createQuadrature( derivativeTimeFunction, quadratureSettings_, independentVariable );

    double time = quadrature->getQuadrature( );

    return time;
}

//! Compute current time from azimuth angle. - Done
double ExposinsShaping::computeCurrentTimeFromAzimuthAngle( const double currentAzimuthAngle )
{
    // Define the derivative of the time function w.r.t the azimuth angle theta.
    std::function< double( const double ) > derivativeTimeFunction = [ = ] ( const double currentAzimuthAngle )
    {
        return computeBaseTimeOfFlightFunction( currentAzimuthAngle );
    };

    // Create numerical quadrature from quadrature settings.
    std::shared_ptr< numerical_quadrature::NumericalQuadrature< double, double > > quadrature =
            numerical_quadrature::createQuadrature( derivativeTimeFunction, quadratureSettings_, currentAzimuthAngle );

    return quadrature->getQuadrature( );
}




//! Compute base TOF function from azimuth angle. - Done
double ExposinsShaping::computeBaseTimeOfFlightFunction( double currentAzimuthAngle )
{

    double windTravelProduct = windingParameter_*travelledAzimuthAngle_;
    double logRadiiFraction = std::log(initialCylindricalRadius_/finalCylindricalRadius_);

    double dynamicRangeValue = std::sqrt(std::pow((logRadiiFraction + (std::tan(requiredGamma_)/windingParameter_)*std::sin(windTravelProduct))/(1-std::cos(windTravelProduct)),2) + std::pow(std::tan(requiredGamma_),2)/std::pow(windingParameter_,2));
    double dynamicRangeSign  = (logRadiiFraction + (std::tan(requiredGamma_)/windingParameter_)*std::sin(windTravelProduct));
    double dynamicRangeParameter = std::copysign(dynamicRangeValue,dynamicRangeSign); // possible problems with zero

    double phaseAngle = std::acos(std::tan(requiredGamma_)/(windingParameter_*dynamicRangeParameter));

    double scalingFactor = initialCylindricalRadius_/(std::exp(dynamicRangeParameter*std::sin(phaseAngle))); // Unit is meter [m] if multiplied with physical_constants::ASTRONOMICAL_UNIT*


    double cosineValue = std::cos(windingParameter_*currentAzimuthAngle + phaseAngle);
    double sineValue   = std::sin(windingParameter_*currentAzimuthAngle + phaseAngle);
    double flightAngle = std::atan(dynamicRangeParameter*windingParameter_*cosineValue);
    double radiusValue = scalingFactor*std::exp(dynamicRangeParameter*sineValue);

    return std::sqrt(std::pow(radiusValue,3)*
                     (std::pow(std::tan(flightAngle),2) + dynamicRangeParameter*std::pow(windingParameter_,2)*sineValue + 1)
                     /centralBodyGravitationalParameter_);

}
//! Compute normalised TOF
double ExposinsShaping::computeNormalizedTimeOfFlight()
{
    // Define the derivative of the time function w.r.t the azimuth angle theta.
    std::function< double( const double ) > derivativeTimeFunction = [ = ] ( const double currentAzimuthAngle )
    {
        return computeBaseTimeOfFlightFunction( currentAzimuthAngle );
    };

    // Create numerical quadrature from quadrature settings.
    std::shared_ptr< numerical_quadrature::NumericalQuadrature< double, double > > quadrature =
            numerical_quadrature::createQuadrature( derivativeTimeFunction, quadratureSettings_, travelledAzimuthAngle_ );

    return quadrature->getQuadrature( );
}





//! Iterate to match the required time of flight
void ExposinsShaping::iterateToMatchRequiredTimeOfFlight( std::shared_ptr< root_finders::RootFinderSettings > rootFinderSettings,
                                                           const double lowerBound,
                                                           const double upperBound,
                                                           const double initialGuess )
{

    // Define the structure updating the time of flight from the free coefficient value, while still satisfying the boundary conditions.
    std::function< void ( const double ) > resetRequiredGammaFunction = std::bind( &ExposinsShaping::resetRequiredGamma, this, std::placeholders::_1 );
    std::function< double ( ) >  computeTOFfunction = std::bind( &ExposinsShaping::computeNormalizedTimeOfFlight, this );
    std::function< double ( ) > getRequiredTOFfunction = std::bind( &ExposinsShaping::getNormalizedRequiredTimeOfFlight, this );

    std::shared_ptr< basic_mathematics::Function< double, double > > timeOfFlightFunction =
            std::make_shared< ExposinsShaping::TimeOfFlightFunction >( resetRequiredGammaFunction, computeTOFfunction, getRequiredTOFfunction );

    // Create root finder from root finder settings.
    std::shared_ptr< root_finders::RootFinderCore< double > > rootFinder = root_finders::createRootFinder( rootFinderSettings, lowerBound, upperBound, initialGuess );

    // Iterate to find the free coefficient value that matches the required time of flight.
    double updatedRequiredGamma = rootFinder->execute( timeOfFlightFunction, initialGuess );

}

//! Update shape equation coefficients
void ExposinsShaping::updateShapeCoefficients()
{

    double windTravelProduct = windingParameter_*travelledAzimuthAngle_;
    double logRadiiFraction = std::log(initialCylindricalRadius_/finalCylindricalRadius_);

    double dynamicRangeValue = std::sqrt(std::pow((logRadiiFraction + (std::tan(requiredGamma_)/windingParameter_)*std::sin(windTravelProduct))/(1-std::cos(windTravelProduct)),2) + std::pow(std::tan(requiredGamma_),2)/std::pow(windingParameter_,2));
    double dynamicRangeSign  = (logRadiiFraction + (std::tan(requiredGamma_)/windingParameter_)*std::sin(windTravelProduct));


    dynamicRangeParameter_ = std::copysign(dynamicRangeValue,dynamicRangeSign); // possible problems with zero

    phaseAngle_ = std::acos(std::tan(requiredGamma_)/(windingParameter_*dynamicRangeParameter_));

    scalingFactor_ = initialCylindricalRadius_/(std::exp(dynamicRangeParameter_*std::sin(phaseAngle_))); // Unit is meter [m] if multiplied with physical_constants::ASTRONOMICAL_UNIT*

    coefficientsRadialDistanceFunction_[0] = scalingFactor_;
    coefficientsRadialDistanceFunction_[1] = dynamicRangeParameter_;
    coefficientsRadialDistanceFunction_[2] = windingParameter_;
    coefficientsRadialDistanceFunction_[3] = phaseAngle_;

    radialDistanceCompositeFunction_ = std::make_shared< CompositeRadialFunctionExposinsShaping >( coefficientsRadialDistanceFunction_ );

}





//! Computed travelled angle theta_f
void ExposinsShaping::computeTravelledAzimuthAngle()
{

    double dotProductPositions            = initialPositionUnit_.dot(finalPositionUnit_);
    Eigen::Vector3d crossProductPositions = initialPositionUnit_.cross(finalPositionUnit_);

    travelledAzimuthAngle_ = std::atan2(crossProductPositions.norm(),dotProductPositions);

    Eigen::Vector3d planetVelocity =  initialState_.segment(3,3);

    Eigen::Vector3d planetRotationAxis = initialPositionUnit_.cross(planetVelocity.normalized());

    dotProductPositions   = axisOfRotation_.dot(planetRotationAxis);
    if (dotProductPositions < 0)
    {
    // issues
        // what if travel angle is negative?
        // zAxis is harcoded, change to be rotation of planet?

    axisOfRotation_ = -axisOfRotation_ ;
    travelledAzimuthAngle_ = 2.0 * mathematical_constants::PI  - travelledAzimuthAngle_;
    }

    travelledAzimuthAngle_ = travelledAzimuthAngle_ + 2.0*mathematical_constants::PI*numberOfRevolutions_;
}



//! Compute first derivative Azimuth Exposins
double  ExposinsShaping::computeFirstDerivativeAzimuthAngleWrtTime(const double currentAzimuthAngle )
{


    double windTravelProduct = windingParameter_*travelledAzimuthAngle_;
    double logRadiiFraction = std::log(initialCylindricalRadius_/finalCylindricalRadius_);

    double dynamicRangeValue = std::sqrt(std::pow((logRadiiFraction + (std::tan(requiredGamma_)/windingParameter_)*std::sin(windTravelProduct))/(1-std::cos(windTravelProduct)),2) + std::pow(std::tan(requiredGamma_),2)/std::pow(windingParameter_,2));
    double dynamicRangeSign  = (logRadiiFraction + (std::tan(requiredGamma_)/windingParameter_)*std::sin(windTravelProduct));
    double dynamicRangeParameter = std::copysign(dynamicRangeValue,dynamicRangeSign); // possible problems with zero

    double phaseAngle = std::acos(std::tan(requiredGamma_)/(windingParameter_*dynamicRangeParameter));

    double scalingFactor = initialCylindricalRadius_/(std::exp(dynamicRangeParameter*std::sin(phaseAngle))); // Unit is meter [m] if multiplied with physical_constants::ASTRONOMICAL_UNIT*

    double tanGamma = dynamicRangeParameter*windingParameter_*std::cos(windingParameter_*currentAzimuthAngle + phaseAngle);
    double k1k22s   = dynamicRangeParameter*std::pow(windingParameter_,2)*std::sin(windingParameter_*currentAzimuthAngle + phaseAngle);
    double currentRadius = scalingFactor*std::exp(dynamicRangeParameter*std::sin(windingParameter_*currentAzimuthAngle + phaseAngle)) * physical_constants::ASTRONOMICAL_UNIT;
    double centralBodyGravitationalParameter = centralBodyGravitationalParameter_ / std::pow( physical_constants::JULIAN_YEAR, 2.0 )
            * std::pow( physical_constants::ASTRONOMICAL_UNIT, 3.0 );

    // return dimensional value [rad/s]
    return std::sqrt((centralBodyGravitationalParameter/std::pow(currentRadius,3))*1/(std::pow(tanGamma,2) + k1k22s + 1));

}


//! Compute second derivative of the azimuth angle. - not used
double ExposinsShaping::computeSecondDerivativeAzimuthAngleWrtTime( const double currentAzimuthAngle )
{
//    // Compute first derivative azimuth angle w.r.t. time.
//    double firstDerivativeAzimuthAngle = computeFirstDerivativeAzimuthAngleWrtTime( currentAzimuthAngle );

//    // Compute scalar function of the time equation, and its derivative w.r.t. azimuth angle.
//    double scalarFunctionTimeEquation = computeScalarFunctionTimeEquation( currentAzimuthAngle );
//    double derivativeScalarFunctionTimeEquation = computeDerivativeScalarFunctionTimeEquation( currentAzimuthAngle );

//    // Compute radial distance, and its derivative w.r.t. azimuth angle.
//    double radialDistance = radialDistanceCompositeFunction_->evaluateCompositeFunction( currentAzimuthAngle );
//    double firstDerivativeRadialDistance = radialDistanceCompositeFunction_->evaluateCompositeFunctionFirstDerivative( currentAzimuthAngle );

//    return - std::pow( firstDerivativeAzimuthAngle, 2.0 )
//            * ( derivativeScalarFunctionTimeEquation / ( 2.0 * scalarFunctionTimeEquation ) + firstDerivativeRadialDistance / radialDistance );
    return 0;
}




//! Compute current cartesian state.
Eigen::Vector6d ExposinsShaping::computeNormalizedStateVector( const double currentAzimuthAngle )
{
    // Initialisation
    double currentNormalizedRadius      = radialDistanceCompositeFunction_->evaluateCompositeFunction( currentAzimuthAngle ); // need normalised radius

    // Position
    Eigen::Vector3d currentPositionUnit   = initialPositionUnit_*std::cos(currentAzimuthAngle) + axisOfRotation_.cross(initialPositionUnit_)*std::sin(currentAzimuthAngle) + axisOfRotation_*(axisOfRotation_.dot(initialPositionUnit_))*(1-std::cos(currentAzimuthAngle));
    Eigen::Vector3d currentPosition       = currentPositionUnit*currentNormalizedRadius;

    // Velocity
    double firstDerivativeAzimuthAngle = computeFirstDerivativeAzimuthAngleWrtTime( currentAzimuthAngle );
    double tanGamma                    = dynamicRangeParameter_*windingParameter_*std::cos(windingParameter_*currentAzimuthAngle + phaseAngle_);

    double radialVelocityNorm  = currentNormalizedRadius*firstDerivativeAzimuthAngle*tanGamma;
    double angularVelocityNorm = currentNormalizedRadius*firstDerivativeAzimuthAngle;

    Eigen::Vector3d currentRadialVelocity  = radialVelocityNorm*currentPositionUnit;
    Eigen::Vector3d currentAngularVelocity = angularVelocityNorm*axisOfRotation_.cross(currentPositionUnit);
    Eigen::Vector3d currentVelocity        = currentRadialVelocity+currentAngularVelocity;

    // State
    Eigen::Vector6d normalisedStateVector;
    normalisedStateVector.segment( 0, 3 ) = currentPosition;
    normalisedStateVector.segment( 3, 3 ) = currentVelocity;

    return  normalisedStateVector;
}

//! Compute current cartesian state.
Eigen::Vector6d ExposinsShaping::computeCurrentStateVector( const double currentAzimuthAngle )
{
    Eigen::Vector6d dimensionalStateVector;
    dimensionalStateVector.segment( 0, 3 ) = computeNormalizedStateVector( currentAzimuthAngle ).segment( 0, 3 )
            * physical_constants::ASTRONOMICAL_UNIT;
    dimensionalStateVector.segment( 3, 3 ) = computeNormalizedStateVector( currentAzimuthAngle ).segment( 3, 3 )
            * physical_constants::ASTRONOMICAL_UNIT; // no julian year needed as computed velocity is already per second, physical_constants::JULIAN_YEAR;

    return dimensionalStateVector;
}


//! Compute current thrust acceleration in cartesian coordinates. - TO DO
Eigen::Vector3d ExposinsShaping::computeNormalizedThrustAccelerationVector( const double currentAzimuthAngle )
{
    Eigen::Vector3d cartesianAcceleration;
    Eigen::Vector3d currentVelocity;
    // magnitude can still be negative
    double inPlaneThrustMagnitudeNorm = computeCurrentThrustAccelerationMagnitude(currentAzimuthAngle)/(physical_constants::ASTRONOMICAL_UNIT/ std::pow( physical_constants::JULIAN_YEAR, 2.0 ));

    currentVelocity = computeNormalizedStateVector( currentAzimuthAngle ).segment( 3, 3 );

    cartesianAcceleration = currentVelocity.normalized()*inPlaneThrustMagnitudeNorm;

    return cartesianAcceleration;
}

//! Compute magnitude cartesian acceleration.
double  ExposinsShaping::computeCurrentThrustAccelerationMagnitude( double currentAzimuthAngle )
{
    // not normalised
    double windTravelProduct = windingParameter_*travelledAzimuthAngle_;
    double logRadiiFraction = std::log(initialCylindricalRadius_/finalCylindricalRadius_);

    double dynamicRangeValue = std::sqrt(std::pow((logRadiiFraction + (std::tan(requiredGamma_)/windingParameter_)*std::sin(windTravelProduct))/(1-std::cos(windTravelProduct)),2) + std::pow(std::tan(requiredGamma_),2)/std::pow(windingParameter_,2));
    double dynamicRangeSign  = (logRadiiFraction + (std::tan(requiredGamma_)/windingParameter_)*std::sin(windTravelProduct));
    double dynamicRangeParameter = std::copysign(dynamicRangeValue,dynamicRangeSign); // possible problems with zero

    double phaseAngle = std::acos(std::tan(requiredGamma_)/(windingParameter_*dynamicRangeParameter));

    double scalingFactor = initialCylindricalRadius_/(std::exp(dynamicRangeParameter*std::sin(phaseAngle))); // Unit is meter [m] if multiplied with physical_constants::ASTRONOMICAL_UNIT*

    double flightAngle = std::atan(dynamicRangeParameter*windingParameter_* std::cos(windingParameter_*currentAzimuthAngle+phaseAngle)); // diff cosine vamue

    double tanGamma = dynamicRangeParameter*windingParameter_*std::cos(windingParameter_*currentAzimuthAngle + phaseAngle);
    double k1k22s   = dynamicRangeParameter*std::pow(windingParameter_,2)*std::sin(windingParameter_*currentAzimuthAngle + phaseAngle);

    double currentRadius = scalingFactor*std::exp(dynamicRangeParameter*std::sin(windingParameter_*currentAzimuthAngle + phaseAngle)) * physical_constants::ASTRONOMICAL_UNIT;
    double centralBodyGravitationalParameter = centralBodyGravitationalParameter_ / std::pow( physical_constants::JULIAN_YEAR, 2.0 )
            * std::pow( physical_constants::ASTRONOMICAL_UNIT, 3.0 );

    double C1 = tanGamma/(2*std::cos(flightAngle));
    double C2 = 1;
    double C3 = std::pow(tanGamma,2) + k1k22s + 1;
    double C4 = std::pow(windingParameter_,2)*(1-2*dynamicRangeParameter*std::sin(windingParameter_*currentAzimuthAngle + phaseAngle));
    double C5 = std::pow(C3,2);
    double C6 = centralBodyGravitationalParameter/std::pow(currentRadius,2);

    // return dimensional value [m/s^2]
    return (C1*(C2/C3 - C4/C5))*C6;
}

//! Compute direction thrust acceleration in cartesian coordinates.
Eigen::Vector3d ExposinsShaping::computeCurrentThrustAccelerationDirection( double currentAzimuthAngle )
{
    return computeNormalizedThrustAccelerationVector( currentAzimuthAngle ).normalized();
}

//! Compute Normalized thrust acceleration magnitude - Exposins
double  ExposinsShaping::computeCurrentNormalizedThrustAccelerationMagnitude( double currentAzimuthAngle )
{
    return computeCurrentThrustAccelerationMagnitude( currentAzimuthAngle ) / (physical_constants::ASTRONOMICAL_UNIT / std::pow( physical_constants::JULIAN_YEAR, 2.0 ));
}






//! Compute velocity difference.
double ExposinsShaping::computeDeltaVBoundaries()
{

    Eigen::Vector3d initialVelocity = initialState_.segment( 3, 3 )/ physical_constants::JULIAN_YEAR * physical_constants::ASTRONOMICAL_UNIT;
    Eigen::Vector3d finalVelocity   = finalState_.segment( 3, 3 )/ physical_constants::JULIAN_YEAR * physical_constants::ASTRONOMICAL_UNIT;


    Eigen::Vector3d deltaV1 = initialVelocity - computeCurrentStateVector(0.0).segment(3,3);
    Eigen::Vector3d deltaV2 = finalVelocity   - computeCurrentStateVector(travelledAzimuthAngle_).segment(3,3);

    return deltaV1.norm() + deltaV2.norm();

}

//! Compute final deltaV.
double ExposinsShaping::computeDeltaV( )
{
    // Define the derivative of the deltaV, ie thrust acceleration function, as a function of the azimuth angle.
    std::function< double( const double ) > derivativeFunctionDeltaV = [ = ] ( const double currentAzimuthAngle ){


        return std::abs(computeCurrentThrustAccelerationMagnitude( currentAzimuthAngle ))/
                computeFirstDerivativeAzimuthAngleWrtTime( currentAzimuthAngle );


    };

    // Define numerical quadrature from quadratrure settings.
    std::shared_ptr< numerical_quadrature::NumericalQuadrature< double, double > > quadrature =
            numerical_quadrature::createQuadrature( derivativeFunctionDeltaV, quadratureSettings_, travelledAzimuthAngle_ );

    // Return dimensional deltaV [m/s]
    return quadrature->getQuadrature( );
}

//! Compute final deltaV - Hohmann.
double ExposinsShaping::computeDeltaVHohmann( double centralBodyGravitationalParameter )
{
    double deltaV1 = std::sqrt(centralBodyGravitationalParameter/(initialCylindricalRadius_*physical_constants::ASTRONOMICAL_UNIT))*(std::sqrt((2*finalCylindricalRadius_)/(initialCylindricalRadius_+finalCylindricalRadius_)) - 1);
    double deltaV2 = std::sqrt(centralBodyGravitationalParameter/(finalCylindricalRadius_*physical_constants::ASTRONOMICAL_UNIT))*(1 - std::sqrt((2*initialCylindricalRadius_)/(initialCylindricalRadius_+finalCylindricalRadius_)));

    return deltaV1 + deltaV2;

}




////! Get low-thrust acceleration model from shaping method.
//std::shared_ptr< propulsion::ThrustAcceleration > ExposinsShaping::getLowThrustAccelerationModel(
////        simulation_setup::NamedBodyMap& bodyMap,
////        const std::string& bodyToPropagate,
//        std::function< double( const double ) > specificImpulseFunction,
//        std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > > interpolatorPolarAngleFromTime )
//{

//    std::shared_ptr< simulation_setup::Body > vehicle = bodyMap_[ bodyToPropagate_ ];

//    // Define thrust magnitude function from the shaped trajectory.
//    std::function< double( const double ) > thrustMagnitudeFunction = [ = ]( const double currentTime )
//    {

//        // Compute current azimuth angle.
//        double currentAzimuthAngle = interpolatorPolarAngleFromTime->interpolate( currentTime );

//        // Compute current acceleration.
//        double currentAcceleration = computeCurrentThrustAccelerationMagnitude( currentAzimuthAngle ) * physical_constants::ASTRONOMICAL_UNIT
//                / std::pow( physical_constants::JULIAN_YEAR, 2.0 );

//        // Compute current mass of the vehicle.
//        double currentMass = vehicle->getBodyMass();

//        // Compute and return magnitude of the low-thrust force.
//        return currentAcceleration * currentMass;
//    };

//    // Define thrust magnitude settings from thrust magnitude function.
//    std::shared_ptr< simulation_setup::FromFunctionThrustMagnitudeSettings > thrustMagnitudeSettings
//            = std::make_shared< simulation_setup::FromFunctionThrustMagnitudeSettings >(
//                thrustMagnitudeFunction, specificImpulseFunction );


//    // Define thrust direction function from the shaped trajectory.
//    std::function< Eigen::Vector3d( const double ) > thrustDirectionFunction = [ = ]( const double currentTime )
//    {
//        // Compute current azimuth angle.
//        double currentAzimuthAngle = interpolatorPolarAngleFromTime->interpolate( currentTime );

//        // Compute current direction of the acceleration vector.
//        Eigen::Vector3d currentAccelerationDirection = computeCurrentThrustAccelerationDirection( currentAzimuthAngle );

//        // Return direction of the low-thrust acceleration.
//        return currentAccelerationDirection;
//    };

//    // Define thrust direction settings from the direction of thrust acceleration retrieved from the shaping method.
//    std::shared_ptr< simulation_setup::CustomThrustDirectionSettings > thrustDirectionSettings =
//            std::make_shared< simulation_setup::CustomThrustDirectionSettings >( thrustDirectionFunction );

//    // Define thrust acceleration settings.
//    std::shared_ptr< simulation_setup::ThrustAccelerationSettings > thrustAccelerationSettings =
//            std::make_shared< simulation_setup::ThrustAccelerationSettings >(
//                thrustDirectionSettings, thrustMagnitudeSettings );

//    // Create low thrust acceleration model.
//    std::shared_ptr< propulsion::ThrustAcceleration > lowThrustAccelerationModel = createThrustAcceleratioModel(
//                thrustAccelerationSettings, bodyMap_, bodyToPropagate_ );

//    return lowThrustAccelerationModel;
//}


void ExposinsShaping::computeShapedTrajectoryAndFullPropagation(
//        simulation_setup::NamedBodyMap& bodyMap,
        std::function< double( const double ) > specificImpulseFunction,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
                std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > >& propagatorSettings,
        std::map< double, Eigen::VectorXd >& fullPropagationResults,
        std::map< double, Eigen::VectorXd >& shapingMethodResults,
        std::map< double, Eigen::VectorXd >& dependentVariablesHistory,
        const bool isMassPropagated ){

    fullPropagationResults.clear();
    shapingMethodResults.clear();
    dependentVariablesHistory.clear();

    std::string bodyToPropagate = propagatorSettings.first->bodiesToIntegrate_[ 0 ];

    // Retrieve initial step size.
    double initialStepSize = integratorSettings->initialTimeStep_;

    // Vector of azimuth angles at which the time should be computed.
    Eigen::VectorXd azimuthAnglesToComputeAssociatedEpochs =
            Eigen::VectorXd::LinSpaced( std::ceil( computeNormalizedTimeOfFlight() * physical_constants::JULIAN_YEAR / initialStepSize ),
                                        initialAzimuthAngle_, finalAzimuthAngle_ );

    std::map< double, double > dataToInterpolate;
    for ( int i = 0 ; i < azimuthAnglesToComputeAssociatedEpochs.size() ; i++ )
    {
        dataToInterpolate[ computeCurrentTimeFromAzimuthAngle( azimuthAnglesToComputeAssociatedEpochs[ i ] ) * physical_constants::JULIAN_YEAR ]
                = azimuthAnglesToComputeAssociatedEpochs[ i ];
    }

    // Create interpolator.
    std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings =
            std::make_shared< interpolators::LagrangeInterpolatorSettings >( 10 );

    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > > interpolator =
            interpolators::createOneDimensionalInterpolator( dataToInterpolate, interpolatorSettings );

    // Compute halved time of flight.
    double halvedTimeOfFlight = computeNormalizedTimeOfFlight() / 2.0;

    // Compute azimuth angle at half of the time of flight.
    double azimuthAngleAtHalvedTimeOfFlight = interpolator->interpolate( halvedTimeOfFlight * physical_constants::JULIAN_YEAR );

    // Compute state at half of the time of flight.
    Eigen::Vector6d initialStateAtHalvedTimeOfFlight = computeNormalizedStateVector( azimuthAngleAtHalvedTimeOfFlight );
    initialStateAtHalvedTimeOfFlight.segment( 0, 3 ) *= physical_constants::ASTRONOMICAL_UNIT;
    initialStateAtHalvedTimeOfFlight.segment( 3, 3 ) *= physical_constants::ASTRONOMICAL_UNIT / physical_constants::JULIAN_YEAR;

    // Create low thrust acceleration model.
    std::shared_ptr< propulsion::ThrustAcceleration > lowThrustAccelerationModel =
            getLowThrustAccelerationModel( /*bodyMap, propagatorSettings.first->bodiesToIntegrate_[ 0 ],*/ specificImpulseFunction/*,
                                           std::bind( &ExposinsShaping::computeCurrentThrustAccelerationDirection, this, std::placeholders::_1 ),
                                           std::bind( &ExposinsShaping::computeCurrentThrustAccelerationMagnitude, this, std::placeholders::_1 )*/ /*, interpolator*/ );

    basic_astrodynamics::AccelerationMap accelerationMap = propagators::getAccelerationMapFromPropagatorSettings(
                std::dynamic_pointer_cast< propagators::SingleArcPropagatorSettings< double > >( propagatorSettings.first ) );

    accelerationMap[ propagatorSettings.first->bodiesToIntegrate_[ 0 ] ][ propagatorSettings.first->bodiesToIntegrate_[ 0 ] ]
            .push_back( lowThrustAccelerationModel );


    // Create complete propagation settings (backward and forward propagations).
    std::pair< std::shared_ptr< propagators::PropagatorSettings< double > >,
            std::shared_ptr< propagators::PropagatorSettings< double > > > completePropagatorSettings;


    // Define translational state propagation settings
    std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
            std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > > translationalStatePropagatorSettings;

    // Define backward translational state propagation settings.
    translationalStatePropagatorSettings.first = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
                        ( propagatorSettings.first->centralBodies_, accelerationMap, propagatorSettings.first->bodiesToIntegrate_,
                          initialStateAtHalvedTimeOfFlight, propagatorSettings.first->getTerminationSettings(),
                          propagatorSettings.first->propagator_, propagatorSettings.first->getDependentVariablesToSave() );

    // Define forward translational state propagation settings.
    translationalStatePropagatorSettings.second = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
                        ( propagatorSettings.second->centralBodies_, accelerationMap, propagatorSettings.second->bodiesToIntegrate_,
                          initialStateAtHalvedTimeOfFlight, propagatorSettings.second->getTerminationSettings(),
                          propagatorSettings.second->propagator_, propagatorSettings.second->getDependentVariablesToSave() );


    // If translational state and mass are propagated concurrently.
    if ( isMassPropagated )
    {
        // Create mass rate models
        std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModels;
        massRateModels[ bodyToPropagate ] = createMassRateModel( bodyToPropagate, std::make_shared< simulation_setup::FromThrustMassModelSettings >( 1 ),
                                                           bodyMap_, accelerationMap );


        // Propagate mass until half of the time of flight.
        std::shared_ptr< propagators::PropagatorSettings< double > > massPropagatorSettingsToHalvedTimeOfFlight =
                std::make_shared< propagators::MassPropagatorSettings< double > >( std::vector< std::string >{ bodyToPropagate }, massRateModels,
                    ( Eigen::Vector1d() << bodyMap_[ bodyToPropagate ]->getBodyMass() ).finished(),
                    std::make_shared< propagators::PropagationTimeTerminationSettings >( halvedTimeOfFlight * physical_constants::JULIAN_YEAR, true ) );

        integratorSettings->initialTime_ = 0.0;

        // Create dynamics simulation object.
        propagators::SingleArcDynamicsSimulator< double, double > dynamicsSimulator(
                    bodyMap_, integratorSettings, massPropagatorSettingsToHalvedTimeOfFlight, true, false, false );

        // Propagate spacecraft mass until half of the time of flight.
        std::map< double, Eigen::VectorXd > propagatedMass = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
        double massAtHalvedTimeOfFlight = propagatedMass.rbegin()->second[ 0 ];

        // Create settings for propagating the mass of the vehicle.
        std::pair< std::shared_ptr< propagators::MassPropagatorSettings< double > >,
                std::shared_ptr< propagators::MassPropagatorSettings< double > > > massPropagatorSettings;

        // Define backward mass propagation settings.
        massPropagatorSettings.first = std::make_shared< propagators::MassPropagatorSettings< double > >(
                    std::vector< std::string >{ bodyToPropagate }, massRateModels, ( Eigen::Matrix< double, 1, 1 >( ) << massAtHalvedTimeOfFlight ).finished( ),
                                                                      propagatorSettings.first->getTerminationSettings() );

        // Define forward mass propagation settings.
        massPropagatorSettings.second = std::make_shared< propagators::MassPropagatorSettings< double > >(
                    std::vector< std::string >{ bodyToPropagate }, massRateModels, ( Eigen::Matrix< double, 1, 1 >( ) << massAtHalvedTimeOfFlight ).finished( ),
                                                                      propagatorSettings.second->getTerminationSettings() );


        // Create list of propagation settings.
        std::pair< std::vector< std::shared_ptr< propagators::SingleArcPropagatorSettings< double > > >,
                std::vector< std::shared_ptr< propagators::SingleArcPropagatorSettings< double > > > > propagatorSettingsVector;

        // Backward propagator settings vector.
        propagatorSettingsVector.first.push_back( translationalStatePropagatorSettings.first );
        propagatorSettingsVector.first.push_back( massPropagatorSettings.first );

        // Forward propagator settings vector.
        propagatorSettingsVector.second.push_back( translationalStatePropagatorSettings.second );
        propagatorSettingsVector.second.push_back( massPropagatorSettings.second );


        // Backward hybrid propagation settings.
        completePropagatorSettings.first = std::make_shared< propagators::MultiTypePropagatorSettings< double > >( propagatorSettingsVector.first,
                    propagatorSettings.first->getTerminationSettings(), propagatorSettings.first->getDependentVariablesToSave() );

        // Forward hybrid propagation settings.
        completePropagatorSettings.second = std::make_shared< propagators::MultiTypePropagatorSettings< double > >( propagatorSettingsVector.second,
                    propagatorSettings.second->getTerminationSettings(), propagatorSettings.second->getDependentVariablesToSave() );


    }

    // If only translational state is propagated.
    else
    {
        // Backward hybrid propagation settings.
        completePropagatorSettings.first = translationalStatePropagatorSettings.first;

        // Forward hybrid propagation settings.
        completePropagatorSettings.second =  translationalStatePropagatorSettings.second;
    }


    // Define forward propagator settings variables.
    integratorSettings->initialTime_ = halvedTimeOfFlight * physical_constants::JULIAN_YEAR;

//    // Define propagation settings
//    std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > propagatorSettingsForwardPropagation;
//    std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > propagatorSettingsBackwardPropagation;

//    // Define forward propagation settings.
//    propagatorSettingsForwardPropagation = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
//                        ( propagatorSettings.second->centralBodies_ /*centralBodies*/, accelerationMap, propagatorSettings.second->bodiesToIntegrate_ /*bodiesToPropagate*/,
//                          initialStateAtHalvedTimeOfFlight, propagatorSettings.second->getTerminationSettings() /*terminationSettings.second*/,
//                          propagatorSettings.second->propagator_ /*propagatorType*/, propagatorSettings.second->getDependentVariablesToSave() /*dependentVariablesToSave*/ );

//    // Define backward propagation settings.
//    propagatorSettingsBackwardPropagation = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
//                        ( propagatorSettings.first->centralBodies_ /*centralBodies*/, accelerationMap, propagatorSettings.first->bodiesToIntegrate_ /*bodiesToPropagate*/,
//                          initialStateAtHalvedTimeOfFlight, propagatorSettings.first->getTerminationSettings() /*terminationSettings.first*/,
//                          propagatorSettings.first->propagator_ /*propagatorType*/, propagatorSettings.first->getDependentVariablesToSave() /*dependentVariablesToSave*/ );

    // Perform forward propagation.
    propagators::SingleArcDynamicsSimulator< > dynamicsSimulatorIntegrationForwards( bodyMap_, integratorSettings, completePropagatorSettings.second );
    std::map< double, Eigen::VectorXd > stateHistoryFullProblemForwardPropagation = dynamicsSimulatorIntegrationForwards.getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > dependentVariableHistoryForwardPropagation = dynamicsSimulatorIntegrationForwards.getDependentVariableHistory( );

    // Compute and save full propagation and shaping method results along the forward propagation direction.
    for( std::map< double, Eigen::VectorXd >::iterator itr = stateHistoryFullProblemForwardPropagation.begin( );
         itr != stateHistoryFullProblemForwardPropagation.end( ); itr++ )
    {
        double currentAzimuthAngle = interpolator->interpolate( itr->first );

        Eigen::Vector6d currentNormalisedState = computeNormalizedStateVector( currentAzimuthAngle );
        Eigen::Vector6d currentState;
        currentState.segment( 0, 3 ) = currentNormalisedState.segment( 0, 3 ) * physical_constants::ASTRONOMICAL_UNIT;
        currentState.segment( 3, 3 ) = currentNormalisedState.segment( 3, 3 ) * physical_constants::ASTRONOMICAL_UNIT / physical_constants::JULIAN_YEAR;
        shapingMethodResults[ itr->first ] = currentState;
        fullPropagationResults[ itr->first ] = itr->second;
        dependentVariablesHistory[ itr->first ] = dependentVariableHistoryForwardPropagation[ itr->first ];
    }


    // Define backward propagator settings variables.
    integratorSettings->initialTimeStep_ = - integratorSettings->initialTimeStep_;
    integratorSettings->initialTime_ = halvedTimeOfFlight * physical_constants::JULIAN_YEAR;

    // Perform the backward propagation.
    propagators::SingleArcDynamicsSimulator< > dynamicsSimulatorIntegrationBackwards( bodyMap_, integratorSettings, completePropagatorSettings.first );
    std::map< double, Eigen::VectorXd > stateHistoryFullProblemBackwardPropagation = dynamicsSimulatorIntegrationBackwards.getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > dependentVariableHistoryBackwardPropagation = dynamicsSimulatorIntegrationBackwards.getDependentVariableHistory( );

    // Compute and save full propagation and shaping method results along the backward propagation direction
    for( std::map< double, Eigen::VectorXd >::iterator itr = stateHistoryFullProblemBackwardPropagation.begin( );
         itr != stateHistoryFullProblemBackwardPropagation.end( ); itr++ )
    {
        double currentAzimuthAngle = interpolator->interpolate( itr->first );

        Eigen::Vector6d currentNormalisedState = computeNormalizedStateVector( currentAzimuthAngle );
        Eigen::Vector6d currentState;
        currentState.segment( 0, 3 ) = currentNormalisedState.segment( 0, 3 ) * physical_constants::ASTRONOMICAL_UNIT;
        currentState.segment( 3, 3 ) = currentNormalisedState.segment( 3, 3 ) * physical_constants::ASTRONOMICAL_UNIT / physical_constants::JULIAN_YEAR;
        shapingMethodResults[ itr->first ] = currentState;
        fullPropagationResults[ itr->first ] = itr->second;
        dependentVariablesHistory[ itr->first ] = dependentVariableHistoryBackwardPropagation[ itr->first ];
    }

    // Reset initial integrator settings.
    integratorSettings->initialTimeStep_ = - integratorSettings->initialTimeStep_;

}


} // namespace shape_based_methods
} // namespace tudat
