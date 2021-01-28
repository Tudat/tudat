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
#include "invPolyShaping.h"

namespace tudat
{
namespace shape_based_methods
{


//! Constructur for InvPoly shaping.
InvPolyShaping::InvPolyShaping( Eigen::Vector6d initialState,
                                    Eigen::Vector6d finalState,
                                    double requiredTimeOfFlight,
                                    int numberOfRevolutions,
                                    simulation_setup::NamedBodyMap& bodyMap,
                                    const std::string bodyToPropagate,
                                    const std::string centralBody,
                                    double initialValueFreeCoefficient,
                                    std::shared_ptr< root_finders::RootFinderSettings >& rootFinderSettings,
                                    const double lowerBoundFreeCoefficient,
                                    const double upperBoundFreeCoefficient,
                                    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings ):
    ShapeBasedMethodLeg( initialState, finalState, requiredTimeOfFlight, bodyMap, bodyToPropagate, centralBody, integratorSettings ),
    numberOfRevolutions_( numberOfRevolutions ),
    bodyMap_( bodyMap ), bodyToPropagate_( bodyToPropagate ), centralBody_( centralBody ),
    initialValueFreeCoefficient_( initialValueFreeCoefficient ),
    rootFinderSettings_( rootFinderSettings ),
    lowerBoundFreeCoefficient_( lowerBoundFreeCoefficient ),
    upperBoundFreeCoefficient_( upperBoundFreeCoefficient ),
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


    // Compute initial state in spherical coordinates.
    initialStateSphericalCoordinates_ = coordinate_conversions::convertCartesianToSphericalState( initialState_ );

    // Compute final state in spherical coordinates.
    finalStateSphericalCoordinates_ = coordinate_conversions::convertCartesianToSphericalState( finalState_ );

    // Compute initial state in cylindrical coordinates.
    initialStateCylindricalCoordinates_ = coordinate_conversions::convertCartesianToCylindricalState( initialState_ );

    // Compute final state in cylindrical coordinates.
    finalStateCylindricalCoordinates_ = coordinate_conversions::convertCartesianToCylindricalState( finalState_ );

    if ( initialStateSphericalCoordinates_( 1 ) < 0.0 )
    {
        initialStateSphericalCoordinates_( 1 ) += 2.0 * mathematical_constants::PI;
    }
    if ( finalStateSphericalCoordinates_( 1 ) < 0.0 )
    {
        finalStateSphericalCoordinates_( 1 ) += 2.0 * mathematical_constants::PI;
    }

    // Retrieve the initial value of the azimuth angle.
    initialAzimuthAngle_ = initialStateSphericalCoordinates_[ 1 ];

    // Compute final value of the azimuth angle.
    if ( ( finalStateSphericalCoordinates_( 1 ) - initialStateSphericalCoordinates_( 1 ) ) < 0.0 )
    {
        finalAzimuthAngle_ = finalStateSphericalCoordinates_( 1 ) + 2.0 * mathematical_constants::PI * ( numberOfRevolutions_ + 1.0 );
    }
    else
    {
        finalAzimuthAngle_ = finalStateSphericalCoordinates_( 1 ) + 2.0 * mathematical_constants::PI * ( numberOfRevolutions_ );
    }



    // Initialise coefficients for radial distance and elevation angle functions. -
    coefficientsRadialDistanceFunction_.resize( 7 );
    for ( int i = 0 ; i < 7 ; i++ )
    {
        coefficientsRadialDistanceFunction_[ i ] = 1.0;
    }
    coefficientsZFunction_.resize( 4 );
    for ( int i = 0 ; i < 4 ; i++ )
    {
        coefficientsZFunction_[ i ] = 1.0;
    }

    // Define coefficients for radial distance and elevation angle composite functions.
    radialDistanceCompositeFunction_ = std::make_shared< CompositeRadialFunctionInvPolyShaping >( coefficientsRadialDistanceFunction_ );
    zCompositeFunction_ = std::make_shared< CompositeZFunctionInvPolyShaping >( coefficientsZFunction_ );


    //A Travelled azimuth angle
    initialPositionUnit_   = initialState_.segment( 0, 3 ).normalized();
    finalPositionUnit_     = finalState_.segment( 0, 3 ).normalized();
    axisOfRotation_        = initialPositionUnit_.cross(finalPositionUnit_).normalized(); // needs to be unit vector
    travelledAzimuthAngle_ = finalAzimuthAngle_ - initialAzimuthAngle_;

    // Define settings for numerical quadrature, to be used to compute time of flight and final deltaV.  TO DO add correctnumber for quadrature + start angle
    quadratureSettings_ = std::make_shared< numerical_quadrature::GaussianQuadratureSettings < double > >( 0, 64 );

    // Compute bounds for d coefficient / TOF coefficient
    Eigen::Vector2d dCoefficientBounds = computeDCoefficientBounds();

    // Compute TOF for d parameter bounds
    coefficientsRadialDistanceFunction_[3] = dCoefficientBounds[0];
    satisfyBoundaryConditions();
    TOFdUpperLimit_ = computeNormalizedTimeOfFlight();


    coefficientsRadialDistanceFunction_[3] = dCoefficientBounds[0];
    satisfyBoundaryConditions();
    double TOFlower = computeNormalizedTimeOfFlight();
    coefficientsRadialDistanceFunction_[3] = dCoefficientBounds[1];
    satisfyBoundaryConditions();
    double TOFupper = computeNormalizedTimeOfFlight();

    // Check if requested TOF is feasible
    if ((requiredTimeOfFlight_ - TOFlower) * (requiredTimeOfFlight_ - TOFupper) <= 0)
    {
        // Iterate on the free coefficient value until the time of flight matches its required value.
        coefficientsRadialDistanceFunction_[3] = 0.0;
        iterateToMatchRequiredTimeOfFlight( rootFinderSettings_, dCoefficientBounds[0], dCoefficientBounds[1], 0.0);
        infeasibleTOF_ = false;

//        // Retrieve initial step size.
//        double initialStepSize = integratorSettings->initialTimeStep_;

//        // Vector of azimuth angles at which the time should be computed.
//        Eigen::VectorXd azimuthAnglesToComputeAssociatedEpochs =
//                Eigen::VectorXd::LinSpaced( std::ceil( computeNormalizedTimeOfFlight() * physical_constants::JULIAN_YEAR / initialStepSize ),
//                                            initialAzimuthAngle_, finalAzimuthAngle_ );

//        std::map< double, double > dataToInterpolate;
//        for ( int i = 0 ; i < azimuthAnglesToComputeAssociatedEpochs.size() ; i++ )
//        {
//            dataToInterpolate[ computeCurrentTimeFromAzimuthAngle( azimuthAnglesToComputeAssociatedEpochs[ i ] ) * physical_constants::JULIAN_YEAR ]
//                    = azimuthAnglesToComputeAssociatedEpochs[ i ];
//        }

//        // Create interpolator. - can give errors
//        std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings =
//                std::make_shared< interpolators::LagrangeInterpolatorSettings >( 10 );

//        interpolator_ = interpolators::createOneDimensionalInterpolator( dataToInterpolate, interpolatorSettings );

    }
    else
    {
        infeasibleTOF_ = true;
    }

}


//! Convert time to independent variable. - TO DO
double InvPolyShaping::convertTimeToIndependentVariable( const double time )
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

//! Convert independent variable to time. - TO DO
double InvPolyShaping::convertIndependentVariableToTime( const double independentVariable )
{
    // Define the derivative of the time function w.r.t the azimuth angle theta.
    std::function< double( const double ) > derivativeTimeFunction = [ = ] ( const double currentAzimuthAngle )
    {
        return computeBaseTimeOfFlightFunction( currentAzimuthAngle );
    };

    // Create numerical quadrature from quadrature settings.
    std::shared_ptr< numerical_quadrature::NumericalQuadrature< double, double > > quadrature =
            numerical_quadrature::createQuadrature( derivativeTimeFunction, quadratureSettings_, independentVariable );

    return quadrature->getQuadrature( );
}

////! Computed travelled angle theta_f
//void InvPolyShaping::computeTravelledAzimuthAngle()
//{

//    double dotProductPositions            = initialPositionUnit_.dot(finalPositionUnit_);
//    Eigen::Vector3d crossProductPositions = initialPositionUnit_.cross(finalPositionUnit_);

//    travelledAzimuthAngle_ = std::atan2(crossProductPositions.norm(),dotProductPositions);

//    Eigen::Vector3d planetVelocity =  initialState_.segment(3,3);

//    Eigen::Vector3d planetRotationAxis = initialPositionUnit_.cross(planetVelocity.normalized());

//    dotProductPositions   = axisOfRotation_.dot(planetRotationAxis);
//    if (dotProductPositions < 0)
//    {
//    // issues
//        // what if travel angle is negative?
//        // zAxis is harcoded, change to be rotation of planet?
//        // what if planet arrival and departura are retrograde/prograde

//    axisOfRotation_ = -axisOfRotation_ ;
//    travelledAzimuthAngle_ = 2.0 * mathematical_constants::PI  - travelledAzimuthAngle_;
////    std::cout << "angleCheck zAxis [deg] " << dotProductPositions << std::endl;
//    }
////    std::cout << "axisOfRotation_ " << axisOfRotation_ << std::endl;
////    std::cout << "travelledAzimuthAngle_ [deg] " << travelledAzimuthAngle_*180/mathematical_constants::PI << std::endl;
//    // is umber of revolutions done correctly?
//    travelledAzimuthAngle_ = travelledAzimuthAngle_ + 2.0*mathematical_constants::PI*numberOfRevolutions_;
////    std::cout << "travelledAzimuthAngle_ [deg] " << travelledAzimuthAngle_*180/mathematical_constants::PI << std::endl;


//}

//! Computed bounds for d coefficient
Eigen::Vector2d InvPolyShaping::computeDCoefficientBounds( )
{

    // Compute boundary values
    // 2D velocity
    Eigen::Vector3d initialPositionUnit = initialState_.segment( 0, 3 ).normalized();
    initialPositionUnit[2] = 0.0;
    initialPositionUnit = initialPositionUnit.normalized();

    Eigen::Vector3d initialVelocityUnit = initialState_.segment( 3, 3 ).normalized();
    initialVelocityUnit[2] = 0.0;
    initialVelocityUnit = initialVelocityUnit.normalized();

    Eigen::Vector3d finalPositionUnit   = finalState_.segment( 0, 3 ).normalized();
    finalPositionUnit[2] = 0.0;
    finalPositionUnit = finalPositionUnit.normalized();

    Eigen::Vector3d finalVelocityUnit   = finalState_.segment( 3, 3 ).normalized();
    finalVelocityUnit[2] = 0.0;
    finalVelocityUnit = finalVelocityUnit.normalized();

    // initial flight angle
    double initialDotProductPositions            = initialPositionUnit.dot(initialVelocityUnit);
    Eigen::Vector3d initialCrossProductPositions = initialPositionUnit.cross(initialVelocityUnit);
    double initialComplementFlightAngle          = std::atan2(initialCrossProductPositions.norm(),initialDotProductPositions);

    double initialFlightAngle     = mathematical_constants::PI/2 - initialComplementFlightAngle;

    // final flight angle
    double finalDotProductPositions            = finalPositionUnit.dot(finalVelocityUnit);
    Eigen::Vector3d finalCrossProductPositions = finalPositionUnit.cross(finalVelocityUnit);
    double finalComplementFlightAngle          = std::atan2(finalCrossProductPositions.norm(),finalDotProductPositions);

    double finalFlightAngle       = mathematical_constants::PI/2 - finalComplementFlightAngle;

    // radiis
    double initialRadius          = initialStateCylindricalCoordinates_[0];
    double finalRadius            = finalStateCylindricalCoordinates_[0];

    // Angular velocities
    double initialAngularVelocity = initialState_.segment( 3, 2 ).norm()* std::cos(initialFlightAngle)/initialRadius;
    double finalAngularVelocity   = finalState_.segment( 3, 2 ).norm()  * std::cos(finalFlightAngle)/finalRadius; // /physical_constants::JULIAN_YEAR

    double a = 1.0/initialRadius;
    double b = -std::tan(initialFlightAngle)/initialRadius;
    double c = (1.0/(2.0*initialRadius))*((centralBodyGravitationalParameter_)/(std::pow(initialRadius,3.0)*std::pow(initialAngularVelocity,2.0)) - 1);

    double theta_f1 = travelledAzimuthAngle_;
    double theta_f2 = std::pow(travelledAzimuthAngle_,2.0);
    double theta_f3 = std::pow(travelledAzimuthAngle_,3.0);
    double theta_f4 = std::pow(travelledAzimuthAngle_,4.0);


    Eigen::MatrixXd matrixThetaF = Eigen::MatrixXd::Zero( 3, 3 );
    matrixThetaF( 0, 0 ) =  30.0 * theta_f2;
    matrixThetaF( 0, 1 ) = -10.0 * theta_f3;
    matrixThetaF( 0, 2 ) =         theta_f4;
    matrixThetaF( 1, 0 ) = -48.0 * theta_f1;
    matrixThetaF( 1, 1 ) =  18.0 * theta_f2;
    matrixThetaF( 1, 2 ) = -2.0  * theta_f3;
    matrixThetaF( 2, 0 ) =  20.0;
    matrixThetaF( 2, 1 ) = -8.0  * theta_f1;
    matrixThetaF( 2, 2 ) =         theta_f2;

    Eigen::Vector3d vectorABC;
    vectorABC[0] = 1/finalRadius                                                                                      - (a + b*theta_f1 + c*theta_f2);
    vectorABC[1] = -std::tan(finalFlightAngle)/finalRadius                                                            - (b + 2*c*theta_f1);
    vectorABC[2] = centralBodyGravitationalParameter_/(std::pow(finalRadius, 4.0)*std::pow(finalAngularVelocity,2.0)) - (1/finalRadius + 2*c);

    Eigen::Vector3d vectorD;
    vectorD[0] =     -theta_f3;
    vectorD[1] = -3.0*theta_f2;
    vectorD[2] = -6.0*theta_f1;


    double currentAzimuthAngle;
    double stepsize = travelledAzimuthAngle_/500;
    double dLowerBound =  -100.0;
    double dUpperBound = 100.0;


    for ( currentAzimuthAngle = 0 ; currentAzimuthAngle <= travelledAzimuthAngle_ ; currentAzimuthAngle = currentAzimuthAngle + stepsize )
    {
        // PA
        double PA;
        PA = a + b*currentAzimuthAngle + c*(std::pow(currentAzimuthAngle,2.0) + 2);

        // PC
        Eigen::Vector3d vectorEFG;
        vectorEFG[0] = std::pow(currentAzimuthAngle, 4.0) + 12.0*std::pow(currentAzimuthAngle, 2.0);
        vectorEFG[1] = std::pow(currentAzimuthAngle, 5.0) + 20.0*std::pow(currentAzimuthAngle, 3.0);
        vectorEFG[2] = std::pow(currentAzimuthAngle, 6.0) + 30.0*std::pow(currentAzimuthAngle, 4.0);

        double PC;
        PC = (1.0/(2.0*std::pow(travelledAzimuthAngle_,6.0)))*vectorEFG.transpose()*matrixThetaF*vectorABC;

        // PD
        double PD;
        PD = (1.0/(2.0*std::pow(travelledAzimuthAngle_,6.0)))*vectorEFG.transpose()*matrixThetaF*vectorD;

        // PE
        double PE;
        PE = std::pow(currentAzimuthAngle, 3.0) + 6.0*currentAzimuthAngle + PD;

        // dcomparer
        double d = (-PA - PC)/PE;

        if ( PE < 0 )
        {
            if (d < dUpperBound)
            {
                dUpperBound = d - 0.001*std::abs(d);
            }
        }
        else
        {
            if (d > dLowerBound)
            {
                dLowerBound = d + 0.001*std::abs(d); // to ensure the boundary value itself does not lead to TOF = nan
            }
        }
    }

    return ( Eigen::Vector2d() <<
             dLowerBound,
             dUpperBound).finished();

}

//! Satisfy boundary conditions
void InvPolyShaping::satisfyBoundaryConditions( )
{

    // Compute boundary values
    // 2D velocity
    Eigen::Vector3d initialPositionUnit = initialState_.segment( 0, 3 ).normalized();
    initialPositionUnit[2] = 0.0;
    initialPositionUnit = initialPositionUnit.normalized();

    Eigen::Vector3d initialVelocityUnit = initialState_.segment( 3, 3 ).normalized();
    initialVelocityUnit[2] = 0.0;
    initialVelocityUnit = initialVelocityUnit.normalized();

    Eigen::Vector3d finalPositionUnit   = finalState_.segment( 0, 3 ).normalized();
    finalPositionUnit[2] = 0.0;
    finalPositionUnit = finalPositionUnit.normalized();

    Eigen::Vector3d finalVelocityUnit   = finalState_.segment( 3, 3 ).normalized();
    finalVelocityUnit[2] = 0.0;
    finalVelocityUnit = finalVelocityUnit.normalized();

    // initial flight angle
    double initialDotProductPositions            = initialPositionUnit.dot(initialVelocityUnit);
    Eigen::Vector3d initialCrossProductPositions = initialPositionUnit.cross(initialVelocityUnit);
    double initialComplementFlightAngle          = std::atan2(initialCrossProductPositions.norm(),initialDotProductPositions);

    double initialFlightAngle     = mathematical_constants::PI/2 - initialComplementFlightAngle;

    // final flight angle
    double finalDotProductPositions            = finalPositionUnit.dot(finalVelocityUnit);
    Eigen::Vector3d finalCrossProductPositions = finalPositionUnit.cross(finalVelocityUnit);
    double finalComplementFlightAngle          = std::atan2(finalCrossProductPositions.norm(),finalDotProductPositions);

    double finalFlightAngle       = mathematical_constants::PI/2 - finalComplementFlightAngle;

    // radiis
    double initialRadius          = initialStateCylindricalCoordinates_[0];
    double finalRadius            = finalStateCylindricalCoordinates_[0];

    // angular velocities
    double initialAngularVelocity = initialState_.segment( 3, 2 ).norm()* std::cos(initialFlightAngle)/initialRadius;
    double finalAngularVelocity   = finalState_.segment( 3, 2 ).norm()  * std::cos(finalFlightAngle)/finalRadius; // /physical_constants::JULIAN_YEAR

    double theta_f1 = travelledAzimuthAngle_;
    double theta_f2 = std::pow(travelledAzimuthAngle_,2.0);
    double theta_f3 = std::pow(travelledAzimuthAngle_,3.0);
    double theta_f4 = std::pow(travelledAzimuthAngle_,4.0);


    // Define first three components a,b,c // use of initialSphericalcoords - adjust for cylindrical coords
    // 2D case
    double a = 1.0/initialRadius;
    double b = -std::tan(initialFlightAngle)/initialRadius;
    double c = (1.0/(2.0*initialRadius))*((centralBodyGravitationalParameter_)/(std::pow(initialRadius,3.0)*std::pow(initialAngularVelocity,2.0)) - 1);
    double d = coefficientsRadialDistanceFunction_[3];// make sure other functions use it correctly

    Eigen::MatrixXd matrixA = Eigen::MatrixXd::Zero( 3, 3 );

    matrixA( 0, 0 ) =  30.0 * theta_f2;
    matrixA( 0, 1 ) = -10.0 * theta_f3;
    matrixA( 0, 2 ) =         theta_f4;
    matrixA( 1, 0 ) = -48.0 * theta_f1;
    matrixA( 1, 1 ) =  18.0 * theta_f2;
    matrixA( 1, 2 ) = -2.0  * theta_f3;
    matrixA( 2, 0 ) =  20.0;
    matrixA( 2, 1 ) = -8.0  * theta_f1;
    matrixA( 2, 2 ) =         theta_f2;

    Eigen::Vector3d vectorB;

    vectorB[0] = 1/finalRadius                                                                                      - (a + b*theta_f1 + c*theta_f2 + d*theta_f3);
    vectorB[1] = -std::tan(finalFlightAngle)/finalRadius                                                            - (b + 2*c*theta_f1 + 3*d*theta_f2);
    vectorB[2] = centralBodyGravitationalParameter_/(std::pow(finalRadius, 4.0)*std::pow(finalAngularVelocity,2.0)) - (1/finalRadius + 2*c + 6*d*theta_f1);

    Eigen::Vector3d efgCoefficients = (1.0/(2.0*std::pow(travelledAzimuthAngle_,6.0)))*matrixA*vectorB;

    coefficientsRadialDistanceFunction_[0] = a;
    coefficientsRadialDistanceFunction_[1] = b;
    coefficientsRadialDistanceFunction_[2] = c;
    coefficientsRadialDistanceFunction_[3] = d;
    coefficientsRadialDistanceFunction_[4] = efgCoefficients[0];
    coefficientsRadialDistanceFunction_[5] = efgCoefficients[1];
    coefficientsRadialDistanceFunction_[6] = efgCoefficients[2];

    // 3D case

    double az = initialState_[2];
    double bz = initialState_[5]/initialAngularVelocity;
    double q  = 7.0; // power of Z function

    Eigen::MatrixXd matrixC = Eigen::MatrixXd::Zero( 2, 2 );
    matrixC( 0, 0 ) =  q*travelledAzimuthAngle_ ;
    matrixC( 0, 1 ) = -std::pow(travelledAzimuthAngle_, 2.0);
    matrixC( 1, 0 ) = -(q-1.0);
    matrixC( 1, 1 ) =  travelledAzimuthAngle_;

    Eigen::Vector2d vectorD;

    vectorD[0] = finalState_[2] - az - bz*travelledAzimuthAngle_;
    vectorD[1] = (finalState_[5]/finalAngularVelocity) - bz;

    Eigen::Vector2d cdzCoefficients = (1.0/std::pow(travelledAzimuthAngle_,q))*matrixC*vectorD;

    coefficientsZFunction_[0] = az;
    coefficientsZFunction_[1] = bz;
    coefficientsZFunction_[2] = cdzCoefficients[0];
    coefficientsZFunction_[3] = cdzCoefficients[1];


    radialDistanceCompositeFunction_->resetCompositeFunctionCoefficients( coefficientsRadialDistanceFunction_ );
    zCompositeFunction_->resetCompositeFunctionCoefficients( coefficientsZFunction_ );

}


//! Compute current time from azimuth angle.
double InvPolyShaping::computeCurrentTimeFromAzimuthAngle( const double currentAzimuthAngle )
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

//! Base Time Of Flight function (integrand)
double InvPolyShaping::computeBaseTimeOfFlightFunction( double currentAzimuthAngle )
{
    double radialDistance = radialDistanceCompositeFunction_->evaluateCompositeFunction( currentAzimuthAngle );

    return std::sqrt( (std::pow(radialDistance,4.0)/centralBodyGravitationalParameter_)*
           (  (1.0/radialDistance)
            + 2.0 *coefficientsRadialDistanceFunction_[2]
            + 6.0 *coefficientsRadialDistanceFunction_[3]*currentAzimuthAngle
            + 12.0*coefficientsRadialDistanceFunction_[4]*std::pow(currentAzimuthAngle,2.0)
            + 20.0*coefficientsRadialDistanceFunction_[5]*std::pow(currentAzimuthAngle,3.0)
            + 30.0*coefficientsRadialDistanceFunction_[6]*std::pow(currentAzimuthAngle,4.0)));

}

//! Compute Time Of Flight (numerical integral)
double InvPolyShaping::computeNormalizedTimeOfFlight()
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
void InvPolyShaping::iterateToMatchRequiredTimeOfFlight( std::shared_ptr< root_finders::RootFinderSettings > rootFinderSettings,
                                                           const double lowerBound,
                                                           const double upperBound,
                                                           const double initialGuess )
{

    // Define the structure updating the time of flight from the free coefficient value, while still satisfying the boundary conditions.
    std::function< void ( const double ) > resetFreeCoefficientFunction = std::bind( &InvPolyShaping::resetValueFreeCoefficient, this, std::placeholders::_1 );
    std::function< void( ) > satisfyBoundaryConditionsFunction = std::bind( &InvPolyShaping::satisfyBoundaryConditions, this );
    std::function< double ( ) >  computeTOFfunction = std::bind( &InvPolyShaping::computeNormalizedTimeOfFlight, this );
    std::function< double ( ) > getRequiredTOFfunction = std::bind( &InvPolyShaping::getNormalizedRequiredTimeOfFlight, this );

    std::shared_ptr< basic_mathematics::Function< double, double > > timeOfFlightFunction =
            std::make_shared< InvPolyShaping::TimeOfFlightFunction >( resetFreeCoefficientFunction, satisfyBoundaryConditionsFunction, computeTOFfunction, getRequiredTOFfunction );

    // Create root finder from root finder settings.
    std::shared_ptr< root_finders::RootFinderCore< double > > rootFinder = root_finders::createRootFinder( rootFinderSettings, lowerBound, upperBound, initialGuess );

    // Iterate to find the free coefficient value that matches the required time of flight.
    double updatedFreeCoefficient = rootFinder->execute( timeOfFlightFunction, initialGuess );

}

//! Compute current spherical position. - TO DO - NOT USED
Eigen::Vector3d InvPolyShaping::computePositionVectorInSphericalCoordinates( const double currentAzimuthAngle )
{

    return ( Eigen::Vector3d() <<
             radialDistanceCompositeFunction_->evaluateCompositeFunction( currentAzimuthAngle ),
             currentAzimuthAngle,
             0).finished();

}

//! Compute current derivative of the azimuth angle. - check units
double InvPolyShaping::computeFirstDerivativeAzimuthAngleWrtTime( const double currentAzimuthAngle )
{
    double radialDistance = radialDistanceCompositeFunction_->evaluateCompositeFunction( currentAzimuthAngle );


    return std::sqrt( (centralBodyGravitationalParameter_/std::pow(radialDistance,4.0))*
           1.0/ (  1.0/radialDistance
              + 2.0 *coefficientsRadialDistanceFunction_[2]
              + 6.0 *coefficientsRadialDistanceFunction_[3]*currentAzimuthAngle
              + 12.0*coefficientsRadialDistanceFunction_[4]*std::pow(currentAzimuthAngle,2.0)
              + 20.0*coefficientsRadialDistanceFunction_[5]*std::pow(currentAzimuthAngle,3.0)
              + 30.0*coefficientsRadialDistanceFunction_[6]*std::pow(currentAzimuthAngle,4.0)));

}

//! Compute second derivative of the azimuth angle. - TO DO - check units
double InvPolyShaping::computeSecondDerivativeAzimuthAngleWrtTime( const double currentAzimuthAngle )
{
    double currentRadius             = radialDistanceCompositeFunction_->evaluateCompositeFunction(currentAzimuthAngle);
    double currentTanFlightPathAngle = computeTanFlightPathAngle(currentAzimuthAngle);

    double c = coefficientsRadialDistanceFunction_[2];
    double d = coefficientsRadialDistanceFunction_[3];
    double e = coefficientsRadialDistanceFunction_[4];
    double f = coefficientsRadialDistanceFunction_[5];
    double g = coefficientsRadialDistanceFunction_[6];

    double A1 = -centralBodyGravitationalParameter_/(2.0*std::pow(currentRadius, 4.0));
    double A2 = 4.0*currentTanFlightPathAngle;
    double A3 = (1.0/currentRadius) + 2.0*c + 6.0*d*currentAzimuthAngle + 12.0*e*std::pow(currentAzimuthAngle, 2.0) + 20.0*f*std::pow(currentAzimuthAngle, 3.0) + 30.0*g*std::pow(currentAzimuthAngle, 4.0);
    double A4 = 6.0*d + 24.0*e*currentAzimuthAngle + 60.0*f*std::pow(currentAzimuthAngle, 2.0) + 120.0*g*std::pow(currentAzimuthAngle, 3.0) - currentTanFlightPathAngle/currentRadius;
    double A5 = std::pow(A3, 2.0);

    return A1*(A2/A3 + A4/A5);

}


//! Compute current velocity in spherical coordinates. - TO DO - NOT USED
Eigen::Vector3d InvPolyShaping::computeVelocityVectorInSphericalCoordinates(const double currentAzimuthAngle )
{
    // Compute first derivative of the azimuth angle w.r.t. time.
    double derivativeAzimuthAngle = computeFirstDerivativeAzimuthAngleWrtTime( currentAzimuthAngle );

    // Compute and return current velocity vector in spherical coordinates.
    return derivativeAzimuthAngle * computeCurrentVelocityParametrizedByAzimuthAngle( currentAzimuthAngle );
}


//! Compute current state vector in spherical coordinates. - TO DO - NOT USED
Eigen::Vector6d InvPolyShaping::computeStateVectorInSphericalCoordinates( const double currentAzimuthAngle )
{
    Eigen::Vector6d currentSphericalState;
    currentSphericalState.segment( 0, 3 ) = computePositionVectorInSphericalCoordinates( currentAzimuthAngle );
    currentSphericalState.segment( 3, 3 ) = computeVelocityVectorInSphericalCoordinates( currentAzimuthAngle );

    return currentSphericalState;
}

//! Compute tangents of flight path angle
double InvPolyShaping::computeTanFlightPathAngle( const double currentAzimuthAngle )
{
    double currentNormalizedRadius      = radialDistanceCompositeFunction_->evaluateCompositeFunction( currentAzimuthAngle );
    return -currentNormalizedRadius*(coefficientsRadialDistanceFunction_[1]
            + 2.0*coefficientsRadialDistanceFunction_[2]*currentAzimuthAngle
            + 3.0*coefficientsRadialDistanceFunction_[3]*std::pow(currentAzimuthAngle, 2.0)
            + 4.0*coefficientsRadialDistanceFunction_[4]*std::pow(currentAzimuthAngle, 3.0)
            + 5.0*coefficientsRadialDistanceFunction_[5]*std::pow(currentAzimuthAngle, 4.0)
            + 6.0*coefficientsRadialDistanceFunction_[6]*std::pow(currentAzimuthAngle, 5.0));
}

//! Compute current cartesian state. - TO DO - check formulas for velocities
Eigen::Vector6d InvPolyShaping::computeNormalizedStateVector( const double currentAzimuthAngle )
{
    // Initialisation
    double currentNormalizedRadius  = radialDistanceCompositeFunction_->evaluateCompositeFunction( currentAzimuthAngle ); // need normalised radius
    double currentNormalizedZ       = zCompositeFunction_->evaluateCompositeFunction( currentAzimuthAngle ); // need normalised radius

    // Position
    Eigen::Vector3d currentPositionUnit   = initialPositionUnit_*std::cos(currentAzimuthAngle) + axisOfRotation_.cross(initialPositionUnit_)*std::sin(currentAzimuthAngle) + axisOfRotation_*(axisOfRotation_.dot(initialPositionUnit_))*(1-std::cos(currentAzimuthAngle));
    Eigen::Vector3d currentPosition       = currentPositionUnit*currentNormalizedRadius;

    // Velocity
    double firstDerivativeAzimuthAngle  = computeFirstDerivativeAzimuthAngleWrtTime( currentAzimuthAngle );
    double currentTanGamma              = computeTanFlightPathAngle( currentAzimuthAngle );

    double radialVelocityNorm  = currentNormalizedRadius*firstDerivativeAzimuthAngle*currentTanGamma;
    double angularVelocityNorm = currentNormalizedRadius*firstDerivativeAzimuthAngle;

    Eigen::Vector3d currentRadialVelocity  = radialVelocityNorm*currentPositionUnit;
    Eigen::Vector3d currentAngularVelocity = angularVelocityNorm*axisOfRotation_.cross(currentPositionUnit);
    Eigen::Vector3d currentVelocity        = currentRadialVelocity+currentAngularVelocity;

    Eigen::Vector6d normalisedStateVector;
    normalisedStateVector.segment( 0, 3 ) = currentPosition;
    normalisedStateVector.segment( 3, 3 ) = currentVelocity;


    Eigen::Vector6d currentCylindricalStateVector;

    // r,theta,z,Vr,Vtheta,,Vz
    currentCylindricalStateVector[0] = currentNormalizedRadius;
    currentCylindricalStateVector[1] = currentAzimuthAngle + initialAzimuthAngle_;
    currentCylindricalStateVector[2] = currentNormalizedZ;
    currentCylindricalStateVector[3] = radialDistanceCompositeFunction_->evaluateCompositeFunctionFirstDerivative(currentAzimuthAngle)*firstDerivativeAzimuthAngle;
    currentCylindricalStateVector[4] = radialDistanceCompositeFunction_->evaluateCompositeFunction(currentAzimuthAngle)*firstDerivativeAzimuthAngle;
    currentCylindricalStateVector[5] = zCompositeFunction_->evaluateCompositeFunctionFirstDerivative(currentAzimuthAngle)*firstDerivativeAzimuthAngle;

    normalisedStateVector = coordinate_conversions::convertCylindricalToCartesianState(currentCylindricalStateVector);

    return normalisedStateVector;
}


//! Compute current cartesian state. - TO DO - check units
Eigen::Vector6d InvPolyShaping::computeCurrentStateVector( const double currentAzimuthAngle )
{
    Eigen::Vector6d dimensionalStateVector;
    dimensionalStateVector.segment( 0, 3 ) = computeNormalizedStateVector( currentAzimuthAngle ).segment( 0, 3 )
            * physical_constants::ASTRONOMICAL_UNIT;
    dimensionalStateVector.segment( 3, 3 ) = computeNormalizedStateVector( currentAzimuthAngle ).segment( 3, 3 )
            * physical_constants::ASTRONOMICAL_UNIT/ physical_constants::JULIAN_YEAR; // no julian year needed as computed velocity is already per second, physical_constants::JULIAN_YEAR;

    return dimensionalStateVector;
}


//! Compute current velocity in spherical coordinates parametrized by azimuth angle theta. - TO DO - NOT USED
Eigen::Vector3d InvPolyShaping::computeCurrentVelocityParametrizedByAzimuthAngle( const double currentAzimuthAngle )
{

    // Retrieve current radial distance and elevation angle, as well as their derivatives w.r.t. azimuth angle.
    double radialDistance = radialDistanceCompositeFunction_->evaluateCompositeFunction( currentAzimuthAngle );
    double derivativeRadialDistance = radialDistanceCompositeFunction_->evaluateCompositeFunctionFirstDerivative( currentAzimuthAngle );

    // Compute and return velocity vector parametrized by azimuth angle.
    return ( Eigen::Vector3d() << derivativeRadialDistance,
             radialDistance * std::cos( 0 ),
             1.0 ).finished();
}


//! Compute current acceleration in spherical coordinates parametrized by azimuth angle theta. - TO DO
Eigen::Vector3d InvPolyShaping::computeCurrentAccelerationParametrizedByAzimuthAngle( const double currentAzimuthAngle )
{
    // Retrieve spherical coordinates and their derivatives w.r.t. to the azimuth angle.



    double radialDistance                = radialDistanceCompositeFunction_->evaluateCompositeFunction( currentAzimuthAngle );
    double firstDerivativeAzimuthAngle   = computeFirstDerivativeAzimuthAngleWrtTime( currentAzimuthAngle );
    double tanFlightPathAngle            = computeTanFlightPathAngle( currentAzimuthAngle );

    double flightPathAngle = std::atan(tanFlightPathAngle);

    double A1 = -centralBodyGravitationalParameter_/(2*std::pow(currentAzimuthAngle, 3.0)*std::cos(flightPathAngle));
    double A2 =   6.0*coefficientsRadialDistanceFunction_[3]
                + 24.0*coefficientsRadialDistanceFunction_[4]*currentAzimuthAngle
                + 60.0*coefficientsRadialDistanceFunction_[5]*std::pow(currentAzimuthAngle, 2.0)
                + 120.0*coefficientsRadialDistanceFunction_[6]*std::pow(currentAzimuthAngle, 3.0)
                - tanFlightPathAngle/radialDistance;

    double A3 =   1/radialDistance
                + 2.0*coefficientsRadialDistanceFunction_[2]
                + 6.0*coefficientsRadialDistanceFunction_[3]*currentAzimuthAngle
                + 12.0*coefficientsRadialDistanceFunction_[4]*std::pow(currentAzimuthAngle, 2.0)
                + 20.0*coefficientsRadialDistanceFunction_[5]*std::pow(currentAzimuthAngle, 3.0)
                + 30.0*coefficientsRadialDistanceFunction_[6]*std::pow(currentAzimuthAngle, 4.0);

    double thrustAccelerationMagnitude = A1*(A2/std::pow(A3, 2.0));


    // Compute and return acceleration vector parametrized by the azimuth angle theta.
    Eigen::Vector3d accelerationParametrizedByAzimuthAngle;
    accelerationParametrizedByAzimuthAngle[ 0 ] = 1.0;
    accelerationParametrizedByAzimuthAngle[ 1 ] = 1.0;
    accelerationParametrizedByAzimuthAngle[ 2 ] = 1.0;

    return accelerationParametrizedByAzimuthAngle;

}

//! Compute thrust acceleration vector in spherical coordinates. - TO DO - NOT USED ??
Eigen::Vector3d InvPolyShaping::computeThrustAccelerationInSphericalCoordinates( const double currentAzimuthAngle )
{
    // Compute current radial distance.
    double radialDistance = radialDistanceCompositeFunction_->evaluateCompositeFunction( currentAzimuthAngle );

    // Compute first and second derivatives of the azimuth angle w.r.t. time.
    double firstDerivativeAzimuthAngleWrtTime = computeFirstDerivativeAzimuthAngleWrtTime( currentAzimuthAngle );
    double secondDerivativeAzimuthAngleWrtTime = computeSecondDerivativeAzimuthAngleWrtTime( currentAzimuthAngle );

    // Compute velocity vector parametrized by azimuth angle theta.
    Eigen::Vector3d velocityParametrizedByAzimuthAngle = computeCurrentVelocityParametrizedByAzimuthAngle( currentAzimuthAngle );

    // Compute acceleration vector parametrized by azimuth angle theta.
    Eigen::Vector3d accelerationParametrizedByAzimuthAngle = computeCurrentAccelerationParametrizedByAzimuthAngle( currentAzimuthAngle );


    // Compute and return the current thrust acceleration vector in spherical coordinates.
    return std::pow( firstDerivativeAzimuthAngleWrtTime, 2.0 ) * accelerationParametrizedByAzimuthAngle
            + secondDerivativeAzimuthAngleWrtTime * velocityParametrizedByAzimuthAngle
            + centralBodyGravitationalParameter_ / std::pow( radialDistance, 3.0 ) * ( Eigen::Vector3d() << radialDistance, 0.0, 0.0 ).finished();

}

//! Compute current thrust acceleration in cartesian coordinates. - TO DO
Eigen::Vector3d InvPolyShaping::computeNormalizedThrustAccelerationVector( const double currentAzimuthAngle )
{
    Eigen::Vector6d normalizedStateVector;
    Eigen::Vector3d cartesianAcceleration;
    Eigen::Vector2d inPlaneThrustAcceleration;
    Eigen::Vector2d normalisedInPlaneThrustAcceleration;


    double inPlaneThrustMagnitude    = computeCurrentInPlaneThrustAccelerationMagnitude( currentAzimuthAngle );
    double outOfPlaneThrustMagnitude = computeCurrentOutOfPlaneThrustAccelerationMagnitude( currentAzimuthAngle );

    normalizedStateVector    = computeNormalizedStateVector(currentAzimuthAngle);

    inPlaneThrustAcceleration(0) = normalizedStateVector(3);
    inPlaneThrustAcceleration(1) = normalizedStateVector(4);

    normalisedInPlaneThrustAcceleration = inPlaneThrustAcceleration.normalized();

    cartesianAcceleration(0) = normalisedInPlaneThrustAcceleration(0)*inPlaneThrustMagnitude;
    cartesianAcceleration(1) = normalisedInPlaneThrustAcceleration(1)*inPlaneThrustMagnitude;
    cartesianAcceleration(2) = outOfPlaneThrustMagnitude;

    return cartesianAcceleration;
}


//! Compute in-plane magnitude cartesian acceleration. - TO DO - TO CHECK
double  InvPolyShaping::computeCurrentInPlaneThrustAccelerationMagnitude( double currentAzimuthAngle )
{
    double currentRadius             = radialDistanceCompositeFunction_->evaluateCompositeFunction(currentAzimuthAngle);
    double currentTanFlightPathAngle = computeTanFlightPathAngle(currentAzimuthAngle);
    double currentFlightPathAngle    = std::atan(currentTanFlightPathAngle);

    double c = coefficientsRadialDistanceFunction_[2];
    double d = coefficientsRadialDistanceFunction_[3];
    double e = coefficientsRadialDistanceFunction_[4];
    double f = coefficientsRadialDistanceFunction_[5];
    double g = coefficientsRadialDistanceFunction_[6];

    double A1 = -centralBodyGravitationalParameter_/(2.0*std::pow(currentRadius,3.0)*std::cos(currentFlightPathAngle));
    double A2 = 6.0*d + 24.0*e*currentAzimuthAngle + 60.0*f*std::pow(currentAzimuthAngle, 2.0) + 120.0*g*std::pow(currentAzimuthAngle, 3.0) - currentTanFlightPathAngle/currentRadius;
    double A3 = (1.0/currentRadius) + 2.0*c + 6.0*d*currentAzimuthAngle + 12.0*e*std::pow(currentAzimuthAngle, 2.0) + 20.0*f*std::pow(currentAzimuthAngle, 3.0) + 30.0*g*std::pow(currentAzimuthAngle, 4.0);

    return A1*(A2/std::pow(A3, 2.0));
}

//! Compute out-of-plane magnitude cartesian acceleration. - TO DO - TO CHECK
double  InvPolyShaping::computeCurrentOutOfPlaneThrustAccelerationMagnitude( double currentAzimuthAngle )
{
    double currentRadius             = radialDistanceCompositeFunction_->evaluateCompositeFunction(currentAzimuthAngle);

    double currentFirstDerivativeAzimuthAngleWrtTime  = computeFirstDerivativeAzimuthAngleWrtTime(currentAzimuthAngle);
    double currentSecondDerivativeAzimuthAngleWrtTime = computeSecondDerivativeAzimuthAngleWrtTime(currentAzimuthAngle);
    double currentZ = zCompositeFunction_->evaluateCompositeFunction(currentAzimuthAngle);

    double bz = coefficientsZFunction_[1];
    double cz = coefficientsZFunction_[2];
    double dz = coefficientsZFunction_[3];
    double q  = 7.0;// zth power

    double B1 = (centralBodyGravitationalParameter_/std::pow(currentRadius, 3.0))*currentZ;
    double B2 = cz*(q-1.0)*(q-2.0)*std::pow(currentAzimuthAngle,q-3.0) + dz*q*(q-1.0)*std::pow(currentAzimuthAngle,q-2.0);
    double B3 = bz + cz*(q-1.0)*std::pow(currentAzimuthAngle,q-2.0) + dz*q*std::pow(currentAzimuthAngle,q-1.0);

    return B1 - B2*std::pow(currentFirstDerivativeAzimuthAngleWrtTime,2.0) + B3*currentSecondDerivativeAzimuthAngleWrtTime;
}


//! Compute current thrust acceleration magnitude.
double InvPolyShaping::computeNormalizedThrustAccelerationMagnitude( double currentAzimuthAngle )
{
    double Tain = computeCurrentInPlaneThrustAccelerationMagnitude( currentAzimuthAngle );
    double Taz  = computeCurrentOutOfPlaneThrustAccelerationMagnitude( currentAzimuthAngle );

    // cant unitialise here cause this function is used for delta V
    return std::sqrt(std::pow(Tain, 2.0) + std::pow(Taz, 2.0));
}

//! Compute magnitude cartesian acceleration. - TO DO - TO CHECK
double  InvPolyShaping::computeCurrentThrustAccelerationMagnitude( double currentAzimuthAngle )
{
    return computeNormalizedThrustAccelerationMagnitude( currentAzimuthAngle ) * physical_constants::ASTRONOMICAL_UNIT / physical_constants::JULIAN_YEAR / physical_constants::JULIAN_YEAR;;//  *physical_constants::ASTRONOMICAL_UNIT / physical_constants::JULIAN_YEAR / physical_constants::JULIAN_YEAR;
}

//! Compute direction thrust acceleration in cartesian coordinates. - TO DO
Eigen::Vector3d InvPolyShaping::computeCurrentThrustAccelerationDirection( double currentAzimuthAngle )
{
    return computeNormalizedThrustAccelerationVector( currentAzimuthAngle ).normalized();
}

//! Compute final deltaV. - TO DO
double InvPolyShaping::computeDeltaV( )
{
    // Define the derivative of the deltaV, ie thrust acceleration function, as a function of the azimuth angle.
    std::function< double( const double ) > derivativeFunctionDeltaV = [ = ] ( const double currentAzimuthAngle ){
        return std::abs(computeNormalizedThrustAccelerationMagnitude( currentAzimuthAngle ))/
                computeFirstDerivativeAzimuthAngleWrtTime( currentAzimuthAngle );
    };

    // Define numerical quadrature from quadratrure settings.
    std::shared_ptr< numerical_quadrature::NumericalQuadrature< double, double > > quadrature =
            numerical_quadrature::createQuadrature( derivativeFunctionDeltaV, quadratureSettings_, travelledAzimuthAngle_ );

    // Return dimensional deltaV [m/s]
    return quadrature->getQuadrature( ) * physical_constants::ASTRONOMICAL_UNIT / physical_constants::JULIAN_YEAR;
}




////! Get low-thrust acceleration model from shaping method.
//std::shared_ptr< propulsion::ThrustAcceleration > InvPolyShaping::getLowThrustAccelerationModel(
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


void InvPolyShaping::computeShapedTrajectoryAndFullPropagation(
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
                                           std::bind( &InvPolyShaping::computeCurrentThrustAccelerationDirection, this, std::placeholders::_1 ),
                                           std::bind( &InvPolyShaping::computeCurrentThrustAccelerationMagnitude, this, std::placeholders::_1 )*/ /*, interpolator*/ );

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
