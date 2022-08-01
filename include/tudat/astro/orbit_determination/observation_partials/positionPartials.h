/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_POSITIONPARTIALS_H
#define TUDAT_POSITIONPARTIALS_H

#include <vector>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <functional>

#include "tudat/basics/basicTypedefs.h"

#include "tudat/astro/ephemerides/simpleRotationalEphemeris.h"
#include "tudat/astro/orbit_determination/observation_partials/rotationMatrixPartial.h"
#include "tudat/astro/orbit_determination/observation_partials/observationPartial.h"
#include "tudat/astro/observation_models/linkTypeDefs.h"

namespace tudat
{

namespace observation_partials
{

//! Function to calculate the partial of position of a point on body wrt position of that body.
/*!
 *  Function to calculate the partial of position of a point on body wrt position of that body,
 *  with both positions expressed in the same frame.
 *  \return Requested partial (3x3 identity matrix)
 */
Eigen::Matrix3d calculatePartialOfPointPositionWrtBodyPosition( );

//! Function to calculate the partial of position of a point on a body wrt its body-fixed position
/*!
 *  Function to calculate the partial of position of a point on a body (expressed in a non-corotating, non-body fixed frame)
 *  wrt its body-fixed position, e.g. the partial of a ground station's inertial
 *  position wrt its body-fixed position.
 *  \param rotationMatrixToInertialFrame Rotation matrix from body-fixed to inertial frame.
 *  \return Partial of position of a point on a body wrt its body-fixed position
 */
Eigen::Matrix3d calculatePartialOfPointPositionWrtBodyFixedPointPosition(
        const Eigen::Matrix3d& rotationMatrixToInertialFrame );


//! Base class for calculating the partial of an inertial position wrt a parameter.
/*!
 *  Base class for calculating the partial of an inertial position wrt a parameter. A derived class is implemented for
 *  each (type of) estimatable parameter. A separate instance of the class must be made for each distinct state.
 *  Note that partials wrt parameters describing a property of a rotation matrix from a local to the inertial frame is
 *  implemented in the RotationMatrixPartial class, which is
 *  then used by the CartesianStatePartialWrtRotationMatrixParameter derived class of this class
 */
class CartesianStatePartial
{
public:

    CartesianStatePartial( ){ }

    //! Destructor.
    virtual ~CartesianStatePartial( ){ }

    //! Pure virtual base class function for determining partial at current time and body state.
    /*!
     *  Pure virtual base class function for determining partial at current time and body state wrt parameter of specific
     *  implemented derived. class.
     *  \param state Current inertial state of point of which partial is to be calculated by derived class
     *  \param time Current time
     *  \return Partial of point position wrt parameter (with specific parameter determined by derived class implementation).
     */
    virtual Eigen::Matrix< double, 3, Eigen::Dynamic > calculatePartialOfPosition(
            const Eigen::Vector6d& state, const double time ) = 0;

    virtual Eigen::Matrix< double, 3, Eigen::Dynamic > calculatePartialOfVelocity(
            const Eigen::Vector6d& state, const double time ) = 0;
};

//! Class to compute the partial derivative of the Cartesian state of a body w.r.t. to inertial three-dimensional
//! position of this body
class CartesianStatePartialWrtCartesianState: public CartesianStatePartial
{
public:

    //! Constructor
    CartesianStatePartialWrtCartesianState( )
    {
        positionPartial_.setZero( );
        positionPartial_.block( 0, 0, 3, 3 ) = Eigen::Matrix3d::Identity( );
        positionPartial_.block( 0, 3, 3, 3 ) = Eigen::Matrix3d::Zero( );

        velocityPartial_.setZero( );
        velocityPartial_.block( 0, 0, 3, 3 ) = Eigen::Matrix3d::Zero( );
        velocityPartial_.block( 0, 3, 3, 3 ) = Eigen::Matrix3d::Identity( );

    }

    //! Destructor
    ~CartesianStatePartialWrtCartesianState( ){ }

    //! Function for determining partial of position at current time and body state.
    /*!
     *  Function for determining partial of position  at current time and body state wrt three-dimensional state
     *  \param state Current inertial state of point of which partial is to be calculated.
     *  \param time Current time
     *  \return Partial of point state wrt position
     */
    Eigen::Matrix< double, 3, Eigen::Dynamic > calculatePartialOfPosition(
            const Eigen::Vector6d& state,
            const double time )
    {
        return positionPartial_;
    }
    //! Function for determining partial of velocity at current time and body state.
    /*!
     *  Function for determining partial of velocity  at current time and body state wrt three-dimensional state
     *  \param state Current inertial state of point of which partial is to be calculated.
     *  \param time Current time
     *  \return Partial of point state wrt velocity
     */
    Eigen::Matrix< double, 3, Eigen::Dynamic > calculatePartialOfVelocity(
                const Eigen::Vector6d& state, const double time )
    {
        return velocityPartial_;
    }
private:

    //! Partial of current state w.r.t. position
    Eigen::Matrix< double, 3, 6 > positionPartial_;

    //! Partial of current state w.r.t. velocity
    Eigen::Matrix< double, 3, 6 > velocityPartial_;

};



//! Class to compute the partial derivative of the Cartesian state of a body w.r.t. to a property of a rotation
//! matrix to/from a body-fixed frame.
class CartesianStatePartialWrtRotationMatrixParameter: public CartesianStatePartial
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param rotationMatrixPartialObject Object to compute the associated partial of a rotation matrix
     * \param positionFunctionInLocalFrame Function returning the body-fixed position of the point at which the partial
     * is to be computed.
     */
    CartesianStatePartialWrtRotationMatrixParameter(
            const std::shared_ptr< RotationMatrixPartial > rotationMatrixPartialObject,
            const std::function< Eigen::Vector3d( const double ) > positionFunctionInLocalFrame ):
        rotationMatrixPartialObject_( rotationMatrixPartialObject ),
        positionFunctionInLocalFrame_( positionFunctionInLocalFrame ){ }

    //! Destructor
    ~CartesianStatePartialWrtRotationMatrixParameter( ){ }

    //! Function for determining partial of position at current time and body state.
    /*!
     *  Function for determining partial of position  at current time and body state wrt three-dimensional state
     *  \param state Current inertial state of point of which partial is to be calculated.
     *  \param time Current time
     *  \return Partial of point state wrt rotation property
     */
    Eigen::Matrix< double, 3, Eigen::Dynamic > calculatePartialOfPosition(
            const Eigen::Vector6d& state,
            const double time )
    {
        return rotationMatrixPartialObject_->calculatePartialOfInertialPositionWrtParameter(
                    time, positionFunctionInLocalFrame_( time ) );
    }

    //! Function for determining partial of velocity at current time and body state.
    /*!
     *  Function for determining partial of velocity  at current time and body state wrt three-dimensional state
     *  \param state Current inertial state of point of which partial is to be calculated.
     *  \param time Current time
     *  \return Partial of point state wrt rotation property
     */
    Eigen::Matrix< double, 3, Eigen::Dynamic > calculatePartialOfVelocity(
            const Eigen::Vector6d& state,
            const double time )
    {
        return rotationMatrixPartialObject_->calculatePartialOfInertialVelocityWrtParameter( time, positionFunctionInLocalFrame_( time ) );
    }

private:

    //! Object to compute the associated partial of a rotation matrix
    std::shared_ptr< RotationMatrixPartial > rotationMatrixPartialObject_;

    //! Function returning the body-fixed position of the point at which the partial is to be computed.
    std::function< Eigen::Vector3d( const double ) > positionFunctionInLocalFrame_;
};

//! Class to compute the partial derivative of the inertial Cartesian state of a point on a body w.r.t. the constant
//! body-fixed position of that point.
class CartesianPartialWrtBodyFixedPosition: public CartesianStatePartial
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param bodyRotationModel Rotation model for body.
     */
    CartesianPartialWrtBodyFixedPosition( const std::shared_ptr< ephemerides::RotationalEphemeris > bodyRotationModel ):
        bodyRotationModel_( bodyRotationModel ){ }

    //! Destructor
    ~CartesianPartialWrtBodyFixedPosition( ){ }

    //! Function for determining partial of position at current time and body state.
    /*!
     *  Function for determining partial of position  at current time and body state wrt three-dimensional body-fixed state
     *  \param state Current inertial state of point of which partial is to be calculated.
     *  \param time Current time
     *  \return Partial of point state wrt body-fixed position.
     */
    Eigen::Matrix< double, 3, Eigen::Dynamic > calculatePartialOfPosition(
            const Eigen::Vector6d& state,
            const double time )
    {
        return Eigen::Matrix3d( bodyRotationModel_->getRotationToBaseFrame( time ) );
    }

    //! Function for determining partial of velocity at current time and body state.
    /*!
     *  Function for determining partial of velocity  at current time and body state wrt three-dimensional state
     *  \param state Current inertial state of point of which partial is to be calculated.
     *  \param time Current time
     *  \return Partial of point state wrt body-fixed position.
     */
    Eigen::Matrix< double, 3, Eigen::Dynamic > calculatePartialOfVelocity(
            const Eigen::Vector6d& state,
            const double time )
    {
        return bodyRotationModel_->getDerivativeOfRotationToBaseFrame( time );
    }

private:

    //! Rotation model for body.
    std::shared_ptr< ephemerides::RotationalEphemeris > bodyRotationModel_;
};

//! Derived class for scaling three-dimensional position partial to position observable partial
/*!
 *  Derived class for scaling three-dimensional position partial to position observable partial. Although the implementation
 *  is trivial for non-relativistic reference frames, it is included in the architecture pending future implementation
 *  of more rigorous reference frames.
 */
class PositionObservationScaling: public PositionPartialScaling
{
public:

    //! Destructor
    ~PositionObservationScaling( ){ }

    //! Update the scaling object to the current times and states (no functionality needed).
    /*!
     *  Update the scaling object to the current times and states (no functionality needed).
     *  \param linkEndStates List of states at each link end during observation.
     *  \param times List of times at each link end during observation.
     *  \param fixedLinkEnd Link end at which observation time is defined, i.e. link end for which associated time
     *  is kept constant when computing observable.
     *  \param currentObservation Value of observation for which partial scaling is to be computed
     */
    void update( const std::vector< Eigen::Vector6d >& linkEndStates,
                 const std::vector< double >& times,
                 const observation_models::LinkEndType fixedLinkEnd,
                 const Eigen::VectorXd currentObservation ){ }

    //! Function to retrieve the scaling factor for specific link end
    /*!
     * Function to retrieve the scaling factor for specific link end
     * \param linkEndType Link end for which scaling factor is to be returned
     * \return Position partial scaling factor at current link end
     */
    Eigen::Matrix< double, 3, 3 > getScalingFactor(
            const observation_models::LinkEndType linkEndType )
    {
        return Eigen::Matrix3d::Identity( );
    }

private:

};

//! Class to compute the partial derivatives of a three-dimensional position observable.
class PositionObervationPartial: public ObservationPartial< 3 >
{

public:

    //! Local typedef for return type (list of partial matrices and associated evaluation times).
    typedef std::vector< std::pair< Eigen::Matrix< double, 3, Eigen::Dynamic >, double > >
    PositionObservationPartialReturnType;

    //! Constructor
    /*!
     * Constructor
     * \param positionObservationScaler Scaling object used for mapping partials of positions to partials of observable
     * \param positionPartialList List of position partial per link end.
     * \param parameterIdentifier Id of parameter for which instance of class computes partial derivatives
     */
    PositionObervationPartial(
            const std::shared_ptr< PositionObservationScaling > positionObservationScaler,
            const std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > >& positionPartialList,
            const estimatable_parameters::EstimatebleParameterIdentifier parameterIdentifier ):
        ObservationPartial< 3 >( parameterIdentifier ), positionObservationScaler_( positionObservationScaler ),
        positionPartialList_( positionPartialList )
    { }

    //! Destructor
    ~PositionObervationPartial( ) { }

    //! Fnuction to calculate the observation partial(s) at required time and state
    /*!
     *  Function to calculate the observation partial(s) at required time and state. State and time
     *  are typically obtained from evaluation of observation model.
     *  \param states Link end stats. Index maps to link end for a given ObsevableType through getLinkEndIndex function.
     *  \param times Link end time.
     *  \param linkEndOfFixedTime Link end that is kept fixed when computing the observable.
     *  \param currentObservation Value of the observation for which the partial is to be computed (default NaN for
     *  compatibility purposes)
     *  \return Vector of pairs containing partial values and associated times.
     */
    virtual PositionObservationPartialReturnType calculatePartial(
            const std::vector< Eigen::Vector6d >& states,
            const std::vector< double >& times,
            const observation_models::LinkEndType linkEndOfFixedTime,
            const Eigen::Vector3d& currentObservation = Eigen::Vector3d::Constant( TUDAT_NAN ) )
    {
        PositionObservationPartialReturnType returnPartial;


        // Iterate over all link ends.
        for( positionPartialIterator_ = positionPartialList_.begin( );
             positionPartialIterator_ != positionPartialList_.end( );
             positionPartialIterator_++ )
        {
            // Retrieve link end time and state
            if( positionPartialIterator_->first == observation_models::observed_body )
            {
                currentState_ = states[ 0 ];
                currentTime_ = times[ 0 ];
            }
            else
            {
                throw std::runtime_error(
                            "Error when calculating position observation partial, invalid link end type requested" );
            }

            // Get position partial and scale with associated term.
            returnPartial.push_back(
                        std::make_pair(
                            positionObservationScaler_->getScalingFactor(
                                positionPartialIterator_->first ) *
                            ( positionPartialIterator_->second->calculatePartialOfPosition(
                                  currentState_, currentTime_ ) ), currentTime_ ) );
        }

        return returnPartial;
    }

protected:

    //!  Scaling object used for mapping partials of positions to partials of observable
    std::shared_ptr< PositionObservationScaling > positionObservationScaler_;

    //! List of position partial per link end.
    std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > > positionPartialList_;

    //! Iterator over list of position partial per link end (predeclared for efficiency).
    std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > >::iterator positionPartialIterator_;


    //! Pre-declared state variable to be used in calculatePartial function.
    Eigen::Vector6d currentState_;

    //! Pre-declared time variable to be used in calculatePartial function.
    double currentTime_;
};


//! Derived class for scaling three-dimensional velocity partial to velocity observable partial
/*!
 *  Derived class for scaling three-dimensional velocity partial to velocity observable partial. Although the implementation
 *  is trivial for non-relativistic reference frames, it is included in teh architecture pending future implementation
 *  of more rigorous reference frames.
 */
class VelocityObservationScaling: public PositionPartialScaling
{
public:

    //! Destructor
    ~VelocityObservationScaling( ){ }

    //! Update the scaling object to the current times and states (no functionality needed).
    /*!
     *  Update the scaling object to the current times and states (no functionality needed).
     *  \param linkEndStates List of states at each link end during observation.
     *  \param times List of times at each link end during observation.
     *  \param fixedLinkEnd Link end at which observation time is defined, i.e. link end for which associated time
     *  is kept constant when computing observable.
     *  \param currentObservation Value of observation for which partial scaling is to be computed
     */
    void update( const std::vector< Eigen::Vector6d >& linkEndStates,
                 const std::vector< double >& times,
                 const observation_models::LinkEndType fixedLinkEnd,
                 const Eigen::VectorXd currentObservation ){ }

    //! Function to retrieve the scaling factor for specific link end
    /*!
     * Function to retrieve the scaling factor for specific link end
     * \param linkEndType Link end for which scaling factor is to be returned
     * \return Velocity partial scaling factor at current link end
     */
    Eigen::Matrix< double, 3, 3 > getScalingFactor(
            const observation_models::LinkEndType linkEndType )
    {
        return Eigen::Matrix3d::Identity( );
    }

private:

};


//! Class to compute the partial derivatives of a three-dimensional velocity observable.
class VelocityObervationPartial: public ObservationPartial< 3 >
{

public:

    //! Local typedef for return type (list of partial matrices and associated evaluation times).
    typedef std::vector< std::pair< Eigen::Matrix< double, 3, Eigen::Dynamic >, double > >
    VelocityObservationPartialReturnType;

    //! Constructor
    /*!
     * Constructor
     * \param velocityObservationScaler Scaling object used for mapping partials of velocitys to partials of observable
     * \param velocityPartialList List of velocity partial per link end.
     * \param parameterIdentifier Id of parameter for which instance of class computes partial derivatives
     */
    VelocityObervationPartial(
            const std::shared_ptr< VelocityObservationScaling > velocityObservationScaler,
            const std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > >& velocityPartialList,
            const estimatable_parameters::EstimatebleParameterIdentifier parameterIdentifier ):
        ObservationPartial< 3 >( parameterIdentifier ), velocityObservationScaler_( velocityObservationScaler ),
        velocityPartialList_( velocityPartialList )
    { }

    //! Destructor
    ~VelocityObervationPartial( ) { }

    //! Fnuction to calculate the observation partial(s) at required time and state
    /*!
     *  Function to calculate the observation partial(s) at required time and state. State and time
     *  are typically obtained from evaluation of observation model.
     *  \param states Link end stats. Index maps to link end for a given ObsevableType through getLinkEndIndex function.
     *  \param times Link end time.
     *  \param linkEndOfFixedTime Link end that is kept fixed when computing the observable.
     *  \param currentObservation Value of the observation for which the partial is to be computed (default NaN for
     *  compatibility purposes)
     *  \return Vector of pairs containing partial values and associated times.
     */
    virtual VelocityObservationPartialReturnType calculatePartial(
            const std::vector< Eigen::Vector6d >& states,
            const std::vector< double >& times,
            const observation_models::LinkEndType linkEndOfFixedTime,
            const Eigen::Vector3d& currentObservation = Eigen::Vector3d::Constant( TUDAT_NAN ) )
    {
        VelocityObservationPartialReturnType returnPartial;


        // Iterate over all link ends.
        for( velocityPartialIterator_ = velocityPartialList_.begin( );
             velocityPartialIterator_ != velocityPartialList_.end( );
             velocityPartialIterator_++ )
        {
            // Retrieve link end time and state
            if( velocityPartialIterator_->first == observation_models::observed_body )
            {
                currentState_ = states[ 0 ];
                currentTime_ = times[ 0 ];
            }
            else
            {
                throw std::runtime_error(
                            "Error when calculating velocity observation partial, invalid link end type requested" );
            }

            // Get velocity partial and scale with associated term.
            returnPartial.push_back(
                        std::make_pair(
                            velocityObservationScaler_->getScalingFactor(
                                velocityPartialIterator_->first ) *
                            ( velocityPartialIterator_->second->calculatePartialOfVelocity(
                                  currentState_, currentTime_ ) ), currentTime_ ) );
        }

        return returnPartial;
    }

protected:

    //!  Scaling object used for mapping partials of velocitys to partials of observable
    std::shared_ptr< VelocityObservationScaling > velocityObservationScaler_;

    //! List of velocity partial per link end.
    std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > > velocityPartialList_;

    //! Iterator over list of velocity partial per link end (predeclared for efficiency).
    std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > >::iterator velocityPartialIterator_;


    //! Pre-declared state variable to be used in calculatePartial function.
    Eigen::Vector6d currentState_;

    //! Pre-declared time variable to be used in calculatePartial function.
    double currentTime_;
};

}

}


#endif // TUDAT_POSITIONPARTIALS_H
