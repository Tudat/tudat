/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


#ifndef TUDAT_ONEWAYDOPPLERPARTIAL_H
#define TUDAT_ONEWAYDOPPLERPARTIAL_H


#include <functional>

#include <Eigen/Core>

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/observation_models/oneWayDopplerObservationModel.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"
#include "tudat/astro/orbit_determination/observation_partials/observationPartial.h"
#include "tudat/astro/orbit_determination/observation_partials/positionPartials.h"
#include "tudat/astro/orbit_determination/observation_partials/lightTimeCorrectionPartial.h"
#include "tudat/astro/observation_models/linkTypeDefs.h"

namespace tudat
{

namespace observation_partials
{

//! Base class to compute the state partial scaling factors for the proper time component of one-way Doppler observables
class OneWayDopplerProperTimeComponentScaling: public PositionPartialScaling
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param linkEndWithPartial Link end for which this object computes partials
     */
    OneWayDopplerProperTimeComponentScaling( const observation_models::LinkEndType linkEndWithPartial ):
        linkEndWithPartial_( linkEndWithPartial ){ }

    //! Destructor
    virtual ~OneWayDopplerProperTimeComponentScaling( ){ }

    //! Function to retrieve the scaling factor for the derivative w.r.t. the position of a given link end
    /*!
     * Function to retrieve the scaling factor for the derivative w.r.t. the position of a given link end
     * \param linkEndType Link end position for which partial scaling is to be retrieved.
     * \return The scaling factor for the derivative w.r.t. the position of a given link end
     */
    virtual Eigen::Matrix< double, 1, 3 > getPositionScalingFactor( const observation_models::LinkEndType linkEndType ) = 0;

    //! Function to retrieve the scaling factor for the derivative w.r.t. the velocity of a given link end
    /*!
     * Function to retrieve the scaling factor for the derivative w.r.t. the velocity of a given link end
     * \param linkEndType Link end velocity for which partial scaling is to be retrieved.
     * \return The scaling factor for the derivative w.r.t. the velocity of a given link end
     */
    virtual Eigen::Matrix< double, 1, 3 > getVelocityScalingFactor( const observation_models::LinkEndType linkEndType ) = 0;

    //! Function to get the direct partial derivative, and associated time, of proper time
    /*!
     *  Function to get the direct partial derivative (e.g. independent of state partials), and associated time, of proper
     *  time, w.r.t. parameter defined by parameterType.
     *  \param parameterType Parameter for which partial is to be checked
     *  \return Direct partial derivative, and associated time, of proper time rate
     */
    virtual std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > getProperTimeParameterPartial(
            const estimatable_parameters::EstimatebleParameterIdentifier parameterType ) = 0;

    //! Function to get the size of the direct dependency of proper time rate on parameter
    /*!
     * Function to get the size of the direct dependency (e.g. independent of state dependency) of proper time rate on parameter.
     * Returns zero for no dependency.
     * \param parameterType Parameter for which dependency is to be checked.
     * \return Number of columns in direct partial derivative of proper time rates w.r.t. parameterType
     */
    virtual int getParameterDependencySize(
            const estimatable_parameters::EstimatebleParameterIdentifier parameterType ) = 0;

protected:

    //! Link end for which this object computes partials
    observation_models::LinkEndType linkEndWithPartial_;

};

//! Class to compute the state partial scaling factors for first-order proper time component of one-way Doppler observables
class OneWayDopplerDirectFirstOrderProperTimeComponentScaling: public OneWayDopplerProperTimeComponentScaling
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param properTimeRateModel Object used to compute proper time rate
     * \param linkEndWithPartial Link end for which this partial is created
     * \param computeStatePartials Boolean to denote whether state partials are to be computed. It is false if the link end for
     * which this object computes the proper time partials is fixed to the perturbing body.
     */
    OneWayDopplerDirectFirstOrderProperTimeComponentScaling(
            const std::shared_ptr< observation_models::DirectFirstOrderDopplerProperTimeRateInterface > properTimeRateModel,
            const observation_models::LinkEndType linkEndWithPartial,
            const bool computeStatePartials ):
        OneWayDopplerProperTimeComponentScaling( linkEndWithPartial ),
        properTimeRateModel_( properTimeRateModel ),
        computeStatePartials_( computeStatePartials ){ }

    //! Update the scaling object to the current times and states
    /*!
     *  Update the scaling object to the current times and states
     *  \param linkEndStates List of states at each link end during observation.
     *  \param times List of times at each link end during observation.
     *  \param fixedLinkEnd Link end at which observation time is defined, i.e. link end for which associated time
     *  is kept constant when computing observable.
     *  \param currentObservation Value of observation (proper time rate) for which partial scaling is to be computed
     */
    void update( const std::vector< Eigen::Vector6d >& linkEndStates,
                 const std::vector< double >& times,
                 const observation_models::LinkEndType fixedLinkEnd,
                 const Eigen::VectorXd currentObservation );

    //! Function to retrieve the scaling factor for the derivative w.r.t. the position of a given link end
    /*!
     * Function to retrieve the scaling factor for the derivative w.r.t. the position of a given link end
     * \param linkEndType Link end position for which partial scaling is to be retrieved.
     * \return The scaling factor for the derivative w.r.t. the position of a given link end
     */
    Eigen::Matrix< double, 1, 3 > getPositionScalingFactor( const observation_models::LinkEndType linkEndType );

    //! Function to retrieve the scaling factor for the derivative w.r.t. the velocity of a given link end
    /*!
     * Function to retrieve the scaling factor for the derivative w.r.t. the velocity of a given link end
     * \param linkEndType Link end velocity for which partial scaling is to be retrieved.
     * \return The scaling factor for the derivative w.r.t. the velocity of a given link end
     */
    Eigen::Matrix< double, 1, 3 > getVelocityScalingFactor( const observation_models::LinkEndType linkEndType );

    //! Function to compute partial of proper time rate w.r.t. equivalence principle violation parameter
    /*!
     * Function to compute partial of proper time rate w.r.t. equivalence principle violation parameter
     * \return Partial of proper time rate w.r.t. equivalence principle violation parameter
     */
    double getEquivalencePrincipleViolationParameterPartial( )
    {
        return -currentGravitationalParameter_ / currentDistance_ * physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT;
    }

    //! Function to get the direct partial derivative, and associated time, of proper time
    /*!
     *  Function to get the direct partial derivative (e.g. independent of state partials), and associated time, of proper
     *  time, w.r.t. parameter defined by parameterType.
     *  \param parameterType Parameter for which partial is to be checked
     *  \return Direct partial derivative, and associated time, of proper time rate
     */
    std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > getProperTimeParameterPartial(
            const estimatable_parameters::EstimatebleParameterIdentifier parameterType );

    //! Function to get the size of the direct dependency of proper time rate on parameter
    /*!
     * Function to get the size of the direct dependency (e.g. independent of state dependency) of proper time rate on parameter.
     * Returns zero for no dependency.
     * \param parameterType Parameter for which dependency is to be checked.
     * \return Number of columns in direct partial derivative of proper time rates w.r.t. parameterType
     */
    int getParameterDependencySize( const estimatable_parameters::EstimatebleParameterIdentifier parameterType );

private:

    //! Object used to compute proper time rate
    std::shared_ptr< observation_models::DirectFirstOrderDopplerProperTimeRateInterface > properTimeRateModel_;

    //! Partial of proper time rate w.r.t. position, as computed by last call to update function.
    Eigen::Matrix< double, 1, 3 > partialWrPosition_;

    //! Partial of proper time rate w.r.t. velocity, as computed by last call to update function.
    Eigen::Matrix< double, 1, 3 > partialWrtVelocity_;

    //! Current time at link end for which this object computes proper time ratre partials
    double currentLinkEndTime_;

    //! Current distance between perturbing body center of mass and computation point.
    double currentDistance_;

    //! Current value of gravitational parameter of central body.
    double currentGravitationalParameter_;

    //! Boolean to denote whether state partials are to be computed
    /*!
     *  Boolean to denote whether state partials are to be computed. It is false if the link end for whicj this object computes
     *  the proper time partials is fixed to the perturbing body.
     */
    bool computeStatePartials_;
};

//! Derived class for scaling three-dimensional position partial to one-way doppler observable partial
/*!
 *  Derived class for scaling three-dimensional position partial to one-way doppler observable partial. Implementation is taken
 *  from Moyer(2000) and is separately implemented for fixed receiver and transmitter.
 */
class OneWayDopplerScaling: public PositionPartialScaling
{
public:

    //! Destructor
    /*!
     * Destructor
     * \param transmitterAccelerationFunction Function returning the Cartesian acceleration of the transmitter as a function of
     * time.
     * \param receiverAccelerationFunction Function returning the Cartesian acceleration of the receiver as a function of
     * time.
     * \param transmitterProperTimePartials Object used to compute the contribution of receiver proper time rate to the scaling
     * \param receiverProperTimePartials Object used to compute the contribution of transmitter proper time rate to the scaling
     */
    OneWayDopplerScaling(
            const std::function< Eigen::Vector3d( const double ) > transmitterAccelerationFunction,
            const std::function< Eigen::Vector3d( const double ) > receiverAccelerationFunction,
            const std::shared_ptr< OneWayDopplerProperTimeComponentScaling > transmitterProperTimePartials = nullptr,
            const std::shared_ptr< OneWayDopplerProperTimeComponentScaling > receiverProperTimePartials = nullptr ):
        transmitterAccelerationFunction_( transmitterAccelerationFunction ),
        receiverAccelerationFunction_( receiverAccelerationFunction ),
        transmitterProperTimePartials_( transmitterProperTimePartials ),
        receiverProperTimePartials_( receiverProperTimePartials ){ }

    //! Destructor
    ~OneWayDopplerScaling( ){ }

    //! Update the scaling object to the current times and states
    /*!
     *  Update the scaling object to the current times and states
     *  \param linkEndStates List of states at each link end during observation Index of vector maps to link end for a
     *  given ObsevableType through getLinkEndIndex function.
     *  \param times List of times at each link end during observation.
     *  \param fixedLinkEnd Link end at which observation time is defined, i.e. link end for which associated time
     *  is kept constant when computing observable.
     *  \param currentObservation Value of observation for which partial scaling is to be computed
     */
    void update( const std::vector< Eigen::Vector6d >& linkEndStates,
                 const std::vector< double >& times,
                 const observation_models::LinkEndType fixedLinkEnd,
                 const Eigen::VectorXd currentObservation = Eigen::VectorXd::Constant( 1, TUDAT_NAN ) );


    //! Function to retrieve the position scaling factor for specific link end
    /*!
     * Function to retrieve the position scaling factor for specific link end
     * \param linkEndType Link end for which scaling factor is to be returned
     * \return Position partial scaling factor at current link end
     */
    Eigen::Matrix< double, 1, 3 > getPositionScalingFactor( const observation_models::LinkEndType linkEndType );

    //! Function to retrieve the velocity scaling factor for specific link end
    /*!
     * Function to retrieve the velocity scaling factor for specific link end
     * \param linkEndType Link end for which scaling factor is to be returned
     * \return Velocity partial scaling factor at current link end
     */
    Eigen::Matrix< double, 1, 3 > getVelocityScalingFactor( const observation_models::LinkEndType linkEndType );

    //! Function to get the fixed link end for last computation of update() function.
    /*!
     * Fixed link end for last computation of update() function.
     * \return Function to get the fixed link end for last computation of update() function.
     */
    observation_models::LinkEndType getCurrentLinkEndType( )
    {
        return currentLinkEndType_;
    }

    //! Function to return factor by which light-time correction state partial is to be scaled for  one-way Doppler partial
    /*!
     *  Function to return factor by which light-time correction state partial is to be scaled to be added to one-way
     *  \return Factor by which light-time correction state partial is to be scaled for  one-way Doppler partial
     */
    double getLightTimeCorrectionPartialScaling( )
    {
        return lightTimeEffectPositionScalingFactor_;
    }

    //! Function to return object used to compute the contribution of transmitter proper time rate to the scaling
    /*!
     * Function to return object used to compute the contribution of transmitter proper time rate to the scaling
     * \return Object used to compute the contribution of transmitter proper time rate to the scaling
     */
    std::shared_ptr< OneWayDopplerProperTimeComponentScaling > getTransmitterProperTimePartials( )
    {
        return transmitterProperTimePartials_;
    }

    //! Function to return object used to compute the contribution of receiver proper time rate to the scaling
    /*!
     * Function to return object used to compute the contribution of receiver proper time rate to the scaling
     * \return Object used to compute the contribution of receiver proper time rate to the scaling
     */
    std::shared_ptr< OneWayDopplerProperTimeComponentScaling > getReceiverProperTimePartials( )
    {
        return receiverProperTimePartials_;
    }

    //! Function to get the size of the direct dependency of proper time rate on parameter
    /*!
     * Function to get the size of the direct dependency (e.g. independent of state dependency) of proper time rate on parameter.
     * Returns zero for no dependency. Requires that the transmitter/receiver have the same size of dependency, or that
     * one of them has dependency size zero (no dependency).
     * \param parameterType Parameter for which dependency is to be checked.
     * \return Number of columns in direct partial derivative of proper time rates w.r.t. parameterType
     */
    int getProperTimeParameterDependencySize(
            const estimatable_parameters::EstimatebleParameterIdentifier parameterType );

    //! Function to get the direct partial derivatives, and associated times, of proper time components of Doppler partials
    /*!
     *  Function to get the direct partial derivatives (e.g. independent of state partials), and associated times, of proper
     *  time components of Doppler partials, w.r.t. parameter defined by parameterType.
     *  \param parameterType Parameter for which partial is to be checked
     *  \return Direct partial derivatives, and associated times, of proper time components of Doppler partials
     */
    std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > getProperTimeParameterPartial(
            const estimatable_parameters::EstimatebleParameterIdentifier parameterType  );

private:

    //! Computed position scaling factor, for relative position vector (transmitter to receiver)
    Eigen::Matrix< double, 1, 3 > positionScalingFactor_;

    //! Factor by which light time correction state partial is to be scaled to be added to one-way Doppler partial
    /*!
     *  Factor by which light time correction state partial is to be scaled to be added to one-way Doppler partial. ALso forms
     *  part of the total partial w.r.t. the relative position.
     */
    double lightTimeEffectPositionScalingFactor_;

    //! Computed scaling factor for receiver velocity partials.
    Eigen::Matrix< double, 1, 3 > receiverVelocityScalingFactor_;

    //! Computed scaling factor for transmitter velocity partials.
    Eigen::Matrix< double, 1, 3 > transmitterVelocityScalingFactor_;

    //! Fixed link end for last computation of update() function.
    observation_models::LinkEndType currentLinkEndType_;

    //! Function returning the Cartesian acceleration of the transmitter as a function of time.
    std::function< Eigen::Vector3d( const double ) > transmitterAccelerationFunction_;

    //! Function returning the Cartesian acceleration of the receiver as a function of time.
    std::function< Eigen::Vector3d( const double ) > receiverAccelerationFunction_;

    //! Object used to compute the contribution of receiver proper time rate to the scaling
    std::shared_ptr< OneWayDopplerProperTimeComponentScaling > transmitterProperTimePartials_;

    //! Object used to compute the contribution of transmitter proper time rate to the scaling
    std::shared_ptr< OneWayDopplerProperTimeComponentScaling > receiverProperTimePartials_;
};

//! Function to computed the derivative of the unit vector from transmitter to receiver w.r.t. the observation time
/*!
 * Function to computed the derivative of the unit vector from transmitter to receiver w.r.t. the observation time
 * \param vectorToReceiver Vector from transmitter to receiver (unnormalized).
 * \param unitVectorToReceiver Vector from transmitter to receiver (normalized).
 * \param linkEndDistance Distance between transmitter and receiver
 * \param linkEndVelocity Velocity of link end at which the time is varied
 * \return Derivative of the unit vector from transmitter to receiver w.r.t. the observation time
 */
Eigen::Vector3d computePartialOfUnitVectorWrtLinkEndTime(
        const Eigen::Vector3d& vectorToReceiver,
        const Eigen::Vector3d& unitVectorToReceiver,
        const double linkEndDistance,
        const Eigen::Vector3d linkEndVelocity );


//! Function to computed the derivative of velocity component along line-of-sight vector w.r.t. the observation time
/*!
 * Function to computed the derivative of velocity component along line-of-sight vector w.r.t. the observation time
 * \param vectorToReceiver Vector from transmitter to receiver (unnormalized).
 * \param projectedLinkEndVelocity Velocity vector of link end, projected along line-of-sight vector
 * \param variableLinkEndVelocity Velocity of link end at which the time is varied
 * \param projectedLinkEndAcceleration Acceleration vector of link end, projected along line-of-sight vector
 * \param linkEndIsReceiver Boolean denoting whether the link end at which partial is taken at receiver (if true) or transmitter
 * (if false)
 * \param projectedLinkEndIsVariableLinkEnd Boolean denoting whether the link end at which partial is the reference link end
 * for the observation (if false) or not (if true)
 * \return Derivative of velocity component along line-of-sight vector w.r.t. the observation time
 */
double computePartialOfProjectedLinkEndVelocityWrtAssociatedTime(
        const Eigen::Vector3d& vectorToReceiver,
        const Eigen::Vector3d& projectedLinkEndVelocity,
        const Eigen::Vector3d& variableLinkEndVelocity,
        const Eigen::Vector3d& projectedLinkEndAcceleration,
        const bool linkEndIsReceiver,
        const bool projectedLinkEndIsVariableLinkEnd = true );

//! Class to compute the partial derivatives of a one-way doppler observation partial.
class OneWayDopplerPartial: public ObservationPartial< 1 >
{

public:

    typedef std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > OneWayDopplerPartialReturnType;
    typedef std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > SingleOneWayDopplerPartialReturnType;

    //! Constructor
    /*!
     * Constructor
     * \param oneWayDopplerScaler Scaling object used for mapping partials of positions to partials of observable
     * \param positionPartialList List of position partials per link end.
     * \param parameterIdentifier Id of parameter for which instance of class computes partial derivatives.
     * \param lighTimeCorrectionPartials List if light-time correction partial objects.
     */
    OneWayDopplerPartial(
            const std::shared_ptr< OneWayDopplerScaling > oneWayDopplerScaler,
            const std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > >& positionPartialList,
            const estimatable_parameters::EstimatebleParameterIdentifier parameterIdentifier,
            const std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > >&
            lighTimeCorrectionPartials =
            std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > >( ) );

    //! Destructor.
    ~OneWayDopplerPartial( ) { }

    //! Function to calculate the observation partial(s) at required time and state
    /*!
     *  Function to calculate the observation partial(s) at required time and state. State and time
     *  are typically obtained from evaluation of observation model.
     *  \param states Link end states. Index maps to link end for a given ObsevableType through getLinkEndIndex function.
     *  \param times Link end time.
     *  \param linkEndOfFixedTime Link end that is kept fixed when computing the observable.
     *  \param currentObservation Value of the observation for which the partial is to be computed (default NaN for
     *  compatibility purposes).
     *  \return Vector of pairs containing partial values and associated times.
     */
    OneWayDopplerPartialReturnType calculatePartial(
            const std::vector< Eigen::Vector6d >& states,
            const std::vector< double >& times,
            const observation_models::LinkEndType linkEndOfFixedTime,
            const Eigen::Vector1d& currentObservation = Eigen::Vector1d::Constant( TUDAT_NAN ) );

    //! Function to get scaling object used for mapping partials of positions to partials of observable
    /*!
     * Function to get scaling object used for mapping partials of positions to partials of observable
     * \return
     */
    std::shared_ptr< OneWayDopplerScaling > getOneWayDopplerScaler( )
    {
        return oneWayDopplerScaler_;
    }

    //! Function to get the number of light-time correction partial functions.
    /*!
     * Number of light-time correction partial functions.
     * \return Number of light-time correction partial functions.
     */
    int getNumberOfLighTimeCorrectionPartialsFunctions( )
    {
        return lighTimeCorrectionPartialsFunctions_.size( );
    }

protected:

    //! Scaling object used for mapping partials of positions to partials of observable
    std::shared_ptr< OneWayDopplerScaling > oneWayDopplerScaler_;

    //! List of position partials per link end.
    std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > > positionPartialList_;

    //! Iterator over list of position partials per link end.
    std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > >::iterator positionPartialIterator_;

    //! List of light-time correction partial functions.
    std::vector< std::function< SingleOneWayDopplerPartialReturnType(
            const std::vector< Eigen::Vector6d >&, const std::vector< double >& ) > >
    lighTimeCorrectionPartialsFunctions_;

    //! List of light-time correction partial objects.
    std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > > lighTimeCorrectionPartials_;

    //! Pre-declared state variable to be used in calculatePartial function.
    Eigen::Vector6d currentState_;

    //! Pre-declared time variable to be used in calculatePartial function.
    double currentTime_;

    //! Boolean to denote whether partials of proper time w.r.t. parameters should be added to final partial.
    bool addProperTimeParameterPartials_;

};

}

}

#endif // TUDAT_ONEWAYDOPPLERPARTIAL_H
