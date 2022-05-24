/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_ESTIMATABLEPARAMETERSETTINGS_H
#define TUDAT_ESTIMATABLEPARAMETERSETTINGS_H

#include "tudat/astro/basic_astro/empiricalAcceleration.h"
#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/astro/observation_models/linkTypeDefs.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"

namespace tudat
{

namespace estimatable_parameters
{

//! Class for providing settings for parameters to estimate
/*!
 *  Class for providing settings for parameters to estimate. This class is a functional base class for parameters that
 *  require no information in addition to their type.
 *  This class can be used for the easy setup of parameter objects (see createEstimatableParameters.h), but users may also
 *  chose to do so manually. (Derived) Class members are all public, for ease of access and modification.
 */
class EstimatableParameterSettings
{
public:

    //! Constructor.
    /*!
     *  Constructor, takes parameter type and body of which it is a property.
     *  \param associatedBody Body of which parameter is a property.
     *  \param parameterType Type of parameter.
     *  \param pointOnBodyId Reference point on body associated with parameter (empty by default).
     */
    EstimatableParameterSettings( const std::string associatedBody ,
                                  const EstimatebleParametersEnum parameterType,
                                  const std::string pointOnBodyId = "" ):
        parameterType_( std::make_pair( parameterType, std::make_pair( associatedBody, pointOnBodyId ) ) ){ }

    //! Virtual destructor
    /*!
     *  Virtual destructor
     */
    virtual ~EstimatableParameterSettings( ){ }


    //! Identifier for parameter.
    /*!
     *  Identifier for parameter, contains type of parameter and body of which parameter is a property.
     */
    EstimatebleParameterIdentifier parameterType_;

};


//! Function to retrieve full list of degree/order of spherical harmonic coeficients for given range of degrees and orders
/*!
 * Function to retrieve full list of degree/order of spherical harmonic coeficients for given range of degrees and orders
 * \param minimumDegree Minimum degree of field
 * \param minimumOrder Minimum order of field
 * \param maximumDegree Maximum degree of field
 * \param maximumOrder Maximum order of field
 * \return List of paird containing (degree,order) of field.
 */
inline std::vector< std::pair< int, int > > getSphericalHarmonicBlockIndices(
        const int minimumDegree,
        const int minimumOrder,
        const int maximumDegree,
        const int maximumOrder )
{
    std::vector< std::pair< int, int > > blockIndices;
    for( int i = minimumDegree; i <= maximumDegree; i++ )
    {
        for( int j = minimumOrder; ( ( j <= i ) && ( j <= maximumOrder ) ); j++ )
        {
            blockIndices.push_back( std::make_pair( i, j ) );
        }
    }
    return blockIndices;
}

//! Class for providing settings spherical harmonic coefficient(s) parameter
class SphericalHarmonicEstimatableParameterSettings: public EstimatableParameterSettings
{
public:

    //! Constructor
    /*!
     * Constructor, used to set user-defined list of degrees and orders.
     * \param blockIndices List of degrees and orders that are to estimated (first and second of each entry are
     * degree and order).
     * \param associatedBody Body for which coefficients are to be estimated.
     * \param parameterType Type of parameter that is to be estimated (must be spherical_harmonics_cosine_coefficient_block
     * of spherical_harmonics_sine_coefficient_block).
     */
    SphericalHarmonicEstimatableParameterSettings( const std::vector< std::pair< int, int > > blockIndices,
                                                   const std::string associatedBody,
                                                   const EstimatebleParametersEnum parameterType ):
        EstimatableParameterSettings( associatedBody, parameterType ), blockIndices_( blockIndices ),
        minimumDegree_( -1 ), minimumOrder_( -1 ), maximumDegree_( -1 ), maximumOrder_( -1 )
    {
        if( ( parameterType != spherical_harmonics_cosine_coefficient_block ) &&
                ( parameterType != spherical_harmonics_sine_coefficient_block ) )
        {

            throw std::runtime_error(
                        "Error when making spherical harmonic parameter settings, input parameter type is inconsistent." );
        }
    }

    //! Constructor
    /*!
     * Constructor, used to set a full block of degrees and orders that are to be estimated
     * \param minimumDegree Minimum degree of field that is to be estimated.
     * \param minimumOrder Minimum order of field that is to be estimated.
     * \param maximumDegree Maximum degree of field that is to be estimated.
     * \param maximumOrder Maximum order of field that is to be estimated.
     * \param associatedBody Body for which  coefficients are to be estimated.
     * \param parameterType Type of parameter that is to be estimated (must be spherical_harmonics_cosine_coefficient_block
     * of spherical_harmonics_sine_coefficient_block).
     */
    SphericalHarmonicEstimatableParameterSettings( const int minimumDegree,
                                                   const int minimumOrder,
                                                   const int maximumDegree,
                                                   const int maximumOrder,
                                                   const std::string associatedBody,
                                                   const EstimatebleParametersEnum parameterType ):
        EstimatableParameterSettings( associatedBody, parameterType ),
        minimumDegree_( minimumDegree ), minimumOrder_( minimumOrder ), maximumDegree_( maximumDegree ),
        maximumOrder_( maximumOrder )
    {
        if( ( parameterType != spherical_harmonics_cosine_coefficient_block ) &&
                ( parameterType != spherical_harmonics_sine_coefficient_block ) )
        {
            throw std::runtime_error( "Error when making spherical harmonic parameter settings, input parameter type is inconsistent." );
        }

        blockIndices_ = getSphericalHarmonicBlockIndices( minimumDegree, minimumOrder, maximumDegree, maximumOrder );
    }

    //! List of degrees and orders that are to estimated (first and second of each entry are degree and order.
    std::vector< std::pair< int, int > > blockIndices_;

    //!  Minimum degree of field that is to be estimated.
    int minimumDegree_;

    //! Minimum order of field that is to be estimated.
    int minimumOrder_;

    //! Maximum degree of field that is to be estimated.
    int maximumDegree_;

    //! Maximum order of field that is to be estimated.
    int maximumOrder_;
};

//! Class to define settings for estimation of constant observation biases (absolute or relative)
class ConstantObservationBiasEstimatableParameterSettings: public EstimatableParameterSettings
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param linkEnds Observation link ends for which the bias is to be estimated.
     * \param observableType Observable type for which the bias is to be estimated.
     * \param isBiasAdditive True if bias is absolute, false if it is relative
     */
    ConstantObservationBiasEstimatableParameterSettings(
            const observation_models::LinkEnds& linkEnds,
            const observation_models::ObservableType observableType,
            const bool isBiasAdditive ):
        EstimatableParameterSettings(
            linkEnds.begin( )->second.first,
            isBiasAdditive ? constant_additive_observation_bias : constant_relative_observation_bias,
            linkEnds.begin( )->second.second ), linkEnds_( linkEnds ), observableType_( observableType ){ }

    //! Destructor
    ~ConstantObservationBiasEstimatableParameterSettings( ){ }

    //! Observation link ends for which the bias is to be estimated.
    observation_models::LinkEnds linkEnds_;

    //! Observable type for which the bias is to be estimated.
    observation_models::ObservableType observableType_;
};

//! Class to define settings for estimation of arc-wise constant observation biases (absolute or relative)
class ArcWiseConstantObservationBiasEstimatableParameterSettings: public EstimatableParameterSettings
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param linkEnds Observation link ends for which the bias is to be estimated.
     * \param observableType Observable type for which the bias is to be estimated.
     * \param arcStartTimes Start times for arcs in which biases are defined
     * \param linkEndForTime Link end index from which the 'current time' is determined
     * \param isBiasAdditive True if bias is absolute, false if it is relative
     */
    ArcWiseConstantObservationBiasEstimatableParameterSettings(
            const observation_models::LinkEnds& linkEnds,
            const observation_models::ObservableType observableType,
            const std::vector< double > arcStartTimes,
            const observation_models::LinkEndType linkEndForTime,
            const bool isBiasAdditive ):
        EstimatableParameterSettings(
            linkEnds.begin( )->second.first,
            isBiasAdditive ? arcwise_constant_additive_observation_bias : arcwise_constant_relative_observation_bias,
            linkEnds.begin( )->second.second ), linkEnds_( linkEnds ), observableType_( observableType ),
        arcStartTimes_( arcStartTimes ), linkEndForTime_( linkEndForTime ){ }

    //! Destructor
    ~ArcWiseConstantObservationBiasEstimatableParameterSettings( ){ }

    //! Observation link ends for which the bias is to be estimated.
    observation_models::LinkEnds linkEnds_;

    //! Observable type for which the bias is to be estimated.
    observation_models::ObservableType observableType_;

    //! Start times for arcs in which biases are defined
    std::vector< double > arcStartTimes_;

    //! Link end index from which the 'current time' is determined
    observation_models::LinkEndType linkEndForTime_;

};

//! Class to define settings for estimating an initial translational state.
template< typename InitialStateParameterType = double >
class InitialTranslationalStateEstimatableParameterSettings: public EstimatableParameterSettings
{
public:

    //! Constructor, sets initial value of translational state.
    /*!
     * Constructor, sets initial value of translational state.
     * \param associatedBody Body for which initial state is to be estimated.
     * \param initialStateValue Current value of initial state (w.r.t. centralBody)
     * \param centralBody Body w.r.t. which the initial state is to be estimated.
     * \param frameOrientation Orientation of the frame in which the state is defined.
     */
    InitialTranslationalStateEstimatableParameterSettings(
            const std::string& associatedBody,
            const Eigen::Matrix< InitialStateParameterType, 6, 1 > initialStateValue,
            const std::string& centralBody = "SSB", const std::string& frameOrientation = "ECLIPJ2000" ):
        EstimatableParameterSettings( associatedBody, initial_body_state ), initialTime_( TUDAT_NAN ),
        initialStateValue_( initialStateValue ),
        centralBody_( centralBody ), frameOrientation_( frameOrientation ){ }

    //! Constructor, without initial value of translational state.
    /*!
     * Constructor, without initial value of translational state. Current initial state is retrieved from environment
     * (ephemeris objects) during creation of parameter object.
     * \param associatedBody Body for which initial state is to be estimated.
     * \param initialTime Time at which initial state is defined.
     * \param centralBody Body w.r.t. which the initial state is to be estimated.
     * \param frameOrientation Orientation of the frame in which the state is defined.
     */
    InitialTranslationalStateEstimatableParameterSettings(
            const std::string& associatedBody,
            const double initialTime,
            const std::string& centralBody = "SSB", const std::string& frameOrientation = "ECLIPJ2000" ):
        EstimatableParameterSettings( associatedBody, initial_body_state ), initialTime_( initialTime ),
        centralBody_( centralBody ), frameOrientation_( frameOrientation ){ }

    //! Time at which initial state is defined (NaN for user-defined initial state value).
    double initialTime_;

    //! Current value of initial state (w.r.t. centralBody), set manually by used.
    Eigen::Matrix< InitialStateParameterType, 6, 1 > initialStateValue_;

    //! Body w.r.t. which the initial state is to be estimated.
    std::string centralBody_;

    //! Orientation of the frame in which the state is defined.
    std::string frameOrientation_;

};

//! Class to define settings for estimating an arcwise initial translational state.
template< typename InitialStateParameterType >
class ArcWiseInitialTranslationalStateEstimatableParameterSettings: public EstimatableParameterSettings
{
public:


    //! Constructor, sets initial value of translational state and a single central body.
    /*!
     * Constructor, sets initial value of translational state and a single central body
     * \param associatedBody Body for which initial state is to be estimated.
     * \param initialStateValue Current value of initial arc states (concatenated in same order as arcs)
     * \param arcStartTimes Start times for separate arcs
     * \param centralBody Body w.r.t. which the initial state is to be estimated.
     * \param frameOrientation Orientation of the frame in which the state is defined.
     */
    ArcWiseInitialTranslationalStateEstimatableParameterSettings(
            const std::string& associatedBody,
            const Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > initialStateValue,
            const std::vector< double >& arcStartTimes,
            const std::string& centralBody = "SSB", const std::string& frameOrientation = "ECLIPJ2000" ):
        EstimatableParameterSettings( associatedBody, arc_wise_initial_body_state ), initialStateValue_( initialStateValue ),
        arcStartTimes_( arcStartTimes ), frameOrientation_( frameOrientation ),
        isStateSet_( 1 )
    {
        for( unsigned int i = 0; i < arcStartTimes.size( ); i++ )
        {
            centralBodies_.push_back( centralBody );
        }
    }

    //! Constructor, sets initial value of translational state and an arc-wise variable central body.
    /*!
     * Constructor, sets initial value of translational state and an arc-wise variable central body
     * \param associatedBody Body for which initial state is to be estimated.
     * \param initialStateValue Current value of initial arc states (concatenated in same order as arcs)
     * \param arcStartTimes Start times for separate arcs
     * \param centralBodies List of central bodies (per arc) w.r.t. which the initial state is to be estimated.
     * \param frameOrientation Orientation of the frame in which the state is defined.
     */
    ArcWiseInitialTranslationalStateEstimatableParameterSettings(
            const std::string& associatedBody,
            const Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > initialStateValue,
            const std::vector< double >& arcStartTimes,
            const std::vector< std::string > centralBodies, const std::string& frameOrientation = "ECLIPJ2000" ):
        EstimatableParameterSettings( associatedBody, arc_wise_initial_body_state ), initialStateValue_( initialStateValue ),
        arcStartTimes_( arcStartTimes ), centralBodies_( centralBodies ), frameOrientation_( frameOrientation ),
        isStateSet_( 1 ){ }

    //! Constructor, without initial value of translational state, for a single central body
    /*!
     * Constructor, without initial value of translational state. Current initial state is retrieved from environment
     * (ephemeris objects) during creation of parameter object.
     * \param associatedBody Body for which initial state is to be estimated.
     * \param arcStartTimes Start times for separate arcs
     * \param centralBody Body w.r.t. which the initial state is to be estimated.
     * \param frameOrientation Orientation of the frame in which the state is defined.
     */
    ArcWiseInitialTranslationalStateEstimatableParameterSettings(
            const std::string& associatedBody,
            const std::vector< double >& arcStartTimes,
            const std::string& centralBody = "SSB", const std::string& frameOrientation = "ECLIPJ2000" ):
        EstimatableParameterSettings( associatedBody, arc_wise_initial_body_state ),
        arcStartTimes_( arcStartTimes ), frameOrientation_( frameOrientation ),
        isStateSet_( 0 )
    {
        for( unsigned int i = 0; i < arcStartTimes.size( ); i++ )
        {
            centralBodies_.push_back( centralBody );
        }
    }

    //! Constructor, without initial value of translational state, for an arc-wise variable central body
    /*!
     * Constructor, without initial value of translational state, for an arc-wise variable central body.
     * Current initial state is retrieved from environment
     * (ephemeris objects) during creation of parameter object.
     * \param associatedBody Body for which initial state is to be estimated.
     * \param arcStartTimes Start times for separate arcs
     * \param centralBodies List of central bodies (per arc) w.r.t. which the initial state is to be estimated.
     * \param frameOrientation Orientation of the frame in which the state is defined.
     */
    ArcWiseInitialTranslationalStateEstimatableParameterSettings(
            const std::string& associatedBody,
            const std::vector< double >& arcStartTimes,
            const std::vector< std::string > centralBodies, const std::string& frameOrientation = "ECLIPJ2000" ):
        EstimatableParameterSettings( associatedBody, arc_wise_initial_body_state ),
        arcStartTimes_( arcStartTimes ), centralBodies_( centralBodies ), frameOrientation_( frameOrientation ),
        isStateSet_( 0 ){ }

    //! Current value of initial arc states (concatenated in same order as arcs)
    Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > initialStateValue_;

    //! Start times for separate arcs
    std::vector< double > arcStartTimes_;

    //! List of bodies (arc-wise) w.r.t. which the initial state is to be estimated.
    std::vector< std::string > centralBodies_;

    //!Orientation of the frame in which the state is defined.
    std::string frameOrientation_;

    //! Boolean to denote whether initial states are set, or if they need to be computed
    bool isStateSet_;

};

//! Class to define settings for estimating an initial rotational state.
template< typename InitialStateParameterType = double >
class InitialRotationalStateEstimatableParameterSettings: public EstimatableParameterSettings
{
public:

    //! Constructor, sets initial value of rotational state.
    /*!
     * Constructor, sets initial value of rotational state.
     * \param associatedBody Body for which initial state is to be estimated.
     * \param initialStateValue Current value of initial state (w.r.t. baseOrientation)
     * \param baseOrientation Orientation w.r.t. which the initial state is to be estimated.
     * \param frameOrientation Orientation of the frame in which the state is defined.
     */
    InitialRotationalStateEstimatableParameterSettings(
            const std::string& associatedBody,
            const Eigen::Matrix< InitialStateParameterType, 7, 1 > initialStateValue,
            const std::string& baseOrientation = "SSB" ):
        EstimatableParameterSettings( associatedBody, initial_rotational_body_state ), initialTime_( TUDAT_NAN ),
        initialStateValue_( initialStateValue ),
        baseOrientation_( baseOrientation ){ }

    //! Constructor, without initial value of rotational state.
    /*!
     * Constructor, without initial value of rotational state. Current initial state is retrieved from environment
     * (ephemeris objects) during creation of parameter object.
     * \param associatedBody Body for which initial state is to be estimated.
     * \param initialTime Time at which initial state is defined.
     * \param baseOrientation Orientation w.r.t. which the initial state is to be estimated.
     * \param frameOrientation Orientation of the frame in which the state is defined.
     */
    InitialRotationalStateEstimatableParameterSettings(
            const std::string& associatedBody,
            const double initialTime,
            const std::string& baseOrientation = "SSB"):
        EstimatableParameterSettings( associatedBody, initial_rotational_body_state ), initialTime_( initialTime ),
        baseOrientation_( baseOrientation ){ }

    //! Time at which initial state is defined (NaN for user-defined initial state value).
    double initialTime_;

    //! Current value of initial state (w.r.t. baseOrientation), set manually by used.
    Eigen::Matrix< InitialStateParameterType, 7, 1 > initialStateValue_;

    //! Orientation w.r.t. which the initial state is to be estimated.
    std::string baseOrientation_;

};

//! Class to define settings for estimating an initial rotational state.
template< typename InitialStateParameterType = double >
class InitialMassEstimatableParameterSettings: public EstimatableParameterSettings
{
public:

    InitialMassEstimatableParameterSettings(
            const std::string& associatedBody,
            const double initialStateValue ):
        EstimatableParameterSettings( associatedBody, initial_mass_state ),
        initialStateValue_( initialStateValue ){ }

    double initialStateValue_;
};


//! Class to define settings for estimating time-independent empirical acceleration components
class EmpiricalAccelerationEstimatableParameterSettings: public EstimatableParameterSettings
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param associatedBody Name of body undergoing acceleration
     * \param centralBody Name of central body
     * \param componentsToEstimate List of components of empirical acceleration that are to be estimated.
     */
    EmpiricalAccelerationEstimatableParameterSettings(
            const std::string associatedBody,
            const std::string centralBody,
            const std::map< basic_astrodynamics::EmpiricalAccelerationComponents,
            std::vector< basic_astrodynamics::EmpiricalAccelerationFunctionalShapes > > componentsToEstimate ):
        EstimatableParameterSettings( associatedBody, empirical_acceleration_coefficients, centralBody ),
        componentsToEstimate_( componentsToEstimate ){ }

    //!  List of components of empirical acceleration that are to be estimated.
    std::map< basic_astrodynamics::EmpiricalAccelerationComponents,
    std::vector< basic_astrodynamics::EmpiricalAccelerationFunctionalShapes > > componentsToEstimate_;

};

//! Class to define settings for estimating time-dependent (arcwise constant) empirical acceleration components
class ArcWiseEmpiricalAccelerationEstimatableParameterSettings: public EstimatableParameterSettings
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param associatedBody Name of body undergoing acceleration
     * \param centralBody Name of central body
     * \param componentsToEstimate List of components of empirical acceleration that are to be estimated.
     * \param arcStartTimeList List of times at which empirical acceleration arcs are to start
     */
    ArcWiseEmpiricalAccelerationEstimatableParameterSettings(
            const std::string associatedBody,
            const std::string centralBody,
            const std::map< basic_astrodynamics::EmpiricalAccelerationComponents,
            std::vector< basic_astrodynamics::EmpiricalAccelerationFunctionalShapes > > componentsToEstimate,
            const std::vector< double > arcStartTimeList ):
        EstimatableParameterSettings( associatedBody, arc_wise_empirical_acceleration_coefficients, centralBody ),
        componentsToEstimate_( componentsToEstimate ), arcStartTimeList_( arcStartTimeList ){ }

    //! List of components of empirical acceleration that are to be estimated.
    std::map< basic_astrodynamics::EmpiricalAccelerationComponents,
    std::vector< basic_astrodynamics::EmpiricalAccelerationFunctionalShapes > > componentsToEstimate_;

    //! List of times at which empirical acceleration arcs are to start
    std::vector< double > arcStartTimeList_;


};

//! Class to define settings for estimating time-dependent (arcwise constant) empirical acceleration components
class ArcWiseRadiationPressureCoefficientEstimatableParameterSettings: public EstimatableParameterSettings
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param associatedBody Name of body undergoing acceleration
     * \param arcStartTimeList List of times at which radiation pressure coefficient arcs are to start
     */
    ArcWiseRadiationPressureCoefficientEstimatableParameterSettings(
            const std::string associatedBody,
            const std::vector< double > arcStartTimeList ):
        EstimatableParameterSettings( associatedBody, arc_wise_radiation_pressure_coefficient ),
        arcStartTimeList_( arcStartTimeList ){ }

    //! List of times at which radiation pressure coefficient arcs are to start
    std::vector< double > arcStartTimeList_;


};

//! Class to define settings for estimating time-dependent (arcwise constant) drag coefficients
class ArcWiseDragCoefficientEstimatableParameterSettings: public EstimatableParameterSettings
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param associatedBody Name of body undergoing acceleration
     * \param arcStartTimeList List of times at which drag coefficient arcs are to start
     */
    ArcWiseDragCoefficientEstimatableParameterSettings(
            const std::string associatedBody,
            const std::vector< double > arcStartTimeList):
        EstimatableParameterSettings( associatedBody, arc_wise_constant_drag_coefficient ),
        arcStartTimeList_( arcStartTimeList ){ }

    //! List of times at which drag coefficient arcs are to start
    std::vector< double > arcStartTimeList_;


};

//! Class to define settings for estimating a Tidal Love number (k_{n}) at a single degree that is constant for all orders
/*!
 *  Class to define settings for estimating a Tidal Love number (k_{n}) at a single degree that is constant for all orders.
 *  Either a real or a complex Love number may be estimated (represented by entries of a VectorXd).
 *  The constructor argument representing the deforming body/bodies must correspond exactly to the deforming bodies in a
 *  BasicSolidBodyTideGravityFieldVariations member object of the deformed body. Alternatively, if only one
 *  BasicSolidBodyTideGravityFieldVariations object is present, the deforming body list may be left empty.
 */
class FullDegreeTidalLoveNumberEstimatableParameterSettings: public EstimatableParameterSettings
{
public:

    //! Constructor for a single deforming body
    /*!
     * Constructor for a single deforming body
     * \param associatedBody Deformed body
     * \param degree Degree of Love number that is to be estimated
     * \param deformingBody Name of body causing tidal deformation
     * \param useComplexValue True if the complex Love number is estimated, false if only the real part is considered
     */
    FullDegreeTidalLoveNumberEstimatableParameterSettings(  const std::string& associatedBody,
                                                            const int degree,
                                                            const std::string deformingBody,
                                                            const bool useComplexValue = 0 ):
        EstimatableParameterSettings( associatedBody, full_degree_tidal_love_number ), degree_( degree ),
        useComplexValue_( useComplexValue )
    {
        if( deformingBody != "" )
        {
            deformingBodies_.push_back( deformingBody );
        }
    }

    //! Constructor for a list of deforming bodyies
    /*!
     * Constructor for a list of deforming bodyies
     * \param associatedBody Deformed body
     * \param degree Degree of Love number that is to be estimated
     * \param deformingBodies Names of bodies causing tidal deformation
     * \param useComplexValue True if the complex Love number is estimated, false if only the real part is considered
     */
    FullDegreeTidalLoveNumberEstimatableParameterSettings(  const std::string& associatedBody,
                                                            const int degree ,
                                                            const std::vector< std::string >& deformingBodies,
                                                            const bool useComplexValue = 0 ):
        EstimatableParameterSettings( associatedBody, full_degree_tidal_love_number ), degree_( degree ),
        deformingBodies_( deformingBodies ), useComplexValue_( useComplexValue ){ }

    //! Degree of Love number that is to be estimated
    int degree_;

    //! Names of bodies causing tidal deformation
    std::vector< std::string > deformingBodies_;

    //! True if the complex Love number is estimated, false if only the real part is considered
    bool useComplexValue_;

};

//! Class to define settings for estimating a set of Tidal Love number (k_{n,m}) at a single degree.
/*!
 *  Class to define settings for estimating a set of Tidal Love number (k_{n,m}) at a single degree and a set of orders at this
 *  degree. The estimation will provide separate Love numbers for each order
 *  Either a real or a complex Love number may be estimated (represented by entries of a VectorXd).
 *  The constructor argument representing the deforming body/bodies must correspond exactly to the deforming bodies in a
 *  BasicSolidBodyTideGravityFieldVariations member object of the deformed body. Alternatively, if only one
 *  BasicSolidBodyTideGravityFieldVariations object is present, the deforming body list may be left empty.
 */
class SingleDegreeVariableTidalLoveNumberEstimatableParameterSettings: public EstimatableParameterSettings
{
public:

    //! Constructor for a list of deforming bodyies
    /*!
     * Constructor for a list of deforming bodyies
     * \param associatedBody Deformed body
     * \param degree Degree of Love number that is to be estimated
     * \param orders List of orders at which Love numbers are to be estimated.
     * \param deformingBody Names of body causing tidal deformation
     * \param useComplexValue True if the complex Love number is estimated, false if only the real part is considered
     */
    SingleDegreeVariableTidalLoveNumberEstimatableParameterSettings(  const std::string associatedBody,
                                                                      const int degree,
                                                                      const std::vector< int > orders,
                                                                      const std::string& deformingBody,
                                                                      const bool useComplexValue = 0 ):
        EstimatableParameterSettings( associatedBody, single_degree_variable_tidal_love_number ), degree_( degree ),
        orders_( orders ), useComplexValue_( useComplexValue )
    {
        if( deformingBody != "" )
        {
            deformingBodies_.push_back( deformingBody );
        }
    }

    //! Constructor for a list of deforming bodyies
    /*!
     * Constructor for a list of deforming bodyies
     * \param associatedBody Deformed body
     * \param degree Degree of Love number that is to be estimated
     * \param orders List of orders at which Love numbers are to be estimated.
     * \param deformingBodies Names of bodies causing tidal deformation
     * \param useComplexValue True if the complex Love number is estimated, false if only the real part is considered
     */
    SingleDegreeVariableTidalLoveNumberEstimatableParameterSettings(  const std::string associatedBody,
                                                                      const int degree,
                                                                      const std::vector< int > orders,
                                                                      const std::vector< std::string >& deformingBodies,
                                                                      const bool useComplexValue = 0 ):
        EstimatableParameterSettings( associatedBody, single_degree_variable_tidal_love_number ), degree_( degree ),
        orders_( orders ), deformingBodies_( deformingBodies ), useComplexValue_( useComplexValue ){ }

    //! Degree of Love number that is to be estimated
    int degree_;

    //! List of orders at which Love numbers are to be estimated.
    const std::vector< int > orders_;

    //! Names of bodies causing tidal deformation
    std::vector< std::string > deformingBodies_;

    //! True if the complex Love number is estimated, false if only the real part is considered
    bool useComplexValue_;

};

//! Class to define settings for estimating the tidal time lag of a direct tidal acceleration model
/*!
 *  Class to define settings for estimating the tidal time lag of a direct tidal acceleration model, it links to one or more
 *  objects of type DirectTidalDissipationAcceleration. The user can provide a list of bodies cause deformation, and the
 *  associated DirectTidalDissipationAcceleration objects will be used. If the list of bodies causing deformation is left empty,
 *  all DirectTidalDissipationAcceleration objects for the given body undergoing deformation will be used
 */
class DirectTidalTimeLagEstimatableParameterSettings: public EstimatableParameterSettings
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param associatedBody Body being deformed
     * \param deformingBody Body causing deformed
     */
    DirectTidalTimeLagEstimatableParameterSettings( const std::string& associatedBody,
                                                    const std::string& deformingBody ):
        EstimatableParameterSettings( associatedBody, direct_dissipation_tidal_time_lag )
    {
        if( deformingBody != "" )
        {
            deformingBodies_.push_back( deformingBody );
        }
    }

    //! Constructor
    /*!
     * Constructor
     * \param associatedBody Body being deformed
     * \param deformingBodies Names of bodies causing tidal deformation
     */
    DirectTidalTimeLagEstimatableParameterSettings( const std::string& associatedBody,
                                                    const std::vector< std::string >& deformingBodies ):
        EstimatableParameterSettings( associatedBody, direct_dissipation_tidal_time_lag ),
        deformingBodies_( deformingBodies ){ }


    //! Names of bodies causing tidal deformation
    std::vector< std::string > deformingBodies_;

};


inline std::shared_ptr< EstimatableParameterSettings > gravitationalParameter( const std::string bodyName )
{
    return std::make_shared< EstimatableParameterSettings >( bodyName, gravitational_parameter );
}

inline std::shared_ptr< EstimatableParameterSettings > constantDragCoefficient( const std::string bodyName )
{
    return std::make_shared< EstimatableParameterSettings >( bodyName, constant_drag_coefficient );
}

inline std::shared_ptr< EstimatableParameterSettings > radiationPressureCoefficient( const std::string bodyName )
{
    return std::make_shared< EstimatableParameterSettings >( bodyName, radiation_pressure_coefficient );
}

inline std::shared_ptr< EstimatableParameterSettings > sphericalHarmonicsCosineBlock(
        const std::string bodyName,
        const int minimumDegree,
        const int minimumOrder,
        const int maximumDegree,
        const int maximumOrder )
{
    return std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                minimumDegree, minimumOrder, maximumDegree, maximumOrder,
                bodyName, spherical_harmonics_cosine_coefficient_block );
}

inline std::shared_ptr< EstimatableParameterSettings > sphericalHarmonicsCosineBlock(
        const std::string bodyName,
        const std::vector< std::pair< int, int > > blockIndices )
{
    return std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                blockIndices,  bodyName, spherical_harmonics_cosine_coefficient_block );
}

inline std::shared_ptr< EstimatableParameterSettings > sphericalHarmonicsSineBlock(
        const std::string bodyName,
        const int minimumDegree,
        const int minimumOrder,
        const int maximumDegree,
        const int maximumOrder )
{
    return std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                minimumDegree, minimumOrder, maximumDegree, maximumOrder,
                bodyName, spherical_harmonics_sine_coefficient_block );
}

inline std::shared_ptr< EstimatableParameterSettings > sphericalHarmonicsSineBlock(
        const std::string bodyName,
        const std::vector< std::pair< int, int > > blockIndices )
{
    return std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                blockIndices,  bodyName, spherical_harmonics_sine_coefficient_block );
}


inline std::shared_ptr< EstimatableParameterSettings > arcwiseRadiationPressureCoefficient(
        std::string associatedBody,
        const std::vector< double > arcStartTimeList )
{
    return std::make_shared< ArcWiseRadiationPressureCoefficientEstimatableParameterSettings >(
                associatedBody, arcStartTimeList );
}

inline std::shared_ptr< EstimatableParameterSettings > arcwiseDragCoefficient(
        std::string associatedBody,
        const std::vector< double > arcStartTimeList )
{
    return std::make_shared< ArcWiseDragCoefficientEstimatableParameterSettings >(
                associatedBody, arcStartTimeList );
}

inline std::shared_ptr< EstimatableParameterSettings > constantRotationRate(
        std::string bodyName)
{
    return std::make_shared< EstimatableParameterSettings >( bodyName, constant_rotation_rate );
}

inline std::shared_ptr< EstimatableParameterSettings > rotationPolePosition(
        std::string bodyName)
{
    return std::make_shared< EstimatableParameterSettings >( bodyName, rotation_pole_position );
}

inline std::shared_ptr< EstimatableParameterSettings > observationBias(
        const observation_models::LinkEnds& linkEnds,
        const observation_models::ObservableType observableType )
{
    return std::make_shared< ConstantObservationBiasEstimatableParameterSettings >(
                linkEnds, observableType, true );
}

inline std::shared_ptr< EstimatableParameterSettings > relativeObservationBias(
        const observation_models::LinkEnds& linkEnds,
        const observation_models::ObservableType observableType )
{
    return std::make_shared< ConstantObservationBiasEstimatableParameterSettings >(
                linkEnds, observableType, false );
}

inline std::shared_ptr< EstimatableParameterSettings > arcwiseObservationBias(
        const observation_models::LinkEnds& linkEnds,
        const observation_models::ObservableType observableType,
        const std::vector< double > arcStartTimes,
        const observation_models::LinkEndType linkEndForTime = observation_models::receiver )
{
    return std::make_shared< ArcWiseConstantObservationBiasEstimatableParameterSettings >(
                linkEnds, observableType, arcStartTimes, linkEndForTime, true );
}

inline std::shared_ptr< EstimatableParameterSettings > arcwiseRelativeObservationBias(
        const observation_models::LinkEnds& linkEnds,
        const observation_models::ObservableType observableType,
        const std::vector< double > arcStartTimes,
        const observation_models::LinkEndType linkEndForTime = observation_models::receiver )
{
    return std::make_shared< ArcWiseConstantObservationBiasEstimatableParameterSettings >(
                linkEnds, observableType, arcStartTimes, linkEndForTime, false );
}


inline std::shared_ptr< EstimatableParameterSettings > constantEmpiricalAccelerationMagnitudes(
        const std::string associatedBody,
        const std::string centralBody )
{
    std::map< basic_astrodynamics::EmpiricalAccelerationComponents,
            std::vector< basic_astrodynamics::EmpiricalAccelerationFunctionalShapes > > componentsToEstimate;
    componentsToEstimate[ basic_astrodynamics::radial_empirical_acceleration_component ].push_back(
            basic_astrodynamics::constant_empirical );
    componentsToEstimate[ basic_astrodynamics::along_track_empirical_acceleration_component ].push_back(
            basic_astrodynamics::constant_empirical );
    componentsToEstimate[ basic_astrodynamics::across_track_empirical_acceleration_component ].push_back(
            basic_astrodynamics::constant_empirical );

    return std::make_shared< EmpiricalAccelerationEstimatableParameterSettings >(
                associatedBody, centralBody, componentsToEstimate );
}


inline std::shared_ptr< EstimatableParameterSettings > empiricalAccelerationMagnitudes(
        const std::string associatedBody,
        const std::string centralBody,
        const std::map< basic_astrodynamics::EmpiricalAccelerationComponents,
        std::vector< basic_astrodynamics::EmpiricalAccelerationFunctionalShapes > > componentsToEstimate  )
{
    return std::make_shared< EmpiricalAccelerationEstimatableParameterSettings >(
                associatedBody, centralBody, componentsToEstimate );
}

inline std::shared_ptr< EstimatableParameterSettings > constantArcWiseEmpiricalAccelerationMagnitudes(
        const std::string associatedBody,
        const std::string centralBody,
        const std::vector< double > arcStartTimes )
{
    std::map< basic_astrodynamics::EmpiricalAccelerationComponents,
            std::vector< basic_astrodynamics::EmpiricalAccelerationFunctionalShapes > > componentsToEstimate;
    componentsToEstimate[ basic_astrodynamics::radial_empirical_acceleration_component ].push_back(
            basic_astrodynamics::constant_empirical );
    componentsToEstimate[ basic_astrodynamics::along_track_empirical_acceleration_component ].push_back(
            basic_astrodynamics::constant_empirical );
    componentsToEstimate[ basic_astrodynamics::across_track_empirical_acceleration_component ].push_back(
            basic_astrodynamics::constant_empirical );

    return std::make_shared< ArcWiseEmpiricalAccelerationEstimatableParameterSettings >(
                associatedBody, centralBody, componentsToEstimate, arcStartTimes );
}


inline std::shared_ptr< EstimatableParameterSettings > arcWiseEmpiricalAccelerationMagnitudes(
        const std::string associatedBody,
        const std::string centralBody,
        const std::map< basic_astrodynamics::EmpiricalAccelerationComponents,
        std::vector< basic_astrodynamics::EmpiricalAccelerationFunctionalShapes > > componentsToEstimate,
        const std::vector< double > arcStartTimes )
{
    return std::make_shared< ArcWiseEmpiricalAccelerationEstimatableParameterSettings >(
                associatedBody, centralBody, componentsToEstimate, arcStartTimes );
}


inline std::shared_ptr< EstimatableParameterSettings > ppnParameterGamma( )
{
    return std::make_shared< EstimatableParameterSettings >( "", ppn_parameter_gamma );
}

inline std::shared_ptr< EstimatableParameterSettings > ppnParameterBeta( )
{
    return std::make_shared< EstimatableParameterSettings >( "", ppn_parameter_beta );
}

inline std::shared_ptr< EstimatableParameterSettings > groundStationPosition(
        const std::string& body,
        const std::string& groundStationName )
{
    return std::make_shared< EstimatableParameterSettings >( body, ground_station_position, groundStationName );
}

inline std::shared_ptr< EstimatableParameterSettings > directTidalDissipationLagTime(
        const std::string& body,
        const std::vector< std::string >& deformingBodies )
{
    return std::make_shared< DirectTidalTimeLagEstimatableParameterSettings >(
                body, deformingBodies);
}


inline std::shared_ptr< EstimatableParameterSettings > directTidalDissipationLagTime(
        const std::string& body,
        const std::string& deformingBody)
{
    return directTidalDissipationLagTime( body, std::vector< std::string >( { deformingBody } ) );
}

inline std::shared_ptr< EstimatableParameterSettings > meanMomentOfInertia(
        const std::string& body )
{
    return std::make_shared< EstimatableParameterSettings >( body, mean_moment_of_inertia );
}



inline std::shared_ptr< EstimatableParameterSettings > orderInvariantKLoveNumber(
        const std::string& associatedBody,
        const int degree,
        const std::string deformingBody,
        const bool useComplexValue = 0 )
{
    return std::make_shared< FullDegreeTidalLoveNumberEstimatableParameterSettings >(
                associatedBody, degree, deformingBody, useComplexValue );
}


inline std::shared_ptr< EstimatableParameterSettings > orderInvariantKLoveNumber(
        const std::string& associatedBody,
        const int degree,
        const std::vector< std::string >& deformingBodies,
        const bool useComplexValue = 0 )
{
    return std::make_shared< FullDegreeTidalLoveNumberEstimatableParameterSettings >(
                associatedBody, degree, deformingBodies, useComplexValue );
}


inline std::shared_ptr< EstimatableParameterSettings > orderInvariantKLoveNumber(
        const std::string& associatedBody,
        const int degree,
        const bool useComplexValue = 0 )
{
    return std::make_shared< FullDegreeTidalLoveNumberEstimatableParameterSettings >(
                associatedBody, degree, std::vector< std::string >( ), useComplexValue );
}

inline std::shared_ptr< EstimatableParameterSettings > orderVaryingKLoveNumber(
        const std::string& associatedBody,
        const int degree,
        const std::vector< int >& orders,
        const std::string deformingBody,
        const bool useComplexValue = 0 )
{
    return std::make_shared< SingleDegreeVariableTidalLoveNumberEstimatableParameterSettings >(
                associatedBody, degree, orders, deformingBody, useComplexValue );
}


inline std::shared_ptr< EstimatableParameterSettings > orderVaryingKLoveNumber(
        const std::string& associatedBody,
        const int degree,
        const std::vector< int >& orders,
        const std::vector< std::string >& deformingBodies,
        const bool useComplexValue = 0 )
{
    return std::make_shared< SingleDegreeVariableTidalLoveNumberEstimatableParameterSettings >(
                associatedBody, degree, orders, deformingBodies, useComplexValue );
}


inline std::shared_ptr< EstimatableParameterSettings > orderVaryingKLoveNumber(
        const std::string& associatedBody,
        const int degree,
        const std::vector< int >& orders,
        const bool useComplexValue = 0 )
{
    return std::make_shared< SingleDegreeVariableTidalLoveNumberEstimatableParameterSettings >(
                associatedBody, degree, orders, std::vector< std::string >( ), useComplexValue );
}


inline std::shared_ptr< EstimatableParameterSettings > coreFactor(
        const std::string& associatedBody )
{
    return std::make_shared< EstimatableParameterSettings >(
                associatedBody, core_factor );
}

inline std::shared_ptr< EstimatableParameterSettings > freeCoreNutationRate(
        const std::string& associatedBody )
{
    return std::make_shared< EstimatableParameterSettings >(
                associatedBody, free_core_nutation_rate );
}

inline std::shared_ptr< EstimatableParameterSettings > periodicSpinVariations(
        const std::string& associatedBody )
{
    return std::make_shared< EstimatableParameterSettings >(
                associatedBody, periodic_spin_variation );
}

inline std::shared_ptr< EstimatableParameterSettings > polarMotionAmplitudes(
        const std::string& associatedBody )
{
    return std::make_shared< EstimatableParameterSettings >(
                associatedBody, polar_motion_amplitude );
}

inline std::shared_ptr< EstimatableParameterSettings > quasiImpulsiveShots(
        const std::string& associatedBody )
{
    return std::make_shared< EstimatableParameterSettings >(
                associatedBody, desaturation_delta_v_values );
}

} // namespace estimatable_parameters

} // namespace tudat

#endif // TUDAT_ESTIMATABLEPARAMETERSETTINGS_H
