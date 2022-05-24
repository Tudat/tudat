/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_OBSERVATIONBIAS_H
#define TUDAT_OBSERVATIONBIAS_H

#include <vector>
#include <iostream>

#include <memory>
#include <boost/make_shared.hpp>
#include <functional>

#include <Eigen/Core>

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/basics/basicTypedefs.h"
#include "tudat/astro/observation_models/linkTypeDefs.h"
#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/math/interpolators/lookupScheme.h"


namespace tudat
{

namespace observation_models
{

//! Enum listing types of observation biases that are availabe
enum ObservationBiasTypes
{
    multiple_observation_biases,
    constant_absolute_bias,
    constant_relative_bias,
    arc_wise_constant_absolute_bias,
    arc_wise_constant_relative_bias,
};

//! Base class (non-functional) for describing observation biases
/*!
 * Base class (non-functional) for describing observation biases. In this context, an observation bias denotes any deviation
 * from the ideal observable between two reference points.
 */
template< int ObservationSize = 1 >
class ObservationBias
{
public:

    //! Constructor
    ObservationBias( ){ }

    //! Destructor
    virtual ~ObservationBias( ){ }

    //! Pure virtual function to retrieve the observation bias.
    /*!
     * Pure virtual function to retrieve the observation bias as a function of the observation times and states (which are
     * typically computed by the ObservationModel)
     * \param linkEndTimes List of times at each link end during observation.
     * \param linkEndStates List of states at each link end during observation.
     * \param currentObservableValue  Unbiased value of the observable (default NAN for compatibility purposes with original
     * version of code).
     * \return Observation bias at given times and states.
     */
    virtual Eigen::Matrix< double, ObservationSize, 1 > getObservationBias(
            const std::vector< double >& linkEndTimes,
            const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
            const Eigen::Matrix< double, ObservationSize, 1 >& currentObservableValue =
            Eigen::Matrix< double, ObservationSize, 1 >::Constant( TUDAT_NAN ) ) = 0;

    //! Function to return the size of the associated observation
    /*!
     * Function to return the size of the associated observation
     * \return Size of the associated observation
     */
    int getObservationSize( )
    {
        return ObservationSize;
    }
};

//! Class for a constant absolute observation bias of a given size
/*!
 *  Class for a constant absolute observation bias of a given size. For unbiases observation h and bias A, the biased observation
 *  is computed as h + A
 */
template< int ObservationSize = 1 >
class ConstantObservationBias: public ObservationBias< ObservationSize >
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param observationBias Constant (entry-wise) observation bias.
     */
    ConstantObservationBias( const Eigen::Matrix< double, ObservationSize, 1 > observationBias ):
        observationBias_( observationBias ){ }

    //! Destructor
    ~ConstantObservationBias( ){ }

    //! Function to retrieve the constant observation bias.
    /*!
     * Function to retrieve the constant observation bias.
     * \param linkEndTimes List of times at each link end during observation (unused).
     * \param linkEndStates List of states at each link end during observation (unused).
     * \param currentObservableValue  Unbiased value of the observable (unused and default NAN).
     * \return Constant observation bias.
     */
    Eigen::Matrix< double, ObservationSize, 1 > getObservationBias(
            const std::vector< double >& linkEndTimes,
            const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
            const Eigen::Matrix< double, ObservationSize, 1 >& currentObservableValue =
            ( Eigen::Matrix< double, ObservationSize, 1 >( ) << TUDAT_NAN ).finished( ) )
    {
        return observationBias_;
    }


    //! Function retrieve the constant (entry-wise) absolute observation bias.
    /*!
     * Function retrieve the constant (entry-wise) absolute observation bias.
     * \return The constant (entry-wise) absolute observation bias.
     */
    Eigen::Matrix< double, ObservationSize, 1 > getConstantObservationBias( )
    {
        return observationBias_;
    }

    //! Function to reset the constant (entry-wise) absolute observation bias.
    /*!
     * Function to reset the constant (entry-wise) absolute observation bias.
     * \param observationBias The new constant (entry-wise) absolute observation bias.
     */
    void resetConstantObservationBias( const Eigen::Matrix< double, ObservationSize, 1 >& observationBias )
    {
        observationBias_ = observationBias;
    }

    //! Function retrieve the constant (entry-wise) absolute observation bias as a variable-size vector.
    /*!
     * Function retrieve the constant (entry-wise) absolute observation bias as a variable-size vector
     * \return The constant (entry-wise) absolute observation bias.
     */
    Eigen::VectorXd getTemplateFreeConstantObservationBias( )
    {
        return observationBias_;
    }

    //! Function to reset the constant (entry-wise) absolute observation bias with variable-size input.
    /*!
     * Function to reset the constant (entry-wise) absolute observation bias with variable-size input. Input VectorXd size
     * must match ObservationSize class template parameter.
     * \param observationBias The new constant (entry-wise) absolute observation bias.
     */
    void resetConstantObservationBiasTemplateFree( const Eigen::VectorXd& observationBias )
    {
        if( observationBias.rows( ) == ObservationSize )
        {
            observationBias_ = observationBias;
        }
        else
        {
            throw std::runtime_error( "Error when resetting constant bias, size is inconsistent" );
        }
    }


private:

    //! Constant (entry-wise) observation bias.
    Eigen::Matrix< double, ObservationSize, 1 > observationBias_;

};

//! Class for an arc-wise constant absolute observation bias of a given size
/*!
 *  Class for an  arc-wise constant absolute observation bias of a given size. For unbiases observation h and bias A,
 *  the biased observation is computed as h + A. The bias A is provided per arc, with the arc start times provided to the
 *  class constructor.
 */
template< int ObservationSize = 1 >
class ConstantArcWiseObservationBias: public ObservationBias< ObservationSize >
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param arcStartTimes Start times for arcs in which biases (observationBiases) are used
     * \param observationBiases Absolute biases, constant per arc
     * \param linkEndIndexForTime Link end index from which the 'current time' is determined (e.g. entry from linkEndTimes used
     * in getObservationBias function.
     */
    ConstantArcWiseObservationBias(
            const std::vector< double >& arcStartTimes,
            const std::vector< Eigen::Matrix< double, ObservationSize, 1 > >& observationBiases,
            const int linkEndIndexForTime ):
        arcStartTimes_( arcStartTimes ), observationBiases_( observationBiases ), linkEndIndexForTime_( linkEndIndexForTime )
    {
        resetArcStartTimes( arcStartTimes_, observationBiases_ );
    }

    //! Destructor
    ~ConstantArcWiseObservationBias( ){ }

    //! Function to retrieve the observation bias, determining the current arc from linkEndTimes.
    /*!
     * Function to retrieve the constant observation bias, determining the current arc from linkEndTimes.
     * \param linkEndTimes List of times at each link end during observation
     * \param linkEndStates List of states at each link end during observation (unused).
     * \param currentObservableValue  Unbiased value of the observable (unused and default NAN).
     * \return Constant observation bias.
     */
    Eigen::Matrix< double, ObservationSize, 1 > getObservationBias(
            const std::vector< double >& linkEndTimes,
            const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
            const Eigen::Matrix< double, ObservationSize, 1 >& currentObservableValue =
            ( Eigen::Matrix< double, ObservationSize, 1 >( ) << TUDAT_NAN ).finished( ) )
    {
        return observationBiases_.at(
                    lookupScheme_->findNearestLowerNeighbour( linkEndTimes.at( linkEndIndexForTime_ ) ) );
    }

    //! Function retrieve the list of arc-wise constant absolute observation biases as a variable-size vector.
    /*!
     * Function retrieve the  list of absolute observation biases as a variable-size vector, conatenated first by
     * entry, and then by arc.
     * \return List of absolute observation biases as a variable-size vector, conatenated first by
     * entry, and then by arc.
     */
    std::vector< Eigen::VectorXd > getTemplateFreeConstantObservationBias( )
    {
        std::vector< Eigen::VectorXd > templateFreeObservationBiases;
        for( unsigned int i = 0; i < observationBiases_.size( ); i++ )
        {
            templateFreeObservationBiases.push_back( observationBiases_.at( i ) );
        }
        return templateFreeObservationBiases;
    }

    //! Function reset the list of absolute observation biases
    /*!
     *  Function reset the list of absolute observation biases
     *  \param observationBiases The new list of absolute observation biases as a variable-size vector, with the bias for arc i
     *  in index i of the input vector
     */
    void resetConstantObservationBiasTemplateFree( const std::vector< Eigen::VectorXd >& observationBiases )
    {
        if( observationBiases_.size( ) == observationBiases.size( ) )
        {
            for( unsigned int i = 0; i < observationBiases.size( ); i++ )
            {
                if( ! ( observationBiases.at( i ).rows( ) == ObservationSize ) )
                {
                    throw std::runtime_error( "Error when resetting arc-wise constant bias, single entry size is inconsistent" );
                }
                else
                {
                    observationBiases_[ i ] = observationBiases.at( i );
                }
            }
        }
        else
        {
            throw std::runtime_error( "Error when resetting arc-wise constant bias, size is inconsistent" );
        }
    }

    //! Function to reset the arc start times for the bias (and associated biases)
    /*!
     * Function to reset the arc start times for the bias (and associated biases)
     * \param arcStartTimes Start times for arcs in which biases (observationBiases) are used
     * \param observationBiases Absolute biases, constant per arc
     */
    void resetArcStartTimes( const std::vector< double >& arcStartTimes,
                             const std::vector< Eigen::Matrix< double, ObservationSize, 1 > >& observationBiases )
    {
        arcStartTimes_ = arcStartTimes;
        observationBiases_ = observationBiases;

        if( arcStartTimes_.size( ) != observationBiases_.size( ) )
        {
            throw std::runtime_error( "Error when creating constant arc-wise biases, input is inconsistent" );
        }

        // Create current arc lookup scheme
        std::vector< double > lookupSchemeTimes = arcStartTimes_;
        lookupSchemeTimes.push_back( std::numeric_limits< double >::max( ) );
        lookupScheme_ = std::make_shared< interpolators::HuntingAlgorithmLookupScheme< double > >(
                    lookupSchemeTimes );

    }
    //! Function to retrieve start times for arcs in which biases (observationBiases) are used
    /*!
     * Function to retrieve start times for arcs in which biases (observationBiases) are used
     * \return Start times for arcs in which biases (observationBiases) are used
     */
    std::vector< double > getArcStartTimes( )
    {
        return arcStartTimes_;
    }

    //! Function to retrieve link end index from which the 'current time' is determined
    /*!
     * Function to retrieve link end index from which the 'current time' is determined
     * \return Link end index from which the 'current time' is determined
     */
    int getLinkEndIndexForTime( )
    {
        return linkEndIndexForTime_;
    }

    //! Function to retrieve object used to determine the index from observationBiases_ to be used, based on the current time.
    /*!
     * Function to retrieve object used to determine the index from observationBiases_ to be used, based on the current time.
     * \return Object used to determine the index from observationBiases_ to be used, based on the current time.
     */
    std::shared_ptr< interpolators::LookUpScheme< double > > getLookupScheme( )
    {
        return lookupScheme_;
    }

private:

    //! Start times for arcs in which biases (observationBiases) are used
    std::vector< double > arcStartTimes_;

    //! Absolute biases, constant per arc
    std::vector< Eigen::Matrix< double, ObservationSize, 1 > > observationBiases_;

    //! Link end index from which the 'current time' is determined (e.g. entry from linkEndTimes used in getObservationBias
    //! function.
    int linkEndIndexForTime_;

    //! Object used to determine the index from observationBiases_ to be used, based on the current time.
    std::shared_ptr< interpolators::LookUpScheme< double > > lookupScheme_;
};

//! Class for a constant relative observation bias of a given size
/*!
 *  Class for a constant relative observation bias of a given size. For unbiases observation h and bias A, the biased observation
 *  is computed as h .* A, where .* is the component-wise multiplication.
 */
template< int ObservationSize = 1 >
class ConstantRelativeObservationBias: public ObservationBias< ObservationSize >
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param relativeObservationBias Constant (entry-wise) observation bias.
     */
    ConstantRelativeObservationBias( const Eigen::Matrix< double, ObservationSize, 1 > relativeObservationBias ):
        relativeObservationBias_( relativeObservationBias ){ }

    //! Destructor
    ~ConstantRelativeObservationBias( ){ }

    //! Function to retrieve the constant relative observation bias.
    /*!
     * Function to retrieve the constant relative observation bias.
     * \param linkEndTimes List of times at each link end during observation (unused).
     * \param linkEndStates List of states at each link end during observation (unused).
     * \param currentObservableValue Unbiased value of the observable
     * \return Relative observation bias.
     */
    Eigen::Matrix< double, ObservationSize, 1 > getObservationBias(
            const std::vector< double >& linkEndTimes,
            const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
            const Eigen::Matrix< double, ObservationSize, 1 >& currentObservableValue )
    {
        return relativeObservationBias_.cwiseProduct( currentObservableValue );
    }

    //! Function retrieve the constant (entry-wise) relative observation bias.
    /*!
     * Function retrieve the constant (entry-wise) relative observation bias.
     * \return The constant (entry-wise) relative observation bias.
     */
    Eigen::Matrix< double, ObservationSize, 1 > getConstantObservationBias( )
    {
        return relativeObservationBias_;
    }

    //! Function to reset the constant (entry-wise) relative observation bias.
    /*!
     * Function to reset the constant (entry-wise) relative observation bias.
     * \param relativeObservationBias The new constant (entry-wise) relative observation bias.
     */
    void resetConstantObservationBias( const Eigen::Matrix< double, ObservationSize, 1 >& relativeObservationBias )
    {
        relativeObservationBias_ = relativeObservationBias;
    }

    //! Function retrieve the constant (entry-wise) relative observation bias as a variable-size vector.
    /*!
     * Function retrieve the constant (entry-wise) relative observation bias as a variable-size vector
     * \return The constant (entry-wise) relative observation bias.
     */
    Eigen::VectorXd getTemplateFreeConstantObservationBias( )
    {
        return relativeObservationBias_;
    }

    //! Function to reset the constant (entry-wise) relative observation bias with variable-size input.
    /*!
     * Function to reset the constant (entry-wise) relative observation bias with variable-size input. Input VectorXd size
     * must match ObservationSize class template parameter.
     * \param relativeObservationBias The new constant (entry-wise) relative observation bias.
     */
    void resetConstantObservationBiasTemplateFree( const Eigen::VectorXd& relativeObservationBias )
    {
        if( relativeObservationBias.rows( ) == ObservationSize )
        {
            relativeObservationBias_ = relativeObservationBias;
        }
        else
        {
            throw std::runtime_error( "Error when resetting constant bias, size is inconsistent" );
        }
    }

private:

    //! Constant (entry-wise) relative observation bias.
    Eigen::Matrix< double, ObservationSize, 1 > relativeObservationBias_;

};


//! Class for an arc-wise constant relative observation bias of a given size
/*!
 *  Class for an  arc-wise constant relative observation bias of a given size. For unbiases observation h and bias A,
 *  the biased observation is computed as h .* A, where .* is the component-wise multiplication.The bias A is provided per arc,
 *  with the arc start times provided to the class constructor.
 */
template< int ObservationSize = 1 >
class ConstantRelativeArcWiseObservationBias: public ObservationBias< ObservationSize >
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param arcStartTimes Start times for arcs in which biases (observationBiases) are used
     * \param observationBiases Relative biases, constant per arc
     * \param linkEndIndexForTime Link end index from which the 'current time' is determined (e.g. entry from linkEndTimes used
     * in getObservationBias function.
     */
    ConstantRelativeArcWiseObservationBias(
            const std::vector< double >& arcStartTimes,
            const std::vector< Eigen::Matrix< double, ObservationSize, 1 > >& observationBiases,
            const int linkEndIndexForTime ):
        arcStartTimes_( arcStartTimes ), observationBiases_( observationBiases ), linkEndIndexForTime_( linkEndIndexForTime )
    {
        if( arcStartTimes_.size( ) != observationBiases_.size( ) )
        {
            throw std::runtime_error( "Error when creating constant arc-wise biases, input is inconsistent" );
        }

        // Create current arc lookup scheme
        std::vector< double > lookupSchemeTimes = arcStartTimes_;
        lookupSchemeTimes.push_back( std::numeric_limits< double >::max( ) );
        lookupScheme_ = std::make_shared< interpolators::HuntingAlgorithmLookupScheme< double > >(
                    lookupSchemeTimes );
    }

    //! Destructor
    ~ConstantRelativeArcWiseObservationBias( ){ }

    //! Function to retrieve the constant observation bias, determining the current arc from linkEndTimes.
    /*!
     * Function to retrieve the constant observation bias, determining the current arc from linkEndTimes.
     * \param linkEndTimes List of times at each link end during observation
     * \param linkEndStates List of states at each link end during observation (unused).
     * \param currentObservableValue  Unbiased value of the observable (unused and default NAN).
     * \return Constant observation bias.
     */
    Eigen::Matrix< double, ObservationSize, 1 > getObservationBias(
            const std::vector< double >& linkEndTimes,
            const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
            const Eigen::Matrix< double, ObservationSize, 1 >& currentObservableValue =
            ( Eigen::Matrix< double, ObservationSize, 1 >( ) << TUDAT_NAN ).finished( ) )
    {
        return observationBiases_.at( lookupScheme_->findNearestLowerNeighbour( linkEndTimes.at( linkEndIndexForTime_ ) ) ).
                cwiseProduct( currentObservableValue );
    }

    //! Function retrieve the constant (entry-wise) relative observation bias as a variable-size vector.
    /*!
         * Function retrieve the constant (entry-wise) relative observation bias as a variable-size vector
         * \return The constant (entry-wise) relative observation bias.
         */
    std::vector< Eigen::VectorXd > getTemplateFreeConstantObservationBias( )
    {
        std::vector< Eigen::VectorXd > templateFreeObservationBiases;
        for( unsigned int i = 0; i < observationBiases_.size( ); i++ )
        {
            templateFreeObservationBiases.push_back( observationBiases_.at( i ) );
        }
        return templateFreeObservationBiases;
    }

    //! Function to reset the constant (entry-wise) relative observation bias with variable-size input.
    /*!
     *  Function to reset the constant (entry-wise) relative observation bias with variable-size input. Input VectorXd size
     *  must match ObservationSize class template parameter.
     *  \param observationBiases The new constant arc-wise list of (entry-wise) relative observation bias, with the bias for arc i
     *  in index i of the input vector
     */
    void resetConstantObservationBiasTemplateFree( const std::vector< Eigen::VectorXd >& observationBiases )
    {
        if( observationBiases_.size( ) == observationBiases.size( ) )
        {
            for( unsigned int i = 0; i < observationBiases.size( ); i++ )
            {
                if( ! ( observationBiases.at( i ).rows( ) == ObservationSize ) )
                {
                    throw std::runtime_error( "Error when resetting arc-wise constant bias, single entry size is inconsistent" );
                }
                else
                {
                    observationBiases_[ i ] = observationBiases.at( i );
                }
            }
        }
        else
        {
            throw std::runtime_error( "Error when resetting arc-wise constant bias, size is inconsistent" );
        }
    }

    //! Function to retrieve start times for arcs in which biases (observationBiases) are used
    /*!
     * Function to retrieve start times for arcs in which biases (observationBiases) are used
     * \return Start times for arcs in which biases (observationBiases) are used
     */
    std::vector< double > getArcStartTimes( )
    {
        return arcStartTimes_;
    }

    //! Function to retrieve link end index from which the 'current time' is determined
    /*!
     * Function to retrieve link end index from which the 'current time' is determined
     * \return Link end index from which the 'current time' is determined
     */
    int getLinkEndIndexForTime( )
    {
        return linkEndIndexForTime_;
    }

    //! Function to retrieve object used to determine the index from observationBiases_ to be used, based on the current time.
    /*!
     * Function to retrieve object used to determine the index from observationBiases_ to be used, based on the current time.
     * \return Object used to determine the index from observationBiases_ to be used, based on the current time.
     */
    std::shared_ptr< interpolators::LookUpScheme< double > > getLookupScheme( )
    {
        return lookupScheme_;
    }

private:

    //! Start times for arcs in which biases (observationBiases) are used
    std::vector< double > arcStartTimes_;

    //! Absolute biases, constant per arc
    std::vector< Eigen::Matrix< double, ObservationSize, 1 > > observationBiases_;

    //! Link end index from which the 'current time' is determined (e.g. entry from linkEndTimes used in getObservationBias
    //! function.
    int linkEndIndexForTime_;

    //! Object used to determine the index from observationBiases_ to be used, based on the current time.
    std::shared_ptr< interpolators::LookUpScheme< double > > lookupScheme_;
};

//! Class for combining multiple observation bias models into a single bias value
/*!
 *  Class for combining multiple observation bias models into a single bias value. This class computes a list of biases,
 *  all based on the nominal, unbiased, observation and sums them up to form the total observation bias.
 */
template< int ObservationSize = 1 >
class MultiTypeObservationBias: public ObservationBias< ObservationSize >
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param biasList List of bias objects that are to be combined.
     */
    MultiTypeObservationBias( const std::vector< std::shared_ptr< ObservationBias< ObservationSize > > > biasList ):
        biasList_( biasList ){ }

    //! Destructor
    ~MultiTypeObservationBias( ){ }

    //! Function to retrieve the total observation bias.
    /*!
     * Function to retrieve the total observation bias.
     * \param linkEndTimes List of times at each link end during observation.
     * \param linkEndStates List of states at each link end during observation.
     * \param currentObservableValue  Unbiased value of the observable.
     * \return Total observation bias.
     */
    Eigen::Matrix< double, ObservationSize, 1 > getObservationBias(
            const std::vector< double >& linkEndTimes,
            const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
            const Eigen::Matrix< double, ObservationSize, 1 >& currentObservableValue )
    {
        Eigen::Matrix< double, ObservationSize, 1 > totalBias = Eigen::Matrix< double, ObservationSize, 1 >::Zero( );
        for( unsigned int i = 0; i < biasList_.size( ); i++ )
        {
            totalBias += biasList_.at( i )->getObservationBias( linkEndTimes, linkEndStates, currentObservableValue );
        }
        return totalBias;
    }

    //! Function to retrieve the list of bias objects that are to be combined.
    /*!
     * Function to retrieve the list of bias objects that are to be combined.
     * \return The list of bias objects that are to be combined.
     */
    std::vector< std::shared_ptr< ObservationBias< ObservationSize > > > getBiasList( )
    {
        return biasList_;
    }


private:


    //! List of bias objects that are to be combined.
    std::vector< std::shared_ptr< ObservationBias< ObservationSize > > > biasList_;

};

//! Function to retrieve the type of an observation bias
/*!
 *  Function to retrieve the type of an observation bias.
 *  \param biasObject Bias for which the type is to be determined
 *  \return Bias type for biasObject.
 */
template< int ObservationSize >
ObservationBiasTypes getObservationBiasType(
        const std::shared_ptr< ObservationBias< ObservationSize > > biasObject )
{
    ObservationBiasTypes biasType;

    // Check available bias types
    if( std::dynamic_pointer_cast< ConstantObservationBias< ObservationSize > >( biasObject ) != nullptr )
    {
        biasType = constant_absolute_bias;
    }
    else if( std::dynamic_pointer_cast< ConstantArcWiseObservationBias< ObservationSize > >( biasObject ) != nullptr )
    {
        biasType = arc_wise_constant_absolute_bias;
    }
    else if( std::dynamic_pointer_cast< ConstantRelativeObservationBias< ObservationSize > >( biasObject ) != nullptr )
    {
        biasType = constant_relative_bias;
    }
    else if( std::dynamic_pointer_cast< ConstantRelativeArcWiseObservationBias< ObservationSize > >( biasObject ) != nullptr )
    {
        biasType = arc_wise_constant_relative_bias;
    }
    else if( std::dynamic_pointer_cast< MultiTypeObservationBias< ObservationSize > >( biasObject ) != nullptr )
    {
        biasType = multiple_observation_biases;
    }
    else
    {
        std::string errorMessage = "Error, did not recognize observation bias when retrieveing bias type";
        throw std::runtime_error( errorMessage );
    }
    return biasType;
}

} // namespace observation_models

} // namespace tudat

#endif // TUDAT_OBSERVATIONBIAS_H
