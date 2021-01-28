/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_INVPOLY_SHAPING_OPTIMISATION_SETUP_H
#define TUDAT_INVPOLY_SHAPING_OPTIMISATION_SETUP_H

#include <vector>
#include <utility>
#include <limits>

#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include <Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h>
#include <Tudat/Basics/testMacros.h>

#include <Tudat/Mathematics/BasicMathematics/mathematicalConstants.h>

#include "Tudat/Mathematics/NumericalIntegrators/createNumericalIntegrator.h"
#include "Tudat/Mathematics/NumericalQuadrature/createNumericalQuadrature.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include <Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <Tudat/InputOutput/basicInputOutput.h>

#include "Tudat/Astrodynamics/ShapeBasedMethods/invPolyShaping.h"
#include "Tudat/Astrodynamics/ShapeBasedMethods/exposinsShaping.h"
#include "Tudat/Astrodynamics/ShapeBasedMethods/sphericalShaping.h"

#include <random>

#include "pagmo/island.hpp"
#include "pagmo/io.hpp"
#include "pagmo/problem.hpp"
#include <pagmo/rng.hpp>
#include <Eigen/Core>

typedef Eigen::Matrix< double, 6, 1 > StateType;

using namespace tudat::ephemerides;
using namespace tudat::basic_astrodynamics;
using namespace tudat::orbital_element_conversions;
using namespace tudat::basic_mathematics;
using namespace tudat::input_output;
using namespace pagmo;

namespace tudat
{

namespace shape_based_methods
{

//! Test function for a new low-thrust trajectory class in Tudat
struct InvPolyShapingOptimisationProblem
{

    InvPolyShapingOptimisationProblem( ){ }

    InvPolyShapingOptimisationProblem(
            const std::function< Eigen::Vector6d( const double ) >& initialStateFunction,
            const std::function< Eigen::Vector6d( const double, const double ) >& finalStateFunction,
//            double timeOfFlight,
            int numberOfRevolutions,
            simulation_setup::NamedBodyMap bodyMap,
            const std::string bodyToPropagate,
            const std::string centralBody,
            double initialValueFreeCoefficient,
            std::shared_ptr< root_finders::RootFinderSettings >& rootFinderSettings,
            const std::vector< std::vector< double > >& freeCoefficientsBounds,
            const double lowerBoundFreeCoefficient = TUDAT_NAN,
            const double upperBoundFreeCoefficient = TUDAT_NAN
//            std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings = std::shared_ptr< numerical_integrators::IntegratorSettings< > >( )
            ):
        initialStateFunction_( initialStateFunction ),
        finalStateFunction_( finalStateFunction ),
//        timeOfFlight_( timeOfFlight ),
        numberOfRevolutions_( numberOfRevolutions ),
        bodyMap_( bodyMap ), bodyToPropagate_( bodyToPropagate ), centralBody_( centralBody ),
        initialValueFreeCoefficient_( initialValueFreeCoefficient ),
        rootFinderSettings_( rootFinderSettings ),
        problemBounds_( freeCoefficientsBounds ),
        lowerBoundFreeCoefficient_( lowerBoundFreeCoefficient ),
        upperBoundFreeCoefficient_( upperBoundFreeCoefficient )
//        integratorSettings_( integratorSettings )
    {  }

    // Calculates the fitness
    std::vector< double > fitness( const std::vector< double > &x ) const;

    std::pair< std::vector< double >, std::vector< double > > get_bounds() const
    {
        return { problemBounds_[ 0 ], problemBounds_[ 1 ] };
    }

    template <typename Archive>
    void serialize(Archive &ar)
    {
        ar(problemBounds_);
    }

    vector_double::size_type get_nobj() const
    {
        return 1u;
    }

protected:

private:

    const std::function< Eigen::Vector6d( const double ) > initialStateFunction_;

    const std::function< Eigen::Vector6d( const double, const double ) > finalStateFunction_;

//    double timeOfFlight_;

    int numberOfRevolutions_;

    //! Body map object.
    mutable simulation_setup::NamedBodyMap bodyMap_;

    //! Name of the body to be propagated.
    std::string bodyToPropagate_;

    //! Name of the central body.
    std::string centralBody_;

    //! Initial guess for the free coefficient (i.e. coefficient of the second order component of the radial inverse polynomial).
    double initialValueFreeCoefficient_;

    //! Root finder settings, to be used to find the free coefficient value that ensures the time of flight is correct.
    mutable std::shared_ptr< root_finders::RootFinderSettings > rootFinderSettings_;

    const std::vector< std::vector< double > > problemBounds_;

    //! Lower bound for the free coefficient, to be used when trying to match the required time of flight.
    double lowerBoundFreeCoefficient_;

    //! Upper bound for the free coefficient, to be used when trying to match the required time of flight.
    double upperBoundFreeCoefficient_;

//    //! Integrator settings.
//    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings_;



}; // Full optimisation

struct ExposinsShapingOptimisationProblem
{

    ExposinsShapingOptimisationProblem( ){ }

    ExposinsShapingOptimisationProblem(
            const std::function< Eigen::Vector6d( const double ) >& initialStateFunction,
            const std::function< Eigen::Vector6d( const double, const double ) >& finalStateFunction,
//            double timeOfFlight,
            int numberOfRevolutions,
            simulation_setup::NamedBodyMap bodyMap,
            const std::string bodyToPropagate,
            const std::string centralBody,
            double initialValueFreeCoefficient,
            std::shared_ptr< root_finders::RootFinderSettings >& rootFinderSettings,
            const std::vector< std::vector< double > >& freeCoefficientsBounds,
            const double lowerBoundFreeCoefficient = TUDAT_NAN,
            const double upperBoundFreeCoefficient = TUDAT_NAN
//            std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings = std::shared_ptr< numerical_integrators::IntegratorSettings< > >( )
            ):
        initialStateFunction_( initialStateFunction ),
        finalStateFunction_( finalStateFunction ),
//        timeOfFlight_( timeOfFlight ),
        numberOfRevolutions_( numberOfRevolutions ),
        bodyMap_( bodyMap ), bodyToPropagate_( bodyToPropagate ), centralBody_( centralBody ),
        initialValueFreeCoefficient_( initialValueFreeCoefficient ),
        rootFinderSettings_( rootFinderSettings ),
        problemBounds_( freeCoefficientsBounds ),
        lowerBoundFreeCoefficient_( lowerBoundFreeCoefficient ),
        upperBoundFreeCoefficient_( upperBoundFreeCoefficient )
//        integratorSettings_( integratorSettings )
    {  }

    // Calculates the fitness
    std::vector< double > fitness( const std::vector< double > &x ) const;

    std::pair< std::vector< double >, std::vector< double > > get_bounds() const
    {
        return { problemBounds_[ 0 ], problemBounds_[ 1 ] };
    }

    template <typename Archive>
    void serialize(Archive &ar)
    {
        ar(problemBounds_);
    }

    vector_double::size_type get_nobj() const
    {
        return 1u;
    }

protected:

private:

    const std::function< Eigen::Vector6d( const double ) > initialStateFunction_;

    const std::function< Eigen::Vector6d( const double, const double ) > finalStateFunction_;

//    double timeOfFlight_;

    int numberOfRevolutions_;

    //! Body map object.
    mutable simulation_setup::NamedBodyMap bodyMap_;

    //! Name of the body to be propagated.
    std::string bodyToPropagate_;

    //! Name of the central body.
    std::string centralBody_;

    //! Initial guess for the free coefficient (i.e. coefficient of the second order component of the radial inverse polynomial).
    double initialValueFreeCoefficient_;

    //! Root finder settings, to be used to find the free coefficient value that ensures the time of flight is correct.
    mutable std::shared_ptr< root_finders::RootFinderSettings > rootFinderSettings_;

    const std::vector< std::vector< double > > problemBounds_;

    //! Lower bound for the free coefficient, to be used when trying to match the required time of flight.
    double lowerBoundFreeCoefficient_;

    //! Upper bound for the free coefficient, to be used when trying to match the required time of flight.
    double upperBoundFreeCoefficient_;

//    //! Integrator settings.
//    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings_;



}; // exposins optimisation

struct AnvPolyShapingOptimisationProblem
{

    AnvPolyShapingOptimisationProblem( ){ }

    AnvPolyShapingOptimisationProblem(
            const std::function< Eigen::Vector6d( const double ) >& initialStateFunction,
            const std::function< Eigen::Vector6d( const double, const double ) >& finalStateFunction,
//            double timeOfFlight,
            int numberOfRevolutions,
            simulation_setup::NamedBodyMap bodyMap,
            const std::string bodyToPropagate,
            const std::string centralBody,
            double initialValueFreeCoefficient,
            std::shared_ptr< root_finders::RootFinderSettings >& rootFinderSettings,
            const std::vector< std::vector< double > >& freeCoefficientsBounds,
            const double lowerBoundFreeCoefficient = TUDAT_NAN,
            const double upperBoundFreeCoefficient = TUDAT_NAN
//            std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings = std::shared_ptr< numerical_integrators::IntegratorSettings< > >( )
            ):
        initialStateFunction_( initialStateFunction ),
        finalStateFunction_( finalStateFunction ),
//        timeOfFlight_( timeOfFlight ),
        numberOfRevolutions_( numberOfRevolutions ),
        bodyMap_( bodyMap ), bodyToPropagate_( bodyToPropagate ), centralBody_( centralBody ),
        initialValueFreeCoefficient_( initialValueFreeCoefficient ),
        rootFinderSettings_( rootFinderSettings ),
        problemBounds_( freeCoefficientsBounds ),
        lowerBoundFreeCoefficient_( lowerBoundFreeCoefficient ),
        upperBoundFreeCoefficient_( upperBoundFreeCoefficient )
//        integratorSettings_( integratorSettings )
    {  }

    // Calculates the fitness
    std::vector< double > fitness( const std::vector< double > &x ) const;

    std::pair< std::vector< double >, std::vector< double > > get_bounds() const
    {
        return { problemBounds_[ 0 ], problemBounds_[ 1 ] };
    }

    template <typename Archive>
    void serialize(Archive &ar)
    {
        ar(problemBounds_);
    }

    vector_double::size_type get_nobj() const
    {
        return 1u;
    }

protected:

private:

    const std::function< Eigen::Vector6d( const double ) > initialStateFunction_;

    const std::function< Eigen::Vector6d( const double, const double ) > finalStateFunction_;

//    double timeOfFlight_;

    int numberOfRevolutions_;

    //! Body map object.
    mutable simulation_setup::NamedBodyMap bodyMap_;

    //! Name of the body to be propagated.
    std::string bodyToPropagate_;

    //! Name of the central body.
    std::string centralBody_;

    //! Initial guess for the free coefficient (i.e. coefficient of the second order component of the radial inverse polynomial).
    double initialValueFreeCoefficient_;

    //! Root finder settings, to be used to find the free coefficient value that ensures the time of flight is correct.
    mutable std::shared_ptr< root_finders::RootFinderSettings > rootFinderSettings_;

    const std::vector< std::vector< double > > problemBounds_;

    //! Lower bound for the free coefficient, to be used when trying to match the required time of flight.
    double lowerBoundFreeCoefficient_;

    //! Upper bound for the free coefficient, to be used when trying to match the required time of flight.
    double upperBoundFreeCoefficient_;

//    //! Integrator settings.
//    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings_;



}; // invpoly optimisation

struct SphericalShapingOptimisationProblem
{

    SphericalShapingOptimisationProblem( ){ }

    SphericalShapingOptimisationProblem(
            const std::function< Eigen::Vector6d( const double ) >& initialStateFunction,
            const std::function< Eigen::Vector6d( const double, const double ) >& finalStateFunction,
//            double timeOfFlight,
            int numberOfRevolutions,
            simulation_setup::NamedBodyMap bodyMap,
            const std::string bodyToPropagate,
            const std::string centralBody,
            double initialValueFreeCoefficient,
            std::shared_ptr< root_finders::RootFinderSettings >& rootFinderSettings,
            const std::vector< std::vector< double > >& freeCoefficientsBounds,
            const double lowerBoundFreeCoefficient = TUDAT_NAN,
            const double upperBoundFreeCoefficient = TUDAT_NAN
//            std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings = std::shared_ptr< numerical_integrators::IntegratorSettings< > >( )
            ):
        initialStateFunction_( initialStateFunction ),
        finalStateFunction_( finalStateFunction ),
//        timeOfFlight_( timeOfFlight ),
        numberOfRevolutions_( numberOfRevolutions ),
        bodyMap_( bodyMap ), bodyToPropagate_( bodyToPropagate ), centralBody_( centralBody ),
        initialValueFreeCoefficient_( initialValueFreeCoefficient ),
        rootFinderSettings_( rootFinderSettings ),
        problemBounds_( freeCoefficientsBounds ),
        lowerBoundFreeCoefficient_( lowerBoundFreeCoefficient ),
        upperBoundFreeCoefficient_( upperBoundFreeCoefficient )
//        integratorSettings_( integratorSettings )
    {  }

    // Calculates the fitness
    std::vector< double > fitness( const std::vector< double > &x ) const;

    std::pair< std::vector< double >, std::vector< double > > get_bounds() const
    {
        return { problemBounds_[ 0 ], problemBounds_[ 1 ] };
    }

    template <typename Archive>
    void serialize(Archive &ar)
    {
        ar(problemBounds_);
    }

    vector_double::size_type get_nobj() const
    {
        return 1u;
    }

protected:

private:

    const std::function< Eigen::Vector6d( const double ) > initialStateFunction_;

    const std::function< Eigen::Vector6d( const double, const double ) > finalStateFunction_;

//    double timeOfFlight_;

    int numberOfRevolutions_;

    //! Body map object.
    mutable simulation_setup::NamedBodyMap bodyMap_;

    //! Name of the body to be propagated.
    std::string bodyToPropagate_;

    //! Name of the central body.
    std::string centralBody_;

    //! Initial guess for the free coefficient (i.e. coefficient of the second order component of the radial inverse polynomial).
    double initialValueFreeCoefficient_;

    //! Root finder settings, to be used to find the free coefficient value that ensures the time of flight is correct.
    mutable std::shared_ptr< root_finders::RootFinderSettings > rootFinderSettings_;

    const std::vector< std::vector< double > > problemBounds_;

    //! Lower bound for the free coefficient, to be used when trying to match the required time of flight.
    double lowerBoundFreeCoefficient_;

    //! Upper bound for the free coefficient, to be used when trying to match the required time of flight.
    double upperBoundFreeCoefficient_;

//    //! Integrator settings.
//    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings_;



}; // spherical optimisation

struct HodographicShapingOptimisationProblem
{

    HodographicShapingOptimisationProblem( ){ }

    HodographicShapingOptimisationProblem(
            const std::function< Eigen::Vector6d( const double ) >& initialStateFunction,
            const std::function< Eigen::Vector6d( const double, const double ) >& finalStateFunction,
//            double timeOfFlight,
            int numberOfRevolutions,
            simulation_setup::NamedBodyMap bodyMap,
            const std::string bodyToPropagate,
            const std::string centralBody,
            double initialValueFreeCoefficient,
            std::shared_ptr< root_finders::RootFinderSettings >& rootFinderSettings,
            const std::vector< std::vector< double > >& freeCoefficientsBounds,
            const double lowerBoundFreeCoefficient = TUDAT_NAN,
            const double upperBoundFreeCoefficient = TUDAT_NAN
//            std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings = std::shared_ptr< numerical_integrators::IntegratorSettings< > >( )
            ):
        initialStateFunction_( initialStateFunction ),
        finalStateFunction_( finalStateFunction ),
//        timeOfFlight_( timeOfFlight ),
        numberOfRevolutions_( numberOfRevolutions ),
        bodyMap_( bodyMap ), bodyToPropagate_( bodyToPropagate ), centralBody_( centralBody ),
        initialValueFreeCoefficient_( initialValueFreeCoefficient ),
        rootFinderSettings_( rootFinderSettings ),
        problemBounds_( freeCoefficientsBounds ),
        lowerBoundFreeCoefficient_( lowerBoundFreeCoefficient ),
        upperBoundFreeCoefficient_( upperBoundFreeCoefficient )
//        integratorSettings_( integratorSettings )
    {  }

    // Calculates the fitness
    std::vector< double > fitness( const std::vector< double > &x ) const;

    std::pair< std::vector< double >, std::vector< double > > get_bounds() const
    {
        return { problemBounds_[ 0 ], problemBounds_[ 1 ] };
    }

    template <typename Archive>
    void serialize(Archive &ar)
    {
        ar(problemBounds_);
    }

    vector_double::size_type get_nobj() const
    {
        return 1u;
    }

protected:

private:

    const std::function< Eigen::Vector6d( const double ) > initialStateFunction_;

    const std::function< Eigen::Vector6d( const double, const double ) > finalStateFunction_;

//    double timeOfFlight_;

    int numberOfRevolutions_;

    //! Body map object.
    mutable simulation_setup::NamedBodyMap bodyMap_;

    //! Name of the body to be propagated.
    std::string bodyToPropagate_;

    //! Name of the central body.
    std::string centralBody_;

    //! Initial guess for the free coefficient (i.e. coefficient of the second order component of the radial inverse polynomial).
    double initialValueFreeCoefficient_;

    //! Root finder settings, to be used to find the free coefficient value that ensures the time of flight is correct.
    mutable std::shared_ptr< root_finders::RootFinderSettings > rootFinderSettings_;

    const std::vector< std::vector< double > > problemBounds_;

    //! Lower bound for the free coefficient, to be used when trying to match the required time of flight.
    double lowerBoundFreeCoefficient_;

    //! Upper bound for the free coefficient, to be used when trying to match the required time of flight.
    double upperBoundFreeCoefficient_;

//    //! Integrator settings.
//    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings_;



}; // spherical optimisation




} // namespace shape_based_methods
} // namespace tudat

#endif // TUDAT_INVPOLY_SHAPING_OPTIMISATION_SETUP_H
