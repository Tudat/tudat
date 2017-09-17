/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include <sstream>
#include <stdexcept>

#include <boost/exception/all.hpp>
#include <boost/math/special_functions/factorials.hpp>

#include "Tudat/Mathematics/BasicMathematics/legendrePolynomials.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

namespace tudat
{
namespace basic_mathematics
{



//! Default constructor, initializes cache object with 0 maximum degree and order.
LegendreCache::LegendreCache( const bool useGeodesyNormalization )
{
    useGeodesyNormalization_  = useGeodesyNormalization;

    if( useGeodesyNormalization_ )
    {
        legendrePolynomialFunction_ = geodesyNormalizedLegendrePolynomialFunction;
    }
    else
    {
        legendrePolynomialFunction_ = regularLegendrePolynomialFunction;
    }

    resetMaximumDegreeAndOrder( 0, 0 );
    computeSecondDerivatives_ = 0;

}

//! Constructor
LegendreCache::LegendreCache( const int maximumDegree, const int maximumOrder, const bool useGeodesyNormalization  )
{
    useGeodesyNormalization_  = useGeodesyNormalization;

    if( useGeodesyNormalization_ )
    {
        legendrePolynomialFunction_ = geodesyNormalizedLegendrePolynomialFunction;
    }
    else
    {
        legendrePolynomialFunction_ = regularLegendrePolynomialFunction;
    }

    resetMaximumDegreeAndOrder( maximumDegree, maximumOrder );
    computeSecondDerivatives_ = 0;
}

//! Get Legendre polynomial from cache when possible, and from direct computation otherwise.
void LegendreCache::update( const double polynomialParameter  )
{
    // Check if cache needs update
    if( !( polynomialParameter == currentPolynomialParameter_ ) )
    {
        currentPolynomialParameter_ = polynomialParameter;

        // Set complement of argument (assuming it to be sine of latitude) cosine of latitude is always positive.
        currentPolynomialParameterComplement_ = std::sqrt( 1.0 - polynomialParameter * polynomialParameter );

        LegendreCache& thisReference = *this;

        int jMax = -1;
        for( int i = 0; i <= maximumDegree_; i++ )
        {
            jMax = std::min( i, maximumOrder_ );
            for( int j = 0; j <= jMax ; j++ )
            {
                // Compute legendre polynomial
                legendreValues_[ i * ( maximumOrder_ + 1 ) + j ] = legendrePolynomialFunction_( i, j, thisReference );

                if( j != 0 )
                {
                    // Compute legendre polynomial derivative
                    if( useGeodesyNormalization_ )
                    {
                        legendreDerivatives_[ i * ( maximumOrder_ + 1 ) + ( j - 1 ) ] =
                                computeGeodesyLegendrePolynomialDerivative(
                                    i, j - 1, currentPolynomialParameter_,
                                    legendreValues_[ i * ( maximumOrder_ + 1 ) + ( j - 1 ) ],
                                legendreValues_[ i * ( maximumOrder_ + 1 ) + j ],
                                derivativeNormalizations_[ i * ( maximumOrder_ + 1 ) + ( j - 1 ) ] );
                    }
                    else
                    {
                        legendreDerivatives_[ i * ( maximumOrder_ + 1 ) + ( j - 1 ) ] =
                                computeLegendrePolynomialDerivative(
                                    j - 1, currentPolynomialParameter_,
                                    legendreValues_[ i * ( maximumOrder_ + 1 ) + ( j - 1 ) ],
                                legendreValues_[ i * ( maximumOrder_ + 1 ) + j ] );
                    }

                }
            }

            // Compute legendre polynomial derivative for i = j  (if needed)
            if( jMax == i )
            {
                if( useGeodesyNormalization_ )
                {
                    legendreDerivatives_[ i * ( maximumOrder_ + 1 ) +  jMax ] =
                            computeGeodesyLegendrePolynomialDerivative(
                                i, jMax, currentPolynomialParameter_,
                                legendreValues_[ i * ( maximumOrder_ + 1 ) + jMax ], 0.0,
                            derivativeNormalizations_[ i * ( maximumOrder_ + 1 ) + jMax ] );
                }
                else
                {
                    legendreDerivatives_[ i * ( maximumOrder_ + 1 ) + jMax ] =
                            computeLegendrePolynomialDerivative(
                                jMax, currentPolynomialParameter_,
                                legendreValues_[ i * ( maximumOrder_ + 1 ) + jMax ], 0.0 );
                }
            }
        }

        // Compute second derivatives of Legendre polynomials if needed
        if( computeSecondDerivatives_ )
        {
            for( int i = 0; i <= maximumDegree_; i++ )
            {
                jMax = std::min( i, maximumOrder_ );
                for( int j = 0; j <= jMax ; j++ )
                {
                    if( j != 0 )
                    {
                        // Compute legendre polynomial second derivatives
                        if( useGeodesyNormalization_ )
                        {
                            legendreSecondDerivatives_[ i * ( maximumOrder_ + 1 ) + ( j - 1 ) ] =
                                    computeGeodesyLegendrePolynomialSecondDerivative(
                                        i, j - 1, currentPolynomialParameter_,
                                        legendreValues_[ i * ( maximumOrder_ + 1 ) + ( j - 1 ) ],
                                    legendreValues_[ i * ( maximumOrder_ + 1 ) + j ],
                                    legendreDerivatives_[ i * ( maximumOrder_ + 1 ) + ( j - 1 ) ],
                                    legendreDerivatives_[ i * ( maximumOrder_ + 1 ) + j ],
                                    derivativeNormalizations_[ i * ( maximumOrder_ + 1 ) + ( j - 1 ) ] );
                        }
                        else
                        {
                            legendreSecondDerivatives_[ i * ( maximumOrder_ + 1 ) + ( j - 1 ) ] =
                                    computeGeodesyLegendrePolynomialSecondDerivative(
                                        i, j - 1, currentPolynomialParameter_,
                                        legendreValues_[ i * ( maximumOrder_ + 1 ) + ( j - 1 ) ],
                                    legendreValues_[ i * ( maximumOrder_ + 1 ) + j ],
                                    legendreDerivatives_[ i * ( maximumOrder_ + 1 ) + ( j - 1 ) ],
                                    legendreDerivatives_[ i * ( maximumOrder_ + 1 ) + j ], 1.0 );
                        }

                    }
                }
                // Compute legendre polynomial second derivative for i = j  (if needed)
                if( jMax == i )
                {
                    if( useGeodesyNormalization_ )
                    {
                        legendreSecondDerivatives_[ i * ( maximumOrder_ + 1 ) +  jMax ] =
                                computeGeodesyLegendrePolynomialSecondDerivative(
                                    i, jMax, currentPolynomialParameter_,
                                    legendreValues_[ i * ( maximumOrder_ + 1 ) + jMax ], 0.0,
                                legendreDerivatives_[ i * ( maximumOrder_ + 1 ) + jMax ], 0.0,
                                derivativeNormalizations_[ i * ( maximumOrder_ + 1 ) + jMax ] );
                    }
                    else
                    {
                        legendreSecondDerivatives_[ i * ( maximumOrder_ + 1 ) +  jMax ] =
                                computeGeodesyLegendrePolynomialSecondDerivative(
                                    i, jMax, currentPolynomialParameter_,
                                    legendreValues_[ i * ( maximumOrder_ + 1 ) + jMax ], 0.0,
                                legendreDerivatives_[ i * ( maximumOrder_ + 1 ) + jMax ], 0.0,
                                1.0 );
                    }
                }

            }
        }
    }
}

//! Update maximum degree and order of cache
void LegendreCache::resetMaximumDegreeAndOrder( const int maximumDegree, const int maximumOrder )
{
    maximumDegree_ = maximumDegree;
    maximumOrder_ = maximumOrder;

    if( maximumOrder_ > maximumDegree_ )
    {
        maximumOrder_ = maximumDegree_;
    }

    legendreValues_.resize( ( maximumDegree_ + 1 ) * ( maximumOrder_ + 1 ) );
    legendreDerivatives_.resize( ( maximumDegree_ + 1 ) * ( maximumOrder_ + 1 ) );
    legendreSecondDerivatives_.resize( ( maximumDegree_ + 1 ) * ( maximumOrder_ + 1 ) );

    derivativeNormalizations_.resize( ( maximumDegree_ + 1 ) * ( maximumOrder_ + 1 ) );

    for( int i = 0; i <= maximumDegree_; i++ )
    {
        for( int j = 0; ( ( j <= i ) && ( j <= maximumOrder_ ) ) ; j++ )
        {
            // Compute normalization correction factor.
            derivativeNormalizations_[ i * ( maximumOrder_ + 1 ) + j ] = std::sqrt(
                        ( static_cast< double >( i + j + 1 ) )
                        * ( static_cast< double >( i - j ) ) );

            // If order is zero apply multiplication factor.
            if ( j == 0 )
            {
                derivativeNormalizations_[ i * ( maximumOrder_ + 1 ) + j ] *= std::sqrt( 0.5 );
            }
        }
    }

    currentPolynomialParameter_ = TUDAT_NAN;
    currentPolynomialParameterComplement_ = TUDAT_NAN;
}


//! Get Legendre polynomial value from the cache.
double LegendreCache::getLegendrePolynomial(
        const int degree, const int order )
{
    if( degree > maximumDegree_ || order > maximumOrder_ )
    {
        std::string errorMessage = "Error when requesting legendre cache, maximum degree or order exceeded " +
                boost::lexical_cast< std::string >( degree ) + " " +
                boost::lexical_cast< std::string >( maximumDegree_ ) + " " +
                boost::lexical_cast< std::string >( order ) + " " +
                boost::lexical_cast< std::string >( maximumOrder_ );
        throw std::runtime_error( errorMessage );
        return TUDAT_NAN;
    }
    else if( order > degree )
    {
        return 0.0;
    }
    else
    {
        return legendreValues_[ degree * ( maximumOrder_ + 1  ) + order ];
    };
}

//! Get first derivative of Legendre polynomial value from the cache.
double LegendreCache::getLegendrePolynomialDerivative(
        const int degree, const int order )
{
    if( degree > ( maximumDegree_ ) || order > maximumOrder_ )
    {
        std::string errorMessage = "Error when requesting legendre cache first derivatives, maximum degree or order exceeded " +
                boost::lexical_cast< std::string >( degree ) + " " +
                boost::lexical_cast< std::string >( maximumDegree_ ) + " " +
                boost::lexical_cast< std::string >( order ) + " " +
                boost::lexical_cast< std::string >( maximumOrder_ );
        throw std::runtime_error( errorMessage );
        return TUDAT_NAN;
    }
    else if( order > degree )
    {
        return 0.0;
    }
    else
    {
        return legendreDerivatives_[ degree * ( maximumOrder_ + 1  ) + order ];
    };
}

//! Get second derivative of Legendre polynomial value from the cache.
double LegendreCache::getLegendrePolynomialSecondDerivative(
        const int degree, const int order )
{
    if( degree > ( maximumDegree_  ) || order > maximumOrder_ )
    {
        std::string errorMessage = "Error when requesting legendre cache second derivatives, maximum degree or order exceeded " +
                boost::lexical_cast< std::string >( degree ) + " " +
                boost::lexical_cast< std::string >( maximumDegree_ ) + " " +
                boost::lexical_cast< std::string >( order ) + " " +
                boost::lexical_cast< std::string >( maximumOrder_ );
        throw std::runtime_error( errorMessage );
        return TUDAT_NAN;
    }
    else if( computeSecondDerivatives_ == 0 )
    {
        throw std::runtime_error( "Error when requesting legendre cache second derivatives, no computations performed" );
    }
    else if( order > degree )
    {
        return 0.0;
    }
    else
    {
        return legendreSecondDerivatives_[ degree * ( maximumOrder_ + 1  ) + order ];
    };
}

//! Compute unnormalized associated Legendre polynomial.
double computeLegendrePolynomialFromCache( const int degree,
                                           const int order,
                                           LegendreCache& legendreCache )
{
    if( legendreCache.getUseGeodesyNormalization( ) )
    {
        throw std::runtime_error( "Error when computing Legendre polynomial, input uses normalization" );
    }

    // If degree or order is negative...
    if ( degree < 0 || order < 0 )
    {
        // Set error message.
        std::stringstream errorMessage;
        errorMessage << "Error: the Legendre polynomial of = " << degree << " and order = "
                     << order << " is undefined." << std::endl;

        // Throw a run-time error.
        throw std::runtime_error( errorMessage.str( ) );
    }

    // Else if order is greater than degree...
    else if ( order > degree && degree >= 0 )
    {
        // Return zero.
        return 0.0;
    }

    // Else if order and degree are lower than 2...
    else if ( degree <= 1 && order <= 1 )
    {
        // Compute polynomial explicitly.
        return computeLegendrePolynomialExplicit( degree, order, legendreCache.getCurrentPolynomialParameter( ) );
    }

    // Else if degree and order are sectoral...
    else if ( degree == order )
    {
        // Obtain polynomial of degree one and order one.
        const double degreeOneOrderOnePolynomial = legendreCache.getLegendrePolynomial(
                    1, 1 );

        // Obtain prior sectoral polynomial.
        const double priorSectoralPolynomial = legendreCache.getLegendrePolynomial(
                    degree - 1, order - 1 );

        // Compute polynomial.
        return computeLegendrePolynomialDiagonal(
                    degree, degreeOneOrderOnePolynomial, priorSectoralPolynomial );
    }

    // Else degree and order are zonal/tessoral...
    else
    {
        // Obtain prior degree polynomial.
        const double oneDegreePriorPolynomial = legendreCache.getLegendrePolynomial(
                    degree - 1, order );

        // Obtain two degrees prior polynomial.
        const double twoDegreesPriorPolynomial = legendreCache.getLegendrePolynomial(
                    degree - 2, order );

        // Compute polynomial.
        return computeLegendrePolynomialVertical( degree,
                                                  order,
                                                  legendreCache.getCurrentPolynomialParameter( ),
                                                  oneDegreePriorPolynomial,
                                                  twoDegreesPriorPolynomial );
    }
}


double computeLegendrePolynomial( const int degree,
                                  const int order,
                                  const double legendreParameter )
{
    LegendreCache legendreCache( degree, order, 0 );
    legendreCache.update( legendreParameter );
    return computeLegendrePolynomialFromCache( degree, order, legendreCache );
}


//! Compute geodesy-normalized associated Legendre polynomial.
double computeGeodesyLegendrePolynomialFromCache( const int degree,
                                                  const int order,
                                                  LegendreCache& geodesyLegendreCache )
{

    if( !geodesyLegendreCache.getUseGeodesyNormalization( ) )
    {
        throw std::runtime_error( "Error when computing Legendre polynomial, input uses no normalization" );
    }

    // If degree or order is negative...
    if ( degree < 0 || order < 0 )
    {
        // Set error message.
        std::stringstream errorMessage;
        errorMessage << "Error: the Legendre polynomial of = " << degree << " and order = "
                     << order << " is undefined." << std::endl;

        // Throw a run-time error.
        throw std::runtime_error( errorMessage.str( ) );
    }

    // Else if order is greater than degree...
    else if ( order > degree && degree >= 0 )
    {
        // Return zero.
        return 0.0;
    }

    // Else if order and degree are lower than 2...
    else if ( degree <= 1 && order <= 1 )
    {
        // Compute polynomial explicitly.
        return computeGeodesyLegendrePolynomialExplicit( degree, order, geodesyLegendreCache.getCurrentPolynomialParameter( ) );
    }

    // Else if degree and order are sectoral...
    else if ( degree == order )
    {
        // Obtain polynomial of degree one and order one.
        double degreeOneOrderOnePolynomial = geodesyLegendreCache.getLegendrePolynomial(
                    1, 1 );

        // Obtain prior sectoral polynomial.
        double priorSectoralPolynomial = geodesyLegendreCache.getLegendrePolynomial(
                    degree - 1, order - 1 );

        // Compute polynomial.
        return computeGeodesyLegendrePolynomialDiagonal(
                    degree, degreeOneOrderOnePolynomial, priorSectoralPolynomial );
    }

    // Else degree and order are zonal/tessoral...
    else
    {
        // Obtain prior degree polynomial.
        double oneDegreePriorPolynomial = geodesyLegendreCache.getLegendrePolynomial(
                    degree - 1, order );

        // Obtain two degrees prior polynomial.
        double twoDegreesPriorPolynomial = geodesyLegendreCache.getLegendrePolynomial(
                    degree - 2, order );

        // Compute polynomial.
        return computeGeodesyLegendrePolynomialVertical( degree,
                                                         order,
                                                         geodesyLegendreCache.getCurrentPolynomialParameter( ),
                                                         oneDegreePriorPolynomial,
                                                         twoDegreesPriorPolynomial );
    }
}

//! Compute geodesy-normalized associated Legendre polynomial.
double computeGeodesyLegendrePolynomial( const int degree,
                                         const int order,
                                         const double legendreParameter )
{
    LegendreCache legendreCache( degree, order, 1 );
    legendreCache.update( legendreParameter );
    return computeGeodesyLegendrePolynomialFromCache( degree, order, legendreCache );
}

//! Compute derivative of unnormalized Legendre polynomial.
double computeLegendrePolynomialDerivative( const int order,
                                            const double polynomialParameter,
                                            const double currentLegendrePolynomial,
                                            const double incrementedLegendrePolynomial )
{
    // Return polynomial derivative.
    return incrementedLegendrePolynomial
            / std::sqrt( 1.0 - polynomialParameter * polynomialParameter )
            - static_cast< double >( order ) * polynomialParameter
            / ( 1.0 - polynomialParameter * polynomialParameter )
            * currentLegendrePolynomial;
}

//! Compute derivative of geodesy-normalized Legendre polynomial.
double computeGeodesyLegendrePolynomialDerivative( const int degree,
                                                   const int order,
                                                   const double polynomialParameter,
                                                   const double currentLegendrePolynomial,
                                                   const double incrementedLegendrePolynomial,
                                                   const double normalizationCorrection )
{
    // Return polynomial derivative.
    return normalizationCorrection * incrementedLegendrePolynomial
            / std::sqrt( 1.0 - polynomialParameter * polynomialParameter )
            - static_cast< double >( order ) * polynomialParameter
            / ( 1.0 - polynomialParameter * polynomialParameter )
            * currentLegendrePolynomial;
}

//! Compute derivative of geodesy-normalized Legendre polynomial.
double computeGeodesyLegendrePolynomialDerivative( const int degree,
                                                   const int order,
                                                   const double polynomialParameter,
                                                   const double currentLegendrePolynomial,
                                                   const double incrementedLegendrePolynomial )
{
    // Compute normalization correction factor.
    double normalizationCorrection = std::sqrt(
                ( static_cast< double >( degree + order + 1 ) )
                * ( static_cast< double >( degree - order ) ) );

    // If order is zero apply multiplication factor.
    if ( order == 0 )
    {
        normalizationCorrection *= std::sqrt( 0.5 );
    }

    // Return polynomial derivative.
    return computeGeodesyLegendrePolynomialDerivative(
                degree, order, polynomialParameter, currentLegendrePolynomial,
                incrementedLegendrePolynomial, normalizationCorrection );
}

//! Compute second derivative of geodesy-normalized associated Legendre polynomial.
double computeGeodesyLegendrePolynomialSecondDerivative( const int degree,
                                                         const int order,
                                                         const double polynomialParameter,
                                                         const double currentLegendrePolynomial,
                                                         const double incrementedLegendrePolynomial,
                                                         const double currentLegendrePolynomialDerivative,
                                                         const double incrementedLegendrePolynomialDerivative,
                                                         const double normalizationCorrection )
{
    double polynomialParameterSquare = polynomialParameter * polynomialParameter;

    // Return polynomial derivative.
    return normalizationCorrection * (
                incrementedLegendrePolynomialDerivative / std::sqrt( 1.0 - polynomialParameterSquare ) +
                polynomialParameter * std::pow( 1.0 - polynomialParameterSquare, -1.5 ) * incrementedLegendrePolynomial ) -
            static_cast< double >( order ) *
            ( polynomialParameter / ( 1.0 - polynomialParameter * polynomialParameter ) * currentLegendrePolynomialDerivative +
              ( 1.0 + polynomialParameterSquare ) / ( ( 1.0 - polynomialParameterSquare ) * ( 1.0 - polynomialParameterSquare ) ) * currentLegendrePolynomial );
}


//! Compute low degree/order unnormalized Legendre polynomial explicitly.
double computeLegendrePolynomialExplicit( const int degree,
                                          const int order,
                                          const double polynomialParameter )
{
    // Check which order is required for Legendre polynomial.
    switch( degree )
    {
    case 0:
        switch( order )
        {
        case 0:
            return 1.0;
        default:
        {
            std::string errorMessage = "Error, explicit legendre polynomial not possible for " +
                    boost::lexical_cast< std::string >( degree ) + ", " +
                    boost::lexical_cast< std::string >( order );
            throw std::runtime_error( errorMessage );
        }
        }
        break;
    case 1:
        switch( order )
        {
        case 0:
            return polynomialParameter;
        case 1:
            return std::sqrt( 1 - polynomialParameter * polynomialParameter );
        default:
        {
            std::string errorMessage = "Error, explicit legendre polynomial not possible for " +
                    boost::lexical_cast< std::string >( degree ) + ", " +
                    boost::lexical_cast< std::string >( order );
            throw std::runtime_error( errorMessage );
        }
        }
        break;
    case 2:
        switch( order )
        {
        case 0:
            return 0.5 * ( 3.0 * polynomialParameter * polynomialParameter - 1.0 );
        case 1:
            return 3.0 * polynomialParameter
                    * std::sqrt( 1.0 - polynomialParameter * polynomialParameter );
        case 2:
            return 3.0 * ( 1.0 - polynomialParameter * polynomialParameter );
        default:
        {
            std::string errorMessage = "Error, explicit legendre polynomial not possible for " +
                    boost::lexical_cast< std::string >( degree ) + ", " +
                    boost::lexical_cast< std::string >( order );
            throw std::runtime_error( errorMessage );
        }
        }
        break;
    case 3:
        switch( order )
        {
        case 0:
            return 0.5 * polynomialParameter
                    * ( 5.0 * polynomialParameter * polynomialParameter - 3.0 );
        case 1:
            return 1.5 * ( 5.0 * polynomialParameter * polynomialParameter - 1.0 )
                    * std::sqrt( 1.0 - polynomialParameter * polynomialParameter );
        case 2:
            return 15.0 * polynomialParameter * ( 1.0 - polynomialParameter * polynomialParameter );
        case 3:
            return 15.0 * ( 1.0 - polynomialParameter * polynomialParameter )
                    * std::sqrt( 1.0 - polynomialParameter * polynomialParameter );
        default:
        {
            std::string errorMessage = "Error, explicit legendre polynomial not possible for " +
                    boost::lexical_cast< std::string >( degree ) + ", " +
                    boost::lexical_cast< std::string >( order );
            throw std::runtime_error( errorMessage );
        }
        }
        break;
    case 4:
        switch( order )
        {
        case 0:
            return ( 35.0 * polynomialParameter * polynomialParameter
                     * polynomialParameter * polynomialParameter
                     - 30.0 * polynomialParameter * polynomialParameter + 3.0 ) / 8.0;
        case 1:
            return -2.5 * ( 7.0 * polynomialParameter * polynomialParameter * polynomialParameter
                            - 3.0 * polynomialParameter )
                    * std::sqrt( 1.0 - polynomialParameter * polynomialParameter );
        case 2:
            return 15.0 / 2.0 * ( - 1.0 + 7.0 * polynomialParameter * polynomialParameter )
                    * ( 1.0 - polynomialParameter * polynomialParameter );
        case 3:
            return -105.0 * polynomialParameter * ( 1.0 - polynomialParameter * polynomialParameter )
                    * std::sqrt( 1.0 - polynomialParameter * polynomialParameter );
        case 4:
            return 105.0 * ( 1.0 - polynomialParameter * polynomialParameter )
                    * ( 1.0 - polynomialParameter * polynomialParameter );

        default:
        {
            std::string errorMessage = "Error, explicit legendre polynomial not possible for " +
                    boost::lexical_cast< std::string >( degree ) + ", " +
                    boost::lexical_cast< std::string >( order );
            throw std::runtime_error( errorMessage );
        }
        }
        break;
    default:
    {
        std::string errorMessage = "Error, explicit legendre polynomial not possible for " +
                boost::lexical_cast< std::string >( degree ) + ", " +
                boost::lexical_cast< std::string >( order );
        throw std::runtime_error( errorMessage );
    }
    }
    return TUDAT_NAN;
}

//! Compute low degree/order geodesy-normalized Legendre polynomials explicitly.
double computeGeodesyLegendrePolynomialExplicit( const int degree,
                                                 const int order,
                                                 const double polynomialParameter )
{
    // If 0,0 term is requested return Legendre polynomial value.
    if ( degree == 0 && order == 0 )
    {
        return 1.0;
    }

    // Else if 1,0 term is requested return polynomial value.
    else if ( degree == 1 && order == 0 )
    {
        return std::sqrt( 3.0 ) * polynomialParameter;
    }

    // Else if 1,1 term is requested return polynomial value.
    else if ( degree == 1 && order == 1 )
    {
        return std::sqrt( 3.0 - 3.0 * polynomialParameter * polynomialParameter );
    }

    // Else the requested term cannot be computed; throw a run-time error.
    else
    {
        // Set error message.
        std::stringstream errorMessage;
        errorMessage  <<  "Error: computation of Legendre polynomial of = "  <<  degree
                       <<  " and order = "  <<  order  <<  " is not supported."  <<  std::endl;

        // Throw a run-time error.
        throw std::runtime_error( errorMessage.str( ) );
    }
}

//! Compute unnormalized Legendre polynomial through sectoral recursion.
double computeLegendrePolynomialDiagonal( const int degree,
                                          const double degreeOneOrderOnePolynomial,
                                          const double priorSectoralPolynomial )
{
    // Return polynomial.
    return ( 2.0 * static_cast< double >( degree ) - 1.0 )
            * degreeOneOrderOnePolynomial * priorSectoralPolynomial;

}

//! Compute geodesy-normalized Legendre polynomial through sectoral recursion.
double computeGeodesyLegendrePolynomialDiagonal( const int degree,
                                                 const double degreeOneOrderOnePolynomial,
                                                 const double priorSectoralPolynomial )
{
    // Return polynomial.
    return std::sqrt( ( 2.0 * static_cast< double >( degree ) + 1.0 )
                      / ( 6.0 * static_cast< double >( degree ) ) )
            * degreeOneOrderOnePolynomial * priorSectoralPolynomial;
}

//! Compute unnormalized Legendre polynomial through degree recursion.
double computeLegendrePolynomialVertical( const int degree,
                                          const int order,
                                          const double polynomialParameter,
                                          const double oneDegreePriorPolynomial,
                                          const double twoDegreesPriorPolynomial )
{
    // Return polynomial.
    return ( ( 2.0 * static_cast< double >( degree ) - 1.0 ) * polynomialParameter
             * oneDegreePriorPolynomial - ( static_cast< double >( degree + order ) - 1.0 )
             * twoDegreesPriorPolynomial ) / ( static_cast< double >( degree - order ) );
}

//! Compute geodesy-normalized Legendre polynomial through degree recursion.
double computeGeodesyLegendrePolynomialVertical( const int degree,
                                                 const int order,
                                                 const double polynomialParameter,
                                                 const double oneDegreePriorPolynomial,
                                                 const double twoDegreesPriorPolynomial )
{
    // Return polynomial.
    return std::sqrt( ( 2.0 * static_cast< double >( degree ) + 1.0 )
                      / ( ( static_cast< double >( degree + order ) ) * ( static_cast< double >( degree - order ) ) ) )
            * ( std::sqrt( 2.0 * static_cast< double >( degree ) - 1.0 ) * polynomialParameter * oneDegreePriorPolynomial
                - std::sqrt( ( static_cast< double >( degree + order ) - 1.0 )
                             * ( static_cast< double >( degree - order ) - 1.0 )
                             / ( 2.0 * static_cast< double >( degree ) - 3.0 ) )
                * twoDegreesPriorPolynomial );
}

//! Function to calculate the normalization factor for Legendre polynomials to geodesy-normalized.
double calculateLegendreGeodesyNormalizationFactor( const int degree, const int order )
{

    double deltaFunction = 0.0;
    if( order == 0 )
    {
        deltaFunction = 1.0;
    }

    double factor = std::sqrt(
                boost::math::factorial< double >( static_cast< double >( degree + order ) )
                / ( ( 2.0 - deltaFunction ) * ( 2.0 * static_cast< double >( degree ) + 1.0 )
                    * boost::math::factorial< double >( static_cast< double >( degree - order ) ) ) );
    return 1.0 / factor;
}

} // namespace basic_mathematics
} // namespace tudat