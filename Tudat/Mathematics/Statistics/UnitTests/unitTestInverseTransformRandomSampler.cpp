/*    Copyright (c) 2015, R. Hoogendoorn
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120522    R.Hoogendoorn          First creation of code.
 *
 *    References
 *
 *
 *    Notes
 *
 */

//// ================= Unit Test Function ================= //
/// Author  :   R. Hoogendoorn
/// Date    :   03-05-2016

#define BOOST_TEST_MAIN

#include <iostream> // cout sometimes needs this
#include <vector>
#include <limits>

#include <Eigen/Core>

#include <boost/test/unit_test.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/timer/timer.hpp> // Computation time

#include <boost/math/distributions/normal.hpp>

#include "Tudat/Basics/testMacros.h"

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include <Tudat/InputOutput/basicInputOutput.h>
#include <Tudat/Mathematics/Statistics/boostProbabilityDistributions.h>
#include <Tudat/Mathematics/Statistics/basicStatistics.h>

#include <Tudat/Mathematics/Statistics/inverseTransformRandomSampler.h>

//#include <SIP/Statistical_Methods/RandomInitialStateSampling.h>
//#include <SIP/Statistical_Methods/RandomSampling.h>
//#include <SIP/Tools/Tools.h>

class GaussianDistributionNPeaks
{
public:

    //! Constructor.
    GaussianDistributionNPeaks( Eigen::VectorXd Mean, Eigen::VectorXd StandardDeviation )
    {
        using namespace tudat::statistics;
        mean_ = Mean;
        standardDeviation_ = StandardDeviation;
        numberOfDistributions = mean_.rows();

        std::vector< double > parameters(2);

        for( int i = 0 ; i < numberOfDistributions ; i++ )
        {
            parameters[0] = mean_(i);
            parameters[1] = standardDeviation_(i);

            distributionPointers.push_back(
                        createBoostRandomVariable( normal_boost_distribution, parameters ) );
        }
    }

    //! Get probability density.
    double getProbabilityDensity( const double& x )
    {
        double probabilityDensity = 0.0 ;
        for( int i = 0 ; i < numberOfDistributions ; i++ )
        {
            probabilityDensity += distributionPointers[i]->evaluatePdf( x );
        }
        return probabilityDensity / ( static_cast< double >( numberOfDistributions ) );
    }

    //! Get cumulative probability.
    double getCumulativeProbability( const double& x )
    {
        double cumulativeProbability = 0.0 ;
        for( int i = 0 ; i < numberOfDistributions ; i++ )
        {
            cumulativeProbability += distributionPointers[i]->evaluateCdf( x );
        }
        return cumulativeProbability / ( static_cast< double >( numberOfDistributions ) );
    }

    //! Set the properties (mean and standard deviation) of the Gaussian distribution
    void setProperties( const Eigen::VectorXd& mean, Eigen::VectorXd standardDeviation )
    {
        mean_ = mean;
        standardDeviation_ = standardDeviation;
    }

protected:
private:
    //! Mean value of the distribution.
    Eigen::VectorXd mean_ ;

    //! Standard deviation value of the distribution.
    Eigen::VectorXd standardDeviation_ ;

    //! Vector with distribution pointers
    std::vector< boost::shared_ptr<
        tudat::statistics::ContinuousProbabilityDistribution< double > > > distributionPointers;

    //! Number of distributions that are added to model the total
    int numberOfDistributions;
};

BOOST_AUTO_TEST_SUITE( test_Inverse_Transform_Random_Sampler )

BOOST_AUTO_TEST_CASE( test_Inverse_Random_Sampler )
{
    using namespace tudat::statistics;
    int numberOfSamples = 2E4;

    std::vector< double > parameters;
    parameters.push_back( 1.0 );
    parameters.push_back( 2.0 );

    boost::shared_ptr< ContinuousProbabilityDistribution< double > > distributionPointer =
            createBoostRandomVariable(
                normal_boost_distribution, parameters );

    boost::function< double( const double& ) > cdfFunction =
            boost::bind( &ContinuousProbabilityDistribution<double>::evaluateCdf, distributionPointer, _1 ) ;

    InverseTransformRandomSampler randomSampler(
                cdfFunction , 100 , -15.0 , 15.0 , 0.0 , 1E-8) ;

    // Generate samples
    std::vector< double > samples(0);
    for( int i = 0 ; i < numberOfSamples ; i++ ){
        samples.push_back( randomSampler.generateRandomSample() );
    }

    double mean = computeSampleMean( samples );
    double variance = computeSampleVariance( samples ) ;

    BOOST_CHECK_CLOSE_FRACTION( mean , 1.0 , 3E-2 );
    BOOST_CHECK_CLOSE_FRACTION( variance , 4.0 , 2E-3 );
}

BOOST_AUTO_TEST_CASE( test_Inverse_Random_Sampler_2peaks )
{
    using namespace tudat::statistics;
    int numberOfSamples = 1E6;

    Eigen::VectorXd mean(2);
    Eigen::VectorXd standardDeviation(2);
    mean << -1.0 , 1.0;
    standardDeviation << 0.5 , 0.8;

    boost::shared_ptr< GaussianDistributionNPeaks > distributionPointer =
            boost::make_shared< GaussianDistributionNPeaks >( mean , standardDeviation );

    boost::function< double( const double& ) > cdfFunction =
            boost::bind( &GaussianDistributionNPeaks::getCumulativeProbability, distributionPointer, _1 ) ;

    InverseTransformRandomSampler randomSampler(
                cdfFunction , 100 , -15.0 , 15.0 , 0.0 , 1E-9) ;

    std::vector< double > samples(0);
    for( int i = 0 ; i < numberOfSamples ; i++ ){
        samples.push_back( randomSampler.generateRandomSample() );
    }

    double sampleMean = computeSampleMean( samples );
    double sampleVariance = computeSampleVariance( samples ) ;

    // Test with values from Matlab with numerical quadrature
    BOOST_CHECK_SMALL( std::fabs( sampleMean - -1.125766146969909e-13 ) , 1.0E-3 );
    BOOST_CHECK_SMALL( std::fabs( sampleVariance - 1.444999999999195 ) , 1.0E-3 );
}

BOOST_AUTO_TEST_CASE( test_Inverse_Random_Sampler_3peaks )
{
    using namespace tudat::statistics;
    int numberOfSamples = 1E6;

    Eigen::VectorXd mean(3);
    Eigen::VectorXd standardDeviation(3);
    mean << -1.0 , 1.0 , 3.5;
    standardDeviation << 0.5 , 0.8 , 0.3;

    boost::shared_ptr< GaussianDistributionNPeaks > distributionPointer =
            boost::make_shared< GaussianDistributionNPeaks >( mean , standardDeviation );

    boost::function< double( const double& ) > cdfFunction =
            boost::bind( &GaussianDistributionNPeaks::getCumulativeProbability, distributionPointer, _1 ) ;

    InverseTransformRandomSampler randomSampler(
                cdfFunction , 100 , -15.0 , 15.0 , 0.0 , 1E-9) ;

    std::vector< double > samples(0);
    for( int i = 0 ; i < numberOfSamples ; i++ ){
        samples.push_back( randomSampler.generateRandomSample() );
    }

    double sampleMean = computeSampleMean( samples );
    double sampleVariance = computeSampleVariance( samples ) ;

    // Test with values from Matlab with numerical quadrature
    BOOST_CHECK_SMALL( std::fabs( sampleMean - 1.166666666666669 ) , 3.0E-3 );
    BOOST_CHECK_SMALL( std::fabs( sampleVariance - 3.715555555555550 ) , 3.0E-3 );

}

BOOST_AUTO_TEST_SUITE_END( )
