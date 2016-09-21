/*    Copyright (c) 2016, R. Hoogendoorn
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
 *      160503    R.Hoogendoorn          First creation of code.
 *
 *    References
 *
 *
 *    Notes
 *
 */

//// ================= Inverse Transform Random Sampler ================= //
/// Author  :   R. Hoogendoorn
/// Date    :   19-05-2016

#include <Tudat/Mathematics/Statistics/inverseTransformRandomSampler.h>

namespace tudat
{

namespace statistics
{

//! Default constructor
InverseTransformRandomSampler::InverseTransformRandomSampler( CDFFunctionType cdfFunction , int Seed,
               double LowerBoundRootFinder, double UpperBoundRootFinder,
               double InitialGuess, double RootFinderRelativeTolerance )
{
    // Declare distribution
    cumulativeDistributionFunction_ = cdfFunction ;

    // Construct uniform sampler
    seed_ = Seed;
    RandomGeneratorType randomGenerator( seed_ );
    boost::uniform_real<> uniformDistribution(0.0 , 1.0);
    uniformSampler_ = boost::make_shared< uniformGeneratorType >( randomGenerator, uniformDistribution );

    // Construct rootfinder function
    rootFinderFunctionPointer_ = boost::make_shared< InverseTransformRootFinderFunction >( cumulativeDistributionFunction_ ) ;

    // Construct rootfinder
    rootFinderRelativeTolerance_ = RootFinderRelativeTolerance ;
    lowerBoundRootFinder_ = LowerBoundRootFinder ;
    upperBoundRootFinder_ = UpperBoundRootFinder ;
    initialGuessRootFinder_ = InitialGuess ;

    using namespace tudat::root_finders;
    using namespace tudat::root_finders::termination_conditions;

    unsigned int maximumNumberOfIterations = 1e3;
    Bisection::TerminationFunction terminationConditionFunction =
            boost::bind( &RootRelativeToleranceTerminationCondition< double >::checkTerminationCondition,
                         boost::make_shared< RootRelativeToleranceTerminationCondition< double > >(
                             rootFinderRelativeTolerance_ , maximumNumberOfIterations ), _1, _2, _3, _4, _5 );

    rootFinder_ = boost::make_shared< Bisection >(
                terminationConditionFunction, lowerBoundRootFinder_, upperBoundRootFinder_ );
}

//! Set the CDF
void InverseTransformRandomSampler::setCumulativeDistributionFunction( CDFFunctionType cdfFunction){
    // Change CDF
    cumulativeDistributionFunction_ = cdfFunction ;
    rootFinderFunctionPointer_->setCumulativeDistributionFunction( cumulativeDistributionFunction_ );
}

//! Set the CDF
void InverseTransformRandomSampler::setCumulativeDistributionFunction( CDFFunctionType cdfFunction ,
                                        double LowerBoundRootFinder ,
                                        double UpperBoundRootFinder ,
                                        double InitialGuess ){
    // Change CDF
    cumulativeDistributionFunction_ = cdfFunction ;
    rootFinderFunctionPointer_->setCumulativeDistributionFunction( cumulativeDistributionFunction_ );

    // Change root finder boundaries
    lowerBoundRootFinder_ = LowerBoundRootFinder ;
    upperBoundRootFinder_ = UpperBoundRootFinder ;
    initialGuessRootFinder_ = InitialGuess ;
    rootFinder_->resetBoundaries( lowerBoundRootFinder_ , upperBoundRootFinder_ );
}

//! Generate random sample
double InverseTransformRandomSampler::generateRandomSample(){

    // Generate uniformly distributed random sample
    double uniformSample = (*uniformSampler_)() ;

    // Find (x) root of F(x) - u = 0  ( Invert: F^-1(u) = x )
    rootFinderFunctionPointer_->setUniformSample( uniformSample );
    return rootFinder_->execute( rootFinderFunctionPointer_ , initialGuessRootFinder_ ) ;
}

//! Generate random sample
double InverseTransformRandomSampler::generateRandomSample(double uniformSample){

    // Find (x) root of F(x) - u = 0  ( Invert: F^-1(u) = x )
    rootFinderFunctionPointer_->setUniformSample( uniformSample );
    return rootFinder_->execute( rootFinderFunctionPointer_ , initialGuessRootFinder_ ) ;
}

} // Close Namespace statistics

} // Close Namespace tudat

