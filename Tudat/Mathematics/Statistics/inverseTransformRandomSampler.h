#ifndef TUDAT_INVERSE_TRANSFORM_SAMPLER_H
#define TUDAT_INVERSE_TRANSFORM_SAMPLER_H

//// ================= State Conversion ================= //
/// Author  :   R. Hoogendoorn
/// Date    :   25-03-2016

#include <iostream> // cout sometimes needs this
#include <map>

#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

// Pseudo random sampler
#include <boost/random.hpp> // Random generator

#include <Tudat/Mathematics/BasicMathematics/function.h>
#include <Tudat/Mathematics/BasicMathematics/mathematicalConstants.h>
#include <Tudat/Mathematics/RootFinders/bisection.h>

namespace tudat
{

namespace statistics
{

//! Root Finder Function
/*!
 * Function that is used in the root finder. This function is used to invert F(x).
 * Function : F(x) - u = 0
 */
class InverseTransformRootFinderFunction: public tudat::basic_mathematics::Function<double,double>
{
private:

    //! Typedef for pointer to
    typedef boost::function< double ( const double& ) > CDFFunctionType;

public:

    //! Constructor
    InverseTransformRootFinderFunction( CDFFunctionType cdfFunction ){
        cumulativeDistributionFunction_ = cdfFunction ;
    }

    //! Destructor
    ~InverseTransformRootFinderFunction(){}

    //! Set the CDF
    void setCumulativeDistributionFunction( CDFFunctionType cdfFunction ){
        cumulativeDistributionFunction_ = cdfFunction ;
    }

    //! Set the uniform sample
    void setUniformSample( double UniformSample ){
        uniformSample_ = UniformSample;
    }

    //! Evaluate function ( F(x) = u --> F(x) - u = 0 )
    double evaluate( const double inputValue ){
        return cumulativeDistributionFunction_( inputValue ) - uniformSample_ ;
    }

    //! Compute State Derivative
    //! FUNCTION NOT IMPLEMENTED
    double computeDerivative( const unsigned int order, const double independentVariable )
    {
        return TUDAT_NAN;
    }

    //! Compute definite integral
    //! FUNCTION NOT IMPLEMENTED
    double computeDefiniteIntegral( const unsigned int order, const double lowerBound, const double upperbound )
    {
        return TUDAT_NAN;
    }

protected:

private:

    //! Cumulative distribution function
    CDFFunctionType cumulativeDistributionFunction_;

    //! Uniform sample
    double uniformSample_;
};

//! Inverse Transform Random Sampler
/*!
 * This class uses a Cumulative Distribution Function (CDF) and a rootfinder to generate
 * random samples with the distribution of the CDF.
 */
class InverseTransformRandomSampler
{
private:

    //! Typedef for pointer to
    typedef boost::function< double ( const double& ) > CDFFunctionType;

    //! Typedef for Random generator
    typedef boost::mt19937 RandomGeneratorType; // Mersenne Twister

    //! Typedef for uniform sample generator
    typedef boost::variate_generator< RandomGeneratorType, boost::uniform_real<> > uniformGeneratorType;

    typedef boost::shared_ptr< uniformGeneratorType > uniformGeneratorPointer;

public:

    //! Default constructor
    InverseTransformRandomSampler( CDFFunctionType cdfFunction , int Seed,
                   double LowerBoundRootFinder, double UpperBoundRootFinder,
                   double InitialGuess, double RootFinderRelativeTolerance = 1E-10 );

    //! Set the CDF
    void setCumulativeDistributionFunction( CDFFunctionType cdfFunction ,
                                            double LowerBoundRootFinder ,
                                            double UpperBoundRootFinder ,
                                            double InitialGuess );

    //! Set the CDF
    void setCumulativeDistributionFunction( CDFFunctionType cdfFunction);

    //! Generate random sample
    double generateRandomSample();

    //! Generate random sample
    double generateRandomSample(double uniformSample);

protected:

private:

    //! Cumulative Distribution Function (CDF)
    CDFFunctionType cumulativeDistributionFunction_;

    //! Uniform sampler
    uniformGeneratorPointer uniformSampler_;

    //! Seed of the uniform sampler
    int seed_;

    //! Typedef of the function whose root we have to determine.
    boost::shared_ptr< InverseTransformRootFinderFunction > rootFinderFunctionPointer_;

    //! Root finder pointer
    tudat::root_finders::BisectionPointer rootFinder_;

    //! (Relative?) Accuracy of the rootfinder
    double rootFinderRelativeTolerance_;

    //! Lowerbound of the rootfinder
    double lowerBoundRootFinder_;

    //! Upperbound of the rootfinder
    double upperBoundRootFinder_;

    //! Initial guess of rootfinder
    double initialGuessRootFinder_;
};

typedef boost::shared_ptr< InverseTransformRandomSampler > InverseTransformRandomSamplerPointer;

} // Close Namespace statistics

} // Close Namespace tudat

#endif // TUDAT_INVERSE_TRANSFORM_SAMPLER_H
