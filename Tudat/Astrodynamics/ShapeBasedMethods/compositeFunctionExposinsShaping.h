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

#ifndef COMPOSITEFUNCTIONEXPOSINSSHAPING_H
#define COMPOSITEFUNCTIONEXPOSINSSHAPING_H

#include "Tudat/Astrodynamics/ShapeBasedMethods/baseFunctionsExposinsShaping.h"
#include <math.h>
#include <vector>
#include <Eigen/Dense>
#include <iostream>

namespace tudat
{
namespace shape_based_methods
{

class CompositeFunctionExposinsShaping
{
public:

    //! Constructor.
    CompositeFunctionExposinsShaping( ){ }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~CompositeFunctionExposinsShaping( ) { }

    //! Return coefficients.
    /*!
     * Returns coefficients of composite function.
     */
    Eigen::Vector4d getCompositeFunctionCoefficients()
    {
        return compositeFunctionCoefficients_;
    }

//    std::vector< std::shared_ptr< BaseFunctionExposinsShaping > > getCompositeComponents()
//    {
//        return compositeFunctionComponents_;
//    }

//    std::shared_ptr< BaseFunctionExposinsShaping > getCompositeComponent( int indexComponent )
//    {
//        return compositeFunctionComponents_[ indexComponent ];
//    }

    //! Reset coefficients.
    /*!
     * Resets the coefficients of the composite function.
     */
    void resetCompositeFunctionCoefficients( Eigen::Vector4d compositeFunctionCoefficients );



    double evaluateCompositeFunction( const double independentVariable );

    double evaluateCompositeFunctionFirstDerivative( const double independentVariable );

    double evaluateCompositeFunctionSecondDerivative( const double independentVariable );

    double evaluateCompositeFunctionThirdDerivative( const double independentVariable );

    double getComponentFunctionCurrentValue( int componentIndex, double currentTime );

//    double getComponentFunctionFirstDerivative( int componentIndex, double currentTime );

//    double getComponentFunctionSecondDerivative( int componentIndex, double currentTime );

//    double getComponentFunctionThirdDerivative( int componentIndex, double currentTime );


protected:

private:

    std::vector< std::shared_ptr< BaseFunctionExposinsShaping > > compositeFunctionComponents_;

    //! Vector of coefficients.
    /*!
     * Vector4d containing coefficients corresponding to components of the composite function.
     */
    Eigen::Vector4d compositeFunctionCoefficients_;

};


class CompositeRadialFunctionExposinsShaping : public CompositeFunctionExposinsShaping
{
public:

    //! Constructor.
    CompositeRadialFunctionExposinsShaping( Eigen::Vector4d& compositeFunctionCoefficients ){

        compositeFunctionCoefficients_ = compositeFunctionCoefficients;

    }


    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~CompositeRadialFunctionExposinsShaping( ) { }

    //! Return coefficients.
    /*!
     * Returns coefficients of composite function.
     */
    Eigen::Vector4d getCompositeFunctionCoefficients()
    {
        return compositeFunctionCoefficients_;
    }

    //! Reset coefficients.
    /*!
     * Resets the coefficients of the composite function.
     */
    void resetCompositeFunctionCoefficients( Eigen::Vector4d compositeFunctionCoefficients );

    double evaluateCompositeFunction( const double independentVariable );

//    double evaluateCompositeFunctionFirstDerivative( const double independentVariable );

//    double evaluateCompositeFunctionSecondDerivative( const double independentVariable );

//    double evaluateCompositeFunctionThirdDerivative( const double independentVariable );

    double getComponentFunctionCurrentValue( int componentIndex, double currentTime );

//    double getComponentFunctionFirstDerivative( int componentIndex, double currentTime );

//    double getComponentFunctionSecondDerivative( int componentIndex, double currentTime );

//    double getComponentFunctionThirdDerivative( int componentIndex, double currentTime );


protected:

private:

    std::vector<  std::shared_ptr< BaseFunctionExposinsShaping > > compositeFunctionComponents_;

    //! Vector of coefficients.
    /*!
     * Vector4d containing coefficients corresponding to components of the composite function.
     */
    Eigen::Vector4d compositeFunctionCoefficients_;

};


//// NOT USED
//class CompositeElevationFunctionExposinsShaping : public CompositeFunctionExposinsShaping
//{
//public:

//    //! Constructor.
//    CompositeElevationFunctionExposinsShaping( Eigen::VectorXd& compositeFunctionCoefficients ){

//        compositeFunctionCoefficients_ = compositeFunctionCoefficients;

//        std::vector< std::shared_ptr< BaseFunctionExposinsShaping > > compositeFunctionComponents;
//        compositeFunctionComponents.push_back( createBaseFunctionExposinsShaping( cosineExposinsShaping ) );
//        compositeFunctionComponents.push_back( createBaseFunctionExposinsShaping( powerCosineExposinsShaping ) );
//        compositeFunctionComponents.push_back( createBaseFunctionExposinsShaping( sineExposinsShaping ) );
//        compositeFunctionComponents.push_back( createBaseFunctionExposinsShaping( powerSineExposinsShaping ) );
//        compositeFunctionComponents_ = compositeFunctionComponents;

//    }


//    //! Default destructor.
//    /*!
//     * Default destructor.
//     */
//    ~CompositeElevationFunctionExposinsShaping( ) { }

//    //! Return coefficients.
//    /*!
//     * Returns coefficients of composite function.
//     */
//    Eigen::VectorXd getCompositeFunctionCoefficients()
//    {
//        return compositeFunctionCoefficients_;
//    }

//    //! Reset coefficients.
//    /*!
//     * Resets the coefficients of the composite function.
//     */
//    void resetCompositeFunctionCoefficients( Eigen::VectorXd compositeFunctionCoefficients );

//    double evaluateCompositeFunction( const double independentVariable );

//    double evaluateCompositeFunctionFirstDerivative( const double independentVariable );

//    double evaluateCompositeFunctionSecondDerivative( const double independentVariable );

//    double evaluateCompositeFunctionThirdDerivative( const double independentVariable );

//    double getComponentFunctionCurrentValue( int componentIndex, double currentTime );

//    double getComponentFunctionFirstDerivative( int componentIndex, double currentTime );

//    double getComponentFunctionSecondDerivative( int componentIndex, double currentTime );

//    double getComponentFunctionThirdDerivative( int componentIndex, double currentTime );


//protected:

//private:

//    std::vector<  std::shared_ptr< BaseFunctionExposinsShaping > > compositeFunctionComponents_;

//    //! Vector of coefficients.
//    /*!
//     * VectorXd containing coefficients corresponding to components of the composite function.
//     */
//    Eigen::VectorXd compositeFunctionCoefficients_;

//};


} // namespace shape_based_methods
} // namespace tudat

#endif // COMPOSITEFUNCTIONEXPOSINSSHAPING_H
