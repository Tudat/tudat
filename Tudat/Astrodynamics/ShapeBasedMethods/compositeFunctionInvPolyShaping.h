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

#ifndef COMPOSITEFUNCTIONINVPOLYSHAPING_H
#define COMPOSITEFUNCTIONINVPOLYSHAPING_H

#include "Tudat/Astrodynamics/ShapeBasedMethods/baseFunctionsInvPolyShaping.h"
#include <math.h>
#include <vector>
#include <Eigen/Dense>
#include <iostream>

namespace tudat
{
namespace shape_based_methods
{

class CompositeFunctionInvPolyShaping
{
public:

    //! Constructor.
    CompositeFunctionInvPolyShaping( ){ }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~CompositeFunctionInvPolyShaping( ) { }

    //! Return coefficients.
    /*!
     * Returns coefficients of composite function.
     */
    Eigen::VectorXd getCompositeFunctionCoefficients()
    {
        return compositeFunctionCoefficients_;
    }

    std::vector< std::shared_ptr< BaseFunctionInvPolyShaping > > getCompositeComponents()
    {
        return compositeFunctionComponents_;
    }

    std::shared_ptr< BaseFunctionInvPolyShaping > getCompositeComponent( int indexComponent )
    {
        return compositeFunctionComponents_[ indexComponent ];
    }

    //! Reset coefficients.
    /*!
     * Resets the coefficients of the composite function.
     */
    void resetCompositeFunctionCoefficients( Eigen::VectorXd compositeFunctionCoefficients );



    double evaluateCompositeFunction( const double independentVariable );

    double evaluateCompositeFunctionFirstDerivative( const double independentVariable );

    double evaluateCompositeFunctionSecondDerivative( const double independentVariable );

    double evaluateCompositeFunctionThirdDerivative( const double independentVariable );

    double getComponentFunctionCurrentValue( int componentIndex, double currentTime );

    double getComponentFunctionFirstDerivative( int componentIndex, double currentTime );

    double getComponentFunctionSecondDerivative( int componentIndex, double currentTime );

    double getComponentFunctionThirdDerivative( int componentIndex, double currentTime );


protected:

private:

    std::vector< std::shared_ptr< BaseFunctionInvPolyShaping > > compositeFunctionComponents_;

    //! Vector of coefficients.
    /*!
     * VectorXd containing coefficients corresponding to components of the composite function.
     */
    Eigen::VectorXd compositeFunctionCoefficients_;

};

class CompositeRadialFunctionInvPolyShaping : public CompositeFunctionInvPolyShaping
{
public:

    //! Constructor.
    CompositeRadialFunctionInvPolyShaping( Eigen::VectorXd& compositeFunctionCoefficients ){

        compositeFunctionCoefficients_ = compositeFunctionCoefficients;

        std::vector< std::shared_ptr< BaseFunctionInvPolyShaping > > compositeFunctionComponents;
        compositeFunctionComponents.push_back( createBaseFunctionInvPolyShaping( constantInvPolyShaping ) );
        compositeFunctionComponents.push_back( createBaseFunctionInvPolyShaping( linearInvPolyShaping ) );
        compositeFunctionComponents.push_back( createBaseFunctionInvPolyShaping( squaredInvPolyShaping ) );
        compositeFunctionComponents.push_back( createBaseFunctionInvPolyShaping( cubedInvPolyShaping ) );
        compositeFunctionComponents.push_back( createBaseFunctionInvPolyShaping( quarticInvPolyShaping ) );
        compositeFunctionComponents.push_back( createBaseFunctionInvPolyShaping( quinticInvPolyShaping ) );
        compositeFunctionComponents.push_back( createBaseFunctionInvPolyShaping( sexticInvPolyShaping ) );
        compositeFunctionComponents_ = compositeFunctionComponents;

    }


    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~CompositeRadialFunctionInvPolyShaping( ) { }

    //! Return coefficients.
    /*!
     * Returns coefficients of composite function.
     */
    Eigen::VectorXd getCompositeFunctionCoefficients()
    {
        return compositeFunctionCoefficients_;
    }

    //! Reset coefficients.
    /*!
     * Resets the coefficients of the composite function.
     */
    void resetCompositeFunctionCoefficients( Eigen::VectorXd compositeFunctionCoefficients );

    double evaluateCompositeFunction( const double independentVariable );

    double evaluateCompositeFunctionFirstDerivative( const double independentVariable );

    double evaluateCompositeFunctionSecondDerivative( const double independentVariable );

    double evaluateCompositeFunctionThirdDerivative( const double independentVariable );

    double getComponentFunctionCurrentValue( int componentIndex, double currentTime );

    double getComponentFunctionFirstDerivative( int componentIndex, double currentTime );

    double getComponentFunctionSecondDerivative( int componentIndex, double currentTime );

    double getComponentFunctionThirdDerivative( int componentIndex, double currentTime );


protected:

private:

    std::vector<  std::shared_ptr< BaseFunctionInvPolyShaping > > compositeFunctionComponents_;

    //! Vector of coefficients.
    /*!
     * VectorXd containing coefficients corresponding to components of the composite function.
     */
    Eigen::VectorXd compositeFunctionCoefficients_;

};


class CompositeZFunctionInvPolyShaping : public CompositeFunctionInvPolyShaping
{
public:

    //! Constructor.
    CompositeZFunctionInvPolyShaping( Eigen::VectorXd& compositeFunctionCoefficients ){

        compositeFunctionCoefficients_ = compositeFunctionCoefficients;

        std::vector< std::shared_ptr< BaseFunctionInvPolyShaping > > compositeFunctionComponents;
        compositeFunctionComponents.push_back( createBaseFunctionInvPolyShaping( constantInvPolyShaping ) );
        compositeFunctionComponents.push_back( createBaseFunctionInvPolyShaping( linearInvPolyShaping ) );
        compositeFunctionComponents.push_back( createBaseFunctionInvPolyShaping( sexticInvPolyShaping ) );
        compositeFunctionComponents.push_back( createBaseFunctionInvPolyShaping( septicInvPolyShaping ) );
//        compositeFunctionComponents.push_back( createBaseFunctionInvPolyShaping( squaredInvPolyShaping ) ); // q = 3
//        compositeFunctionComponents.push_back( createBaseFunctionInvPolyShaping( cubedInvPolyShaping ) );   // q = 3
        compositeFunctionComponents_ = compositeFunctionComponents;

    }


    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~CompositeZFunctionInvPolyShaping( ) { }

    //! Return coefficients.
    /*!
     * Returns coefficients of composite function.
     */
    Eigen::VectorXd getCompositeFunctionCoefficients()
    {
        return compositeFunctionCoefficients_;
    }

    //! Reset coefficients.
    /*!
     * Resets the coefficients of the composite function.
     */
    void resetCompositeFunctionCoefficients( Eigen::VectorXd compositeFunctionCoefficients );

    double evaluateCompositeFunction( const double independentVariable );

    double evaluateCompositeFunctionFirstDerivative( const double independentVariable );

    double evaluateCompositeFunctionSecondDerivative( const double independentVariable );

    double evaluateCompositeFunctionThirdDerivative( const double independentVariable );

    double getComponentFunctionCurrentValue( int componentIndex, double currentTime );

    double getComponentFunctionFirstDerivative( int componentIndex, double currentTime );

    double getComponentFunctionSecondDerivative( int componentIndex, double currentTime );

    double getComponentFunctionThirdDerivative( int componentIndex, double currentTime );


protected:

private:

    std::vector<  std::shared_ptr< BaseFunctionInvPolyShaping > > compositeFunctionComponents_;

    //! Vector of coefficients.
    /*!
     * VectorXd containing coefficients corresponding to components of the composite function.
     */
    Eigen::VectorXd compositeFunctionCoefficients_;

};


} // namespace shape_based_methods
} // namespace tudat

#endif // COMPOSITEFUNCTIONINVPOLYSHAPING_H
