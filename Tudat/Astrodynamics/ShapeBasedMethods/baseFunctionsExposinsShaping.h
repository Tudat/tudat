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

#ifndef BASEFUNCTIONSEXPOSINS_H
#define BASEFUNCTIONSEXPOSINS_H

#include <math.h>
#include <boost/make_shared.hpp>

namespace tudat
{
namespace shape_based_methods
{

enum baseFunctionExposinsShapingType
{
    constantExposinsShaping,
    linearExposinsShaping,
    squaredExposinsShaping,
    cosineExposinsShaping,
    powerCosineExposinsShaping,
    sineExposinsShaping,
    powerSineExposinsShaping
};

//! Base class for shaping functions to be used in the exposins shaping method.
class BaseFunctionExposinsShaping
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    BaseFunctionExposinsShaping( ){ }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~BaseFunctionExposinsShaping( ) { }

    virtual double evaluateFunction( double independentVariable ) = 0;

    virtual double evaluateFirstDerivative( double independentVariable ) = 0;

    virtual double evaluateSecondDerivative( double independentVariable ) = 0;

    virtual double evaluateThirdDerivative( double independentVariable ) = 0;

protected:

private:

};


//! Constant function.
class ConstantFunctionExposinsShaping : public BaseFunctionExposinsShaping
{
public:
    //! Default constructor.
    /*!
     * Default constructor.
     */
    ConstantFunctionExposinsShaping( ){}

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~ConstantFunctionExposinsShaping( ) { }

    double evaluateFunction( double independentVariable )
    {
        return 1.0;
    }

    double evaluateFirstDerivative( double independentVariable )
    {
        return 0.0;
    }

    double evaluateSecondDerivative( double independentVariable )
    {
        return 0.0;
    }

    double evaluateThirdDerivative( double independentVariable )
    {
        return 0.0;
    }

protected:

private:

};


//! Linear function.
class LinearFunctionExposinsShaping : public BaseFunctionExposinsShaping
{
public:
    //! Default constructor.
    /*!
     * Default constructor.
     */
    LinearFunctionExposinsShaping( ) { }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~LinearFunctionExposinsShaping( ) { }


    double evaluateFunction( double independentVariable )
    {
        return independentVariable;
    }

    double evaluateFirstDerivative( double independentVariable )
    {
        return 1.0;
    }

    double evaluateSecondDerivative( double independentVariable )
    {
        return 0.0;
    }

    double evaluateThirdDerivative( double independentVariable )
    {
        return 0.0;
    }


protected:

private:

};


//! x^2 function.
class SquaredFunctionExposinsShaping : public BaseFunctionExposinsShaping
{
public:
    //! Default constructor.
    /*!
     * Default constructor.
     */
    SquaredFunctionExposinsShaping( ) { }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~SquaredFunctionExposinsShaping( ) { }

    double evaluateFunction( double independentVariable ){
        return std::pow( independentVariable, 2.0 );
    }

    double evaluateFirstDerivative( double independentVariable )
    {
        return 2.0 * independentVariable;
    }

    double evaluateSecondDerivative( double independentVariable ){
        return 2.0;
    }

    double evaluateThirdDerivative( double independentVariable )
    {
        return 0.0;
    }


protected:

private:

};


//! cosine function.
class CosineFunctionExposinsShaping : public BaseFunctionExposinsShaping
{
public:
    //! Default constructor.
    /*!
     * Default constructor.
     */
    CosineFunctionExposinsShaping( ) { }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~CosineFunctionExposinsShaping( ) { }


    double evaluateFunction( double independentVariable )
    {
        return std::cos( independentVariable );
    }

    double evaluateFirstDerivative( double independentVariable )
    {
        return - std::sin( independentVariable );
    }

    double evaluateSecondDerivative( double independentVariable )
    {
        return - std::cos( independentVariable );
    }

    double evaluateThirdDerivative( double independentVariable )
    {
        return std::sin( independentVariable );
    }

protected:

private:

};


//! Power cosine function.
class PowerCosineFunctionExposinsShaping : public BaseFunctionExposinsShaping
{
public:
    //! Default constructor.
    /*!
     * Default constructor.
     */
    PowerCosineFunctionExposinsShaping( ) { }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~PowerCosineFunctionExposinsShaping( ) { }


    double evaluateFunction( double independentVariable ){
        return independentVariable * std::cos( independentVariable );
    }

    double evaluateFirstDerivative( double independentVariable )
    {
        return std::cos( independentVariable ) - independentVariable * std::sin( independentVariable );
    }

    double evaluateSecondDerivative( double independentVariable )
    {
        return - 2.0 * std::sin( independentVariable ) - independentVariable * std::cos( independentVariable );
    }

    double evaluateThirdDerivative( double independentVariable )
    {
        return - 3.0 * std::cos( independentVariable ) + independentVariable * std::sin( independentVariable );
    }

protected:

private:

};


//! Sine function.
class SineFunctionExposinsShaping : public BaseFunctionExposinsShaping
{
public:
    //! Default constructor.
    /*!
     * Default constructor.
     */
    SineFunctionExposinsShaping( ) { }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~SineFunctionExposinsShaping( ) { }


    double evaluateFunction( double independentVariable )
    {
        return std::sin( independentVariable );
    }

    double evaluateFirstDerivative( double independentVariable )
    {
        return std::cos( independentVariable );
    }

    double evaluateSecondDerivative( double independentVariable )
    {
        return - std::sin( independentVariable );
    }

    double evaluateThirdDerivative( double independentVariable )
    {
        return - std::cos( independentVariable );
    }

protected:

private:

};


//! Power times sine function.
class PowerSineFunctionExposinsShaping : public BaseFunctionExposinsShaping
{
public:
    //! Default constructor.
    /*!
     * Default constructor.
     */
    PowerSineFunctionExposinsShaping( ) { }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~PowerSineFunctionExposinsShaping( ) { }


    double evaluateFunction( double independentVariable )
    {
        return independentVariable * std::sin( independentVariable );
    }

    double evaluateFirstDerivative( double independentVariable ){
        return std::sin( independentVariable ) + independentVariable * std::cos( independentVariable );
    }

    double evaluateSecondDerivative( double independentVariable )
    {
        return 2.0 * std::cos( independentVariable ) - independentVariable * std::sin( independentVariable );
    }

    double evaluateThirdDerivative( double independentVariable )
    {
        return - 3.0 * std::sin( independentVariable ) - independentVariable * std::cos( independentVariable );
    }

protected:

private:

};

std::shared_ptr< BaseFunctionExposinsShaping > createBaseFunctionExposinsShaping(
        const baseFunctionExposinsShapingType baseFunctionType );


} // namespace shape_based_methods
} // namespace tudat

#endif // BASEFUNCTIONSEXPOSINS_H
