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

#ifndef BASEFUNCTIONSINVPOLY_H
#define BASEFUNCTIONSINVPOLY_H

#include <math.h>
#include <boost/make_shared.hpp>

namespace tudat
{
namespace shape_based_methods
{

enum baseFunctionInvPolyShapingType
{
    constantInvPolyShaping,
    linearInvPolyShaping,
    squaredInvPolyShaping,
    cubedInvPolyShaping,
    quarticInvPolyShaping,
    quinticInvPolyShaping,
    sexticInvPolyShaping,
    septicInvPolyShaping
};

//! Base class for shaping functions to be used in the InvPoly shaping method.
class BaseFunctionInvPolyShaping
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    BaseFunctionInvPolyShaping( ){ }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~BaseFunctionInvPolyShaping( ) { }

    virtual double evaluateFunction( double independentVariable ) = 0;

    virtual double evaluateFirstDerivative( double independentVariable ) = 0;

    virtual double evaluateSecondDerivative( double independentVariable ) = 0;

    virtual double evaluateThirdDerivative( double independentVariable ) = 0;

protected:

private:

};


//! Constant function.
class ConstantFunctionInvPolyShaping : public BaseFunctionInvPolyShaping
{
public:
    //! Default constructor.
    /*!
     * Default constructor.
     */
    ConstantFunctionInvPolyShaping( ){}

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~ConstantFunctionInvPolyShaping( ) { }

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
class LinearFunctionInvPolyShaping : public BaseFunctionInvPolyShaping
{
public:
    //! Default constructor.
    /*!
     * Default constructor.
     */
    LinearFunctionInvPolyShaping( ) { }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~LinearFunctionInvPolyShaping( ) { }


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
class SquaredFunctionInvPolyShaping : public BaseFunctionInvPolyShaping
{
public:
    //! Default constructor.
    /*!
     * Default constructor.
     */
    SquaredFunctionInvPolyShaping( ) { }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~SquaredFunctionInvPolyShaping( ) { }

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


//! x^3 function.
class CubedFunctionInvPolyShaping : public BaseFunctionInvPolyShaping
{
public:
    //! Default constructor.
    /*!
     * Default constructor.
     */
    CubedFunctionInvPolyShaping( ) { }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~CubedFunctionInvPolyShaping( ) { }


    double evaluateFunction( double independentVariable )
    {
        return std::pow( independentVariable, 3.0 );
    }

    double evaluateFirstDerivative( double independentVariable )
    {
        return 3.0 * std::pow( independentVariable, 2.0 );
    }

    double evaluateSecondDerivative( double independentVariable )
    {
        return 3.0 * 2.0 * independentVariable;
    }

    double evaluateThirdDerivative( double independentVariable )
    {
        return 3.0 * 2.0;
    }

protected:

private:

};


//! x^4 function.
class QuarticFunctionInvPolyShaping : public BaseFunctionInvPolyShaping
{
public:
    //! Default constructor.
    /*!
     * Default constructor.
     */
    QuarticFunctionInvPolyShaping( ) { }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~QuarticFunctionInvPolyShaping( ) { }


    double evaluateFunction( double independentVariable )
    {
        return std::pow( independentVariable, 4.0 );
    }

    double evaluateFirstDerivative( double independentVariable )
    {
        return 4.0 * std::pow( independentVariable, 3.0 );
    }

    double evaluateSecondDerivative( double independentVariable )
    {
        return 4.0 * 3.0 * std::pow( independentVariable, 2.0 );
    }

    double evaluateThirdDerivative( double independentVariable )
    {
        return 4.0 * 3.0 * 2.0 * independentVariable;
    }

protected:

private:

};


//! x^5 function.
class QuinticFunctionInvPolyShaping : public BaseFunctionInvPolyShaping
{
public:
    //! Default constructor.
    /*!
     * Default constructor.
     */
    QuinticFunctionInvPolyShaping( ) { }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~QuinticFunctionInvPolyShaping( ) { }


    double evaluateFunction( double independentVariable )
    {
        return std::pow( independentVariable, 5.0 );
    }

    double evaluateFirstDerivative( double independentVariable )
    {
        return 5.0 * std::pow( independentVariable, 4.0 );
    }

    double evaluateSecondDerivative( double independentVariable )
    {
        return 5.0 * 4.0 * std::pow( independentVariable, 3.0 );
    }

    double evaluateThirdDerivative( double independentVariable )
    {
        return 5.0 * 4.0 * 3.0 * std::pow( independentVariable, 2.0 );
    }

protected:

private:

};


//! x^6 function.
class SexticFunctionInvPolyShaping : public BaseFunctionInvPolyShaping
{
public:
    //! Default constructor.
    /*!
     * Default constructor.
     */
    SexticFunctionInvPolyShaping( ) { }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~SexticFunctionInvPolyShaping( ) { }


    double evaluateFunction( double independentVariable )
    {
        return std::pow( independentVariable, 6.0 );
    }

    double evaluateFirstDerivative( double independentVariable )
    {
        return 6.0 * std::pow( independentVariable, 5.0 );
    }

    double evaluateSecondDerivative( double independentVariable )
    {
        return 6.0 * 5.0 * std::pow( independentVariable, 4.0 );
    }

    double evaluateThirdDerivative( double independentVariable )
    {
        return 6.0 * 5.0 * 4.0 * std::pow( independentVariable, 3.0 );
    }

protected:

private:

};

//! x^6 function.
class SepticFunctionInvPolyShaping : public BaseFunctionInvPolyShaping
{
public:
    //! Default constructor.
    /*!
     * Default constructor.
     */
    SepticFunctionInvPolyShaping( ) { }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~SepticFunctionInvPolyShaping( ) { }


    double evaluateFunction( double independentVariable )
    {
        return std::pow( independentVariable, 7.0 );
    }

    double evaluateFirstDerivative( double independentVariable )
    {
        return 7.0 * std::pow( independentVariable, 6.0 );
    }

    double evaluateSecondDerivative( double independentVariable )
    {
        return 7.0 * 6.0 * std::pow( independentVariable, 5.0 );
    }

    double evaluateThirdDerivative( double independentVariable )
    {
        return 7.0 * 6.0 * 5.0 * std::pow( independentVariable, 4.0 );
    }

protected:

private:

};

std::shared_ptr< BaseFunctionInvPolyShaping > createBaseFunctionInvPolyShaping(
        const baseFunctionInvPolyShapingType baseFunctionType );


} // namespace shape_based_methods
} // namespace tudat

#endif // BASEFUNCTIONSINVPOLY_H
