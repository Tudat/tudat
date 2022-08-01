/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/basic_astro/stateRepresentationConversions.h"
#include "tudat/astro/basic_astro/sphericalBodyShapeModel.h"
#include "tudat/astro/basic_astro/oblateSpheroidBodyShapeModel.h"
#include "tudat/math/basic/coordinateConversions.h"

namespace tudat
{

namespace coordinate_conversions
{

//! Function to convert a position from one representation to another
Eigen::Vector3d convertPositionElements(
        const Eigen::Vector3d& originalElements,
        const PositionElementTypes originalElementType,
        const PositionElementTypes convertedElementType,
        const std::shared_ptr< basic_astrodynamics::BodyShapeModel > shapeModel,
        const double tolerance )
{
    Eigen::Vector3d convertedElements = Eigen::Vector3d::Constant( TUDAT_NAN );
    if( originalElementType == convertedElementType )
    {
        convertedElements = originalElements;
    }
    else
    {
        // Retrieve flattening and equatorial radius if required.
        double flattening = TUDAT_NAN, equatorialRadius = TUDAT_NAN;
        if( ( originalElementType == geodetic_position ) || ( convertedElementType == geodetic_position ) )
        {
            if( shapeModel == nullptr )
            {
                throw std::runtime_error( "Error when converting to/from geodetic position, no shape model is provided" );
            }
            else if( std::dynamic_pointer_cast< basic_astrodynamics::SphericalBodyShapeModel >( shapeModel ) != nullptr )
            {
                flattening = 0.0;
                equatorialRadius = std::dynamic_pointer_cast< basic_astrodynamics::SphericalBodyShapeModel >(
                            shapeModel )->getAverageRadius( );
            }
            else if( std::dynamic_pointer_cast< basic_astrodynamics::OblateSpheroidBodyShapeModel >( shapeModel ) != nullptr )
            {
                std::shared_ptr< basic_astrodynamics::OblateSpheroidBodyShapeModel > oblateSphereModel =
                        std::dynamic_pointer_cast< basic_astrodynamics::OblateSpheroidBodyShapeModel >( shapeModel );
                flattening = oblateSphereModel->getFlattening( );
                equatorialRadius = oblateSphereModel->getEquatorialRadius( );
            }
            else
            {
                throw std::runtime_error( "Error when converting to/from geodetic position, shape model not recognized" );
            }
        }
        // Check original type
        switch( originalElementType )
        {
        case cartesian_position:
        {
            // Convert Cartesian position to requested type.
            switch( convertedElementType )
            {
            case spherical_position:
            {
                convertedElements = convertCartesianToSpherical(
                            originalElements );
                convertedElements( 1 ) = mathematical_constants::PI / 2.0 - convertedElements( 1 );
                break;
            }
            case geodetic_position:
            {
                convertedElements = convertCartesianToGeodeticCoordinates(
                            originalElements, equatorialRadius, flattening, tolerance );
                break;
            }
            default:
                throw std::runtime_error( "Error when converting from cartesian position, target element type not recognized" );
            }
            break;
        }
        case spherical_position:
        {
            // Convert spherical position to requested type.
            switch( convertedElementType )
            {
            case cartesian_position:
            {
                Eigen::Vector3d inputElements = originalElements;
                inputElements( 1 ) = mathematical_constants::PI / 2.0 - inputElements( 1 );

                convertedElements = convertSphericalToCartesian(
                             inputElements );
                break;
            }
            case geodetic_position:
            {
                Eigen::Vector3d inputElements = originalElements;
                inputElements( 1 ) = mathematical_constants::PI / 2.0 - inputElements( 1 );

                 Eigen::Vector3d intermediateElements = convertSphericalToCartesian(
                             inputElements );
                 convertedElements = convertCartesianToGeodeticCoordinates(
                             intermediateElements, equatorialRadius, flattening, tolerance );

                 break;
            }
            default:
                throw std::runtime_error( "Error when converting from spherical position, target element type not recognized" );

            }
            break;
        }
        case geodetic_position:
        {
            // Convert geodetic position to requested type.
            switch( convertedElementType )
            {
            case cartesian_position:
            {
                convertedElements = convertGeodeticToCartesianCoordinates(
                            originalElements, equatorialRadius, flattening );
                break;
            }
            case spherical_position:
            {
                Eigen::Vector3d intermediateElements = convertGeodeticToCartesianCoordinates(
                            originalElements, equatorialRadius, flattening );
                convertedElements = convertCartesianToSpherical(
                            intermediateElements );
                convertedElements( 1 ) = mathematical_constants::PI / 2.0 - convertedElements( 1 );

                break;
            }
            default:
                throw std::runtime_error( "Error when converting from geodetic position, target element type not recognized" );
            }
            break;
        }
        default:
            throw std::runtime_error( "Error when converting position elements, base element type not recognized" );
        }
    }
    return convertedElements;

}

}

}

