/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      The Unified State Model: Derivation and Applications in astro and Navigation,
 *          Vivek Vittaldevs; M.Sc. Thesis TU Delft (2010). Available on repository.tudelft.nl
 *
 */

#ifndef TUDAT_UNIFIED_STATE_MODEL_QUATERNION_ELEMENT_CONVERSIONS_H
#define TUDAT_UNIFIED_STATE_MODEL_QUATERNION_ELEMENT_CONVERSIONS_H

#include "tudat/basics/basicTypedefs.h"

namespace tudat
{

namespace orbital_element_conversions
{

//! Convert Keplerian to unified state model elements with quaternions.
/*!
 * Converts Keplerian to unified state model elements with quaternions.
 * \param keplerianElements Vector containing Keplerian elements. Order of elements is important!
 *         keplerianElements( 0 ) = semi-major axis,                                            [m]
 *         keplerianElements( 1 ) = eccentricity,                                               [-]
 *         keplerianElements( 2 ) = inclination (in range [0,PI]),                            [rad]
 *         keplerianElements( 3 ) = argument of periapsis,                                    [rad]
 *         keplerianElements( 4 ) = longitude of ascending node,                              [rad]
 *         keplerianElements( 5 ) = true anomaly.                                             [rad]
 * \param centralBodyGravitationalParameter Gravitational parameter of central body.      [m^3/s^2]
 * \return convertedUnifiedStateModelElements Converted state in unified state model elements with quaternions.
 *         The order of elements is fixed
 *         convertedUnifiedStateModelElements( 0 ) = C hodograph element,                     [m/s]
 *         convertedUnifiedStateModelElements( 1 ) = Rf1 hodograph element,                   [m/s]
 *         convertedUnifiedStateModelElements( 2 ) = Rf1 hodograph element,                   [m/s]
 *         convertedUnifiedStateModelElements( 3 ) = epsilon1 quaternion element,               [-]
 *         convertedUnifiedStateModelElements( 4 ) = epsilon2 quaternion element,               [-]
 *         convertedUnifiedStateModelElements( 5 ) = epsilon3 quaternion element,               [-]
 *         convertedUnifiedStateModelElements( 6 ) = eta quaternion element.                    [-]
 */
Eigen::Matrix< double, 7, 1 > convertKeplerianToUnifiedStateModelQuaternionsElements(
        const Eigen::Matrix< double, 6, 1 >& keplerianElements,
        const double centralBodyGravitationalParameter );

//! Convert unified state model elements with quaternions to Keplerian elements.
/*!
 * Converts unified state model elements with quaternions to Keplerian elements.
 * \param unifiedStateModelElements Vector containing unified state model elements with quaternions.
 *         Order of elements is important!
 *         unifiedStateModelElements( 0 ) = C hodograph element,                              [m/s]
 *         unifiedStateModelElements( 1 ) = Rf1 hodograph element,                            [m/s]
 *         unifiedStateModelElements( 2 ) = Rf1 hodograph element,                            [m/s]
 *         unifiedStateModelElements( 3 ) = epsilon1 quaternion element,                        [-]
 *         unifiedStateModelElements( 4 ) = epsilon2 quaternion element,                        [-]
 *         unifiedStateModelElements( 5 ) = epsilon3 quaternion element,                        [-]
 *         unifiedStateModelElements( 6 ) = eta quaternion element.                             [-]
 * \param centralBodyGravitationalParameter Gravitational parameter of central body.      [m^3/s^2]
 * \param forceQuaternionNormalization Boolean denoting whether the normalization of the quaternion
 * vector has to be forced.
 * \return convertedKeplerianElements Converted state in Keplerian elements. The order of elements is fixed!
 *         convertedKeplerianElements( 0 ) = semi-major axis,                                   [m]
 *         convertedKeplerianElements( 1 ) = eccentricity,                                      [-]
 *         convertedKeplerianElements( 2 ) = inclination (in range [0,PI]),                   [rad]
 *         convertedKeplerianElements( 3 ) = argument of periapsis,                           [rad]
 *         convertedKeplerianElements( 4 ) = longitude of ascending node,                     [rad]
 *         convertedKeplerianElements( 5 ) = true anomaly.                                    [rad]
 */
Eigen::Matrix< double, 6, 1 > convertUnifiedStateModelQuaternionsToKeplerianElements(
        const Eigen::Matrix< double, 7, 1 >& unifiedStateModelElements,
        const double centralBodyGravitationalParameter,
        const bool forceQuaternionNormalization = false );

//! Convert Cartesian elements to unified state model elements with quaternions.
/*!
 * Converts Cartesian to unified state model elements with quaternions.
 * \param cartesianElements Converted state in Cartesian elements. The order of elements is fixed!
 *         cartesianElements( 0 ) = x-position coordinate,                                      [m]
 *         cartesianElements( 1 ) = y-position coordinate,                                      [m]
 *         cartesianElements( 2 ) = z-position coordinate,                                      [m]
 *         cartesianElements( 3 ) = x-velocity coordinate,                                    [m/s]
 *         cartesianElements( 4 ) = y-velocity coordinate,                                    [m/s]
 *         cartesianElements( 5 ) = z-velocity coordinate.                                    [m/s]
 * \param centralBodyGravitationalParameter Gravitational parameter of central body.      [m^3/s^2]
 * \return convertedUnifiedStateModelElements Converted state in unified state model elements with quaternions.
 *         The order of elements is fixed
 *         convertedUnifiedStateModelElements( 0 ) = C hodograph element,                     [m/s]
 *         convertedUnifiedStateModelElements( 1 ) = Rf1 hodograph element,                   [m/s]
 *         convertedUnifiedStateModelElements( 2 ) = Rf1 hodograph element,                   [m/s]
 *         convertedUnifiedStateModelElements( 3 ) = epsilon1 quaternion element,               [-]
 *         convertedUnifiedStateModelElements( 4 ) = epsilon2 quaternion element,               [-]
 *         convertedUnifiedStateModelElements( 5 ) = epsilon3 quaternion element,               [-]
 *         convertedUnifiedStateModelElements( 6 ) = eta quaternion element.                    [-]
 */
Eigen::Matrix< double, 7, 1 > convertCartesianToUnifiedStateModelQuaternionsElements(
        const Eigen::Matrix< double, 6, 1 >& cartesianElements,
        const double centralBodyGravitationalParameter );

//! Convert unified state model elements with quaternions to Cartesian elements.
/*!
 * Converts unified state model elements with quaternions to Cartesian elements.
 * \param unifiedStateModelElements Vector containing unified state model elements with quaternions.
 *        Order of elements is important!
 *         unifiedStateModelElements( 0 ) = C hodograph element,                              [m/s]
 *         unifiedStateModelElements( 1 ) = Rf1 hodograph element,                            [m/s]
 *         unifiedStateModelElements( 2 ) = Rf1 hodograph element,                            [m/s]
 *         unifiedStateModelElements( 3 ) = epsilon1 quaternion element,                        [-]
 *         unifiedStateModelElements( 4 ) = epsilon2 quaternion element,                        [-]
 *         unifiedStateModelElements( 5 ) = epsilon3 quaternion element,                        [-]
 *         unifiedStateModelElements( 6 ) = eta quaternion element.                             [-]
 * \param centralBodyGravitationalParameter Gravitational parameter of central body.      [m^3/s^2]
 * \param forceQuaternionNormalization Boolean denoting whether the normalization of the quaternion
 * vector has to be forced.
 * \return convertedCartesianElements Converted state in Cartesian elements. The order of elements is fixed!
 *         convertedCartesianElements( 0 ) = x-position coordinate,                            [m]
 *         convertedCartesianElements( 1 ) = y-position coordinate,                            [m]
 *         convertedCartesianElements( 2 ) = z-position coordinate,                            [m]
 *         convertedCartesianElements( 3 ) = x-velocity coordinate,                          [m/s]
 *         convertedCartesianElements( 4 ) = y-velocity coordinate,                          [m/s]
 *         convertedCartesianElements( 5 ) = z-velocity coordinate.                          [m/s]
*/
Eigen::Matrix< double, 6, 1 > convertUnifiedStateModelQuaternionsToCartesianElements(
        const Eigen::Matrix< double, 7, 1 >& unifiedStateModelElements,
        const double centralBodyGravitationalParameter,
        const bool forceQuaternionNormalization = false );

} // namespace orbital_element_conversions

} // close tudat

#endif // TUDAT_UNIFIED_STATE_MODEL_QUATERNION_ELEMENT_CONVERSIONS_H
