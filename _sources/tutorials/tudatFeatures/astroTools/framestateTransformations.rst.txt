.. _tudatFeaturesFrameStateTransformations:

Frame/State Transformations
===========================

State type conversions
~~~~~~~~~~~~~~~~~~~~~~
Depending on your application, you will be using any of a number of state (position and velocity) representations. In Tudat, conversions involving the following state representations are available:

- Cartesian elements.
- Keplerian elements.
- Spherical-orbital elements.
- Modified Equinoctial elements.
- Unified State Model elements.

For each of these element types, conversions to/from Cartesian elements are available. Converting between two element types, where neither is Cartesian, will typically involve first transforming to Cartesian elements, and then transforming to your output state type. For each state type, the physical meaning of each of the elements is defined in the :literal:`statevectorIndices.h` file. In this file, you see for instance:

.. code-block:: cpp

    //! Keplerian elements indices.
    enum KeplerianElementIndices
    {
    semiMajorAxisIndex = 0,
    eccentricityIndex = 1,
    inclinationIndex = 2,
    argumentOfPeriapsisIndex = 3,
    longitudeOfAscendingNodeIndex = 4,
    trueAnomalyIndex = 5,
    semiLatusRectumIndex = 0
    };

This indicates that, for instance the eccentricity is index 1 and the true anomaly is index 5. As a result, you can use the following to retrieve the eccentricity from a vector if Kepler elements:

.. code-block:: cpp

    //! Keplerian elements indices.
    Eigen::Vector6d keplerElements = ....
    double currentEccentricity = keplerElements( orbital_element_conversions::eccentricityIndex );

or alternatively:

.. code-block:: cpp

    //! Keplerian elements indices.
    Eigen::Vector6d keplerElements = ....
    double currentEccentricity = keplerElements( 1 );

which yields the exact same result.

In the definition of the ``KeplerianElementIndices`` enum, you can see something peculiar: both ``semiMajorAxisIndex`` and ``semiLatusRectumIndex`` are defined as index 0. The latter option is only applicable when the orbit is parabolic (when the eccentricity is 1.0). That is, if the orbit is parabolic, element 0 does not represent the semi-major axis (as it is not defined) but the semi-latus rectum. Below, we list the details of the implementation of each of these state types in Tudat:

Kepler elements
***************
The Kepler elements are the standard orbital elements used in classical celestial mechanics, with the element indices shown above.
Converting to/from Cartesian state requires an additional piece of information in addition to the state itself: the gravitational parameter of the body w.r.t. the Keplerian elements are defined. Conversion to/from Cartesian elements is done as

.. code-block:: cpp

    Eigen::Vector6d cartesianState = ....
    double centralBodyGravitationalParameter = ...
    Eigen::Vector6d keplerianState = convertCartesianToKeplerianElements( cartesianState, centralBodyGravitationalParameter );

Similarly, the inverse operation is done as:

.. code-block:: cpp

    Eigen::Vector6d keplerianState = ....
    double centralBodyGravitationalParameter = ...
    Eigen::Vector6d cartesianState = convertKeplerianToCartesianElements( keplerianState, centralBodyGravitationalParameter );

In the definition of the state elements, you will notice that element 5 is the true anomaly, not the eccentric or mean anomaly. Tudat also contains functions to convert to these alternative anomalies. Converting between true and eccentric anomaly is done as follows:

.. code-block:: cpp

    double trueAnomaly = ....
    double eccentricity = ...
    double eccentricAnomaly = convertTrueAnomalyToEccentricAnomaly( trueAnomaly, eccentricity );

or directly from the orbital elements:

.. code-block:: cpp

    Eigen::Vector6d keplerianState = ...
    double eccentricAnomaly = convertTrueAnomalyToEccentricAnomaly( keplerianState( trueAnomalyIndex ), keplerianState( eccentricityIndex ) );

Note that this function automatically identifies whether the orbit is elliptical or hyperbolic, and computes the associated eccentric anomaly. The function for the inverse operation is ``convertEccentricAnomalyToTrueAnomaly``. Similarly, Tudat contains functions to convert from eccentric to mean anomaly (automatically checking whether the orbit is elliptical or hyperbolic):

.. code-block:: cpp

    double trueAnomaly = ....
    double eccentricity = ...
    double eccentricAnomaly = convertTrueAnomalyToEccentricAnomaly( trueAnomaly, eccentricity );
    double meanAnomaly = convertEccentricAnomalyToMeanAnomaly( eccentricAnomaly, eccentricity );

The inverse operation, mean to eccentric anomaly, is done separately for hyperbolic and elliptical orbits, through the functions ``convertMeanAnomalyToEccentricAnomaly`` for elliptical and ``convertMeanAnomalyToHyperbolicEccentricAnomaly`` for hyperbolic orbits. In general, you will use them as follows:

.. code-block:: cpp

    double meanAnomaly = ....
    double eccentricity = ...
    double eccentricAnomaly = convertMeanAnomalyToEccentricAnomaly( eccentricity, meanAnomaly );

However, this conversion involves the solution of an implicit algebraic equation, for which a root finder is used. Root finders are discussed in more detail here. When calling the function as in the above example, a :class:`RootFinder` is created internally. However, in some cases you may want to specify your own root finder, as well as a first initial guess for the eccentric anomaly (which the root finder uses at its first iteration). When doing so, you create a pointer to a root finder object and pass it to the conversion function as follows:

.. code-block:: cpp

    double meanAnomaly = ....
    double eccentricity = ...
    double initialGuess = ...
    boost::shared_ptr< RootFinder > rootFinder = ...
    double eccentricAnomaly = convertMeanAnomalyToEccentricAnomaly( eccentricity, meanAnomaly, false, initialGuess, rootFinder );

where the argument ``false`` indicates that the user-specified initial guess is to be used. If you want to use a custom-defined root finder, but not an initial guess, use the following:

.. code-block:: cpp

    double meanAnomaly = ....
    double eccentricity = ...
    boost::shared_ptr< RootFinder > rootFinder = ...
    double eccentricAnomaly = convertMeanAnomalyToEccentricAnomaly( eccentricity, meanAnomaly, true, TUDAT_NAN, rootFinder );

Spherical-orbital elements
**************************
The spherical elements are typically used to denote the conditions in atmospheric flight. In most applications, they will be used to denote the state in a body-fixed frame. The details of the physical meaning of the elements is discussed here. The element indices in Tudat are the following:

.. code-block:: cpp

    //! Spherical orbital state element indices
    enum SphericalOrbitalStateElementIndices
    {
        radiusIndex = 0,
        latitudeIndex = 1,
        longitudeIndex = 2,
        speedIndex = 3,
        flightPathIndex = 4,
        headingAngleIndex = 5
    };

The spherical elements consist of 6 entries, with no additional information required for the conversion to/from Cartesian elements. The conversion from Cartesian to spherical elements is performed as:

.. code-block:: cpp

    Eigen::Vector6d cartesianState = ....
    Eigen::Vector6d sphericalState = convertCartesianToSphericalOrbitalState( cartesianState );

Similarly, the inverse operation is done as:

.. code-block:: cpp

    Eigen::Vector6d sphericalState = ....
    Eigen::Vector6d cartesianState = convertSphericalOrbitalToCartesianState( sphericalState );

Modified Equinoctial elements
*****************************
The modified equinoctial elements are typically used for orbits with eccentricities near 0 or 1 and/or inclinations near 0 or :math:`\pi`. The element indices in Tudat are the following:

.. code-block:: cpp

   //! Modified equinoctial element vector indices.
   enum ModifiedEquinoctialElementVectorIndices
   {
       semiParameterIndex = 0,
       fElementIndex = 1,
       gElementIndex = 2,
       hElementIndex = 3,
       kElementIndex = 4,
       trueLongitudeIndex = 5 
   };

The modified equinoctial elements consists of 6 elements. The conversion to/from Cartesian elements requires the gravitation parameter of the body w.r.t. which the Modified Equinoctial elements are defined. Furthermore, a :literal:`bool` is used to indicate whether the singularity of this element set occurs for inclinations of 0 or :math:`\pi`. The conversion from Cartesian elements is done as:

.. code-block:: cpp

    Eigen::Vector6d cartesianState = ....
    double centralBodyGravitationalParameter = ...
    Eigen::Vector6d cartesianState = convertCartesianToModifiedEquinoctial( cartesianState, centralBodyGravitationalParameter, flipSingularityToZeroInclination );

.. Note:: :literal:`flipSingularityToZeroInlination` is optional for this conversion. If left empty, an overloaded function will determine whether this value is true or false based on the inclination of the orbit. 

Similarly, the inverse operation is done as:

.. code-block:: cpp

    Eigen::Vector6d modifiedEquinoctialElements = ....
    double centralBodyGravitationalParameter = ...
    Eigen::Vector6d cartesianState = convertModifiedEquinoctialToCartesian( modifiedEquinoctialElements, centralBodyGravitationalParameter, flipSingularityToZeroInclination );


 
Unified State Model elements
****************************
The element indices for the Unified State Model elements in Tudat are the following:

.. code-block:: cpp

   //! Unified State Model indices.
   enum UnifiedStateModelElementIndices
   {
       CHodographIndex = 0,
       Rf1HodographIndex = 1,
       Rf2HodographIndex = 2,
       epsilon1QuaternionIndex = 3,
       epsilon2QuaternionIndex = 4,
       epsilon3QuaternionIndex = 5,
       etaQuaternionIndex = 6
   };

The unified state model elements consists of 7 elements. Only the conversion to/from Keplerian elements is implemented. The conversion requires the gravitation parameter of the body w.r.t. which the Unified State Model elements are defined. The conversion from Keplerian elements is done as:

.. code-block:: cpp

    Eigen::Vector6d keplerianElements = ....
    double centralBodyGravitationalParameter = ...
    Eigen::Vector6d cartesianState = convertKeplerianToUnifiedStateModelElements( keplerianElements, centralBodyGravitationalParameter );

Similarly, the inverse operation is done as:

.. code-block:: cpp

    Eigen::Vector6d unifiedStateModelElements = ....
    double centralBodyGravitationalParameter = ...
    Eigen::Vector6d cartesianState = convertUnifiedStateElementsToKeplerian( unifiedStateModelElements, centralBodyGravitationalParameter );


Frame transformations
~~~~~~~~~~~~~~~~~~~~~
Every state, regardless of its representation is expressed with a particular origin and orientation. This is most easy to understand for Cartesian elements, where the origin represents the (0,0,0) position, and the orientation defines the direction of the x-, y- and z-axes. Below, we discuss how to perform these operations in Tudat.

.. warning:: Do not use the getCurrentState or getCurrentRotation functions in the body objects! These functions are used during numerical propagation, and calling them outside of the numerical propagation will generally not lead to meaningful results.

Frame translations
******************
To change the origin of a Cartesian, one can simply add a Cartesian state that represents the difference between the original and the new origin. For instance, when transforming a vector (state of a vehicle) from Earth-centered to Moon-centered (keeping the orientation constant):

.. code-block:: cpp

    Eigen::Vector6d vehicleCartesianStateInEarthCenteredFrame = ....
    Eigen::Vector6d moonCartesianStateInEarthCenteredFrame = ....
    Eigen::Vector6d vehicleCartesianStateInMoonCenteredFrame = vehicleCartesianStateInEarthCenteredFrame + moonCartesianStateInEarthCenteredFrame;

The challenge here, of course, is determining the ``moonCartesianStateInEarthCenteredFrame`` vector. We provide a few ways in which to achieve this. When performing a numerical simulation using a set of body objects, you can use the following (asuming the the ``bodyMap`` contains both an Earth and Moon entry):

.. code-block:: cpp

    NamedBodyMap bodyMap = ...
    double currentTime = ...
    Eigen::Vector6d moonCartesianStateInEarthCenteredFrame = bodyMap.at( "Moon" )->getStateInBaseFrameFromEphemeris( currentTime ) - 
          bodyMap.at( "Earth" )->getStateInBaseFrameFromEphemeris( currentTime );

You can also bypass the body map altogether, and use Spice to obtain the relative state. Note, however, that this will use whichever ``spice`` kernels you have loaded, and may not be consistent with the states you are using the bodies in your simulation.

.. code-block:: cpp

    double currentTime = ...
    std::string frameOrientation = "J2000";
    Eigen::Vector6d moonCartesianStateInEarthCenteredFrame = 
        spice_interface::getBodyCartesianStateAtEpoch( "Moon", "Earth", frameOrientation, "NONE", currentTime

where the ``NONE`` arguments indicates that no light-time corrections are used, and the frame orientation denotes the orientation of the frame in which the relative state is returned.

Frame rotations
***************
Rotating the frame in which a Cartesian state is expressed requires two pieces of information:

    1. The rotation matrix from one frame to the other
    2. The first time derivative of the rotation matrix from one frame to the other

Manually, the state may then be transformed as:

.. code-block:: cpp

    Eigen::Matrix3d rotationToFrame = ...
    Eigen::Matrix3d timeDerivativeOfRotationToFrame  = ...
    Eigen::Vector6d originalState = ...
    Eigen::Vector6d rotatedState;
    rotatedState.segment( 0, 3 ) = rotationToFrame * originalState( 0, 3 );
    rotatedState.segment( 3, 3 ) = rotationToFrame * originalState( 3, 3 ) + timeDerivativeOfRotationToFrame * originalState( 0, 3 );

In many cases, however, your frame rotation will be from the inertial frame to a body-fixed frame. All information required for this is stored in ``RotationalEphemeris`` objects. This object contains a base (inertial) and target (body-fixed) frame and defines the rotation between the two. Assuming that you are using a body map to store your environment, you can transform the state from an inertial to a body-fixed frame as follows, for the example of transforming a vehicle's Cartesian state from an inertial to the body-fixed frame of the Earth:

.. code-block:: cpp

    NamedBodyMap bodyMap = ...
    double currentTime = ...
    Eigen::Vector6d inertialState = ...
    Eigen::Vector6d bodyFixedState = transformStateToTargetFrame( inertialState, currentTime, bodyMap.at( "Earth" )->getRotationalEphemeris( ) );
    
the inverse rotation is done as follows:

.. code-block:: cpp

    NamedBodyMap bodyMap = ...
    double currentTime = ...
    Eigen::Vector6d bodyFixedState = ...
    Eigen::Vector6d inertialState = transformStateToGlobalFrame( bodyFixedState, currentTime, bodyMap.at( "Earth" )->getRotationalEphemeris( ) );


