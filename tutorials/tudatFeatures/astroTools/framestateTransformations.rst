.. _tudatFeaturesFrameStateTransformations:

Frame/State Transformations
===========================

State Type Conversions
~~~~~~~~~~~~~~~~~~~~~~
Depending on your application, you will be using any of a number of translational state (position and velocity) representations. In Tudat, conversions involving the following state representations are available:

- Cartesian elements.
- Keplerian elements.
- Spherical-orbital elements.
- Modified Equinoctial elements.
- Unified State Model elements.

For each of these element types, conversions to/from Cartesian elements are available. Converting between two element types, where neither is Cartesian, will typically involve first transforming to Cartesian elements, and then transforming to your output state type. 

In case you are also working with rotational motion, in Tudat the following representations for attitude are available:

- Quaternions.
- Modified Rodrigues parameters.
- Exponential map. 

Transformation between these elements is done by passing through quaternions first. In fact, this is the default attitude representation in Tudat. 

For each state type, the physical meaning of each of the elements is defined in the :literal:`statevectorIndices.h` file. In this file, you will see for instance:

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

In the definition of the ``KeplerianElementIndices`` enumeration, you can see something peculiar: both ``semiMajorAxisIndex`` and ``semiLatusRectumIndex`` are defined as index 0. The latter option is only applicable when the orbit is parabolic (when the eccentricity is 1.0). That is, if the orbit is parabolic, element 0 does not represent the semi-major axis (as it is not defined) but the semi-latus rectum. Below, we list the details of the implementation of each of these state types in Tudat:

Kepler Elements
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

Spherical-orbital Elements
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

Modified Equinoctial Elements
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

.. Note:: The input :literal:`flipSingularityToZeroInlination` is optional for this conversion. If left empty, an overloaded function will determine whether this value is true or false based on the inclination of the orbit. 

Similarly, the inverse operation is done as:

.. code-block:: cpp

    Eigen::Vector6d modifiedEquinoctialElements = ....
    double centralBodyGravitationalParameter = ...
    Eigen::Vector6d cartesianState = convertModifiedEquinoctialToCartesian( modifiedEquinoctialElements, centralBodyGravitationalParameter, flipSingularityToZeroInclination );

Unified State Model Elements
****************************
Three different versions of the Unified State Model are present in Tudat. They differ based on the coordinates chosen to represent the rotation from local orbital to inertial frame, which can be expressed in quaternions, modified Rodrigues parameters or exponantial map. The element indices for the Unified State Model elements with quaternions (or USM7) are the following:

.. code-block:: cpp

   //! Unified state model with quaternions indices.
   enum UnifiedStateModelQuaternionsElementIndices
   {
       CHodographUSM7Index = 0,
       Rf1HodographUSM7Index = 1,
       Rf2HodographUSM7Index = 2,
       etaUSM7Index = 3,
       epsilon1USM7Index = 4,
       epsilon2USM7Index = 5,
       epsilon3USM7Index = 6
   };

For the Unified State Model elements with modified Rodrigues parameters (or USM6) the indeces are:

.. code-block:: cpp

   //! Unified state model with modified Rodrigues parameters indices.
   enum UnifiedStateModelModifiedRodriguesParametersElementIndices
   {
       CHodographUSM6Index = 0,
       Rf1HodographUSM6Index = 1,
       Rf2HodographUSM6Index = 2,
       sigma1USM6Index = 3,
       sigma2USM6Index = 4,
       sigma3USM6Index = 5,
       shadowFlagUSM6Index = 6
   };

And finally, for the Unified State Model elements with exponential map (or USMEM) they are:

.. code-block:: cpp
   
   //! Unified state model with exponential map indices.
   enum UnifiedStateModelExponentialMapElementIndices
   {
       CHodographUSMEMIndex = 0,
       Rf1HodographUSMEMIndex = 1,
       Rf2HodographUSMEMIndex = 2,
       e1USMEMIndex = 3,
       e2USMEMIndex = 4,
       e3USMEMIndex = 5,
       shadowFlagUSMEMIndex = 6
   };

Regardless of the rotational coordinates chosen, the Unified State Model elements consists of 7 elements. For each Unified State Model representation, conversion to and from Keplerian and Cartesian coordinates is implemented. As an example, the conversion from Keplerian elements for the USM7 elements is shown here:

.. code-block:: cpp

    Eigen::Vector6d keplerianElements = ....
    double centralBodyGravitationalParameter = ...
    Eigen::Vector6d cartesianState = convertKeplerianToUnifiedStateModelElements( keplerianElements, centralBodyGravitationalParameter );

Similarly, the inverse operation is done as:

.. code-block:: cpp

    Eigen::Vector6d unifiedStateModelElements = ....
    double centralBodyGravitationalParameter = ...
    Eigen::Vector6d cartesianState = convertUnifiedStateElementsToKeplerian( unifiedStateModelElements, centralBodyGravitationalParameter );

Quaternions
***********
As mentioned at the beginning of this chapter, quaternions are the default attitude representation in Tudat. Depending on the location in the Tudat framework, you will find a quaternion element expressed as either of the two types below:

.. method:: As an Eigen::Quaterniond

   This method is used mainly in the :class:`Body` class, to express various rotations, such as the rotation to base frame. The advantage of this class, is that it comes with some very useful member functions. For instance, if you have a :literal:`Quaterniond` object, you can directly transform it to a direction cosine matrix (or transformation matrix) by using the method :literal:`.toRotationMatrix( )`. You can find more details on the definition and use of the :literal:`Quaterniond` in the Eigen website. 

   .. warning:: 
      The definition of rotations in Eigen is slightly different than in Tudat. This means that when using Eigen functions such as :literal:`.toRotationMatrix( )` to quaternions defined within Tudat, the result will be the inverse of what is desired. Thus, the tranformation to matrix should always be followed by :literal:`.transpose( )`, to give the correct rotation matrix.

.. method:: As an Eigen::Vector4d

   This method is mainly used, on the other hand, in propagation. By defining the quaternion as a simple four-dimensional vector, its element can be easily extracted and replaced from the rotational state vector (which also includes rotational velocity).

Transformation between the two methods is defined in the linear algebra module of Tudat (see :literal:`linearAlgebra.h`). To transform a quaternion to vector format, one can use:

   .. code-block:: cpp

      Eigen::Quaterniond quaternionAsQuaternion = ...
      Eigen::Vector4d quaternionAsVector = linear_algebra::convertQuaternionToVectorFormat( quaternionAsQuaternion );

and vice-versa as:

   .. code-block:: cpp

      Eigen::Vector4d quaternionAsVector = ...
      Eigen::Quaterniond quaternionAsQuaternion = linear_algebra::convertVectorToQuaternionFormat( quaternionAsVector );

Finally, the quaternion indices are defined as follows:

.. code-block:: cpp
   
   //! Quaternions indices.
   enum QuaternionsElementIndices
   {
       etaQuaternionIndex = 0,
       epsilon1QuaternionIndex = 1,
       epsilon2QuaternionIndex = 2,
       epsilon3QuaternionIndex = 3
   };

Before continuing to the next rotational coordinates, there is one last fact about quaternions that needs to be discussed. This relates to the property that the sum of the square of the quaternion vector equals one, i.e., :math:`\sum_{i=0}^3 q_i^2 = 1`. In fact, during numerical propagation it is possible that this property is violated, due to various sources of differs from one by more than about :math:`5 \times 10^{-15}`. Note that this also applies to the quaternion vector present in the Unified State Model with quaternions, or USM7, state.

Modified Rodrigues Parameters
*****************************
One of the other two supported attitude representations is the modified Rodrigues parameters (MRPs). The indeces for MRPs are defined as follows:

.. code-block:: cpp
   
   //! Modified Rodrigues parameters indices.
   enum ModifiedRodriguesParametersElementIndices
   {
       sigma1ModifiedRodriguesParametersIndex = 0,
       sigma2ModifiedRodriguesParametersIndex = 1,
       sigma3ModifiedRodriguesParametersIndex = 2,
       shadowFlagModifiedRodriguesParametersIndex = 3
   };

Transformation to and from quaternions is achieved with the functions :literal:`convertModifiedRodriguesParametersToQuaternionElements` and :literal:`convertQuaternionsToModifiedRodriguesParameterElements`, respectively, where the only input is the attitude element (in vector format). 

The last element shown in the :literal:`ModifiedRodriguesParametersElementIndices` enumeration is the flag that triggers the shadow modifed Rodrigues parameters (SMRPs). Its use is introduced to avoid the singularity at :math:`\pm 2\pi` radians. If its value is 0, then the elements are MRPs, whereas if it is 1, then they are SMRPs. The use of SMRPs results in slightly different equations of motion and transformations. The switch between MRPs and SMRPs occurs whenever the magnitude of the rotation represented by the MRP vector is larger than :math:`\pi`.

Exponential Map
***************
The final attitude representations is the exponential map (EM). The indeces for EM are defined as follows:

.. code-block:: cpp
   
   //! Exponential map indices.
   enum ExponentialMapElementIndices
   {
       e1ExponentialMapIndex = 0,
       e2ExponentialMapIndex = 1,
       e3ExponentialMapIndex = 2,
       shadowFlagExponentialMapIndex = 3
   };

and transformation to and from quaternions is achieved with the aid of the functions :literal:`convertExponentialMapToQuaternionElements` and :literal:`convertQuaternionsToExponentialMapElements`, respectively. Also for these equations the only input is the attitude element (in vector format).

Similarly to MRPs, the exponential map elements also make use of the shadow flag. In this case, this flag signals whether the shadow exponential map (SEM) is in use. This flag is also introduces to avoid the singularity at :math:`\pm 2\pi` radians, but interestingly, there is no difference between the equations of motion and transformations in terms of EM or SEM. In fact, they are only introduced to make sure that when converting from EM to quaternions, the resulting quaternion sign history is continuous. The switch between EM and SEM occurs whenever the magnitude of the rotation represented by the EM vector is larger than :math:`\pi`.

Frame Transformations
~~~~~~~~~~~~~~~~~~~~~
Every state, regardless of its representation is expressed with a particular origin and orientation. This is most easy to understand for Cartesian elements, where the origin represents the (0,0,0) position, and the orientation defines the direction of the x-, y- and z-axes. Below, we discuss how to perform these operations in Tudat.

.. warning:: Do not use the :literal:`getCurrentState` or :literal:`getCurrentRotation` functions in the body objects! These functions are used during numerical propagation, and calling them outside of the numerical propagation will generally not lead to meaningful results.

Frame Translations
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

Frame Rotations
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


