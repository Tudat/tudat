The Environment During Propagation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Each :class:`Body` object and its constituent members is updated to the current state and time automatically during the numerical propagation. We stress that only those models that are relevant for a given propagation are updated every time step (this is handled automatically, without user intervention). Some time-dependent properties of the body are set in the environment models themselves. Others are updated and stored directly in the :class:`Body` object. Below is a full list of (possibly) time varying environment models, and how to retrieve them from a Body object during propagation.


**The current translational state**

   Retrieved direct from a :class:`Body` object directly, with the :literal:`getState` function, as Cartesian elements. Note that this state is *always* in the global frame origin and orientation (note that separate :literal:`getPosition` and :literal:`getVelocity` functions are also available.

**The current rotational state**

   Retrieved direct from a :class:`Body` object directly. The current orientation is defined by the body's rotation model, and is retrieved as a quaternion. To get the quaternion to transform a vector from inertial (e.g. with global frame orientation) to the body-fixed frame, the :literal:`getCurrentRotationToLocalFrame` can be used. The inverse rotation can be obtained from :literal:`getCurrentRotationToGlobalFrame`. The time-derivative of the orientation is provided in two formulations (with equivalent information content): the angular velocity vector, and the time derivative of the rotation matrix. The angular velocity vector, in inertial and body-fixed coordinates, is obtained from the :literal:`getCurrentAngularVelocityVectorInGlobalFrame` and :literal:`getCurrentAngularVelocityVectorInLocalFrame` functions respectively. The time-derivative of the rotation matrix from inertial to body-fixed frame is given by :literal:`getCurrentRotationMatrixDerivativeToLocalFrame`, while the derivative of the inverse rotation is taken from :literal:`getCurrentRotationMatrixDerivativeToGlobalFrame`

**The current body mass**

   Retrieved direct from a :class:`Body` object directly, with the :literal:`getBodyMass` function.

**Spherical harmonic gravity field coefficients**

   These coefficients may be time variable (set by using one or more :class:`GravityFieldVariationSettings`). The cosine and sine coefficients can be retrieved from a Body object through its gravity field model. A piece of example code on retrieving these coefficients is given below for the case of Earth:

.. code-block:: cpp

    NamedBodyMap bodyMap = ....

    std::shared_ptr< SphericalHarmonicsGravityField > sphericalHarmonicsGravityField = 
        std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( bodyMap.at( "Earth" ) );
    
    if( sphericalHarmonicsGravityField == nullptr )
    {
        throw std::runtime_error( "Error when retrieving spherical harmonic coefficients: Earth does not have a SphericalHarmonicsGravityField" );
    }
    
    Eigen::MatrixXd cosineCoefficients = sphericalHarmonicsGravityField->getCosineCoefficients( );
    Eigen::MatrixXd sineCoefficients = sphericalHarmonicsGravityField->getSineCoefficients( );

Note the check to see if the :literal:`dynamic_pointer_cast` was succesfull. These checks are very helpful in ascertaining the reason for a crashing program.

**Curent flight conditions**

  The :class:`FlightConditions` class, and its derived class :class:`AtmosphericFlightConditions` stores data relating to altitude, flight angles, atmospheric properties, etc. Follow the links for their detailed description. The :class:`FlightConditions` base class pointer is obtained from the :literal:`getFlightConditions` function of a :class:`Body` object.


  

.. note:: As a user, you will typically not access these variables directly. Important examples of cases where users *can* explicitly access them, are custom aerodynamic or thrust guidance.
