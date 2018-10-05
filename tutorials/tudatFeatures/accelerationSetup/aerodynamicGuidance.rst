.. _tudatFeaturesAerodynamicGuidance:

Aerodynamic Guidance
====================
When including aerodynamics in the orbit propagation, it may often be desirable to modify the aerodynamic properties of the vehicle as a function of the current environment. This is the case for ascent and entry for instance. Also, control surface deflections may be needed to modify the vehicle's behaviour. Both options are included in Tudat, and the interfaces are described below.

.. note:: In applications were no aerodynamic guidance is used, the angle of attack, sideslip angle and bank angle are all implicitly set to zero.

Vehicle orientation
~~~~~~~~~~~~~~~~~~~
Although rotational motion can be propagated in Tudat (see :class:`RotationalStatePropagatorSettings`), sometime the user simply wants to impose this rotational motion on the vehicle. 

For aerodynamics, the vehicle orientation is typically described by the angle of attack, angle of sideslip and bank angle, which is also the approach we take in Tudat. There is a single interface in Tudat through which the functions describing these angles are fed to the trajectory propagation:


.. class:: AerodynamicAngleCalculator

   This class is used for the describing the angle of attack, angle of sideslip and bank angle of the vehicle. The :literal:`setOrientationAngleFunctions` member function of this class is used. However, there is a number of manners in which to use this function, either directly or indirectly, the preferred selection of course depending on your application. The options are:

      - Manual definitions of the three angles of a function of time (can also be used to set constant angles).
      - Pre-defined guidance laws (only angle-of-attack trim presently implemented).
      - A user-defined :class:`AerodynamicGuidance` class, which computes the current angles using any number of input/environment variables. This option is typically the preferred option for a realistic entry propagation.

   Below, we will indicate how to use these three options for the example of the Apollo entry example. Currently, the following code is used in the example, defining a constant angle of attack of 30 degree, and 0 sideslip and bank angle:

   .. code-block:: cpp

      bodyMap.at( "Apollo" )->getFlightConditions( )->getAerodynamicAngleCalculator( )->setOrientationAngleFunctions(
                  [ = ]( ){ return 30.0 * mathematical_constants::PI / 180.0; } );

   .. Warning:: We stress that any definition of the aerodynamic guidance must be done after the acceleration models have been created (by :class:`AccelerationSettings`)(but before the trajectory is propagated).

Manual angle functions
**********************
To manually set the body orientation as constant angles, the :literal:`setOrientationAngleFunctions` function is called directly, as is done in the above example. Note that if you need a definition of aerodynamic angles that varies with time, a specific :class:`AerodynamicGuidance` derived class should be implemented or one of the predefined orientation functions should be used. In general, the manual definition of the angles is done as follows:

.. code-block:: cpp
  
   double constantAngleOfAttack = 30.0 * mathematical_constants::PI / 180.0;
   double constantSideslipAngle = 2.0 * mathematical_constants::PI / 180.0;
   double constantBankAngle = -70.0 * mathematical_constants::PI / 180.0;

   bodyMap.at( "Apollo" )->getFlightConditions( )->getAerodynamicAngleCalculator( )->setOrientationAngleFunctions(
               [ = ]( ){ return constantAngleOfAttack; }, [ = ]( ){ return constantSideslipAngle; }, [ = ]( ){ return constantBankAngle; } );

By using this code, the vehicle will fly at 30 degree angle of attack, 2 degree sideslip angle and -70 degree bank angle during the full propagation.

Predefined guidance law
***********************
Presently, only a single predefined guidance law is implemented in Tudat: imposing pitch-trim. When using this setting, the sideslip and bank angles are set to zero, while the angle of attack is chosen such that the pitch moment coefficient is exactly zero. It is set by using:

.. code-block:: cpp

    std::shared_ptr< aerodynamics::TrimOrientationCalculator > trimCalculator = setTrimmedConditions( bodyMap.at( "Apollo" ) );

After calling this function, no additional action is needed from the user. In fact, using the following:

.. code-block:: cpp

    setTrimmedConditions( bodyMap.at( "Apollo" ) );

will work equally well. The :class:`TrimOrientationCalculator` is returned by the function to keep the object through which the computations are performed available to the user.

User-defined aerodynamic orientation
************************************
For a general description of the vehicle orientation, a custom-defined function is typically required, to fit the needs to the mission/simulation under consideration. To facilitate this process, we have defined a virtual base class called :class:`AerodynamicGuidance`.

.. class:: AerodynamicGuidance

   Virtual base class used to facilitate user-defined derived guidance classes.

A user-defined derived class must be defined, through which the orientation is computed at each time step of the propagation. Below, there are several examples of how to implement such a guidance algorithm. In each case, the final binding to the propagation is done as follows:


.. code-block:: cpp

    std::shared_ptr< aerodynamics::AerodynamicGuidance > aerodynamicGuidance =  // Create user-defined guidance object here
    setGuidanceAnglesFunctions( aerodynamicGuidance, bodyMap.at( "Apollo" ) );

An example of the computation of the three aerodynamic angles as a function of time alone can be done by using the following :class:`AerodynamicGuidance` derived class:

.. code-block:: cpp

    class LinearTimeAerodynamicGuidance: public AerodynamicGuidance
    {
        LinearTimeAerodynamicGuidance( 
            const double angleOfAttackRate, const double sideslipAngleRate, const double bankAngleRate,
                const double referenceTime ):
                    angleOfAttackRate_( angleOfAttackRate ), sideslipAngleRate_( sideslipAngleRate ), bankAngleRate_( bankAngleRate ),
                        referenceTime_( referenceTime ){ }

    void updateGuidance( const double currentTime )
    {
        currentAngleOfAttack_ = angleOfAttackRate_ * ( currentTime - referenceTime_ );
        currentAngleOfSideslip_ = sideslipAngleRate_ * ( currentTime - referenceTime_ );
        currentBankAngle_ = bankAngleRate_ * ( currentTime - referenceTime_ );       
    }

    private:

        double angleOfAttackRate_;

        double sideslipAngleRate_;

        double bankAngleRate_;

        double referenceTime_;
    };

Then, the guidance law can be created and set by:

.. code-block:: cpp

    std::shared_ptr< aerodynamics::AerodynamicGuidance > aerodynamicGuidance = std::make_shared< LinearTimeAerodynamicGuidance >( 
        1.0E-4, -2.0E-6, 1.0E-3, 500.0 );
    setGuidanceAnglesFunctions( aerodynamicGuidance, bodyMap.at( "Apollo" ) );

This creates and sets aerodynamic angles that are zero at t=500 s, where the angles of attack, sideslip and bank change by 10 -4, -2*10 -6 and 10 -3 rad/s. Recall that al units in Tudat are SI unless otherwise indicated. The key behind this implementation in the :class:`AerodynamicGuidance` derived class is the following:

   - A definition of a :literal:`void updateGuidance( const double currentTime )` function in the derived class, which is called every time step to compute the current angles as a function of time.
   - The calculation of :literal:`currentAngleOfAttack_`, :literal:`currentAngleOfSideslip_` and :literal:`currentBankAngle_` in this function. Whichever values these variables are set to in the :literal:`updateGuidance` function are the values that will be used during the current time step.

The example of aerodynamic guidance given above is not very representative, of course. In general, you will want to define your body's orientation as a function of its current state/environment, etc. To accomplish this, you can add the body map (or any its contents) as member variables to your :class:`AerodynamicGuidance` derived class. In many cases, the required information will be stored in the :class:`FlightConditions` object. 

.. class:: FlightConditions

   Class which calculates and stores data on altitude, longitude, latitude, etc. for non-atmospheric flight (stored as a member of a :class:`Body` object). Can be used during propagation to retrieve various data, using the functions:

   - :literal:`getCurrentAltitude`  Returns the current altitude w.r.t. the central body shape model

   - :literal:`getCurrentLongitude` Returns the current longitude, in the body-fixed frame of the central body

   - :literal:`getCurrentGeodeticLatitude` Returns the current geodetic latitude, in the body-fixed frame of the central body

   - :literal:`getAerodynamicAngleCalculator` Returns a :literal:`AerodynamicAngleCalculator` object, from which the current geographic latitude, longitude, flight path angle, heading angle, bank angle, sideslip angle and angle of attack can be retrieved.
 

.. class:: AtmosphericFlightConditions

   Class, derived from :class:`FlightConditions`, which stores data on altitude, longitude, latitude, density, airspeed, etc. for atmospheric flight.  Can be used during propagation to retrieve various data. It derives all functions of :class:`FlightConditions`, and also provides the functions:

   - :literal:`getCurrentDensity` Returns the current atmospheric density at the vehicle current state.

   - :literal:`getCurrentFreestreamTemperature` Returns the current atmospheric temperature at the vehicle current state.

   - :literal:`getCurrentDynamicPressure` Returns the current dynamic pressure at the vehicle current state.

   - :literal:`getCurrentPressure` Returns the current static pressure at the vehicle current state.

   - :literal:`getCurrentAirspeed` Returns the current airspeed at the vehicle current state.

   - :literal:`getCurrentSpeedOfSound` Returns the current speed of sound at the vehicle current state.

   - :literal:`getCurrentMachNumber` Returns the current Mach number at the vehicle current state.

   - :literal:`getCurrentAirspeedBasedVelocity` Returns the current velocity vectorof the vehicle w.r.t. the atmosphere, expressed in the frame fixed to the central body.

   - :literal:`getAerodynamicCoefficientInterface` Returns the object of type :class:`AerodynamicCoefficientInterface` responsible for computing and updating the aerodynamic coefficients.

   - :literal:`getAerodynamicCoefficientIndependentVariables` Returns the current list of independent variables of teh aerodynamic coefficients. For instance, if the coefficients depend on Mach number, angle of attack and sideslip angle, this function returns a vector containing these three variables (in order).

An example of an implementation of an aerodynamic guidance class is given and discussed below.

.. code-block:: cpp

    class FlightConditionsBasedAerodynamicGuidance: public AerodynamicGuidance
    {
        FlightConditionsBasedAerodynamicGuidance( 
                const NamedBodyMap& bodyMap,
                const std::string vehicleName )
        { 
            vehicleFlightConditions_ = 
               std::dynamic_pointer_cast< AtmosphericFlightConditions >( 
                  bodyMap.at( vehicleName )->getFlightConditions( ) );
            if( vehicleFlightConditions_ == nullptr )
	    {
   		throw std::runtime_error( "Error in FlightConditionsBasedAerodynamicGuidance, expected AtmosphericFlightConditions" );
	    }
        }

        void updateGuidance( const double currentTime );

    private:

        std::shared_ptr< AtmosphericFlightConditions > vehicleFlightConditions_;
    };

where the :literal:`updateGuidance` function is not defined directly in the :literal:`.h` file, but instead in the :literal:`.cpp` file (below). Note the ``dynamic_pointer_cast`` in the constructor, and the check to see if the :class:`FlightConditions` object is actually of the type :class:`AtmosphericFlightConditions`.

As an example for a guidance scheme, let's consider the simplified (and still not particularly realistic) aerodynamic guidance where:

   - Angle of attack is 35 degrees is altitude is larger than 60 km, angle of attack is 5 degrees at 30 km, and changes linearly between these two values.
   - Sideslip angle is always zero.
   - Bank angle is 80 degrees if mach number is larger than 8.

The implementation of the ``updateGuidance`` functions in the ``.cpp`` file would then read:

.. code-block:: cpp

    void FlightConditionsBasedAerodynamicGuidance::updateGuidance( const double currentTime )
    {
        if( vehicleFlightConditions_->getCurrentAltitude( ) > 60.0E3 )
        {
            currentAngleOfAttack_ = 35.0 * mathematical_constants::PI / 180.0; 
        }
        else if( vehicleFlightConditions_->getCurrentAltitude( ) < 25.0E3 )
        {
            currentAngleOfAttack_ = 5.0 * mathematical_constants::PI / 180.0; 
        }
        else
        {
            currentAngleOfAttack_ = ( 5.0 + 30.0 * ( vehicleFlightConditions_->getCurrentAltitude( ) - 25.0E3 ) / 35.0E3 ) * mathematical_constants::PI / 180.0; 

        }

        currentAngleOfSideslip_ = 0.0;

        if( vehicleFlightConditions_->getCurrentMachNumber( ) < 8 )
        {
            currentBankAngle_ = 80.0 * mathematical_constants::PI / 180.0; 
        }
        else
        {
            currentBankAngle_ = 0.0;
        }
    }

Although this guidance profile is still not very realistic for full numerical simulations, it does show the manner in which the interface is to be set up for a more realistic approach.

Using the environment models
****************************
In computing your aerodynamic guidance commands, you will likely need to use a number of physical quantities from your environment, as is the case with the example above, where the altitude is used. Below, a list is given with the way in which to retrieve some variables that are typical in aerodynamic guidance:

   - **Current conditions at a vehicle's location w.r.t. a central central body:** These are stored in an object of type :class:`FlightConditions`, which may be of the :class:`AtmosphericFlightConditions` derived class, as in the case above (stored in a :class:`Body` object; retrieved by using the :literal:`getFlightConditions` function). In the :class:`FlightConditions` and :class:`AtmosphericFlightConditions` class code, you will see a number of functions called :literal:`getCurrent...`. These functions are key in linking the environment with the guidance. When called from the :class:`AerodynamicGuidance` derived class, the current value of the associated quantity is returned (e.g. :literal:`getCurrentAltitude` returns altitude, :literal:`getCurrentAirspeed` returns airspeed, etc.).

   - **Aerodynamic coefficients:** These often play a particularly important role in the aerodynamic guidance. Whereas the other dependent variables are computed before updating the angles of attack, sideslip and bank, the aerodynamic coefficients are computed as a function of these angles. Therefore, the 'current aerodynamic coefficients' cannot yet be retrieved from the environment when updating the guidance. However, if the angles on which the aerodynamic coefficients depend have already been locally computed (in :literal:`currentAngleOfAttack_`, etc.), they may be used for determination of subsequent angles. Below is an example of aerodynamic coefficients depending on angle of attack, angle of sideslip and Mach number and the bank angle determined as a function of aerodynamic coefficients. The following can then be used inside the :literal:`updateGuidance` function:
   
   .. code-block:: cpp

        // Define aerodynamic coefficient interface/flight conditions (typically retrieved from body map; may also be a member variable)
        std::shared_ptr< aerodynamics::AerodynamicCoefficientInterface > coefficientInterface_ = ...
        std::shared_ptr< aerodynamics::AtmosphericFlightConditions > flightConditions_ = ...

        // Compute angles of attack and sideslip
        currentAngleOfAttack_ = ...
        currentAngleOfSideslip_ = ...

        // Define input to aerodynamic coefficients: take care of order of input (this depends on how the coefficients are created)!
        std::vector< double > currentAerodynamicCoefficientsInput_;
        currentAerodynamicCoefficientsInput_.push_back( currentAngleOfAttack_ );
        currentAerodynamicCoefficientsInput_.push_back( currentAngleOfSideslip_ );
        currentAerodynamicCoefficientsInput_.push_back( flightConditions_->getCurrentMachNumber( ) );

        // Update and retrieve current aerodynamic coefficients
        coefficientInterface_->updateCurrentCoefficients( currentAerodynamicCoefficientsInput_ );
        Eigen::Vector3d currentAerodynamicCoefficients = coefficientInterface_->getCurrentForceCoefficients( );

        // Compute bank angle
        currentBankAngle_ =  some function of currentAerodynamicCoefficients

   Note that the physical meaning of the coefficients may differ, depending on how they are defined in :class:`AerodynamicCoefficientSettings`: if they are defined in the aerodynamic frame (``C_D``, ``C_S``, ``C_L``)  this is how they are returned.

   - **Current vehicle orientation angles:** In particular, the angles used to define the spherical vehicle state: latitude, longitude, flight path angle and heading angle may be needed. These are retrieved from an object of type :class:`AerodynamicAngleCalculator`, which is retrieved from the :class:`FlightConditions` class with the :literal:`getAerodynamicAngleCalculator` function. The :class:`AerodynamicAngleCalculator` class in turn has a function :literal:`getAerodynamicAngle`, which takes a single argument: the type of angle that is to be returned. You can use any of the first four identifiers in the :class:`AerodynamicsReferenceFrameAngles`. In the aerodynamic guidance, DO NOT use this function to retrieve the angle of attack, sideslip or bank. As an example, you can use:

   .. code-block:: cpp
        
      // Define aerodynamic coefficient interface/flight conditions (typically retrieved from body map; may also be a member variable)
      std::shared_ptr< aerodynamics::FlightConditions > flightConditions_ = ...
      double currentFlightPathAngle = flightConditions_->getAerodynamicAngleCalculator( )->getAerodynamicAngle( reference_frames::flight_path_angle );

   - **Body mass:** The mass of the body at the current time is retrieved directly from the :class:`Body` object using the :literal:`getBodyMass( )` function.

Control surface deflections
~~~~~~~~~~~~~~~~~~~~~~~~~~~
For a realistic vehicle entry/ascent trajectory propagation, it will often be necessary to include control surface deflections in the numerical propagation. How to load/define the aerodynamic influence of control surfaces is discussed at the end of this page.

To use the control surface increments, the control surface deflections have to be set, either to a constant value before stating the propagation, or every time step by a user-defined guidance system. In each case, the control surface deflections are stored in a :class:`VehicleSystems` object, which is a member of a :class:`Body` object. The vehicle systems represent a collection of all physical (hardware) properties of a vehicle including the control surface deflections and properties. Presently, the only quantities that are stored for the control surfaces are the current deflection.

.. tip:: If your application requires more extensive functionality, please open an issue requesting this feature on Github).

In either case, a :class:`VehicleSystems` object must be created and stored in the associated :class:`Body` object:

.. code-block:: cpp

    std::shared_ptr< system_models::VehicleSystems > systemsModels = std::make_shared< system_models::VehicleSystems >( );
    bodyMap[ "Vehicle" ]->setVehicleSystems( systemsModels );

The control surface deflections are then set by:

.. code-block:: cpp

    double elevonDeflection = 0.1;
    std::string controlSurfaceId = "Elevon";
    apolloSystems->setCurrentControlSurfaceDeflection( controlSurfaceId, elevonDeflection );

Note that the deflections of multiple control surfaces can be set in exactly the same manner as follows:

.. code-block:: cpp

    apolloSystems->setCurrentControlSurfaceDeflection( "Elevon", 0.1 );
    apolloSystems->setCurrentControlSurfaceDeflection( "Aileron1", -0.15 );
    apolloSystems->setCurrentControlSurfaceDeflection( "Aileron2", 0.15 );

When only using the above, the control surfaces are set to a costant deflection throughout the propagation. This may not be very realistic but can be useful for preliminary analysis.

In general, however, you will want to determine the control surface deflections as a function of your current state, time, etc. The best way to achieve this is by incorporating the control surface deflections into the aerodynamic guidance, in particular into your specific derived class of :class:`AerodynamicGuidance` (see above). As an example, consider the following:

.. code-block:: cpp

    class FlightConditionsBasedAerodynamicExtendedGuidance: public AerodynamicGuidance
    {
        FlightConditionsBasedAerodynamicExtendedGuidance( 
                const NamedBodyMap& bodyMap,
                const std::string vehicleName )
        { 
            vehicleFlightConditions_ = 
               std::dynamic_pointer_cast< AtmosphericFlightConditions >( 
                  bodyMap.at( vehicleName )->getFlightConditions( ) );
            if( vehicleFlightConditions_ == nullptr )
	    {
   		throw std::runtime_error( "Error in FlightConditionsBasedAerodynamicGuidance, expected AtmosphericFlightConditions" );
	    }
            vehicleSystems_ = bodyMap.at( vehicleName )->getVehicleSystems( );
        }

        void updateGuidance( const double currentTime );

    private:

        std::shared_ptr< AtmosphericFlightConditions > vehicleFlightConditions_;

        std::shared_ptr< system_models::VehicleSystems > vehicleSystems_;

    };

Compared to the :class:`FlightConditionsBasedAerodynamicGuidance` class defined above, you can see that it has been extended with the ability to access the :class:`VehicleSystems` member of the associated body. In the implementation of the :literal:`updateGuidance` function, the control surface deflections may now be incorporated as follows:

.. code-block:: cpp

    void FlightConditionsBasedAerodynamicExtendedGuidance::updateGuidance( const double currentTime )
    {
        if( vehicleFlightConditions_->getCurrentAltitude( ) > 60.0E3 )
        {
            currentAngleOfAttack_ = 35.0 * mathematical_constants::PI / 180.0; 
        }
        else if( vehicleFlightConditions_->getCurrentAltitude( ) < 25.0E3 )
        {
            currentAngleOfAttack_ = 5.0 * mathematical_constants::PI / 180.0; 
        }
        else
        {
            currentAngleOfAttack_ = ( 5.0 + 30.0 * ( vehicleFlightConditions_->getCurrentAltitude( ) - 25.0E3 ) / 35.0E3 ) * mathematical_constants::PI / 180.0; 

        }

        currentAngleOfSideslip_ = 0.0;

        if( vehicleFlightConditions_->getCurrentMachNumber( ) < 8 )
        {
            currentBankAngle_ = 80.0 * mathematical_constants::PI / 180.0; 
        }
        else
        {
            currentBankAngle_ = 0.0;
        }

        double elevonDeflection = ( 1.0 + 5.0 * ( vehicleFlightConditions_->getCurrentAltitude( ) - 25.0E3 ) / 35.0E3 ) * mathematical_constants::PI / 180.0; 
        double aileron1Deflection = ( 2.0 + 7.0 * ( vehicleFlightConditions_->getCurrentAltitude( ) - 25.0E3 ) / 35.0E3 ) * mathematical_constants::PI / 180.0; 
        double aileron2Deflection = -( 2.0 + 7.0 * ( vehicleFlightConditions_->getCurrentAltitude( ) - 25.0E3 ) / 35.0E3 ) * mathematical_constants::PI / 180.0; 

        vehicleSystems_->setCurrentControlSurfaceDeflection( "Elevon", elevonDeflection );
        vehicleSystems_->setCurrentControlSurfaceDeflection( "Aileron1", aileron1Deflection );
        vehicleSystems_->setCurrentControlSurfaceDeflection( "Aileron2", aileron2Deflection );

    }

As with the previous examples, the values to which the control surface deflections are set are quite arbitrary and not based on any particularly realistic model. They are defined for illustration purposes only. A key difference between the manner in which the aerodynamic angles and the control surface deflections are handled by the guidance object is that the angles are computed but not set by the object (the angles are retrieved and set in the body model by the :class:`AerodynamicAngleCalculator`). The control surface deflections on the other hand are both computed and set by the guidance object.


