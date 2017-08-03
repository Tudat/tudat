.. _tudatFeaturesAerodynamicGuidance:

Aerodynamic Guidance
====================
When including aerodynamics in the orbit propagation, it may often be desirable to modify the aerodynamic properties of the vehicle as a function of the current environment. This is the case for ascent and entry for instance. Also, control surface deflections may be needed to modify the vehicle's behaviour. Both options are included in Tudat, and the interfaces are described below.

.. note:: In applications were no aerodynamic guidance is used, the angle of attack, sideslip angle and bank angle are all implicitly set to zero.

Vehicle orientation
~~~~~~~~~~~~~~~~~~~
At present, rotational motion cannot yet be propagated using Tudat. Therefore, to define the behaviour of a vehicle moving through an atmosphere (such as a launcher or entry vehicle), the rotational motion needs to be imposed by the user.

For aerodynamics, the vehicle orientation is typically described by the angle of attack, angle of sideslip and bank angle, which is also the approach we take in Tudat. There is a single interface in Tudat through which the functions describing these angles are fed to the trajectory propagation: the :literal:`AerodynamicAngleCalculator`` member function of the :class:`AerodynamicAngleCalculator` class is used. However, there is a number of manners in which to use this function, either directly or indirectly, the preferred selection of course depending on your application. The options are:

    - Manual definitions of the three angles of a function of time (can also be used to set constant angles).
    - Pre-defined guidance laws (only angle-of-attack trim presently implemented).
    - A user-defined :class:`AerodynamicGuidance` class, which computes the current angles using any number of input/environment variables. This option is typically the preferred option for a realistic entry propagation.

Below, we will indicate how to use these three options for the example of the Apollo entry example. Currently, the following code is used in the example, defining a constant angle of attack of 30 degree, and 0 sideslip and bank angle:

.. code-block:: cpp

    bodyMap.at( "Apollo" )->getFlightConditions( )->getAerodynamicAngleCalculator( )->setOrientationAngleFunctions(
                boost::lambda::constant( 30.0 * mathematical_constants::PI / 180.0 ) );

We stress that any definition of the aerodynamic guidance must be done after the acceleration models have been created (but before the trajectory is propagated).

Manual angle functions
**********************
To manually set the body orientation as constanr angles, the :class:`AerodynamicAngleCalculator` function is called directly, as is done in the above example. Note that if you need a definition of aerodynamic angles that varies with time, a specific :class:`AerodynamicGuidance` derived class should be implemented (see below) or one of the predefined orientation functions should be used. In general, the manual definition of the angles is done as follows:

.. code-block:: cpp
  
    double constantAngleOfAttack = 30.0 * mathematical_constants::PI / 180.0;
    double constantSideslipAngle = 2.0 * mathematical_constants::PI / 180.0;
    double constantBankAngle = -70.0 * mathematical_constants::PI / 180.0;

    bodyMap.at( "Apollo" )->getFlightConditions( )->getAerodynamicAngleCalculator( )->setOrientationAngleFunctions(
                boost::lambda::constant( constantAngleOfAttack ), boost::lambda::constant( constantSideslipAngle ), boost::lambda::constant( constantBankAngle ) );

By using this code, the vehicle will fly at 30 degree angle of attack, 2 degree sideslip angle and -70 degree bank angle during the full propagation.

Predefined guidance law
***********************
Presently, only a single predefined guidance law is implemented in Tudat: imposing pitch-trim. When using this setting, the sideslip and bank angles are set to zero, while the angle of attack is chosen such that the pitch moment coefficient is exactly zero. It is set by using:

.. code-block:: cpp

    boost::shared_ptr< aerodynamics::TrimOrientationCalculator > trimCalculator = setTrimmedConditions( bodyMap.at( "Apollo" ) );

After calling this function, not additional action is needed from the user. In fact, using the following:

.. code-block:: cpp

    setTrimmedConditions( bodyMap.at( "Apollo" ) );

will work equally well. The :class:`TrimOrientationCalculator` is returned by the function to keep the object through which the computations are performed available to the user.

User-defined aerodynamic orientation
************************************
For a general description of the vehicle orientation, a custom-defined function is typically required, to fit the needs to the mission/simulation under consideration. To facilitate this process, we have defined a virtual base class called :class:`AerodynamicGuidance`. A user-defined derived class must be defined, through which the orientation is computed at each time step of the propagation. Below, there are several examples of how to implement such a guidance algorithm. In each case, the final binding to the propagation is done as follows:

.. code-block:: cpp

    boost::shared_ptr< aerodynamics::AerodynamicGuidance > aerodynamicGuidance =  // Create user-defined guidance object here
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

    boost::shared_ptr< aerodynamics::AerodynamicGuidance > aerodynamicGuidance = boost::make_shared< LinearTimeAerodynamicGuidance >( 
        1.0E-4, -2.0E-6, 1.0E-3, 500.0 );
    setGuidanceAnglesFunctions( aerodynamicGuidance, bodyMap.at( "Apollo" ) );

This creates and sets aerodynamic angles that are zero at t=500 s, where the angles of attack, sideslip and bank change by 10 -4, -2*10 -6 and 10 -3 rad/s. Recall that al units in Tudat are SI unless otherwise indicated. The key behind this implementation in the :class:`AerodynamicGuidancederived` class is the following:

    - A definition of a :literal:`void updateGuidance( const double currentTime )` function in the derived class, which is called every time step to compute the current angles as a function of time.
    - The calculation of :literal:`currentAngleOfAttack_`, :literal:`currentAngleOfSideslip_` and :literal:`currentBankAngle_` in this function. Whichever values these variables are set to in the :literal:`updateGuidance` function are the values that will be used during the current time step.

The example of aerodynamic guidance given above is not very representative, of course. In general, you will want to define your body's orientation as a function of its current state/environment, etc. To accomplish this, you can add the body map (or any its contents) as member variables to your :class:`AerodynamicGuidance` derived class. In many cases, the required information will be stored in the :class:`FlightConditions` object, which stores data on altitude, density, airspeed, etc. To compute orientation angles from these flight conditions:

.. code-block:: cpp

    class FlightConditionsBasedAerodynamicGuidance: public AerodynamicGuidance
    {
        FlightConditionsBasedAerodynamicGuidance( 
                const NamedBodyMap& bodyMap,
                const std::string vehicleName )
        { 
            vehicleFlightConditions_ = bodyMap.at( vehicleName )->getFlightConditions( );
        }

        void updateGuidance( const double currentTime );

    private:

        boost::shared_ptr< FlightConditions > vehicleFlightConditions_;
    };

where the :literal:`updateGuidance` function is not defined directly in the :literal:`.h` file, but instead in the :literal:`.cpp` file. As an example, let's consider the simplified (and still not particularly realistic) aerodynamic guidance where:

    - Angle of attack is 35 degrees is altitude is larger than 60 km, angle of attack is 5 degrees at 30 km, and changes linearly between these two values.
    - Sideslip angle is always zero.
    - Bank angle is 80 degrees if mach number is larger than 8.

The implementation of the :class:`updateGuidance` functions would then read:

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

    - **Current conditions at a vehicle's location w.r.t. a central central body:** These are stored in an object of type :class:`FlightConditions` (stored in a :class:`Body` object; retrieved by using the :literal:`getFlightConditions` function). In the :class:`FlightConditions` class, you will see a number of functions called :literal:`getCurrent...`. When called from the :class:`AerodynamicGuidance` derived class, the current value of the associated quantity is returned (e.g. :literal:`getCurrentAltitude` returns altitude, :literal:`getCurrentAirspeed` returns airspeed, etc.).

    - **Aerodynamic coefficients:** These often play a particularly important role in the aerodynamic guidance. Whereas the other dependent variables are computed before updating the angles of attack, sideslip and bank, the aerodynamic coefficients are computed as a function of these angles. Therefore, the 'current aerodynamic coefficients' cannot yet be retrieved from the environment when updating the guidance. However, if the angles on which the aerodynamic coefficients depend have already been locally computed (in :literal:`currentAngleOfAttack_`, etc.), they may be used for determination of subsequent angles. Below is an example of aerodynamic coefficients depending on angle of attack, angle of sideslip and Mach number and the bank angle determined as a function of aerodynamic coefficients. The following can then be used inside the :literal:`updateGuidance` function:

    .. code-block:: cpp

        // Define aerodynamic coefficient interface/flight conditions (typically retrieved from body map; may also be a member variable)
        boost::shared_ptr< aerodynamics::AerodynamicCoefficientInterface > coefficientInterface_ = ...
        boost::shared_ptr< aerodynamics::FlightConditions > flightConditions_ = ...

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

    Note that the physical meaning of the coefficients may differ, depending on how they are defined: if they are defined in the aerodynamic frame (``C_D_``, ``C_S_``, ``C_L_``) this is how they are returned.

    - **Current vehicle orientation angles:** In particular, the angles used to define the spherical vehicle state: latitude, longitude, flight path angle and heading angle may be needed. These are retrieved from an object of type :class:`AerodynamicAngleCalculator`, which is retrieved from the :class:`FlightConditions` class with the :literal:`getAerodynamicAngleCalculator` function. The :class:`AerodynamicAngleCalculator` class in turn has a function :literal:`getAerodynamicAngle`, which takes a single argument: the type of angle that is to be returned. You can use any of the first four identifiers in the :class:`AerodynamicsReferenceFrameAngles` (defined in this file). In the aerodynamic guidance, DO NOT use this function to retrieve the angle of attack, sideslip or bank. As an example, you can use:

    .. code-block:: cpp
        
        // Define aerodynamic coefficient interface/flight conditions (typically retrieved from body map; may also be a member variable)
        boost::shared_ptr< aerodynamics::FlightConditions > flightConditions_ = ...
        double currentFlightPathAngle = flightConditions_->getAerodynamicAngleCalculator( )->getAerodynamicAngle( reference_frames::flight_path_angle );

    - **Body mass:** The mass of the body at the current time is retrieve directly from the :class:`Body` object using the :literal:`getBodyMass( )` function.

Control surface deflections
~~~~~~~~~~~~~~~~~~~~~~~~~~~
For a realistic vehicle entry/ascent trajectory propagation, it will often be necessary to include control surface deflections in the numerical propagation see this page to learn how to load/define the aerodynamic influence of control surfaces).

To use the control surface increments, the control surface deflections have to be set, either to a constant value before stating the propagation, or every time step by a user-defined guidance system. In each case, the control surface deflections are stored in a :class:`VehicleSystems` object, which is a member of a :class:`Body` object. The vehicle systems represent a collection of all physical (hardware) properties of a vehicle including the control surface deflections and properties. Presently, the only quantities that are stored for the control surfaces are the current deflection.

.. tip:: If your application requires more extensive functionality, please open an issue requesting this feature on Github).

In either case, a :class:`VehilceSystems` object must be created and stored in the associated :class:`Body` object:

.. code-block:: cpp

    boost::shared_ptr< system_models::VehicleSystems > systemsModels = boost::make_shared< system_models::VehicleSystems >( );
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
            vehicleFlightConditions_ = bodyMap.at( vehicleName )->getFlightConditions( );
            vehicleSystems_ = bodyMap.at( vehicleName )->getVehicleSystems( );
        }

        void updateGuidance( const double currentTime );

    private:

        boost::shared_ptr< FlightConditions > vehicleFlightConditions_;

        boost::shared_ptr< system_models::VehicleSystems > vehicleSystems_;

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

As with the previous examples, the values to which the control surface deflections are set are quite arbitrary and not based on any particularly realistic model. They are defined for illustration purposes only. A key difference between the manner in which the aerodynamic angles and the control surface deflections are handled by the guidance object is that the angles are computed but not set by the object (the angles are retrieved and set in the body model by the :class:`AerodynamicAngleCalculator`). The control surface deflections on teh other hand are both computed and set by the guidance object.

Reading aerodynamic coefficients from Files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
For many simulations/analyses involving atmospheric flight, the aerodynamic coefficients will be provided in tabulated form. If put into the correct file format, these files can be read into Tudat used during the orbit propagation. By specifying the physical meaning of the independent variables of the aerodynamic coefficients, no action on the side of the user is required to update the aerodynamic coefficients to their correct values during propagation. Here, we give an overview and some examples on how to load aerodynamic coefficients from a file.

Loading the coefficient settings
********************************
As a reminder, aerodynamic coefficients are created in Tudat by creating an object of type :class:`AerodynamicCoefficientSettings`:

    .. code-block:: cpp
    
        boost::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings = .....
        bodyMap[ "VehicleName" ]->setAerodynamicCoefficientInterface(
                            createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "VehicleName" ) );

To create an :class:`AerodynamicCoefficientSettings` object from data in files, we provide two functions named :literal:`readTabulatedAerodynamicCoefficientsFromFiles` (in :literal:`createFlightConditions.h`). For one of the functions, only force coefficients are loaded (with moment coefficients set to zero at all times). The other function allows both force and moment coefficients to be loaded.

For the situation where only force coefficients are considered, several pieces information are needed:

    - A list of files for any of the three aerodynamic coefficients (e.g. C_D, C S, C L or C X, C Y, C Z). Note that the behaviour of each coefficient must be provided in a separate file. Note that not every coefficient needs to be defined. If a file is not provided for one of the coefficients, as will often be the case for C S, zeros are assumed at all points in the propagation.
    - The physical meaning of each of the independent variables of the coefficients.
    - The reference area for the aerodynamics. This is not read from the file and must be provided as an input to the function.
    - Two booleans denoting the orientation and direction of the aerodynamic coefficients. For instance C D, C S, C L denote the strength of the aerodynamic force in the aerodynamic reference frame, in a direction opposite to the axes of that frame. The C X, C Y, C Z coefficients, on the other hand, are defined in the body-fixed frame.

As an example, the following can be used to create AerodynamicCoefficientSettings for force coefficients only from a file:

    .. code-block:: cpp
    
        double referenceArea = 50.0; // Define reference area

        // Define physical meaning of independent variables, in this case Mach number and angle of attack
        std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables > independentVariableNames; 
        independentVariableNames.push_back( aerodynamics::mach_number_dependent );
        independentVariableNames.push_back( aerodynamics::angle_of_attack_dependent ); 

        // Define list of files for force coefficients. Entry 0 denotes the x-direction (C ~D~/C ~X~), 1 the y-direction (C ~S~/C ~Y~) and 2 the z-direction (C ~L~/C ~Z~)
        std::map< int, std::string > forceCoefficientFiles; 
        forceCoefficientFiles[ 0 ] = tudat::input_output::getTudatRootPath( ) + "Astrodynamics/Aerodynamics/UnitTests/aurora_CD.txt"; // Set drag coefficient file
        forceCoefficientFiles[ 2 ] = tudat::input_output::getTudatRootPath( ) + "Astrodynamics/Aerodynamics/UnitTests/aurora_CL.txt"; // Set lift coefficient file

        // Define reference frame in which the loaded coefficients are defined.
        bool areCoefficientsInAerodynamicFrame = true;
        bool areCoefficientsInNegativeAxisDirection = true;

        // Load and parse files; create coefficient settings.
        boost::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings = 
            readTabulatedAerodynamicCoefficientsFromFiles( forceCoefficientFiles, referenceArea, independentVariableNames, areCoefficientsInAerodynamicFrame,         areCoefficientsInNegativeAxisDirection );

        // Create and set aerodynamic coefficients
        bodyMap[ "VehicleName" ]->setAerodynamicCoefficientInterface(
                            createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "VehicleName" ) );

Note that in the above, not side force coefficient file (entry 1 for forceCoefficientFiles) is given, so that C S=0 always. Two independent variables have been defines (Mach number and angle of attack). If either the lift or drag coefficient files encounter a different number of independent variables, the program will terminate with an appropriate error message. Also, if the independent variables used for the lift and drag coefficients are not identical, the program is terminated.

.. tip:: Moment coefficients are added in a completely analogous manner (with separate files for the x-, y- and z-components).

Adding control surface influence
********************************
In addition to defining aerodynamic coefficients for the vehicle itself, the influence of control surface deflections on the values of the coefficients will be needed for certain applications. In Tudat, any number of control surfaces may be defined for a vehicle, the deflection of which may be set by your particular guidance model (see here). Loading the aerodynamic coefficient increments of the control surface is done in a manner similar to those of the total vehicle, but:

    - Reference area and reference frame in which the coefficients are defined are not provided. These are implcitily assumed to be equal to those of the aerodynamic coefficients of the 'main body'. If your application requires these quantities to be different for the body and control surface deflections, please open an issue on Github requesting the functionality.
    - Exactly one of the independent variables of the coefficient increments must be a control surface deflection.

Below, an example is given on how to load the aerodynamic coefficient increments:

    .. code-block:: cpp
    
        // Create coefficient settings for body.
        boost::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings = ...

        // Define physical meaning of independent variables for control surface increments, in this case Mach number, angle of attack and control surface deflection
        std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables > controlSurfaceIndependentVariableNames; 
        controlSurfaceIndependentVariableNames.push_back( aerodynamics::mach_number_dependent );
        controlSurfaceIndependentVariableNames.push_back( aerodynamics::angle_of_attack_dependent ); 
        controlSurfaceIndependentVariableNames.push_back( aerodynamics::control_surface_deflection_dependent ); 

        // Define name of control surface
        std::string controlSurfaceName = "Elevon";

        // Define list of files for force coefficients. 
        std::map< int, std::string > controlSurfaceForceCoefficientFiles; 
        controlSurfaceForceCoefficientFiles[ 0 ] = tudat::input_output::getTudatRootPath( ) + "Astrodynamics/Aerodynamics/UnitTests/dCDwTest.txt"; // Set drag coefficient file

        // Add settings for control surface increments to main aerodynamic coefficients
        aerodynamicCoefficientSettings->setControlSurfaceSettings( 
            readTabulatedControlIncrementAerodynamicCoefficientsFromFiles( controlSurfaceForceCoefficientFiles, controlSurfaceIndependentVariableNames, controlSurfaceName ) );

        // Create and set aerodynamic coefficients
        bodyMap[ "VehicleName" ]->setAerodynamicCoefficientInterface(
            createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "VehicleName" ) );

For this example, only the drag coefficient is affected by the control surface deflections.

