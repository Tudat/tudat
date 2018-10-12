.. _tudatFeaturesAerodynamicCoefficients:

Aerodynamic Coefficients
~~~~~~~~~~~~~~~~~~~~~~~~
For a body that experiences aerodynamic forces, the aerodynamic coefficients of that body need to be defined. Aerodynamic coefficients are created in Tudat by creating an object of type :class:`AerodynamicCoefficientSettings`:

.. code-block:: cpp
    
    std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings = .....
    bodyMap[ "VehicleName" ]->setAerodynamicCoefficientInterface(
                        createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "VehicleName" ) );

The aerodynamic coefficients can be implemented in several different ways that are discussed in this section.

Local Incliniation Methods
**************************
One of the options to determine the aerodynamic coefficients during hypersonic flight is to use a hypersonic local inclination method. These methods calculate the aerodynamic coefficients using the orientation of the vehicle with respect to the incoming flow. A selection of these methods are implemented into Tudat and are able to calculate the aerodynamic coefficients for several pre-defined geometries and user defined grids.

The class which can be used is defined below.

.. class:: HypersonicLocalInclinationAnalysis

The settings for this class are constructed as follows:

.. code-block:: cpp
    
    HypersonicLocalInclinationAnalysis(
            const std::vector< std::vector< double > >& dataPointsOfIndependentVariables,
            const std::shared_ptr< SurfaceGeometry > inputVehicleSurface,
            const std::vector< int >& numberOfLines,
            const std::vector< int >& numberOfPoints,
            const std::vector< bool >& invertOrders,
            const std::vector< std::vector< int > >& selectedMethods,
            const double referenceArea,
            const double referenceLength,
            const Eigen::Vector3d& momentReferencePoint )

where:

- :literal:`dataPointsOfIndependentVariables`:
      A 2-dimensional vector of doubles that contains the data points of each independent variable. Index 0 corresponds to the Mach number, index 1 to the angle of attack, and index 2 to the
      sideslip angle.

- :literal:`inputVehicleSurface`:
      A :literal:`SurfaceGeometry` type variable that defines the geometry of the vehicle used for the determination of the aerodynamic coefficients. 

- :literal:`numberOfLines`:
      The number of discretization points in the first independent surface for each subpart of the vehicle defined in :literal:`inputVehicleSurface`. 

- :literal:`numberOfPoints`:
      The number of discretization points in the second independent surface for each subpart of the vehicle defined in :literal:`inputVehicleSurface`. 

- :literal:`invertOrders`:
      A boolean that determines the orientation of the normal vectors of the panels (the discretized units) of each subpart of the vehicle defined in :literal:`inputVehicleSurface`. This can either be
      inward- or outward-facing.

- :literal:`selectedMethods`:
      A two dimensional vector that selects which specific method is used for the local inclination analysis. The first index determines if a compression (0) or a expansion (1) method is used, and
      the second index determines the specific method (a list of available methods is given below).

- :literal:`referenceArea`:
      A :literal:`double` that gives the reference area of the vehicle that is used for the analysis.

- :literal:`referencelength`:
      A :literal:`double` that gives the reference length of the vehicle that is used for the analysis.

- :literal:`momentReferencePoint`:
      A :literal:`Eigen::Vector3d` that gives the location of the moment reference point, which is used as a reference point to calculate the moments acting on the vehicle.


The :literal:`HypersonicLocalInclinationAnalysis` has several methods incorporated in it, that can be selected by the user through the :literal:`selectedMethods` setting. The first index of this vector determines which kind of method is used: a compression or expansion method. The second index then determines the specific method.

For the compression methods, the following are available:

- 0: Newtonian Method.
- 1: Modified Newtonian.
- 2 and 3: not available at this moment.
- 4: tangent-wedge method.
- 5: tangent-cone method.
- 6: modified Dahlem-Buck method.
- 7: VanDyke unified pressure method.
- 8: Smyth Delta Wing method.
- 9: Hankey flat surface method.

The expansion method has the following options:

- 0: Vacuum Pressure coefficient method.
- 1: Zero Pressure function.
- 4: High Mach base pressure method.
- 3 or 5: Prandtl-Meyer method.
- 6: ACM empirical pressure coefficient. 

To get a better idea of how this can be implemented, an example contained in the aerodynamics unit test of Tudat is used. This example describes how the aerodynamic coefficients are generated for the Apollo capsule. A function called :literal:`getApolloCoefficientInterface()` is made which returns a variable of the type: :literal:`HypersonicLocalInclinationAnalysis`.

First the capsule geometry is made using a pre-defined capsule shape:

.. code-block:: cpp
    
    std::shared_ptr< geometric_shapes::Capsule > capsule
            = std::make_shared< geometric_shapes::Capsule >(
               4.694, 1.956, 2.662, -1.0 * 33.0 * PI / 180.0, 0.196 );

There is also an option to define a grid by the user, but this is not explained here. To define the panels on the vehicle, which are used for the calculation of the local inclination, the number of lines and points are easily definable:

.. code-block:: cpp
    
    std::vector< int > numberOfLines;
    std::vector< int > numberOfPoints;
    std::vector< bool > invertOrders;
    numberOfLines.resize( 4 );
    numberOfPoints.resize( 4 );
    invertOrders.resize( 4 );

    // Set number of analysis points.
    numberOfLines[ 0 ] = 31;
    numberOfPoints[ 0 ] = 31;
    numberOfLines[ 1 ] = 31;
    numberOfPoints[ 1 ] = 31;
    numberOfLines[ 2 ] = 31;
    numberOfPoints[ 2 ] = 10;
    numberOfLines[ 3 ] = 11;
    numberOfPoints[ 3 ] = 11;
    invertOrders[ 0 ] = 0;
    invertOrders[ 1 ] = 0;
    invertOrders[ 2 ] = 0;
    invertOrders[ 3 ] = 0;

A moment reference point is then defined. Afterwards, the independent variable points are defined, which are used to calculate the coefficients in the simulation during specific flight conditions. These points can be defined by the user, but there is also an option to get default values:

.. code-block:: cpp
    
    std::vector< std::vector< double > > independentVariableDataPoints;
    independentVariableDataPoints.resize( 3 );
    independentVariableDataPoints[ 0 ] = getDefaultHypersonicLocalInclinationMachPoints( "Full" );
    std::vector< double > angleOfAttackPoints;
    angleOfAttackPoints.resize( 15 );

    for ( int i = 0; i < 15; i++ )
    {
        angleOfAttackPoints[ i ] = static_cast< double >( i - 6 ) * 5.0 * PI / 180.0;
    }

    independentVariableDataPoints[ 1 ] = angleOfAttackPoints;
    independentVariableDataPoints[ 2 ] =
    getDefaultHypersonicLocalInclinationAngleOfSideslipPoints( );

For the Mach number and angle of sideslip, default values are used, whereas for the angle of attack, the values are defined by the user. The final part of the code is the selection of methods. In this case there are several methods used to determine the coefficients. 

.. code-block:: cpp
    
    std::vector< std::vector< int > > selectedMethods;
    selectedMethods.resize( 2 );
    selectedMethods[ 0 ].resize( 4 );
    selectedMethods[ 1 ].resize( 4 );

    selectedMethods[ 0 ][ 0 ] = 1;
    selectedMethods[ 0 ][ 1 ] = 5;
    selectedMethods[ 0 ][ 2 ] = 5;
    selectedMethods[ 0 ][ 3 ] = 1;
    selectedMethods[ 1 ][ 0 ] = 6;
    selectedMethods[ 1 ][ 1 ] = 3;
    selectedMethods[ 1 ][ 2 ] = 3;
    selectedMethods[ 1 ][ 3 ] = 3;


Finally, in the return statement, the local inclination analysis is made, which can be used to define an :literal:`AerodynamicCoefficientInterface`:

.. code-block:: cpp
    
    // Create analysis object and capsule database.
    return std::make_shared< HypersonicLocalInclinationAnalysis >(
                independentVariableDataPoints, capsule, numberOfLines, numberOfPoints,
                invertOrders, selectedMethods, PI * pow( capsule->getMiddleRadius( ), 2.0 ),
                3.9116, momentReference ); 

.. _tudatFeaturesAerodynamicGuidanceReadingAerodynamicCoefficients:
 
Reading aerodynamic coefficients from Files
*******************************************
For many simulations/analyses involving atmospheric flight, the aerodynamic coefficients will be provided in tabulated form. If put into the correct file format, these files can be read into Tudat used during the orbit propagation. By specifying the physical meaning of the independent variables of the aerodynamic coefficients, no action on the side of the user is required to update the aerodynamic coefficients to their correct values during propagation. Here, we give an overview and some examples on how to load aerodynamic coefficients from a file.

As a reminder, aerodynamic coefficients are created in Tudat by creating an object of type :class:`AerodynamicCoefficientSettings`:

.. code-block:: cpp
    
    std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings = .....
    bodyMap[ "VehicleName" ]->setAerodynamicCoefficientInterface(
                        createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "VehicleName" ) );

To create an :class:`AerodynamicCoefficientSettings` object from data in files, we provide two functions named :literal:`readTabulatedAerodynamicCoefficientsFromFiles` (in :literal:`createFlightConditions.h`). For one of the functions, only force coefficients are loaded (with moment coefficients set to zero at all times). The other function allows both force and moment coefficients to be loaded.

For the situation where only force coefficients are considered, several pieces information are needed:

   - A list of files for any of the three aerodynamic coefficients (e.g. C\ :sub:`D`, C\ :sub:`S`, C\ :sub:`L` or C\ :sub:`X`, C\ :sub:`Y`, C\ :sub:`Z`). Note that the behaviour of each coefficient must be provided in a separate file. Note that not every coefficient needs to be defined. If a file is not provided for one of the coefficients, as will often be the case for C\ :sub:`S`, zeros are assumed at all points in the propagation.
   - The physical meaning of each of the independent variables of the coefficients.
   - The reference area for the aerodynamics. This is not read from the file and must be provided as an input to the function.
   - Two booleans denoting the orientation and direction of the aerodynamic coefficients. For instance C\ :sub:`D`, C\ :sub:`S`, C\ :sub:`L` denote the strength of the aerodynamic force in the aerodynamic reference frame, in a direction opposite to the axes of that frame. The C\ :sub:`X`, C\ :sub:`Y`, C\ :sub:`Z` coefficients, on the other hand, are defined in the body-fixed frame.

As an example, the following can be used to create :class:`AerodynamicCoefficientSettings` for force coefficients only from a file:

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
        std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings = 
            readTabulatedAerodynamicCoefficientsFromFiles( forceCoefficientFiles, referenceArea, independentVariableNames, areCoefficientsInAerodynamicFrame,         areCoefficientsInNegativeAxisDirection );

        // Create and set aerodynamic coefficients
        bodyMap[ "VehicleName" ]->setAerodynamicCoefficientInterface(
                            createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "VehicleName" ) );

Note that in the above, no side force coefficient file (entry 1 for forceCoefficientFiles) is given, so that C\ :sub:`S`=0 always. Two independent variables have been defines (Mach number and angle of attack). If either the lift or drag coefficient files encounter a different number of independent variables, the program will terminate with an appropriate error message. Also, if the independent variables used for the lift and drag coefficients are not identical, the program is terminated.

.. tip:: Moment coefficients are added in a completely analogous manner (with separate files for the x-, y- and z-components).

In addition to defining aerodynamic coefficients for the vehicle itself, the influence of control surface deflections on the values of the coefficients will be needed for certain applications. In Tudat, any number of control surfaces may be defined for a vehicle, the deflection of which may be set by your particular :class:`AerodynamicGuidance` derived class. Loading the aerodynamic coefficient increments of the control surface is done in a manner similar to those of the total vehicle, but:

    - Reference area and reference frame in which the coefficients are defined are not provided. These are implcitily assumed to be equal to those of the aerodynamic coefficients of the 'main body'. If your application requires these quantities to be different for the body and control surface deflections, please open an issue on Github requesting the functionality.
    - Exactly one of the independent variables of the coefficient increments must be a control surface deflection.

Below, an example is given on how to load the aerodynamic coefficient increments:

.. code-block:: cpp
    
    // Create coefficient settings for body.
    std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings = ...

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

.. _tudatFeaturesCustomAerodynamicCoefficients:
 
Custom Aerodynamic Coefficient Settings
****************************************
If a user specific aerodynamic coefficient interface is needed, the :literal:`CustomAerodynamicCoefficientInterface` can be used. This class allows a generic aerodynamic and moment coefficient function input, and allows the user to use all the aerodynamic coefficient interface methods. 

.. class:: CustomAerodynamicCoefficientInterface

The constructor for this class looks as follows:

.. code-block:: cpp
    
    CustomAerodynamicCoefficientInterface(
            const std::function< Eigen::Vector3d( const std::vector< double >& ) >
            forceCoefficientFunction,
            const std::function< Eigen::Vector3d( const std::vector< double >& ) >
            momentCoefficientFunction,
            const double referenceLength,
            const double referenceArea,
            const double lateralReferenceLength,
            const Eigen::Vector3d& momentReferencePoint,
            const std::vector< AerodynamicCoefficientsIndependentVariables >
            independentVariableNames,
            const bool areCoefficientsInAerodynamicFrame = true,
            const bool areCoefficientsInNegativeAxisDirection = true )

where:

- :literal:`forceCoefficientFunction`:
      A function outputting the aerodynamic force coefficients as a function of the independent variables, listed in :literal:`independentVariableNames`. This requires a :literal:`std::function
      and :literal:`std::bind`, which is explained in: :ref:`externalUtility`.

- :literal:`momentCoefficientFunction`:
      A function outputting the aerodynamic moment coefficients as a function of the independent variables, listed in :literal:`independentVariableNames`. This requires a :literal:`std::function
      and :literal:`std::bind`, which is explained in: :ref:`externalUtility`. 

- :literal:`referenceLength`:
      The length with which the coefficients are non-dimensionalized (about the z- and x-axis). 

- :literal:`referenceArea`:
      The area with which the coefficients are non-dimensionalized. 

- :literal:`lateralReferenceLength`:
      The length with which the coefficients are non-dimensionalized (about the y-axis). 

- :literal:`momentReferencePoint`:
      The point with respect to the aerodynamic moments are calculated.

- :literal:`independentVariableNames`:
      A vector containing identifiers for the independent variables that are used in the :literal:`forceCoefficientFunction`.

- :literal:`areCoefficientsInAerodynamicFrame`:
      A :literal:`bool` that determines if the aerodynamic coefficients are in the aerodynamic reference frame or not.

- :literal:`areCoefficientsInNegativeAxisDirection`:
      A :literal:`bool` that determines if the aerodynamic coefficients are in the negative axis direction.