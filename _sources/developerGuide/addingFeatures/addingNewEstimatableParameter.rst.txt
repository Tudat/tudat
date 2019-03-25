.. _addingNewEstimatableParameter:

Adding a New Estimatable Parameter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The list of estimatable parameters already available in Tudat is presented in :ref:`parameterEstimationSettings`. However, it is possible to add another parameter to this list of estimatable parameters if needed. This process requires to modify several files located in different directories and  will be described in details based on the examples of the two parameters :literal:`radiation_pressure_coefficient` and :literal:`rotation_pole_position`.

Updating the list of estimatable parameters
*******************************************

First of all, the name of the new estimatable parameter has to be added to the list of the estimatable parameters available in Tudat, in the file :literal:`estimatableParameter.h`:

.. code-block:: cpp

   //! List of parameters that can be estimated by the orbit determination code
   enum EstimatebleParametersEnum
   {
      arc_wise_initial_body_state,
      initial_body_state,
      initial_rotational_body_state,
      gravitational_parameter,
      constant_drag_coefficient,
      radiation_pressure_coefficient,
      arc_wise_radiation_pressure_coefficient,
      spherical_harmonics_cosine_coefficient_block,
      spherical_harmonics_sine_coefficient_block,
      constant_rotation_rate,
      rotation_pole_position,
      constant_additive_observation_bias,
      ...
      ...
   }

In addition to the file :literal:`estimatableParameters.h`, the file :literal:`estimatableParameters.cpp` also has to be modified. In particular, a short description of each estimatable parameter has to be provided in the function :literal:`getParameterTypeString`, as it is done in the following for the parameters :literal:`radiation_pressure_coefficient` and :literal:`rotation_pole_position`:

.. code-block:: cpp
   
   std::string getParameterTypeString( const EstimatebleParametersEnum parameterType )
   {
   ...
   case radiation_pressure_coefficient:
        parameterDescription = "radiation pressure coefficient ";
        break;
   ...
   case rotation_pole_position:
        parameterDescription = "pole position ";
        break;
   ...
   }

The type of the new estimatable parameter must then be specified within the function :literal:`isDoubleParameter` (still inside the file :literal:`estimatableParameters.cpp`). An estimatable parameter can either be a :literal:`double` or a :literal:`Eigen::VectorXd`, depending on whether the parameter is represented by a single or multiple floating-point values. Regarding the two examples which are considered here, the parameter :literal:`radiation_pressure_coefficient` is a :literal:`double` while the parameter :literal:`rotation_pole_position` is a :literal:`Eigen::VectorXd`, whose first element is the right ascension of the rotation pole and the second one its declination.

.. code-block:: cpp
   
   bool isDoubleParameter( const EstimatebleParametersEnum parameterType )
   {
   ...
   case radiation_pressure_coefficient:
        isDoubleParameter = true;
        break;
   ...
   case rotation_pole_position:
        isDoubleParameter = false;
        break;
   ...
   }

Finally, depending on which estimatable parameter is to be added to the list of available parameters, the functions :literal:`isParameterDynamicalPropertyInitialState`, :literal:`isParameterRotationMatrixProperty`, :literal:`isParameterObservationLinkProperty` and :literal:`isParameterTidalProperty` may also need to be modified to take this new parameter into account (again in :literal:`estimatableParameters.cpp`). These function return a boolean which is set to false as default value, but needs to be set to true for specific types of parameters. Specifically, if the parameter represents an initial dynamical state (e.g. initial translational state), a property influencing rotation matrix (e.g. planetary rotation rate), an observation link (e.g. observation bias) or a property influencing tidal gravity field variations (e.g. Love numbers; tidal lags). For our examples, the parameter :literal:`radiation_pressure_coefficient` is related to none of them but the parameter :literal:`rotation_pole_position` is linked to a rotation matrix so that the function :literal:`isParameterRotationMatrixProperty` includes a case switching the boolean to :literal:`true` for this parameter.

.. code-block:: cpp
   
   bool isParameterRotationMatrixProperty( const EstimatebleParametersEnum parameterType )
   {
      bool flag;
      switch( parameterType )
      {
      ...
      case rotation_pole_position:
         flag = true;
         break;
      }
      return flag;
   }


Creating an :literal:`EstimatableParameter` object
**************************************************

Once the new estimatable parameter has been defined and characterised, a corresponding class is to be created to fully describe this parameter. This class is defined in a separate file which has to be created in the directory:

      .../tudatBundle/tudat/Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/

Each estimatable parameter class is defined in a similar way, usually from the base class :literal:`EstimatableParameter` and includes the definition of three functions :literal:`getParameterValues`, :literal:`setParameterValue` and :literal:`getParameterSize`. Considering the parameter :literal:`radiation_pressure_coefficient`, the definition of its associated class :literal:`RadiationPressureCoefficient` is done as follows (in the file :literal:`radiationPressureCoefficient.h`).

.. code-block:: cpp

   class RadiationPressureCoefficient: public EstimatableParameter< double >
   {

   public:
   //! Constructor.
   RadiationPressureCoefficient(std::shared_ptr< electro_magnetism::RadiationPressureInterface > radiationPressureInterface, std::string& associatedBody ):
   EstimatableParameter< double >( radiation_pressure_coefficient, associatedBody ), radiationPressureInterface_( radiationPressureInterface )
   { }

   //! Destructor.
   ~RadiationPressureCoefficient( ) { }

   //! Function to get the current value of the radiation pressure coefficient that is to be estimated.
   double getParameterValue( )
   {
       return radiationPressureInterface_->getRadiationPressureCoefficient( );
   }

   //! Function to reset the value of the radiation pressure coefficient that is to be estimated.
   void setParameterValue( double parameterValue )
   {
       radiationPressureInterface_->resetRadiationPressureCoefficient( parameterValue );
   }

   //! Function to retrieve the size of the parameter.
   int getParameterSize( ){ return 1; }

   protected:

   private:

   //! Object containing the radiation pressure coefficient to be estimated.
   std::shared_ptr< electro_magnetism::RadiationPressureInterface > radiationPressureInterface_;
   };


For the parameter :literal:`rotation_pole_position`, the class :literal:`ConstantRotationalOrientation` is created in a very similar way, in a file :literal:`constantRotationalOrientation.h`, keeping in mind that this parameter is a :literal:`Eigen::VectorXd` whose size is 2 and not a double as it is the case for the radiation pressure coefficient.


To create the parameter class, a settins class and an associated 'create' function is needs. This is done in the file :literal:`createEstimatableParameters.cpp/.h`, either with the function :literal:`createDoubleParameterToEstimate` if the estimatable parameter is of type :literal:`double` or with :literal:`createVectorParameterToEstimate` if it is a :literal:`Eigen::VectorXd` object. These functions take an :class:`EstimatableParameterSettings` object as input. In the create function, the settings of the parameter to be estimated are checked to ensure that it is properly defined and thus make the estimation possible. 

The information defined in the the base class :class:`EstimatableParameterSettings` is:

* Type of estimated parameter (as an :literal:`EstimatebleParametersEnum`)
* The body associated with the parameter
* A secondary string identifier (e.g. ground station name for :literal:`ground_station_position`)

In some cases, this base class is not detailed enough to give access to all the required properties of the estimatable parameter and a specific class has to be defined. As an example, the parameter :literal:`spherical_harmonics_cosine_coefficient_block` requires the definition of the estimatable parameter settings class called :literal:`SphericalHarmonicEstimatableParameterSettings` (defined in the file :literal:`estimatableParameterSettings.h`). This class manages the different combinations of degrees and orders of the spherical harmonic coefficients that have to be estimated.

Regarding the verification of the estimatable parameter settings of :literal:`radiation_pressure_coefficient`, it is checked that only one single radiation pressure interface is defined before creating the parameter and linking it to this radiation pressure interface and to the propagated body. 

.. code-block:: cpp

   std::shared_ptr< EstimatableParameter< double > > createDoubleParameterToEstimate(
        const std::shared_ptr< EstimatableParameterSettings >& doubleParameterName,
        const NamedBodyMap& bodyMap, const basic_astrodynamics::AccelerationMap& accelerationModelMap )
   {
   ...
   case radiation_pressure_coefficient:
        {
            if( currentBody->getRadiationPressureInterfaces( ).size( ) == 0 )
            {
                std::string errorMessage = "Error, no radiation pressure interfaces found in body " +
                        currentBodyName + " when making Cr parameter.";
                throw std::runtime_error( errorMessage );
            }
            else if( currentBody->getRadiationPressureInterfaces( ).size( ) > 1 )
                std::string errorMessage = "Error, multiple radiation pressure interfaces found in body " +
            {
                        currentBodyName + " when making Cr parameter.";
                throw std::runtime_error( errorMessage );
            }
            else
            {
                doubleParameterToEstimate = std::make_shared< RadiationPressureCoefficient >(
                            currentBody->getRadiationPressureInterfaces( ).begin( )->second,
                            currentBodyName );
            }
            break;
        }
    ...
   }

Concerning the parameter :literal:`rotation_pole_position`, it must be verified that the rotation model is a simple rotational ephemeris for which the position of the rotation pole is indeed defined before creating the estimatable parameter.

.. code-block:: cpp

   std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > createVectorParameterToEstimate(
        const std::shared_ptr< EstimatableParameterSettings >& vectorParameterName,
        const NamedBodyMap& bodyMap, const basic_astrodynamics::AccelerationMap& accelerationModelMap )
   {
   ...
   case rotation_pole_position:
            if( std::dynamic_pointer_cast< SimpleRotationalEphemeris >( currentBody->getRotationalEphemeris( ) ) == nullptr )
            {
                std::string errorMessage = "Warning, no simple rotational ephemeris present in body " + currentBodyName +
                        " when making constant rotation orientation parameter";
                throw std::runtime_error( errorMessage );
            }
            else
            {
                vectorParameterToEstimate = std::make_shared< ConstantRotationalOrientation >
                        ( std::dynamic_pointer_cast< ephemerides::SimpleRotationalEphemeris >
                          ( currentBody->getRotationalEphemeris( ) ), currentBodyName );

            }
            break;
     ...
     }

Note that in both these cases, the :class:`EstimatableParameterSettings` class is sufficient to define the parameter fully, and no dedicated settings class is required.

Implementing the partials required for the estimation
*****************************************************

To allow the parameter estimation to be conducted, partials with respect to this parameter have to be implemented. These derivatives have to be implemented to 'tell' the code how the parameters influence the observations and/or dynamics of a body.

For partial derivatives of observations w.r.t. a parameter, the dependency of observation on parameter is typically a kinematic one, in the sense that a different current parameter value directly changes the observation. A typical example of this is the position of a ground station. This is opposed to parameters that influence the dynamics of a body, which influence the state of the spacecraft, which in turn influences the observation. Examples of such parameters are spherical harmonic coefficients and drag coefficients.

Due to the kinematic nature of the observation partials, they are implemented as partials of the Cartesian state of the propagated body with respect to the estimatable parameter (which are then mapped to an observation partial for each observable by a dedicated class). These Cartesian state partials have to be implemented in the file :literal:`createCartesianStatePartial.cpp`, within the functions :literal:`createCartesianStatePartialsWrtParameter` (two functions exist with two different input types, depending on the type of the parameter (:literal:`double` or :literal:`Eigen::VectorXd`) that is to be considered). 

If the estimatable parameter has been identified as being a property of a rotation matrix, then the function :literal:`createCartesianStatePartialsWrtParameter` calls another function named :literal:`createRotationMatrixPartialsWrtParameter` and defined in :literal:`createCartesianStatePartial.cpp` too (again two functions with the same name exist for the two types of estimatable parameters: VectorXd and double parameters). A specific case has to be added within this function for each parameter which is related to a rotation matrix. For the parameter :literal:`rotation_pole_position`, the following lines of code have been added to first check that the rotation model is consistent with the estimatable parameter (here that is a simple rotational model) and then to call a function that returns the required rotation matrix  partials for this particular model.

.. code-block:: cpp

   std::shared_ptr< RotationMatrixPartial > createRotationMatrixPartialsWrtParameter(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameterToEstimate )

  {
     ...
     case estimatable_parameters::rotation_pole_position:
        if( std::dynamic_pointer_cast< ephemerides::SimpleRotationalEphemeris >(
                    currentBody->getRotationalEphemeris( ) ) == nullptr )
        {
            std::string errorMessage = "Warning, body's rotation model is not simple when making position w.r.t. pole position partial";
            throw std::runtime_error( errorMessage );
        }
        // Create rotation matrix partial object
        rotationMatrixPartial = std::make_shared< RotationMatrixPartialWrtPoleOrientation >(
                    std::dynamic_pointer_cast< SimpleRotationalEphemeris>( currentBody->getRotationalEphemeris( ) ) );
        break;

   ...
   }

The class :literal:`RotationMatrixPartialWrtPoleOrientation` used in the code above to define the :literal:`rotationMatrixPartial` variable is created in the file :literal:`rotationMatrixPartial.h` and must contain two internal functions :literal:`calculatePartialOfRotationMatrixToBaseFrameWrtParameter` and :literal:`calculatePartialOfRotationMatrixDerivativeToBaseFrameWrtParameter`. These functions return the rotation matrix and rotation matrix derivative partials respectively, with respect to the estimatable parameter. Specific functions to calculate these partials have to be added to the file :literal:`rotationMatrixPartial.cpp`. Regarding the parameter :literal:`rotation_pole_position`, these functions are :literal:`calculatePartialOfRotationMatrixFromLocalFrameWrtPoleOrientation` and :literal:`calculatePartialOfRotationMatrixDerivativeFromLocalFrameWrtPoleOrientation`, respectively.

.. note::
	When implementing a partial derivative of a rotation matrix, this will *automatically* also include the dynamic influence of a change in rotational parameter on a spherical harmonic acceleration.

So far, we have only considered the case where the estimatable parameter is related to a rotation matrix. However, when it is not the case, the function :literal:`createCartesianStatePartialsWrtParameter` is to be modified in a different way. A specific case has to be created for each parameter that is not a rotation matrix property. If the parameter has a direct impact on the Cartesian state of the propagated body (eg :literal:`ground_station_position`), the partials of the Cartesian state with respect to the parameters must be directly returned by the function.

A last case arises when the estimatable parameter neither is  a rotation matrix property nor has a direct influence on the Cartesian state of the propagated body but is involved in one of the acceleration models. The function :literal:`createCartesianStatePartialsWrtParameter` does not return any partial in that case. Instead, the partial derivatives have to be implemented in existing acceleration partial classes.  The acceleration partials are implemented in the files of the directory:

   .../tudatBundle/tudat/Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/

When adding a new parameter, no new classes need to be defined. Instead, existing class(es) need to be update with teh new dependency of an acceleration on a parameter. For any estimatable parameter related to an acceleration model, one must define a function :literal:`getParameterPartialFunction` in the associated :literal:`AccelerationPartial` class, which is a method of its associated acceleration partial clas. This function provides the partial of the acceleration model with respect to the estimatable parameter. Looking at the parameter :literal:`radiation_pressure_coefficient` in particular, the acceleration partials are defined from the class :literal:`CannonBallRadiationPressurePartial` which is itself derived from the base class :literal:`AccelerationPartial`). This acceleration partials class is defined in the files :literal:`radiationPressureAccelerationPartial.h` and :literal:`radiationPressureAccelerationPartial.cpp`. Only the radiation pressure coefficient can be estimated for this acceleration model so that the function :literal:`getParameterPartialFunction` only contains a single specific case dedicated to this parameter.

.. code-block:: cpp

   //! Function for setting up and retrieving a function returning a partial w.r.t. a double parameter.
   std::pair< std::function< void( Eigen::MatrixXd& ) >, int > CannonBallRadiationPressurePartial::getParameterPartialFunction(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
   {
      std::function< void( Eigen::MatrixXd& ) > partialFunction;
      int numberOfRows = 0;

      // Check if parameter dependency exists.
      if( parameter->getParameterName( ).second.first == acceleratedBody_ )
      {
         switch( parameter->getParameterName( ).first )
         {
            // Set function returning partial w.r.t. radiation pressure coefficient.
            case estimatable_parameters::radiation_pressure_coefficient:

               partialFunction = std::bind( &CannonBallRadiationPressurePartial::wrtRadiationPressureCoefficient, this, std::placeholders::_1 );
               numberOfRows = 1;

               break;
            default:
               break;
          }
       }
       return std::make_pair( partialFunction, numberOfRows );
    }

The function :literal:`wrtRadiationPressureCoefficient` called in the piece of code above is also a method of the class :literal:`CannonBallRadiationPressurePartial` and is defined in the file :literal:`radiationPressureAccelerationPartial.h`. It directly returns the partial derivative of the radiation pressure acceleration model with respect to the radiation pressure coefficient, as follows:

.. code-block:: cpp

   {
   void wrtRadiationPressureCoefficient( Eigen::MatrixXd& partial )
        partial = computePartialOfCannonBallRadiationPressureAccelerationWrtRadiationPressureCoefficient(
                    radiationPressureFunction_( ), areaFunction_( ), acceleratedBodyMassFunction_( ),
                    ( sourceBodyState_( ) - acceleratedBodyState_( ) ).normalized( ) );
   }


