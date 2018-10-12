.. _externalUtilityExamples:

Dynamic Memory: Usage Examples
==============================

This page provides some examples on how to use :class:`shared_ptr`, :class:`make_shared`, :class:`function` and :class:`bind`, which were introduced in :ref:`externalUtility`, at application level and in combination with Tudat. 

Note that the understanding of the two smart pointers is crucial for any Tudat user, whereas the documentation on :class:`function` and :class:`bind` should be consulted only when working on feature development. The example of :class:`bind` also includes the use of the RK4 numerical integrator from Tudat.

Shared Pointer Example
~~~~~~~~~~~~~~~~~~~~~~

The :class:`shared_ptr` is a template class that stores a smart pointer to a dynamically allocated object. Generally, it behaves much like a normal C++ pointer but should be used instead, because it is more safe. Use of this pointer ensures that the object pointed to is destroyed when the last :class:`shared_ptr` to the object is destroyed or reset. Thus, this feature provides automatic memory management and prevents memory leaks. Also, this eliminates the use of explicit delete statements.

This example considers a shared pointer to a :class:`TestBody` class. This "TestBody" class is essentially a stripped-down version of the full ``Body`` class used in Tudat. For this example it is assumed that the code for this class is in your application folder. The code of this small TestBody class is given at the end of this page.

The constructor of the TestBody class is:

   .. code-block:: cpp

      TestBody( const Eigen::Vector6d& aState, const double aTime = 0.0 );

where the input of the constructor is a state vector (in Cartesian elements) and the time, which has the default value 0.0. Assume that the state vector in Cartesian elements is known: initialState.

Include Statements
******************

For this example, the following include statement is required:

   .. code-block:: cpp

      #include <memory>

Generally, this include statement is required when using shared pointers in your code. Furthermore, the code of the "TestBody" class should be included, either above the main function or as a separate header and source file.

Defining the Type
*****************

The :class:`shared_ptr` class is a template class that stores a smart pointer to a dynamically allocated object. Generally, it is of the following form:

   .. code-block:: cpp

      std::shared_ptr< objectType >

where ``std::shared_ptr`` is used to indicate that it is a shared pointer, and ``objectType`` is the type of the object that it points to. The object type is placed in angular brackets, because it is a template parameter.

The shared pointer to an object of the TestBody class can now be declared as:

   .. code-block:: cpp

      std::shared_ptr< TestBody >

The above type can be roughly considered as the 'equivalent' to the following raw pointer-type:

   .. code-block:: cpp

      TestBody*

Typedef
*******

It can be convenient to create a typedef for the :class:`shared_ptr` to a TestBody class as follows:

   .. code-block:: cpp

      typedef std::shared_ptr< TestBody > TestBodyPointer;

Such a typedef can be placed in the header file of the class below the class itself, or in the main function. Within Tudat you will see many typedefs for shared pointers. Creating such a type definition makes it easier to read the code.

Creation
********

In our main function, we now want to create a shared pointer to an object of the TestBody class, we will name the pointer ``exampleTestBodyPointer``. We can use the following code to create the :class:`shared_ptr`:

   .. code-block:: cpp

      std::shared_ptr< TestBody > exampleTestBodyPointer( new TestBody( initialState ) );

where: 

   - ``std::shared_ptr< TestBody >`` is the type: a shared pointer to an object of type :class:`TestBody`.
   - ``exampleTestBodyPointer`` is the name of the object-pointer.
   - ``new`` is used to dynamically allocate the memory to the object.
   - ``TestBody( )`` is the call to the constructor of the :class:`TestBody` class.
   - ``initialState`` is the input to the TestBody class constructor.

Note that the default value for aTime (0.0) is used. Equivalently, we can use the above typedef instead:

   .. code-block:: cpp

      TestBodyPointer exampleTestBodyPointer( new TestBody( initialState ) );

The above line of code can be roughly considered as the :class:`shared_ptr` 'equivalent' to the following raw pointer code:

   .. code-block:: cpp

      TestBody* exampleTestBodyPointer = new TestBody( initialState );

Make Shared Example
~~~~~~~~~~~~~~~~~~~

Boost :class:`make_shared` is a factory function that creates an object of a given type and returns a :class:`shared_ptr` to it. This way the explicit use of ``new`` is avoided.

For this example, again consider the constructor of the TestBody class (see the :class:`shared_ptr` example), and assume that the input state is known: ``initialState``.

Include Statements
******************

For this example the following include statements are required:

   .. code-block:: cpp

      #include <memory>

Creation
********

In our main function, we again want to create a shared pointer to an object of the :class:`TestBody` class, which we will name ``newExampleTestBodyPointer``. Now we will create it using :class:`make_shared`:

   .. code-block:: cpp

      std::shared_ptr< TestBody > newExampleTestBodyPointer = std::make_shared< TestBody >( initialState );

where:

   - ``std::shared_ptr< TestBody >`` is the type: a shared pointer to an object of type :class:`TestBody`.
   - ``newExampleTestBodyPointer`` is the name of the object-pointer.
   - ``std::make_shared< TestBody >`` is the factory function that creates an object of type :class:`TestBody` and returns a shared pointer to it.
   - ``( initialState )`` is the input to the :class:`TestBody` class constructor.

The above line of code can be roughly considered as the :class:`shared_ptr` and :class:`make_shared` 'equivalent' to the following raw pointer code:

   .. code-block:: cpp

      TestBody* asterix = new TestBody( asterixInitialState );

Note that the :class:`make_shared` constructor cannot take more that 10 arguments. If you want to construct a shared pointer to an object that requires more than 10 arguments, use the following code instead:

   .. code-block:: cpp

      std::shared_ptr< objectType > newExampleTestBodyPointer( new objectType( input1, input2, ..., inputN ) ); // Where N > 10.

.. _externalUtilityExamplesFunction:

Function Objects: Usage Examples
================================

Function Example
~~~~~~~~~~~~~~~~

The :class:`function` function is used in Tudat instead of a normal C++ function pointer because it is an easy way to pass a function as an input to another object or function. Specifically, it allows the passing of free functions and class member functions with a single interface. The :class:`functio` type defines both its input types (and order) and its output type. Any function that fits this profile is accepted when it is passed.

For this example, we will show the use of :class:`function` with a free function. The example with a member function is shown in the :class:`bind` example below. We use the Himmelblau function as an example:

   .. code-block:: cpp

      f( x, y ) = ( x 2 + y - 11 )^2 + ( x + y 2 - 7 )^2

For which we have the following C++ free function:

   .. code-block:: cpp

      double himmelblau( double x, double y )
      {
          return ( ( x * x + y - 11.0 ) * ( x * x + y - 11.0 ) +
                    ( x + y * y - 7.0 ) * ( x + y * y - 7.0 ) );
      }

Include Statements
******************

For this example, the following include statement is required:

   .. code-block:: cpp

      #include <functional> 

Generally, when using :class:`function` in your code, this include statement is required.

Defining the Type
*****************

The :class:`function` type that corresponds to the Himmelblau function defined above is:

   .. code-block:: cpp

      std::function< double( double, double ) >;

where the parameters between the angular brackets are the template parameters. These parameters comply to the form:

   .. code-block:: cpp

      < output( input(s) ) >

where :literal:`output` is the output type and :literal:`input(s)` is the type of the input, this can be a list of input arguments. In case of no output, you should write :literal:`void` instead of :literal:`output`, whereas in case of no input, you can leave the :literal:`input(s)` field empty.

Creation
********

In our main function, we want to create a function pointer that complies to this form, which we name ``exampleFunction``. We can use the following code to create the function pointer:

   .. code-block:: cpp

      std::function< double( double, double ) > exampleFunction = &himmelblau;

where:

   - ``std::function< double( double, double ) >`` is the type of exampleFunction, which describes the form of the function (input: double, double; output: double).
   - ``exampleFunction`` is the name of the object.
   - ``&himmelblau`` is the function object that is assigned (note the use of the ampersand :literal:`&`).

Usage
*****

Later-on in the code the input argument to the function object ``exampleFunction`` can be provided using:

   .. code-block:: cpp

      double initialPoint = exampleFunction( 0.0, 0.0 ); // Computes f( 0.0, 0.0 )
      double randomPoint = exampleFunction( -5.0, 3.0 ); // Computes f( -5.0, 3.0 )

Here, the input values (0.0, 0.0 and -5.0, 3.0 respectively) are routed to the function ``himmelblau``, where the output is computed.

Results
*******

Using the code that was discussed in this example, you should be able to reproduce the following results: ::

    x           y           f(x,y)
    3.0         2.0         0.0
   -3.779310   -3.283189    4.24634e-010
   -2.0         3.0         16.0
    0.0         0.0         170.0

.. _exampleUtilityLambdaExpression:

Lambda Expressions
******************

There are cases in which you will want to create a function that always returns the same value, without any reference to a specific implementation. This is achieved by using a so-called lambda expression. The structure of this type of expression is as follows:

   .. code-block:: cpp

      [ captures ]( params ){ body }

   - :literal:`captures`: How the input values should be captured:

      - :literal:`=`: by copy
      - :literal:`&`: by reference
      - :literal:`_` (blank): no capture

   - :literal:`params`: List of the input values, as defined in the function declaration or in the :class:`function` template parameters.

   - :literal:`body`: The body of the function. In this case we will only focus on constant expressions, i.e., they directly return a value.

Applying the above definition to an example, a double-returning function that requires no input can be created as follows:

   .. code-block:: cpp

      std::function< double( ) > testFunction = [ = ]( ){ return 3; };

Now, if at any point in the code, the following is called:

   .. code-block:: cpp

      double testValue = testFunction( );

will result in testValue holding the value 3.0. The same trick can be used for any other type than a double. For instance:

   .. code-block:: cpp

      std::function< Eigen::Vector3d( ) > testVectorFunction = [ = ]( ){ return Eigen::Vector3d::UnitX( ); };

will create a function that always returns the :math:`[1,0,0]^\mathrm{T}` vector. Also, this interface can be used for functions requiring inputs, so that:

   .. code-block:: cpp

      std::function< double( const double ) > testFunction = [ = ]( const double ){ return 3; };
      double inputValue = 20.0;
      double testValue = testFunction( inputValue );

and the call to :literal:`testFunction` will ignore the input and again assign 3.0 to the :literal:`testValue` variable. Thus, the second set of brackets, i.e., the round brackets, are used to specify the type of input the function lambda expression support.

.. tip:: You can find a much more complete description of the construction of lambda expression, `here <https://en.cppreference.com/w/cpp/language/lambda>`_.

Bind Example
~~~~~~~~~~~~

The :class:`bind` function implements a simple and versatile template argument binding mechanism. This feature is used in Tudat because it is a good way to pass function objects. Its use is required for a unified free function and member function interface when using a class method, due to the way C++ calls function works.

When using :class:`bind` in your code, it can take many forms, varying where the class method is called (e.g. outside or inside the class), and varying the use of placeholders. This example considers the use of placeholders and calling from outside or inside a class.

The example described here is that of a Skydiver who starts its fall from an altitude of 25 kilometers. We are interested in its position and velocity after 50 seconds of free fall. The following assumptions are made:

   - The Earth and the Skydiver are modeled as point masses.
   - Variations in Earth's gravity with altitude are ignored.
   - All perturbations are neglected.
   - The following data is also needed for the calculations:
   - The initial altitude is 25 kilometers.
   - The initial velocity is 0 m/s.
   - The gravitational acceleration is constant at 9.81 m/s^2.
   - The Runge-Kutta 4 (RK4) integrator is used with a fixed step size of 1 second.

Include Statements
******************

For this example, the following include statements are required:

   .. code-block:: cpp

      #include <functional>
      #include <Eigen/Core>
      #include <Tudat/Mathematics/NumericalIntegrators/rungeKutta4Integrator.h>

When using :class:`bind` in your code, the first include statement is required. To make use of ``Eigen`` types, the ``Core`` of Eigen is included as well. The last include statement is used for the RK4 integrator. For more information on Eigen, see Eigen. For more information on numerical integrators, see :ref:`tudatFeaturesIntegrators`.

Calling the Class from Outside
******************************

First, an example is given when calling a class method from outside the class. For this we will define the ``Skydiver`` class, with a constructor and a function that computes the state derivative for its equations of motion:

   .. code-block:: cpp

      class Skydiver
      {
      public:

          //! Constructor of the Skydiver class.
          Skydiver( ){ }

          //! Compute the state derivative for given time and state.
          Eigen::Vector2d computeStateDerivative( const double time, const Eigen::Vector2d& state );

      protected:

      private:

      };

The function to compute the state derivative is:

   .. code-block:: cpp

      Eigen::Vector2d Skydiver::computeStateDerivative( const double time, const Eigen::Vector2d& state )
      {
         // Note that the time is not used in this function, but it is required as input for the integrator.
         Eigen::Vector2d stateDerivative = Eigen::Vector2d::Zero( ); // Declare and initialize to zero.
         stateDerivative( 0 ) = state( 1 );  // Velocity
         stateDerivative( 1 ) = -9.81;       // Acceleration

         return stateDerivative;
      }

The Skydiver class and function to compute the derivative can be placed above the main function (or alternatively in a separate header and source file). Inside the main function, first declare the constants that were given above. Declare these yourself:

   .. code-block:: cpp

      initialTime
      initialState // (which is of the type Eigen::Vector2d)
      timeStep
      endTime

Then, our Skydiver testDiver is constructed:

   .. code-block:: cpp

      Skydiver testDiver;

Next, an intermediate step is taken, before defining the integrator itself. The constructor of the RK4 integrator takes a state-derivative function, an initial time and an initial state. The state-derivative function is a :class:`function` type that corresponds to:

   .. code-block:: cpp

      std::function< Eigen::Vector2d ( const double, const Eigen::Vector2d ) >

The state derivative is the class method ``computeStateDerivative`` of the Skydiver class, because this is a class method, :class:`bind` is used to create a function pointer to the state-derivative function:

   .. code-block:: cpp

      std::function< Eigen::Vector2d ( const double, const Eigen::Vector2d ) > stateDerivativeFunction = 
          std::bind( &Skydiver::computeStateDerivative, &testDiver, std::placeholders::_1, std::placeholders::_2 );

where:

   - ``std::function< Eigen::Vector2d ( const double, const Eigen::Vector2d ) >`` is the type of stateDerivativeFunction, which describes the form of the function (input: double, Vector2d; output: Vector2d.
   - ``stateDerivativeFunction`` is the name of the object.
   - ``std::bind( ... )`` is the call to std::bind that binds the class method of the object testDiver to the object.
   - ``&Skydiver::computeStateDerivative`` is the member function of the Skydiver class that is bound (note the use of the ampersand &).
   - ``&testDiver`` is the object that the member function belongs to. Note the use of the ampersand &: this is used to pass a pointer to the object and not the object itself. Note that you can also pass a shared_ptr.
   - ``std::placeholders::_1`` and ``std::placeholders::_2`` specify placeholders. This is used to route the input arguments such that the input of the function can be provided later-on. In this case it specifies that the argument list for the computeStateDerivative function takes two arguments (time, state), and the order in which these arguments are taken.

This ``std::function`` can then be passed to the constructor of the integrator, the integrator is created as follows:

   .. code-block:: cpp

      tudat::mathematics::numerical_integrators::RungeKutta4IntegratorXd integrator(
          stateDerivativeFunction, initialTime, initialState );

where:

   - ``RungeKuttaIntegratorXd`` is the type of the object created, which is the default ``RungeKutta4Integrator``.
   - ``integrator`` is the name of the object.
   - ``stateDerivativeFunction``, ``initialTime`` and ``initialState`` are the input arguments of the constructor of the RK4 integrator.

Note that the previous two code fragments can also be coded as a single statement, which is more complex:

   .. code-block:: cpp

      tudat::mathematics::numerical_integrators::RungeKutta4IntegratorXd integrator(
          std::bind( &Skydiver::computeStateDerivative, &testDiver, std::placeholders::_1, std::placeholders::_2 ),
          initialTime, initialState );

The end state can now be computed using the class method integrateTo( ) of the integrator:

   .. code-block:: cpp

      const Eigen::Vector2d endState = integrator.integrateTo( endTime, timeStep );

Calling the Class from Inside
*****************************

Second, an example is given when calling a class method from inside the class. For this we will modify the Skydiver class by adding a second class method:

   .. code-block:: cpp

      class Skydiver
      {
      public:

          //! Constructor of the Skydiver class.
          Skydiver( ){ }

          //! Compute the state derivative for given time and state.
          Eigen::Vector2d computeStateDerivative( const double time, const Eigen::Vector2d& state );

          //! Compute the final state for given initial conditions, final time and time step-size.
          Eigen::Vector2d computeFinalState( const double initialTime, const Eigen::Vector2d& initialState,
                                             const double endTime, const double timeStep );

      protected: 

      private:

      };

The function to compute the state derivative is the same as above, and the function to compute the final state is:

   .. code-block:: cpp

      Eigen::Vector2d Skydiver::computeFinalState( const double initialTime, const Eigen::Vector2d& initialState,
                                               const double endTime, const double timeStep )
      {
          tudat::mathematics::numerical_integrators::RungeKutta4IntegratorXd integrator(
                      std::bind( &Skydiver::computeStateDerivative, this, 
                      std::placeholders::_1, std::placeholders::_2 ),
                      initialTime, initialState );

          return integrator.integrateTo( endTime, timeStep );
      }

Now, compare this function to the code of the previous example. The integrator is now used inside the class method ``computeFinalState``. This method calls the state derivative function from within the class. Do you notice the difference?

Focus on the ``std::bind`` call. In the previous example, the object that the method belonged to was ``testDiver``. Now, the class method is called from within the class and the object that the method belongs to is this. Did you spot it? Great! That is the difference between calling from outside or inside a class.

As in the previous example, declare the constants that were given in the main function. In the same manner as before, also the ``testDiver`` is declared. Now the end state can be computed using:

   .. code-block:: cpp

      const Eigen::Vector2d endState = testDiver.computeFinalState( initialTime, initialState, endTime, timeStep );

Results
*******

Using the code that was discussed in this example, you should be able to reproduce the following results: ::

      Time [s]        Position [m]        Velocity [m/s]
      10              24509.5            -98.1
      20              23038.0            -196.2
      50              12737.5            -490.5

Note that the following explanation about function calls in C++ is a simplification of reality, and for a rigorous explanation refer to your local computer scientist. When a new instance of a class is created, this does not mean that all function definitions etc. are copied. There is only one copy of a function in memory, regardless of the amount of class instances. When a class method is called, the function is called with a hidden pointer to the class instance; this way the method knows about the class variables. When a class method is called via a function pointer, then the pointer to the class instance needs to be passed explicitly by the called. This is where :class:`bind` comes into the picture. The :class:`bind` function wraps a function into another object, and stores the passed arguments (which in the examples are respectively ``testDiver`` and this). When the function wrapper is eventually called, the function wrapper will call the wrapped function with both the arguments passed at its construction, as those passed to the call. The arguments passed to the function wrapper call are indicated with the ``std::placeholders::_1``, respectively ``std::placeholders::_2`` in the :class:`bind` call.

.. class:: TestBody

   .. code-block:: cpp

      //! The test body class.
      /*!
       * This class serves as an example of how a data repository can be constructed that stores state
       * and time information, which can be used in conjunction with acceleration models, the Cartesian
       * state derivative model, and the composite state derivative model. It should be noted that this
       * class should NOT be used "as is", without consideration for the application at hand. Classes
       * such as this are application-specific, hence unavailable through the Tudat libraries.
       */
      class TestBody
      {
      public:

          // Constructor that takes an input state and time.
          // The default value of the input time is 0.0.
          TestBody( const Eigen::Vector6d& aState, const double aTime = 0.0 )
              : currentState( aState ),
                currentPosition( aState.segment( 0, 3 ) ),
                currentTime( aTime )
          { }

          // Public member function: set the current time and state.
          // This sets the current time, state and position, that are stored internally.
          void setCurrentTimeAndState( const double aTime,
                                       const Eigen::Vector6d& aState )
         {
              currentTime = aTime;
              currentState = aState;
              currentPosition = aState.segment( 0, 3 );
          }

          // Public member function: get the current position.
          // This returns the internally stored current position vector.
          Eigen::Vector3d getCurrentPosition( ) { return currentPosition; }

      protected:

      private:

          // Private member: the current state.
          Eigen::Vector6d currentState;

          // Private member: the current position.
          Eigen::Vector3d currentPosition;

          // Private member: the current time.
          double currentTime;

      };