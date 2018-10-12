.. _walkthroughsUseOfThrustUserDefinedThrustVector:

Use of Thrust: User-defined Thrust Vector
=========================================
In the previous tutorial, we discussed how to include thrust in a simulation, separately specifying the direction and magnitude of the thrust force. Alternatively, you may want to impose the thrust vector (direction and magnitude in some reference frame) as a function of time explicitly. In this tutorial, we will look at which Tudat interfaces you can use to accomplish this. The code for this tutorial is given on Github, and is also located in your tudat bundle at::

   tudatBundle/tudatExampleApplications/satellitePropagatorExamples/SatellitePropagatorExamples/thrustAccelerationFromFileExample.cpp

The main difference of this example w.r.t. the previous one lies in the following steps:

   1. Loading discrete thrust data from a file.
   2. Creating an interpolator that turns the discrete data set into a continuous function.
   3. Creating settings for the thrust acceleration using these data.

Loading the Data
~~~~~~~~~~~~~~~~
Reading data from a file can be done in many different ways. Since the manner in which a file is read depends on the structure of the particular file, there can be many different approaches (and file readers) for different applications. Here, we chose a very basic setup in which the thrust file is set up as follows::

   0       0 0 5
   6068    0 1 5
   6097    1.0 0 5
   6097.5  0.8 0 5
   6098    0.6 0.1 5
   6099    0.1 0.5 5
   12192   0.2 1.0 4.5
   18288   0.3 1.5 4.0
   243575  0.4 2.0 3.0
   3.999e6 1.0 1.0 2.0
   4e6     1.1 5.0 1.0

The first column of this file represents a time, the following three columns represent the x-, y- and z-components of the thrust force vector. To load these data into Tudat, we firstly need to specify the location of the file we wish to load. To do this, we retrieve the string containing the current file (e.g. the :literal:`.cpp` file with the source code of the example) as follows:

.. code-block:: cpp

   std::string cppFilePath( __FILE__ );

Subsequently, we strip of the name of the current file as follows:

.. code-block:: cpp

   std::string cppFolder = cppFilePath.substr( 0 , cppFilePath.find_last_of("/\\")+1 );

Do not worry if the details of these steps escape you. What is important is to know that the :literal:`cppFolder` now contains the location of the current directory. We now load the data into an :literal:`Eigen::MatrixXd` using the Tudat class :class:`FromFileDataMapSettings`, which can be used to extract data into a map format, which can later be used as input to the interpolation settings:

.. code-block:: cpp

   std::shared_ptr< FromFileDataMapSettings< Eigen::Vector3d > > thrustDataSettings =
           std::make_shared< FromFileDataMapSettings< Eigen::Vector3d > >( cppFolder + "testThrustValues.txt" );

Creating the Interpolator
~~~~~~~~~~~~~~~~~~~~~~~~~
We now have a :class:`FromFileDataMapSettings` object containing the time vs. thrust data from the file. To pass this information to the :class:`AccelerationSettings`, we need to turn this discrete data into a continuous function, for which we use an interpolator. Here, we choose to use a linear interpolator. For a list of the various other interpolation options, details of their implementation, and instructions on how to use/create them, go to :ref:`tudatFeaturesInterpolators`. For this example, we use the following code to create an interpolator of the thrust vector:

.. code-block:: cpp

   // Define interpolator settings.
   std::shared_ptr< InterpolatorSettings > thrustInterpolatorSettings =
           std::make_shared< InterpolatorSettings >( linear_interpolator );

   // Create data interpolation settings
   std::shared_ptr< DataInterpolationSettings< double, Eigen::Vector3d > > thrustDataInterpolatorSettings =
           std::make_shared< DataInterpolationSettings< double, Eigen::Vector3d > >(
               thrustDataSettings, thrustInterpolatorSettings );

The first part:

.. code-block:: cpp

   std::shared_ptr< InterpolatorSettings >
           thrustInterpolatorSettingsPointer = std::make_shared< InterpolatorSettings >( linear_interpolator );

creates an :class:`InterpolatorSettings` object that contains the settings for how to create the interpolator. For this application, this means specifying that the interpolator should be of the type :literal:`linear_interpolator`. Note that this setup is very similar to how an environment/acceleration/etc. model is set up. The interpolator is then created by calling:

.. code-block:: cpp

   // Creating settings for thrust force
   std::shared_ptr< DataInterpolationSettings< double, Eigen::Vector3d > > thrustDataInterpolatorSettings =
           std::make_shared< DataInterpolationSettings< double, Eigen::Vector3d > >(
               thrustDataSettings, thrustInterpolatorSettings );

The :class:`DataInterpolationSettings< double, Eigen::Vector3d >` object is created using the :literal:`thrustDataSettings` and :literal:`thrustInterpolatorSettings` defined above. The first parameter of the :class:`DataInterpolationSettings` denotes that the independent variable is a :literal:`double` (time) and the dependent variable is a :literal:`Eigen::Vector3d` (thrust). 

Creating the Thrust Acceleration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The creation of the thrust acceleration is done similarly as in the previous example, by creating an object of type :class:`ThrustAccelerationSettings`, as follows:

.. code-block:: cpp

   double constantSpecificImpulse = 3000.0;

    accelerationsOfVehicle[ "Vehicle" ].push_back( std::make_shared< ThrustAccelerationSettings >(
                                                       thrustDataInterpolatorSettings, constantSpecificImpulse,
                                                       lvlh_thrust_frame, "Earth" ) );

The input to the :class:`ThrustAccelerationSettings`, however, is different from that used in the previous example. In fact, we use a different constructor here, an example of constructor overloading. The input required to the constructor we use here is:

   - The interpolator used to compute the thrust force vector as a function of time.
   - The constant value of specific impulse (could also be a function depending on the current time).
   - The frame type in which the thrust vector is expressed.
   - The reference body for any frame transformation that may be required.

The last two argument define the frame orientation in which the thrust force produced by the :literal:`thrustDataInterpolatorSettings` is expressed. At present, there are two options:

   1. Inertial frame: if this is the case, there is no need to specify a reference body. The interpolated thrust is used directly in the equations of motion, without and transformation.
   2. Local-Vertical Local-Horizontal. This is a satellite-based frame in which the :math:`x`-axis is colinear and in the direction of the velocity vector (relative to the reference body). The :math:`z`-axis is perpendicular to the orbital plane (direction of cross-product of velocity with postion) and the :math:`y`-axis completes the right-handed system.

In this example, we use the second option, basing the thrust direction on the current Earth-centered position of the spacecraft.

The rest of the application, including the definition of the mass propagation, is set up analogously to the previous example, with a single addition: two dependent variables are saved during the propagation, namely the thrust acceleration, and the rotation matrix from LVLH to inertial frame. Note that when saving an acceleration, it is always saved as expressed in the inertial frame. We also save the rotation matrix here, to reconstruct the original thrust profile that we provided, checking the correct implementation.

Results
~~~~~~~
Below, we show the resulting orbit of the spacecraft w.r.t. the Earth. Clearly, the thrust force that we apply has a significant effect, changing the orbital plane and increasing the spacecraft's mean distance from the Earth.

We also show plots of the acceleration (in an inertial frame) and force (in the LVLH frame) due to the thrust. The thrust profile clearly shows the linearly interpolated behaviour from our input data. For the acceleration, the once-per-orbit signature of the transformation is clearly visible.

.. figure:: images/orbitThrustFromFile.png

.. figure:: images/accelerationThrustFromFile.png

The dependent variable history (accelerations) are obtained from the :literal:`getDependentVariableHistory` function inside the :class:`DynamicsSimulator` class. The resulting :literal:`std::map` can be saved as discussed in :ref:`tudatFeaturesInputOutput`. 