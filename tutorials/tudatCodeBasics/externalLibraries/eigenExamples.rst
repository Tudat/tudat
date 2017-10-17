.. _externalEigenExamples:

Eigen: Usage Examples
=====================
This page provides an example on how to use Eigen vector objects, which were introduced on the Eigen page, at application level and in the Tudat library.

**Example description**
    The example described on this page is that of a Moon orbiting spacecraft. We are interested in the course of the gravitational acceleration that the spacecraft is subject to, and the altitude of the spacecraft, during a single orbit. The following assumptions are made:

    - The Moon and the spacecraft are modeled as point masses
    - All perturbations are neglected
    - The following data is also needed for the calculations:
    - The interval of the true anomaly [0,360] degrees is split up using a 1 degree step-size
    - The gravitational parameter of the Moon is 4.9e12 [m^3 s^-2]
    - The radius of the Moon is 1737.1 [km]
    - The orbital elements of the spacecraft are:
    - Semi-major axis = 3000 [km]
    - Eccentricity = 0.2 [-]
    - The inclination, argument of periapsis, and longitude of the ascending node are 0.0 [rad]

**Include statements**
    For this example the following include statements are required:

    .. code-block:: cpp

        #include <Eigen/Core>
        #include <Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
        #include <Tudat/Mathematics/BasicMathematics/mathematicalConstants.h>
        #include <Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h>

    Generally, when using Eigen, the core should be included. To make use of the Vector6d typedef, the Tudat file with linear algebra types (last include statement) also needs to be included. The other include statements are used for:

    - The conversion from Keplerian to Cartesian elements.
    - The mathematical constant "PI".

**Step 1. Declare and initialize constants**
    For this example, two constants are required, these are declared and initialized:

    .. code-block:: cpp

        const double moonGravitationalParameter = 4.9e12;
        const double moonRadius = 1737.1e3;

**Step 2. True anomaly range: VectorXd and LinSpaced**
    Next, the range of true anomaly is created using the LinSpaced method of the Eigen type VectorXd. This is done in the following way:

    .. code-block:: cpp

        Eigen::VectorXd trueAnomalyRange = Eigen::VectorXd::LinSpaced(
                    361, 0, 2*tudat::basic_mathematics::mathematical_constants::PI );

    This statement creates a VectorXd of size 361 on the interval [0,2pi]. We use the mathematical constant "PI" from the Tudat library (see third include statement). Note that these calculations are always in **radians**, where the assumptions above are made in degrees for convenience. We need to access objects in the Eigen namespace via:

    .. code-block:: cpp

        Eigen::

  The objects in the mathematical constants namespace are accessed via:

    .. code-block:: cpp

        tudat::basic_mathematics::mathematical_constants::

**Step 3. State vector with orbital elements of the spacecraft: Vector6d**
    Next, the Keplerian-element state vector of the spacecraft is created. This is done using the Eigen type Vector6d that is in the Tudat library. We start by declaring the vector and we use the method Zero( ) to initialize all elements to 0.0:

    .. code-block:: cpp

        tudat::basic_mathematics::Vector6d keplerianElements = tudat::basic_mathematics::Vector6d::Zero( );

    Now, we set the value of the semi-major axis and the eccentricity:

    .. code-block:: cpp

        keplerianElements( 0 ) = 3.0e6; // Semi-major axis 3000 [km].
        keplerianElements( 1 ) = 0.2; // Eccentricity 0.2 [-].

    Note that the index of the first element is 0, as opposed to the convention in mathematics that the first index is 1.

**Step 4. Creating the output vectors: VectorXd**
    We need to iterate over 361 calculations ( i.e. the interval of the true anomaly), but before we do this we need to allocate the memory for the output vectors in which the position and acceleration at each true anomaly value is to be stored. This is done using the VectorXd type:

    .. code-block:: cpp

        Eigen::VectorXd position = Eigen::VectorXd::Zero( trueAnomalyRange.rows( ) );
        Eigen::VectorXd acceleration = Eigen::VectorXd::Zero( trueAnomalyRange.rows( ) );

    Note that the Eigen::VectorXd is used here for illustrative purposes, as this type will in most cases be handled by a std::vector container, since no linear algebra operations will be performed on it. Note the following details:

    - The method Zero( size ) is used to initialize all elements to 0.0. The size needs to be specified, because it is dynamically allocated.
    - The method rows() returns the number of rows of the true anomaly range vector. This method is used to specify the size of the new vector, thus creating a vector of the same size.

**Step 5. Iterate over the conversion to Cartesian elements, and the acceleration computation**
    Now, we need to iterate over a number of computations on the range of the true anomaly. This is done using a for-loop:

    .. code-block:: cpp

        for( int i = 0; i < trueAnomalyRange.rows(); i++ )
        {
            // Set the value of the true anomaly.
            keplerianElements( 5 ) = trueAnomalyRange( i );

            // Convert from Keplerian to Cartesian elements.
            tudat::basic_mathematics::Vector6d cartesianElements =
                 tudat::orbital_element_conversions::convertKeplerianToCartesianElements(
                      keplerianElements, moonGravitationalParameter );

            // Extract the position vector.
            Eigen::Vector3d positionVector = cartesianElements.segment( 0, 3 );

            // Compute the position.
            position( i ) = positionVector.norm();

            // Compute the absolute acceleration using the simplified equation a = mu / r^2.
            acceleration( i ) = moonGravitationalParameter / ( position( i ) * position( i ) );
        }

    Note the following details:

    - The expression trueAnomalyRange.rows() is used to define the end condition of the for-loop.
    - The first three elements of the Cartesian state vector represent the position. The position is extracted using the method segment( first, size ). This returns a segment of the vector. In this case first=0 and size=3, so we start at element 0 and take the first 3 elements.
    - The norm of the position is computed using the method norm().

**Step 6. Compute the altitude: VectorXd**
    The altitude can be computed using:

    .. code-block:: cpp

        Eigen::VectorXd altitude = position.array() - moonRadius;

    Here, the array() method is used such that element-wise operations are possible.

**Step 7. Extracting information**
    At this point, we have calculated all information that we are interested in. The information can be extracted by writing it to a file, or printing it to the screen. It is possible, for example, to display the minimum or maximum value of a vector using the methods minCoeff() and maxCoeff(), respectively:

    .. code-block:: cpp

        position.minCoeff()
        position.maxCoeff()

**Results**
    Using the code that was discussed in this tutorial, you should be able to reproduce the following results:

    +-------------------+-------------------+-------------------+------------------------+
    |                   | **Position [km]** | **Altitude [km]** | **Acceleration [m/s]** |
    +-------------------+-------------------+-------------------+------------------------+
    | **Minimum value** | 2400              | 662.9             | 0.378086               |
    +-------------------+-------------------+-------------------+------------------------+
    | **Maximum value** | 3600              | 1882.9            | 0.850694               |
    +-------------------+-------------------+-------------------+------------------------+


