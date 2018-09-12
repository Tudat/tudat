.. _tudatFeaturesProbabilityDistributions:

Probability Distributions
=========================
In Tudat, we have a number of options related to the use of probability distributions, as well as the generation of random variables from these probability distributions. Many of these distributions are taken from the Boost libraries, some we have implemented ourselves.

Implementation of a Probability Distribution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To represent a probability distribution with a given probability distribution function (PDF) and associated cumulative distribution function (CDF), we have defined the :class:`ContinuousProbabilityDistribution` base class. From this class, we have two 'branched' of derived classes:

    - One set for which we can compute the inverse CDF (with derived class :class:`InvertibleContinuousProbabilityDistribution`).
    - One set for which we cannot compute the inverse CDF.

Computing the inverse CDF is needed to generate random variables (discussed below), so that random variables can only be drawn from those distributions of the first category. In most cases, a probability distribution of one of the following types will be sufficient, each of which is implemented in boost:

    - Uniform distribution
    - Normal (Gaussian) distribution
    - Exponential distribution
    - Gamma distribution
    - Lognormal distribution
    - Beta distribution

Each of these distributions can be create through the same Tudat interface. Below we given an example of creating a uniform distribution with bounds at -1 and 2.5:

.. code-block:: cpp

    using namespace statistics;

    // Define properties of distribution
    std::vector< double > parameters;
    parameters.push_back( -1.0 );
    parameters.push_back( 2.5 );

    // Create distribution
    boost::shared_ptr< InvertibleContinuousProbabilityDistribution< double > > probabilityDistribution =
            createBoostRandomVariable( uniform_boost_distribution, parameters );

where the ``probabilityDistribution`` object is created by the ``createBoostRandomVariable`` function, which takes two inputs:

    - The type of probability distribution variable, defined by a Tudat-defined enumeration (all the options are listed below).
    - A list of parameters that define the properties of the distribution.

The :class:`InvertibleContinuousProbabilityDistribution` probability distribution comes with three functions that you can used:

    - The :literal:`evaluatePdf` function, giving the PDF at a given point in sample space.
    - The :literal:`evaluateCdf` function, giving the CDF at a given point in sample space.
    - The :literal:`evaluateInverseCdf` function, giving the point in sample space for a given CDF values (which must be between 0 and 1). Note that this is also referred to as a quantile function.

For instance, we can perform the following:

.. code-block:: cpp

        double localPdf = probabilityDistribution->evaluatePdf( 1.6 );
        double localCdf = probabilityDistribution->evaluateCdf( 1.6 );
        double samplePoint = probabilityDistribution->evaluateInverseCdf( 0.8 );

To compute the three quantities listed above.

Creating Probability Distributions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To create one of the 6 distributions listed above, for each of which the invertible CDF can be computed, call the ``createBoostRandomVariable`` with the following input:

   - **Uniform distribution**

      - First argument: ``uniform_boost_distribution``
      - Second argument: vector containing (in order):

         1) Lower bound of interval for distribution 
         2) Upper bound of interval for distribution

   - **Normal distribution**

      - First argument: ``normal_boost_distribution``
      - Second argument: vector containing (in order):

         1) Mean (:math:`\mu`) of distribution 
         2) Standard deviation (:math:`\sigma`) of distribution

   - **Exponential distribution**

      - First argument: ``exponential_boost_distribution``
      - Second argument: vector containing (in order):

         1) :math:`\lambda` parameter of exponential distribution

   - **Gamma distribution**

      - First argument: ``gamma_boost_distribution``
      - Second argument: vector containing (in order):

         1) shape (:math:`k`) parameter of distribution 
         2) scale (:math:`\theta`) parameter of distribution

   - **Lognormal distribution**

      - First argument: ``lognormal_boost_distribution``
      - Second argument: vector containing (in order):

         1) location (:math:`\mu`) parameter of distribution 
         2) scale (:math:`\sigma`) parameter of distribution

   - **Beta distribution**

      - First argument: ``beta_boost_distribution``
      - Second argument: vector containing (in order):

         1) :math:`\alpha` parameter of distribution 
         2) :math:`\beta` parameter of distribution

A number of multivariate distributions are also available in Tudat. These are:

   - Multivariate Gaussian distribution.
   - Gaussian Cupola distrbution.

For more information on these distributions, you are referred to the in-code Doxygen documentation. Note that only a pdf can be evaluated for these distributions. Also, a Kernel density distribution for multi-variate data is available, see the in-code Doxygen documentation. Note that only a pdf and cdf (no inverse cdf) can be evaluated for this distribution.

Generation of Random Numbers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For each of the invertible random variables described above, you can easily create a random number generator which generates variables according to that distribution. For the 6 boost distributions listed above, you can use the same input (type and parameters) as listed above. We provide two interfaces for random variable generation:

    - A class :class:`ContinuousRandomVariableGenerator`, with a ``getRandomVariableValue`` that produces random numbers. 
    - A boost function ``boost::function< double( ) >`` which produces random numbers every time it is called.

We illustrate these two options with some examples below.

.. code-block:: cpp

    using namespace statistics;

    // Define properties of distribution
    std::vector< double > parameters;
    parameters.push_back( 1.0 );
    parameters.push_back( 3.2 );

    // Create distribution object
    double distributionSeed = 42.0;
    boost::shared_ptr< RandomVariableGenerator< double > > randomNumberGenerator = createBoostContinuousRandomVariableGenerator(
            uniform_boost_distribution, parameters, distributionSeed );

    // Create distrubution function
    distributionSeed = 43.0;
    boost::function< double( ) > randomNumberFunction = createBoostContinuousRandomVariableGeneratorFunction(
            uniform_boost_distribution, parameters, distributionSeed );

    // Generate random variables
    std::vector< double > randomVariables1;
    std::vector< double > randomVariables2; 
    for( unsigned int i = 0; i < 1E6; i++ )
    {
        randomVariables1.push_back( randomNumberGenerator->getRandomVariableValue( ) );
        randomVariables2.push_back( randomNumberFunction( ) );
    }

The above creates two vectors give a million doubles, distributed according to a normal distribution with mean of 1 and standard deviation of 3.2.

.. note:: When creating a random variable generator you must also provide a 'seed' value, which initializes the generator. Fixing this value for a given means that the exact same set of random numbers will be generated everytime you run the program.

