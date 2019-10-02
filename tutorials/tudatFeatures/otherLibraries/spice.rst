.. _tudatFeaturesSpice:

External Libraries: SPICE
=========================

The Spice toolkit consists of various data sets, designated as kernels, which are provided by the NASA's Navigation and Ancillary Information Facility (NAIF), at JPL. Those data include ephemerides of spacecraft, planets, asteroids, and comets, as well as some physical and dynamical characteristics of those bodies (e.g. size, shape, rotational state). The Spice library has been linked to Tudat, making it possible for the user to retrieve data from the Spice kernels. In particular, using ephemerides retrieved from Spice is an interesting solution to reduce the computational load with a limited loss of accuracy, instead of propagating the dynamics of all the bodies involved in the simulation. 

.. note:: 
    A detailed documentation of the Spice toolkit can be found here: https://naif.jpl.nasa.gov/naif/documentation.html  

Spice interface
~~~~~~~~~~~~~~~

The Spice interface allows for retrieving various data from the Spice library in Tudat. The available functions are the following ones:

   - :literal:`getBodyCartesianStateAtEpoch`
       Function returning the cartesian state of a body at a given epoch that must be specified.

   - :literal:`getBodyCartesianPositionAtEpoch`
       Function returning the cartesian position of a body at a given epoch that must be specified.

   - :literal:`computeRotationQuaternionBetweenFrames`
       Function returning the quaternion describing the rotation from two different frames.

   - :literal:`computeRotationMatrixDerivativeBetweenFrames`
       Function returning the derivative of the rotation matrix between two different frames.

   - :literal:`getAngularVelocityVectorOfFrameInOriginalFrame`
       Function returning the angular velocity of one given frame with respect to another one.

   - :literal:`computeRotationQuaternionAndRotationMatrixDerivativeBetweenFrames`
       Function returning both the quaternion and rotation matrix derivative between two different frames.

   - :literal:`getBodyProperties`
       Function returning one specific property of a celestial body. This property could be the gravitational parameter (:literal:`GM`), average radius (:literal:`RADII`), right ascension of the pole (:literal:`POLE_RA`), declination of the pole (:literal:`POLE_DEC`), prime meridian (:literal:`PM`), nutation and precession of the right ascension of the pole (:literal:`NUT_PREC_RA`), nutation and precession of the pole declination (:literal:`NUT_PREC_DEC`), nutation and precession of the prime meridian (:literal:`NUT_PREC_PM`), ...

   - :literal:`getBodyGravitationalParameter` 
       Function returning the gravitational parameter of a given celestial body.

   - :literal:`getAverageRadius`
       Function returning the average radius of a given celestial body.

Before retrieving a property from the Spice data, one might first want to ensure that this property is actually available in Spice. A specific function is available in Tudat to this end. Calling :literal:`checkBodyPropertyInKernelPool` with the name of the body and the name of the property you are interested in will return either :literal:`true` if the property is defined for this specific body, or :literal:`false` if it is not.

.. note:: 
   When using the Spice interface, the celestial body for which one wants to retrieve some data is specified by its name (so by a :literal:`std::string` object), in order to make it more readable by the user. However, using the name of a given body as its identifier is only possible when using one of the Spice interface functions. On the contrary, the Spice functions refer to different bodies via their so-called NAIF identifier (NAIF ID). The function :literal:`convertBodyNameToNaifId` which takes as input the name of a body and returns its corresponding NAIF ID (:literal:`int` object) can be used to make the conversion. A list of the NAIF IDs of the most commonly used celestial bodies can also be found here: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/naif_ids.html . NAIF IDs can typically appear in error messages returned by the SPICE interface, so knowing the name of the body corresponding to this NAIF identifier actually makes the Spice output much more readable.

An important point is still to be mentioned: before trying to use any data derived from Spice, the corresponding Spice kernels must be loaded in Tudat. Two different functions are available to do so:

   - :literal:`loadStandardSpiceKernels`
       The function loadStandardSpiceKernels is usually used to import the most commonly used Spice kernels. It includes the loading of the files containing the gravitational parameters of the natural bodies, as well as the files storing the orientation, shape and size of those bodies. The DE430 file (April 2013 version of the Developmental Ephemerides provided by JPL) is also to be loaded. However,it must be stressed out that the original :literal:`de430` file has been restricted to dates ranging from January, 1rst 1975 to Januray, 1rst 2025, for the sake of memory management. The smaller version of the de430 file (:literal:`de430_small`) is the one loaded by the :literal:`loadStandardSpiceKernels` function. Finally, a file containing the leap seconds elapsed from 1990 to 2016 is loaded too. This function can also take as input a :literal:`std::vector< std::string >` object containing a vector of filenames corresponding to the ephemerides kernel files, to be loaded instead of the standard :literal:`de430_small` one. 

   - :literal:`loadSpiceKernelInTudat`
       This function takes as input a :literal:`std::string` filename corresponding to the Spice kernel file to be loaded in Tudat.


Spice ephemerides
~~~~~~~~~~~~~~~~~

One of the main advantages of the use of Spice in Tudat is the use of Spice ephemerides. It indeed allows for directly retrieving the states of some of the celestial bodies from the kernels, without having to propagte their dynamics. This can save a significant fraction of the total computational time, while not deteriorating the accuracy too severely.

Two different classes have been implemented to generate ephemerides based on the Spice data.

.. class:: SpiceEphemeris

   Class for ephemerides based on Spice data. The state of the body is directly retrieved from the Spice kernels, without any need for propagating it.

Such a Spice ephemeris can be created from the following class of ephemeris settings (see :ref:`tudatFeaturesEnvironmentSetUp` for more details): 

.. class:: DirectSpiceEphemerisSettings

   Class of ephemeris settings to set up a ephemeris model directly retrieved from Spice, specifying the frame origin and orientation. Among the available options, it is also possible to choose to apply correcting for stellar and light-time abberations.


The use of a :literal:`SpiceEphemeris` model requires to retrieve Spice data at each time step of the propagation, while it is possible to save some computational time by limiting the number of calls to the Spice data. The :literal:`InterpolatedSpiceEphemerisSettings` class of ephemeris settings can be used to this end. It specifies the settings required to retrieve the body state from Spice kernels at regular time intervals, and then create a :literal:`TabulatedCartesianEphemeris` ephemeris model, by interpolating the scattered Spice data which have just been retrieved (with a 8th Lagrange interpolator by default). In that case, calls to the Spice data occur outside of the propagation loop, before creating the ephemeris model, so that it can save some computational time. However, the time interval at which the Spice data are collected before generating the :literal:`TabulatedCartesianEphemeris` must be wisely chosen, as too sparse Spice data would deteriorate the ephemeris accuracy.

Commonly encountered errors
~~~~~~~~~~~~~~~~~~~~~~~~~~~

When making use of the Spice interface, some errors regularly happen:

  - The Spice functions cannot deal with :literal:`TUDAT_NAN` inputs. This can occur in particular when using variable step-size integrators in combination with Spice data: if the minimum step-size is exceeded and a :literal:`TUDAT_NAN` is returned for the current time, it can then be passed as input when trying to retrieve the body state from Spice kernels. The following error is then returned:

   .. code-block:: cpp
 
      Toolkit version: N0065
 
      SPICE(SPKINSUFFDATA) --
 
      Insufficient ephemeris data has been loaded to compute the state of 1000093
      (TEMPEL 1) relative to 0 (SOLAR SYSTEM BARYCENTER) at the ephemeris epoch
      5879611 B.C. JUN 23 0-:0-:-23.64.
 
      A traceback follows.  The name of the highest level module is first.    
      spkezr_c --> SPKEZR --> SPKEZ --> SPKACS --> SPKAPS --> SPKLTC --> SPKGEO
   
     Oh, by the way:  The SPICELIB error handling actions are USER-TAILORABLE.  You
     can choose whether the Toolkit aborts or continues when errors occur, which
     error messages to output, and where to send the output.  Please read the ERROR
     "Required Reading" file, or see the routines ERRACT, ERRDEV, and ERRPRT.


  - It can also occur that one tries to retrieve non-existing data from Spice. Especially when retrieving a body ephemeris, the time interval over which the ephemeris data are available is limited and it might then be exceeded, leading to the following error message:

    .. code-block:: cpp
 
         Toolkit version: N0065
  
         SPICE(SPKINSUFFDATA) --
 
         Insufficient ephemeris data has been loaded to compute the state of 1000093
         (TEMPEL 1) relative to 0 (SOLAR SYSTEM BARYCENTER) at the ephemeris epoch 2013
         SEP 24 12:00:00.000.
   
         A traceback follows.  The name of the highest level module is first.
         spkezr_c --> SPKEZR --> SPKEZ --> SPKACS --> SPKAPS --> SPKLTC --> SPKGEO
 
         Oh, by the way:  The SPICELIB error handling actions are USER-TAILORABLE.  You
         can choose whether the Toolkit aborts or continues when errors occur, which
         error messages to output, and where to send the output.  Please read the ERROR
         "Required Reading" file, or see the routines ERRACT, ERRDEV, and ERRPRT.



    As mentioned previously, some ephemerides data files have been reduced to avoid loading too large files in Tudat. You might want to extend the time range of the ephemerides files so that it can cover the timescale of your propagation. Modifying the time interval (and more generally the different properties) of an ephemerides file is done via the SPKmerge program. This executable can be downloaded from https://naif.jpl.nasa.gov/naif/utilities.html (first select the platform and architecture corresponding to your own laptop, and then download the :literal:`spkmerge` executable from the list of available utilities). When running :literal:`SPKmerge`, it will ask for the path and name of an input file containing the required characteristics of the ephemerides file you want to create. The new ephemerides file is generated by combining different already-existing Spice kernels, for specific bodies and dates, all defined by the user in the input file. All already-existing Spice kernels imported in Tudat can be found in :literal:`/tudatBundle/tudat/Tudat/External/SpiceInterface/Kernels/` and other kernels can be downloaded from https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/ if necessary). To give an example, the :literal:`de430_small.inp` input file, from which the :literal:`de430_small.bsp` kernel has been generated, is written as follows.

   .. code-block:: cpp

       LEAPSECONDS_KERNEL   = naif0012.tls

       SPK_KERNEL           = de430_small.bsp
       BEGIN_TIME          = 1975 JAN 01
       END_TIME            = 2025 JAN 01
       SOURCE_SPK_KERNEL   = de430.bsp
       SOURCE_SPK_KERNEL   = mar097.bsp
       BODIES = 499
       SOURCE_SPK_KERNEL   = jup310.bsp
       BODIES = 599

   To give more details about the way such an input file is constructed.
      - :literal:`LEAPSECONDS_KERNEL` gives the name of the file containing the leap seconds data.

      - :literal:`SPK_KERNEL` defines the name of the ephemeris file this input file intends to generate.

      - :literal:`BEGIN_TIME` and :literal:`END_TIME` specify the time at which the ephemeris data should begin and end, respectively. This parameter is typically the one to be modified if the time range of the propagation exceeds the current ephemeris datafile.

      - There are several :literal:`SOURCE_SPK_KERNEL` entries here: the first one designates the name of the global ephemeris kernel to be imported. As it is not directly followed by any :literal:`BODIES` entry, it means that the ephemeris data for ALL bodies will be imported from this Spice kernel. The additional :literal:`SOURCE_SPK_KERNEL` entries define the names of other Spice kernels that are to be merged with the first one. However, they are here followed by a :literal:`BODIES` entry, defining for which body the ephemeris data is supposed to be imported. In that particular example, the full :literal:`de430.bsp` kernel is to be imported. The ephemeris data of Mars (599 is the NAIF id corresponding to Mars) from the :literal:`mar097.bsp` kernel and the ephemeris data of Jupiter (NAIF id 599) from the :literal:`jup310.bsp` kernel will then be merged to it.

   Modifying the entries of the input file from which the ephemeris file is generated (or even creating a new input file to generate an additional Spice kernel) can ensure that ephemeris data are available for all the bodies involved in the propagation, over the required timescale. 


