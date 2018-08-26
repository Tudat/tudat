.. _walkthroughsTabulatedAtmosphere:

Tabulated Atmosphere Usage
====================================

In the example application folder you will find an example called :literal:`tabulatedAtmosphereUsage.cpp`, which can be run by selecting the executable :literal:`application_TabulatedAtmosphereUsage`. If you run this simple script as is, its output should look something like this: ::

   Example 1. ----------------------------------------------------------------------- 
   Altitude: 0 km. Density: 7.61788e-05 kg/m^3. Speed of sound: 205.817 m/s.
   Altitude: 50 km. Density: 7.61788e-05 kg/m^3. Speed of sound: 205.817 m/s.
   Altitude: 10000 km. Density: 8.9881e-18 kg/m^3. Speed of sound: 1596.22 m/s.
   Altitude: 100000 km. Density: 8.9881e-18 kg/m^3. Speed of sound: 1596.22 m/s.
   Warning in interpolator, requesting data point outside of boundaries, requested data at 0 but limit values are 50000 and 10000000, taking boundary value instead.
   Warning in interpolator, requesting data point outside of boundaries, requested data at 0 but limit values are 50000 and 10000000, taking boundary value instead.
   Warning in interpolator, requesting data point outside of boundaries, requested data at 0 but limit values are 50000 and 10000000, taking boundary value instead.
   Warning in interpolator, requesting data point outside of boundaries, requested data at 0 but limit values are 50000 and 10000000, taking boundary value instead.
   Warning in interpolator, requesting data point outside of boundaries, requested data at 100000000 but limit values are 50000 and 10000000, taking boundary value instead.
   Warning in interpolator, requesting data point outside of boundaries, requested data at 100000000 but limit values are 50000 and 10000000, taking boundary value instead.
   Warning in interpolator, requesting data point outside of boundaries, requested data at 100000000 but limit values are 50000 and 10000000, taking boundary value instead.
   Warning in interpolator, requesting data point outside of boundaries, requested data at 100000000 but limit values are 50000 and 10000000, taking boundary value instead.

   Example 2. ----------------------------------------------------------------------- 
   Longitude: 360 deg. Latitude: -180 deg. Altitude: 0 km. 
       Density: 0 kg/m^3. Speed of sound: 205.615 m/s.
   Longitude: 360 deg. Latitude: -180 deg. Altitude: 100 km. 
       Density: 0 kg/m^3. Speed of sound: nan m/s.
   Longitude: 360 deg. Latitude: 45 deg. Altitude: 0 km. 
       Density: 0 kg/m^3. Speed of sound: 205.615 m/s.
   Longitude: 360 deg. Latitude: 45 deg. Altitude: 100 km. 
       Density: 0 kg/m^3. Speed of sound: nan m/s.
   Longitude: 135 deg. Latitude: -180 deg. Altitude: 0 km. 
       Density: nan kg/m^3. Speed of sound: 205.615 m/s.
   Longitude: 135 deg. Latitude: -180 deg. Altitude: 100 km. 
       Density: nan kg/m^3. Speed of sound: nan m/s.
   Longitude: 135 deg. Latitude: 45 deg. Altitude: 0 km. 
       Density: 0.02 kg/m^3. Speed of sound: 205.615 m/s.
   Longitude: 135 deg. Latitude: 45 deg. Altitude: 100 km. 
       Density: 1.97143e-07 kg/m^3. Speed of sound: nan m/s.

.. note::
   Note that the location of the warnings may be different in your output. This is fine, since even on different runs on the same computer, these lines may appear in different places. The important thing is that they do show.

In **example 1**, the settings are very basic. The resulting atmosphere object describes how density, pressure, temperature, gas constant, specific heat ratio and molar mass, change over time for the Martian atmosphere, as a function of altitude. Here, the internal interpolator is instructed to use the boundary value, in case the altitude goes beyond its defined domain. As you can see from the results above, the atmosphere object does indeed output the boundary value when the altitude is both above and below the altitude bounds (which are given by the 2 results in the middle). 

A more interesting case is the **second example**. Here, the atmosphere is three-dimensional, depending on longitude, latitude and altitude. In this case :literal:`use_default_value` is used as boundary handling method and a pair of extrapolation values is defined for each combination of dependent and independent variables. For longitude and latitude, these values are quite straightforward:

   - **Longitude**: if out-of-bounds, give 0.0
   - **Latitude**: if out-of-bounds, give NaN

For **altitude**, on the other hand, a more complicated system is set up, where for each dependent value a specific pair is defined. What is interesting to notice is the order in which the out-of-range methods are processed. For instance, in the first test case (when i, j and k all equal 0), all three variables are outside their domain, however the default value of longitude (i.e., the first independent variable as defined by the user) is used as output. This would also be valid in the case that the third independent variable had was using :literal:`throw_exception_at_boundary` method. Thus, only during the last test case, would the simulation fail and give a runtime error (this is the only case where only the third variable is out-of-range).





