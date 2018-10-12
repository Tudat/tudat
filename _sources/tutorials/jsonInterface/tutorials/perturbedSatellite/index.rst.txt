.. _jsonInterface_tutorials_perturbedOrbitPropagation:

.. role:: jsontype
.. role:: jsonkey

Perturbed Earth-orbiting Satellite
==================================

This page describes how to set up the propagation of a perturbed Earth-orbiting satellite using the :literal:`json_interface`.

This example takes the code from :ref:`jsonInterface_tutorials_basicPropagation` as a base upon which we can start to build our new application.

The initial conditions and general simulation settings are thus the same as for the previous example. The only difference will be in the addition of perturbing accelerations. 

Adding Perturbations
~~~~~~~~~~~~~~~~~~~~

The addition of perturbations requires the specifications of a few more characteristics, both for the environment and the spacecraft. 

In this simple application, we would like to reproduce the results of the example in :ref:`walkthroughsPerturbedEarthOrbitingSatellite`. Hence, we would like to model the following accelerations:

   - **Spherical harmonics gravity** of the Earth up to degree and order 5.
   - **Third body perturbations** due to Sun, Moon, Mars and Venus gravity.
   - **Solar radiation pressure** (with occultations by the Earth).
   - **Aerodynamic drag** due to Earth's atmosphere.

Each of these accelerations require some small modifications to the previous code. In the text below, we will analyze what are the missing elements, and how they can be added to the JSON file.

Spherical Harmonics
*******************

To model spherical harmonics (SH), two steps are necessary. The first is the addition of the SH model to the body. This is done by specifying a new key for the object for which the SH model needs to be added. In case you want to add SH to the Earth, as we want to do for this example, you should add the following specification to the :jsonkey:`Earth` key:

   .. code-block:: json

      {
         "Earth": {
            "gravityField": {
               "type": "sphericalHarmonic",
               "model": "ggm02c"
            }
         }
      }

Here, we have specified the :literal:`ggm02c` model, which will be used to load the SH coefficients :math:`C` and :math:`S`. Note that you can find a list of the available spherical harmonics models at :ref:`jsonInterface_keys`, under Body, Gravity Field, Spherical Harmonic Gravity Field from built-in model. 

.. note:: The order in which the keys :jsonkey:`useDefaultSettings` and :jsonkey:`gravityField` are defined inside the :jsonkey:`bodies.Earth` objects is irrelevant: the default settings will always be loaded first, and then the settings specified in the input file (if any) will override the default ones.

Now, we also need to tell the simulation to model the gravitational acceleration. This is done by replacing the :literal:`"pointMassGravity"` value of the :jsonkey:`type` key under the accelerations of Asterix. In this case, we want to add :literal:`"sphericalHarmonicGravity"`. For this type of acceleration, we also need to specify the maximum degree and order to use. This is done by adding two extra keys, namely :jsonkey:`degree` and :jsonkey:`order`, after :jsonkey:`type`. Thus, for now, the JSON code specifying the accelerations acting on the body Asterix will look like this:

   .. code-block:: json

      {
         "Asterix": {
            "Earth": {
               "type": "sphericalHarmonicGravity",
               "maximumDegree": 5,
               "maximumOrder": 5
            }
         }
      }

Third Bodies
************

The addition of third bodies is very straight-forward. Firstly, we need to add these extra bodies to the simulation. Since we do not have any specific request for these bodies, we can simply load the default settings for each of the bodies we need, i.e., for Sun, Moon, Mars and Venus. Therefore, under the key :jsonkey:`bodies` we can add, taking the Sun as example:

   .. code-block:: json

      {
         "Sun": {
            "useDefaultSettings": true
         }
      }

Very simple! Then, we just need to specify the extra accelerations due to central gravitational field of each body. We can do this by adding to the list of accelerations of Asterix, the following:

   .. code-block:: json

      {
         "Asterix": {
            "Sun": [
               {
                 "type": "pointMassGravity"
               }
            ]
         }
      }

where the Sun is taken as example, once more.

Solar Radiation Pressure
************************

Inclusion of solar radiation pressure to the simulation requires three steps:

   - adding the Sun to the simulation (done in the previous subsection)
   - adding radiation pressure characteristics to Asterix
   - adding the acceleration to the list

The first step has already been done to model the third body gravity, thus we can directly focus on the last two. The radiation pressure interface must be inserted under the body settings of Asterix. We will need to specify three things: mass, reference area and the radiation parameters. Hence, under the key :jsonkey:`bodies`, we can include the following to the specifications of the spacecraft:

   .. code-block:: json

      {
         "Asterix": {
            "mass": 400,
            "referenceArea": 4,
            "radiationPressure": {
               "Sun": {
                  "radiationPressureCoefficient": 1.2,
                  "occultingBodies": [ "Earth" ]
               }
            }
         }
      }

For the radiation parameters, one element is compulsory: the radiation pressure coefficient (which in the example at hand is 1.2). The second element shown above, i.e., :jsonkey:`occultingBodies`, can be added if you want to also model the shadow cone due to some other bodies. In this case, we are modeling the shadow due to the Earth, which for a LEO spacecraft, can account for a large part of the orbit. 

The next step is to add the acceleration. This is very similar to what we have done for the previous accelerations. Thus, we will add to the list of accelerations of Asterix, due to the Sun, the following key and value pair:

   .. code-block:: json

      {
         "Asterix": {
            "Sun": [
               {
                 "type": "cannonBallRadiationPressure"
               }
            ]
         }
      }

Aerodynamics
************

Modeling of the aerodynamic acceleration is somewhat similar to the modeling of the solar radiation. Firstly, we need to specify the atmosphere of the planet we are orbiting. However, since we want to keep the default atmosphere of the Earth (the tabulated 1976 Standard Atmosphere Model), we do not need to specify this. 

The second step is the addition of the aerodynamic properties of the spacecraft. In this case, we want to model only a constant drag coefficient, equal to 1.2. We can do this by adding a new key to :jsonkey:`Asterix`, under :jsonkey:`bodies`:

   .. code-block:: json

      {
         "Asterix": {
            "aerodynamics": {
               "forceCoefficients": [ 1.2, 0, 0 ]
            }
         }
      }

Note that to properly model the aerodynamic acceleration, the value of mass and reference area also need to be specified, but these were already added in the previous subsection.

Finally, we need to include the acceleration in the :jsonkey:`propagators` settings. This is another simple task, and can be done by adding the following to the :jsonkey:`accelerations` of Asterix:

   .. code-block:: json

      {
         "Asterix": {
            "Earth": [
               {
                 "type": "aerodynamic"
               }
            ]
         }
      }

Running the Simulation
~~~~~~~~~~~~~~~~~~~~~~

The resulting input file, where all the settings have been combined, should look something like this:

.. literalinclude:: main.json
  :linenos:
  :language: json
  :caption: :class:`main.json`
  :name: main-json-perturbed

After running in Terminal the following command::

   json_inteface main.json

we will get a :class:`stateHistory.txt` file containing the results, and a :class:`fullSettings.json` file containing the list of the full settings, including the default ones.

Results
~~~~~~~

Finally, we can load the results and plot the difference between perturbed and unperturbed state, which is shown below.

.. image:: perturbedSatellite.png

Clearly, the addition of the perturbations causes quite a difference in the position of the spacecraft over one Julian day. 