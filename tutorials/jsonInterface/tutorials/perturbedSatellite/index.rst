.. _jsonInterface_tutorials_perturbedOrbitPropagation:

.. role:: jsontype
.. role:: jsonkey

Perturbed Earth-orbiting Satellite
==================================

This page describes how to set up the propagation of a perturbed Earth-orbiting satellite can be set up using the :literal:`json_interface`.

This example takes the code from :ref:`jsonInterface_tutorials_basicPropagation` as a starting point, upon which we can start to build our new application.

The initial conditions and general simulation settings are thus the same as for the previous example. The only difference will be in the addition of perturbing accelerations. 

Adding Perturbations
~~~~~~~~~~~~~~~~~~~~

The addition of perturbations requires the specifications of a few more characteristics, both for the environment and the spacecraft. 

In this simple application, we would like to reproduce the results of the example in :ref:`walkthroughsPerturbedEarthOrbitingSatellite`. Hence, we would like to model the following accelerations:

   - **Spherical harmonics gravity** of the Earth up to degree and order 5.
   - **Third body perturbations** due to Sun, Moon, Mars and Venus.
   - **Solar radiation pressure**
   - **Aerodynamics** due to Earth's atmosphere

Each of these accelerations require some small modifications to the previous code. In the text below, we will analyze what are the missing elements, and hwo they can be added to the JSON file.

Spherical Harmonics
*******************

Specification of the gravity field of a body is done by adding the following key to the JSON object:

   .. code-block:: json

      {
         "Earth": {
            "gravityField": {
               "type": "sphericalHarmonic",
               "model": "ggm02c"
            },
         }
      }

Here, we have specified to load the :literal:`ggm02c` spherical harmonics model to the Earth.

Third Bodies
************

Solar Radiation Pressure
************************

Aerodynamics
************

Running the Simulation
~~~~~~~~~~~~~~~~~~~~~~

The resulting input file, containing all the settings, will be:

.. literalinclude:: main.json
  :linenos:
  :language: json
  :caption: :class:`main.json`
  :name: main-json-perturbed

After running in Terminal the following command::

   json_inteface main.json

we will get a :class:`stateHistory.txt` file containing the results, and a :class:`fullSettings.json` file containing default settings and with some keys moved to a different place (e.g. :jsonkey:`bodies.Asterix.initialState` to :jsonkey:`propagators[0].initialStates`):

.. literalinclude:: fullSettings.json
  :linenos:
  :language: json
  :caption: :class:`fullSettings.json`
  :name: fullSettings-json

Results
~~~~~~~

Finally, we can load the results and plot the difference between perturbed and unperturbed state, which is shown below.

.. image:: perturbedSatellite.png