.. _tudatFeaturesTimeStateTemplates:

Time and State Templates
===========================
Various classes like :class:`IntegratorSettings` and :class:`TranslationalStatePropagatorSettings` require inputs called :literal:`StateScalarType` or :literal:`TimeType`. These parameters are called template parameters in C++, and they are used in classes and functions to declare a special kind of parameter type (for more information on templates in C++, go `here <http://www.cplusplus.com/doc/oldtutorial/templates/>`_). The templates can be used to define the accuracy of the state and/or the time used in the class, which can have a large influence on the total accuracy of the simulation, as can be seen in the figure at the end of this page. 

:literal:`StateScalarType` is used in classes which contain variables that are used to define the state of the object. This template has a default value of :literal:`double`. The other option for :literal:`StateScalarType` is :literal:`long double`, which, depending on your computer, has a higher accuracy then a :literal:`double` variable. When one of the options is chosen for a specific class, all variables that are of the :literal:`StateScalarType` in that class will have an accuracy which corresponds to the chosen option. 

The :literal:`TimeType` has slightly different options. The default option is the same as for the :literal:`StateScalarType`, i.e. :literal:`double`. However, the second option is a class called :class:`Time`. This class has a special way of calculating time, because of the reduction of the quality of time variables over longer periods of time. For example, when :literal:`double` or :literal:`long double` variables are used to determine the time, the accuracy of these variables will be in the order of (for a time period of 3 years) :math:`10^{-8}` or :math:`10^{-11}`, respectively, which is insufficient for various applications. The :class:`Time` option uses an :literal:`int` to keep track of the number of hours that have passed, and a :literal:`long double` to represent the number of seconds into the present hour. This provides a resulution of < 1 femtosecond, over a range of 2147483647 hours (about 300,000 years), which is more than sufficient for practical applications.

The most important cases where these types can be chosen are:

- In the :class:`integratorSettings` class, where only the :literal:`TimeType` needs to be chosen. If the default value needs to be chosen, :literal:`<>` can be left empty.

.. code-block:: cpp

	IntegratorSettings<TimeType>()

- In the :class:`TranslationalStatePropagatorSetting` class, where only the :literal:`StateScalarType` needs to be chosen. If the default value needs to be chosen, :literal:`<>` can be left empty. 

.. code-block:: cpp

	TranslationalStatePropagatorSettings<StateScalarType>()

- In the :class:`SingleArcDynamicsSimulator` class, where both the :literal:`StateScalarType` and :literal:`TimeType` needs to be chosen. If the default values need to be chosen, :literal:`<>` can be left empty.

.. code-block:: cpp

	SingleArcDynamicsSimulator<StateScalarType, TimeType>

There are other classes where these types need to be used as an input, but these are the most commonly used ones.

.. figure:: images/roundingErrorFigure.png   
   :align: center
