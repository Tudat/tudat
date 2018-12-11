
.. role:: jsontype
.. role:: jsonkey
.. role:: arrow

- :jsontype:`string` :jsonkey:`coefficientsType` (optional). Possible values: :literal:`"constant"`, :literal:`"tabulated"`. Default value: :literal:`"constant"`.
- :jsontype:`number` :jsonkey:`referenceArea` (optional). Reference area with which aerodynamic forces and moments are non-dimensionalized. If not specified, it will be retrieved from the body's :jsonkey:`referenceArea` key. Note that the reference area has to be defined at least once, either at the body level or in the aerodynamics settings. The definition of the reference area is optional only if it has already been defined in the body’s keys
- :jsontype:`number` :jsonkey:`referenceLength` (mandatory if moment coefficients provided). Reference length with which aerodynamic moments are non-dimensionalized.
- :jsontype:`number` :jsonkey:`lateralReferenceLength` (mandatory if moment coefficients provided). Lateral reference length with which aerodynamic moments are non-dimensionalized.
- :jsontype:`number` :jsonkey:`momentReferencePoint` (mandatory if moment coefficients provided). Point w.r.t. which the arm of the moment on a vehicle panel is determined.
- :jsontype:`string[ ]` :jsonkey:`independentVariableNames` (optional). Vector with identifiers of the physical meaning of each independent variable of the aerodynamic coefficients. Possible values: :literal:`machNumber`, :literal:`angleOfAttack`, :literal:`angleOfSideslip`, :literal:`altitude`, :literal:`controlSurfaceDeflection`, :literal:`undefined`.
- :jsontype:`boolean` :jsonkey:`areCoefficientsInAerodynamicFrame` (optional). Boolean to define whether the aerodynamic coefficients are defined in the aerodynamic frame (lift, drag, side force) or in the body frame (typically denoted as Cx, Cy, Cz). Default value: :literal:`FIXME`.
- :jsontype:`boolean` :jsonkey:`areCoefficientsInNegativeAxisDirection` (optional). Boolean to define whether the aerodynamic coefficients are positive along the positive axes of the body or aerodynamic frame. FIXME Note that for (lift, drag, side force), the coefficients are typically defined in negative direction. Default value: :literal:`true`.

.. container:: toggle

	.. container:: header

		:arrow:`Constant Coefficients Aerodynamics`

	.. include:: ConstantCoefficientsAerodynamics.rst

.. container:: toggle

	.. container:: header

		:arrow:`Tabulated Coefficients Aerodynamics`

	.. include:: TabulatedCoefficientsAerodynamics.rst
