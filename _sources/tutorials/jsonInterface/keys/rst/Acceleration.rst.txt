
.. role:: jsontype
.. role:: jsonkey
.. role:: arrow

- :jsontype:`string` :jsonkey:`type` (mandatory). Possible values: :literal:`"undefined"`, :literal:`"pointMassGravity"`, :literal:`"sphericalHarmonicGravity"`, :literal:`"mutualSphericalHarmonicGravity"`, :literal:`"aerodynamic"`, :literal:`"cannonBallRadiationPressure"`, :literal:`"thrust"`, :literal:`"relativisticCorrection"`, :literal:`"empirical"`.

.. container:: toggle

	.. container:: header

		:arrow:`Spherical Harmonic Gravity Acceleration`

	.. include:: SphericalHarmonicGravityAcceleration.rst

.. container:: toggle

	.. container:: header

		:arrow:`Mutual Spherical Harmonic Gravity Acceleration`

	.. include:: MutualSphericalHarmonicGravityAcceleration.rst

.. container:: toggle

	.. container:: header

		:arrow:`Thrust Acceleration`

	.. include:: ThrustAcceleration.rst

.. container:: toggle

	.. container:: header

		:arrow:`Interpolated Thrust Acceleration`

	.. include:: InterpolatedThrustAcceleration.rst

.. container:: toggle

	.. container:: header

		:arrow:`Relativistic Correction Acceleration`

	.. include:: RelativisticCorrectionAcceleration.rst

.. container:: toggle

	.. container:: header

		:arrow:`Empirical Acceleration`

	.. include:: EmpiricalAcceleration.rst
