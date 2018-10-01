
.. role:: jsontype
.. role:: jsonkey
.. role:: arrow

- :jsontype:`number[ 6 ]` :jsonkey:`initialStateInKeplerianElements` (mandatory). Initial Keplerian state.
- :jsontype:`number` :jsonkey:`epochOfInitialState` (mandatory). Initial epoch from which propagation of Kepler orbit is performed.
- :jsontype:`number` :jsonkey:`centralBodyGravitationalParameter` (mandatory). Gravitational parameter of central body about which the Kepler orbit is defined.
- :jsontype:`number` :jsonkey:`rootFinderAbsoluteTolerance` (optional). Convergence tolerance for root finder used to convert mean to eccentric anomaly. Default value: :literal:`200.0 * std::numeric_limits< double >::epsilon( )`.
- :jsontype:`number` :jsonkey:`rootFinderMaximumNumberOfIterations` (optional). Maximum iteration for root finder used to convert mean to eccentric anomaly. Default value: :literal:`1000`.
