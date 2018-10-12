
.. role:: jsontype
.. role:: jsonkey
.. role:: arrow

- :jsontype:`boolean` :jsonkey:`useStandardKernels` (mandatory). Whether the standard kernels should be preloaded.
- :jsontype:`path[ ]` :jsonkey:`alternativeKernels` (optional). Additional kernels to load when using standard kernels. Default value: :literal:`[ ]`.
- :jsontype:`path[ ]` :jsonkey:`kernels` (mandatory if :jsonkey:`useStandardKernels` set to :literal:`false`). Kernels to load when not using standard kernels.
- :jsontype:`boolean` :jsonkey:`preloadEphemeris` (optional). Whether the kernels should be preloaded from :jsonkey:`.initialEpoch` - :jsonkey:`interpolationOffsets[ 0 ]` to :jsonkey:`.finalEpoch` + :jsonkey:`interpolationOffsets[ 1 ]`. Default value: :literal:`true`.
- :jsontype:`number` :jsonkey:`interpolationStep` (optional). Step to use for interpolation of ephemeris. Default value: :literal:`300.0`.
- :jsontype:`number[ 2 ]` :jsonkey:`interpolationOffsets` (optional). Offsets to use for the interpolation. Default value: :literal:`[10*interpolationStep, 10*interpolationStep]`.
