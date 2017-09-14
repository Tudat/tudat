.. role:: arrow

- :literal:`bool` :class:`useStandardKernels` (mandatory) Whether the standard kernels should be preloaded.
- :literal:`path[]` :class:`alternativeKernels` (optional) Additional kernels to load when using standard kernels. Default value: :literal:`[]`.
- :literal:`path[]` :class:`kernels` (mandatory if :literal:`useStandardKernels == true`) Kernels to load when not using standard kernels.
- :literal:`bool` :class:`preloadKernels` (optional) Whether the kernels should be preloaded from :literal:`~/initialEpoch - interpolationOffsets[0]` to :literal:`~/finalEpoch + interpolationOffsets[1]`. Default value: :literal:`true`.
- :literal:`numeric` :class:`interpolationStep` (optional) Step to use for interpolation of ephemeris. Default value: :literal:`300.0`.
- :literal:`numeric[2]` :class:`interpolationOffsets` (optional) Offsets to use for the interpolation. Default value: :literal:`[10*interpolationStep, 10*interpolationStep]`.
