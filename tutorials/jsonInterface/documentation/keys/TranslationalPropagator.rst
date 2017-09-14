.. role:: arrow

- :literal:`string` :class:`type` (optional) Type of the translational propagator to be used. Possible values: :literal:`"cowell"`, :literal:`"encke"`, :literal:`"gaussKeplerian"`, :literal:`"gaussModifiedEquinoctial"`. Default value: :literal:`"cowell"`.
- :literal:`string[]` :class:`centralBodies` (mandatory) Names of the central bodies.
- :literal:`object[]` :class:`accelerations` (mandatory) Map in which each object contains a map of acceleration lists. The keys of the outer map are the names of the bodies undergoing the accelerations, while the keys of the inner maps are the names of the bodies exerting the accelerations. For instance, :literal:`accelerations/satellite/Earth` is read as: accelerations on satellite caused by Earth.
