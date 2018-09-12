The Environment During Propagation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Each :class:`Body` object and its constituent members is updated to the current state and time automatically during the numerical propagation. We stress that only those models that are relevant for a given propagation are updated every time step (this is handled automatically, without user intervention). Most time-dependent properties of the body are set in the environment models themselves. However, a number are updated and stored directly in the :class:`Body` object. These are:

    - The current translational state of the body
    - The current orientation of the body (and its time derivative)
    - The current mass of the body

.. note:: As a user, you will typically not access these variables directly.