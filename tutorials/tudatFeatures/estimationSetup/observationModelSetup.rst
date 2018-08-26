.. _observationModelSetup:

Setting up Observation Models
=============================

.. _observationTypes:

Observation Types
~~~~~~~~~~~~~~~~~

Tudat supports a diverse set of observation types, which are defined in the :literal:`ObservableType` enum. Below, we provide a list of observable types, with a brief description, as well as the set of link end types that they require:

* :class:`one_way_range` Range (in meters) between two link ends. Observable has size 1. Requires :literal:`transmitter` and :literal:`receiver` link ends.
* :class:`angular_position` Right ascension and declination (in radians) of one link end (the transmitter) as observed by other link end (the receiver) as measured in the base frame in which the states of the link ends are defined. Observable has size 2. Requires :literal:`transmitter` and :literal:`receiver` link ends.
* :class:`position_observable` Three-dimensional position of link end (in m) as measured in base the frame in which the state of the link end is defined. Observable has size 3 (*x*, *y* and *z* components). Requires :literal:`observed_body` link end, which directly defines the body of which an observation is taken. Note that this observable is not realized in practice, but is typically used for testing or analysis purposes. It does not use a light-time calculation
* :class:`one_way_doppler` Doppler effect of a signal (dimensionless) between two link ends. Observable has size 1. The value approximately requal to the range-rate between the link ends, divided by *c*. The full model one-way Doppler observable :math:`h_{D(1),AB}` from link end *A* to link end *B* is computed from: 
 
   .. math::
      h_{D(1),AB}=\frac{d\tau_{A}}{dt_{A}}\frac{t_{A}}{dt_{B}}\frac{dt_{B}}{d\tau_{B}}-1
      
   where :math:`t` and :math:`\tau` denote coordinate anr proper time. The one-way Doppler requires :literal:`transmitter` and :literal:`receiver` link ends.
* :class:`two_way_doppler` The combined effect from two one-way Doppler effects: one from link end :math:`A` to link end  :math:`B`, and one from :math:`B` to :math:`C`, computed as:  
 
   .. math::
      h_{D(2),ABC}=(h_{D(1),AB}+1)(h_{D(1),BC}+1)-1
      
   The two-way Doppler requires :literal:`transmitter`, :literal:`reflector1` and :literal:`receiver` link ends (:math:`A`, :math:`B` and :math:`C` in the above example).
* :class:`one_way_differenced_range` Doppler observable as it is typically obtained in interplanetary tracking in so-called 'closed-loop' mode (in m/s) between two link ends. Observable has size 1. The value is computed from the averaged range-rate over some integration time (see below). Requires :literal:`transmitter` and :literal:`receiver` link ends.
* :class:`n_way_range` Accumulated range (in meters) over any number of signal paths, may be used for two-way range (as in SLR or deep space tracking), as well as more exotic situations where more than 2 signal paths are used in generating the observable (as was the case for, for instance, the SELENE mission) Observable has size 1. The value is computed accumulating the light-time (multiplied by *c*) with any retransmission delays that the user defines (see below) Requires :literal:`transmitter`, :literal:`receiver`, as well as :literal:`reflector1`, :literal:`reflector2` ... :literal:`reflectorX` for X+1 signal paths (when only a :literal:`transmitter` and :literal:`receiver` are defined the observation is identical to a :literal:`one_way_range`)

.. _observationSettings:

Observation Settings
~~~~~~~~~~~~~~~~~~~~

The settings for most obsevation model types are provided to an object of the :class:`ObservationSettings` class, which takes the type of observable, the :ref:`lightTimeCorrections` and observation biases :ref:`observationBiases`. For some observation models, however, specific additional information is required, and a dedicated :class:`ObservationSettings` derived class is created. Below the input to these classes is discussed in detail.

.. class:: ObservationSettings

An :literal:`ObservationSettings` is created as follows:
   .. code-block:: cpp

      boost::shared_ptr< ObservationSettings > observationSettings =
            boost::make_shared< ObservationSettings >( 
                observableType, lightTimeCorrectionsList, biasSettings );

   where:

   - :literal:`observableType`

      :literal:`ObservableType` that defines the type of observable. An :literal:`ObservationSettings` object may be created directly for the following obsevable types:
      
         * :literal:`one_way_range`
         * :literal:`n_way_range` (for zero retransmission time)
         * :literal:`angular_position`
         * :literal:`position_observable`
         * :literal:`one_way_doppler` (when using :math:`d\tau/dt=1` for both link ends)
         * :literal:`two_way_doppler` (when using :math:`d\tau/dt=1` for all link ends)


   - :literal:`lightTimeCorrectionsList`
  
      A list of light-time corrections of type :literal:`std::vector< boost::shared_ptr<'LightTimeCorrectionSettings > >` that are used to compute the light-time between the link ends, see :ref:`lightTimeCorrections`. 

   - :literal:`biasSettings`

      An object that defines the settings for an observation bias, of type :literal:`boost::shared_ptr< ObservationBiasSettings >`, see :ref:`observationBiases`. 
      
Note that the light-time correction and bias settings are empty by default, so that:

   .. code-block:: cpp

      boost::shared_ptr< ObservationSettings > observationSettings =
            boost::make_shared< ObservationSettings >( observableType );
         
and:
   
   .. code-block:: cpp

      boost::shared_ptr< ObservationSettings > observationSettings =
            boost::make_shared< ObservationSettings >( 
                observableType, lightTimeCorrectionsList );
            
may be used as well to create an observation model without light-time corrections or biases (in the case of the former) and no biases (in the case of the latter).

Additionally, a second constructor is provided that takes a single light-time correction setting, instead of a list, as its second argument. So, you may substitute the input of type :literal:`std::vector< boost::shared_ptr< LightTimeCorrectionSettings > >` by an input of type :literal:`boost::shared_ptr< LightTimeCorrectionSettings >`, in which case you can set only a single light-time correction. 


.. class:: OneWayDifferencedRangeRateObservationSettings

The :literal:`OneWayDifferencedRangeRateObservationSettings` class is used to define settings for one-way differenced range observables (typically termed one-way Doppler in the deep-space tracking community). An object of this type is created as follows, similar to an :literal:`ObservationSettings` object:

   .. code-block:: cpp

      boost::shared_ptr< OneWayDifferencedRangeRateObservationSettings > observationSettings =
            boost::make_shared< OneWayDifferencedRangeRateObservationSettings >( 
                integrationTimeFunction, lightTimeCorrectionsList, biasSettings );

where no input on the type of observable is given (it is :literal:`one_way_differenced_range` by default). The new input is:

   - :literal:`integrationTimeFunction`

      :literal:`boost::function< double( const double ) >` that returns the integration time of the observable as a function of observation time (function input). In many cases, the integration time will be constant, and you may use a :literal:`boost::lambda`, so for a constant 60 s integration time:
      
     .. code-block:: cpp

      boost::function< double( const double ) >  integrationTimeFunction = 
          boost::lambda::constant( 60.0 );
      
As is the case for the :class:`ObservationSettings` class, the second and third constructor arguments are optional, and the second argument may be either a :literal:`std::vector` of :literal:`boost::shared_ptr< LightTimeCorrectionSettings >`, or a single such object.  

.. class:: NWayRangeObservationSettings

The :literal:`NWayRangeObservationSettings` class is used to define settings for n-way range observables. An object of this type is created as follows:

   .. code-block:: cpp

      boost::shared_ptr< NWayRangeObservationSettings > observationSettings =
            boost::make_shared< NWayRangeObservationSettings >( 
                oneWayRangeObsevationSettings, retransmissionTimesFunction, biasSettings );

where no input on the type of observable is given (it is :literal:`n_way_range` by default). The bias settings input is handled in the same way as for the :class:`ObservationSettings` class (and is again empty by default). The other input arguments are:

   - :literal:`oneWayRangeObsevationSettings`

      A :literal:`std::vector< boost::shared_ptr< ObservationSettings > >` list, that has the observation settings for each leg of the n-way link. Note that the observable type of each of the :class:`ObservationSettings` in this list must be :literal:`one_way_range`.

   - :literal:`retransmissionTimesFunction`

      A :literal:`boost::function< std::vector< double >( const double ) >` list that returns the retransmission time at each of the intermediate link ends. For instance, for a Graz station -> LRO -> Matera station n-way-range observable, there may be some delay between LRO receiving the signal from Graz, and retransmitting the signal to Matera. The :literal:`retransmissionTimesFunction` list returns this delay as a function of the observation time at the retransmitting link end. As was the case for the integration time in the :class:`OneWayDifferencedRangeRateObservationSettings` class, you can use :literal:`boost::lambda` to define constant retransmission delay. When providing an empty :literal:`std::vector`, no retransmission delay is assumed.

.. class:: OneWayDopperObservationSettings

The :literal:`OneWayDopperObservationSettings` class is used to define settings for one-way Doppler observables. Here, the term Doppler is used in the instantaneous sense, and is distinct from what is typically termed Doppler in the deep-space tracking community, which is defined by the :class:`OneWayDifferencedRangeRateObservationSettings` class. An object of the :class:`OneWayDopperObservationSettings` type is created as follows:

   .. code-block:: cpp

      boost::shared_ptr< OneWayDopperObservationSettings > observationSettings =
            boost::make_shared< OneWayDopperObservationSettings >( 
                lightTimeCorrectionsList, transmitterProperTimeRateSettings, receiverProperTimeRateSettings, biasSettings );
                               
where:

   - :literal:`lightTimeCorrectionsList`
  
      A list of light-time corrections of type :literal:`std::vector< boost::shared_ptr<'LightTimeCorrectionSettings > >` that are used to compute the light-time between the link ends, see :ref:`lightTimeCorrections`. 

   - :literal:`transmitterProperTimeRateSettings`

      An object that defines the settings for the proper time rate :math:`d\tau/dt` of the transmitter, of type :literal:`boost::shared_ptr< DopplerProperTimeRateSettings >`, see :ref:`properTimeRates`.

   - :literal:`receiverProperTimeRateSettings`
   
      An object that defines the settings for the proper time rate :math:`d\tau/dt` of the transmitter, of type :literal:`boost::shared_ptr< DopplerProperTimeRateSettings >`, see :ref:`properTimeRates`.

   - :literal:`biasSettings`
   
      An object that defines the settings for an observation bias, of type :literal:`boost::shared_ptr< ObservationBiasSettings >`, see :ref:`observationBiases`. 

.. class:: TwoWayDopperObservationSettings

The :literal:`TwoWayDopperObservationSettings` class is used to define settings for two-way Doppler observables:

   .. code-block:: cpp

      boost::shared_ptr< TwoWayDopperObservationSettings > observationSettings =
            boost::make_shared< TwoWayDopperObservationSettings >( 
                uplinkOneWayDopplerSettings, downlinkOneWayDopplerSettings, biasSettings );
                               
where:

   - :literal:`uplinkOneWayDopplerSettings`

      An object that defines the settings for uplink one-way Doppler, of type :literal:`boost::shared_ptr< OneWayDopperObservationSettings >`
      
   - :literal:`downlinkOneWayDopplerSettings`
   
      An object that defines the settings for uplink one-way Doppler, of type :literal:`boost::shared_ptr< OneWayDopperObservationSettings >`

   - :literal:`biasSettings`
   
      An object that defines the settings for an observation bias, of type :literal:`boost::shared_ptr< ObservationBiasSettings >`, see :ref:`observationBiases`. 

.. _properTimeRates:
     
Proper-time Rate Settings
*************************

When creating a one- or two-way Doppler model, the user can provide settings for the proper time rate model at the link ends, which is handled by the :class:`DopplerProperTimeRateSettings` class. Each type of proper-time rate model has its own dedicated derived class, which are described below.

.. class:: DirectFirstOrderDopplerProperTimeRateSettings

The :literal:`DirectFirstOrderDopplerProperTimeRateSettings` class is derived from  :literal:`DopplerProperTimeRateSettings` and defines the proper time rate due to a single static point mass. In this sense, the term 'static' indicates that the motion of the body causing a change in proper time rate is not incorporated in the model. An object of this type is created as:

   .. code-block:: cpp

      boost::shared_ptr< DirectFirstOrderDopplerProperTimeRateSettings > properTimeRateSettings =
            boost::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( 
                bodyName );
                               
where:

   - :literal:`bodyName`

      An :literal:`std::string` that defines the name of the body that whose mass is causing a proper time rate.

As an example, a one-way Doppler model, where the Earth's mass is used to perturtb the transmitting station's proper time rate, and the Moon's mass the receiver's proper time rate would be by creating :class:`OneWayDopperObservationSettings` as follows:

   .. code-block:: cpp
   
      boost::shared_ptr< DopplerProperTimeRateSettings > transmitterProperTimeRateSettings = 
         boost::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Earth" );
      boost::shared_ptr< DopplerProperTimeRateSettings > receiverProperTimeRateSettings = 
         boost::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Moon" );
      boost::shared_ptr< OneWayDopperObservationSettings > observationSettings =
            boost::make_shared< OneWayDopperObservationSettings >( 
                std::vector< boost::shared_ptr< LightTimeCorrectionSettings > >, transmitterProperTimeRateSettings, receiverProperTimeRateSettings );
                
For a case where no light-time corrections, and no observation bias, is used
                               


.. _lightTimeCorrections:

Light-time Corrections
~~~~~~~~~~~~~~~~~~~~~~

When computing the basic light-time between two link ends, the Euclidean distance between them is computed for a signal travelling at exactly *c* (the speed of light in vacuum). However, various effects must in reality be accounted for to compute the true light time. Tudat currently supports:

* First-order relativistic light-time correction: the correction to the light time of a (set of) stationary point masses, computed up to :math:`c^{-2}` according to general relativity.

As is the case with many other models, a base class settings object is provided :literal:`LightTimeCorrectionSettings`. Specific light-time corrections are defined in derived classes:

.. class:: FirstOrderRelativisticLightTimeCorrectionSettings

   The :class:`FirstOrderRelativisticLightTimeCorrectionSettings` defines settings for a first-order relativistic correction to the light-time, as formulated by *e.g.* Moyer (2000). The class is created by:
   
    .. code-block:: cpp

      boost::shared_ptr< FirstOrderRelativisticLightTimeCorrectionSettings > lightTimeCorrectionSettings =
            boost::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( 
                perturbingBodies );
   
 The input is:

   - :literal:`perturbingBodies`

      A :literal:`std::vector< std::string >` containing the names of the bodies due to which the light-time correction is to be taken into account.
      
      .. note:: One ambiguity in the model is the time at which the states of the perturbing bodies are evaluated. We distinguish two cases: one where the perturbing body contains one of link ends, and one where it does not. In the case where the perturbing body contains a link end (for instance perturbation due to Earth gravity field, with one of the link ends being an Earth-based station), the time at which the Earth's state is evaluated equals the transmission time if Earth acts as transmitter, and reception time if Earth acts as receiver. In other cases, where the perturbing body is not involved in the link ends, its state is evaluated at the midpoint time between transmitter and receiver.

.. _observationBiases:
      
Observation Biases
~~~~~~~~~~~~~~~~~~

Real observations in reality often differ (slightly) from the physical ideal value due to for instance instrumental effects at the link ends. Part of these effects influence the observation in a predictable manner, and can be modelled deterministically. It is these deterministic effects that we collectively term observation biases in Tudat. The stochastic errors are not used as part of the observation model in Tudat, but can be incorporated when simulating observations. 

In Tudat, the following types of biases are currently incorporated, where :math:`h` is the physically ideal observation, and :math:`\tilde{h}` is the biases version of the observable.

   * Relative bias :math:`K_{r}`, which influences the observable as:

   .. math::
      \tilde{h}=h(1+K_{r})
      
      For an observable with size greater than 1, :math:`K_{r}` is a vector and the multiplication is component-wise.

   * Absolute bias :math:`K_{a}`, which influences the observable as: 
   
   .. math::
      \tilde{h}=h+K_{a}
      
      For an observable with size greater than 1, :math:`K_{a}` is a vector and the multiplication is component-wise.

   * A combined bias, which is computed from any number of the above bias types combined. Note that each contribution of the combined bias is computed from the unbiased observable, so when applying both a relative and absolute bias, we get:

   .. math::
      \tilde{h}=h+K_{a}+hK_{r}

      And, crucially:

   .. math::
      \tilde{h}\neq (h+K_{a})(1+K_{r})

As discussed above, the biases are created by passing a :class:`ObservationBiasSettings` (or derived class) object to the :class:`ObservationSettings` (or derived class)  constructor. Each bias type has a dedicated derived class of :class:`ObservationBiasSettings`, which are defined as follows:

.. class:: ConstantObservationBiasSettings

   The :literal:`ConstantObservationBiasSettings` class is used to define settings for an absolute bias :math:`K_{a}` and is created by:

   .. code-block:: cpp

      boost::shared_ptr< ConstantObservationBiasSettings > observationBiasSettings =
            boost::make_shared< ConstantObservationBiasSettings >( 
                observationBias );

   The input is:

   - :literal:`observationBias`

      An :literal:`Eigen::VectorXd` of the same size as the observable (see AAAAA) containing the values of the vector :math:`K_{a}`. 
      
.. class:: ConstantRelativeObservationBiasSettings

   The :literal:`ConstantRelativeObservationBiasSettings` class is used to define settings for an absolute bias :math:`K_{r}` and is created by:

   .. code-block:: cpp

      boost::shared_ptr< ConstantRelativeObservationBiasSettings > observationBiasSettings =
            boost::make_shared< ConstantRelativeObservationBiasSettings >( 
                observationBias );

   The input is:

   - :literal:`observationBias`

      An :literal:`Eigen::VectorXd` of the same size as the observable (see AAAAA) containing the values of the vector :math:`K_{r}`. 
       
.. class:: MultipleObservationBiasSettings

   The :literal:`MultipleObservationBiasSettings` class is used to define settings for a combined bias, that is composed of multiple bias models.

   .. code-block:: cpp

      boost::shared_ptr< MultipleObservationBiasSettings > observationBiasSettings =
            boost::make_shared< MultipleObservationBiasSettings >( 
                biasSettingsList );

   The input is:

   - :literal:`biasSettingsList`

      A :literal:`std::vector< boost::shared_ptr< ObservationBiasSettings > >` list containing the bias models to be applied.








