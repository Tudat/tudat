.. _troubleshootingFrequentMistakes:

Frequently made  coding mistakes - C++
======================================
This section sums op some of the most made mistakes leading to failures either in runtime or during compilation. This page is under construction and will be elaborated in the future.

	- Accessing vector elements out of vector bounds

		For instance, in

		.. code-block:: cpp

			Eigen::Vector2d vectorOfLengthTwo( 2.0, 3.0 );
			double intermediateResult = vectorOfLengthTwo( 2 );

		the code will compile, but the value obtained in ``intermediateResult`` is the value located just next to the ``vectorOfLengthTwo`` in your memory and could be anything. 

	- SPICE: Negative value for BEGIN address

		This occurs when inputting ``NaN``'s into Spice as discussed `here <https://github.com/Tudat/tudat/issues/263>`_. Take a look at where you call Spice inside your application.
