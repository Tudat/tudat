.. _cppTutorials:

Selected C++ Aspects
====================
This page links to tutorials for C++ features commonly used within Tudat. It contains links to existing tutorials online and gives some code examples of its use within Tudat. 

C++ basics
~~~~~~~~~~

For those who are not familiar with the  C++ language, some excellent tutorials can be found at `cplusplus.com <http://www.cplusplus.com>`_. Also, this site contains documentation of the standard C++ library ``std::`` and could be of use when using functions from it.

Iterators
~~~~~~~~~

Iterators are used in Tudat to acces elements of containers. It provides a general method of accesing elements of containers independent of the container type. 

A short description of iterators is written by  Alex Allain on `Cprogramming.com <https://www.cprogramming.com/tutorial/stl/iterators.html>`_. 
A more elaborate tutorial which also explains the different iterator types is written by K. Hong on `BoboToBogo.com <http://www.bogotobogo.com/cplusplus/stl3_iterators.php>`_.

Below are some examples of how iterators are used within Tudat. The first is from :ref:`walkthroughsInnerSolarSystemPropagation`. Here an iterator is used to obtain the propagation history of all propagated bodies from the ``integrationResult`` container as obtained from the :class:`DynamicsSimulator`:

.. code-block:: cpp

   // Retrieve numerically integrated state for each body.
   std::vector< std::map< double, Eigen::VectorXd > > allBodiesPropagationHistory;
   allBodiesPropagationHistory.resize( numberOfNumericalBodies );
   for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = integrationResult.begin( );
        stateIterator != integrationResult.end( ); stateIterator++ )
   {
       for( unsigned int i = 0; i < numberOfNumericalBodies; i++ )
       {
           allBodiesPropagationHistory[ i ][ stateIterator->first ] = stateIterator->second.segment( i * 6, 6 );
       }
   }

Which loops over all the keys (integration time) from the ``integrationResults`` and obtains the state of all numerical bodies, which are all stored inside one single ``Eigen::VectorXd`` during integration, and puts it in seperate ``Eigen::VectorXd``'s.  


Another example is used in the example application described in :ref:`walkthroughsUseOfThrustUserDefinedThrustVector`:

.. code-block:: cpp

    // Manually add thrust force in LVLH frame to output
    for( std::map< double, Eigen::VectorXd >::iterator outputIterator = dependentVariableResult.begin( );
         outputIterator != dependentVariableResult.end( ); outputIterator++ )
    {
        Eigen::Matrix3d currentRotationMatrix =
                getMatrixFromVectorRotationRepresentation( outputIterator->second.segment( 3, 9 ) );
        Eigen::Vector3d currentThrust = outputIterator->second.segment( 0, 3 );
        Eigen::VectorXd newOutput = Eigen::VectorXd( 15 );
        newOutput.segment( 0, 12 ) = outputIterator->second;
        newOutput.segment( 12, 3 ) =
                integrationResult.at( outputIterator->first )( 6 ) *
                ( currentRotationMatrix.transpose( ) * currentThrust );
        dependentVariableResult[ outputIterator->first ] = newOutput;
    }

Here an iterator is used to acces the elements of a map. In this case the ``dependentVariables`` saved during propagation. One creates an iterator from the container used (``std::map``). The loop starts at the first element of the dependentVariable container (in this case from the initial time of the propagation) by using the ``.begin()`` function and ends when the iterator reaches the last element in the ``std::map`` found by using the ``.end()`` function. The value inside the map is accessed by::

   outputIterator->second






