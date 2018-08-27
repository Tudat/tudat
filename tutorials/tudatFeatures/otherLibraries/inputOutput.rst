.. _tudatFeaturesInputOutput:

Input/Output
============
Handling input and output can be an important part of an application, as Tudat is built primarily for numerical orbit propagation/mission analysis, and not for visualization, which most users choose to do in e.g. Matlab. This page provides an overview of the most important Tudat methods for handling input and output, their set-up and how they can be used in an application.

Available Methods
~~~~~~~~~~~~~~~~~
All functions and classes that handle input/output can be found in the :literal:`InputOutput` directory in the Tudat libraries. These functions and classes provide a variety of handlers and can, for example, skip lines, eliminate comments from a file, or store the propagation history to a particular file. In most cases, however, these specific file readers will not need to be called directly by the user, but are instead used when retrieving a default spherical harmonic gravity field or an aerodynamic database.

The input/output elements from Tudat that will be discussed here, can be subdivided into three categories:

- Path-functions.
- File-readers.
- File-writers.

Path-functions
**************
These are functions that return the root-path of the library (or a particular folder). They can be used to define a path relative to the root-path of Tudat when using a file-reader/writer. The following path functions are available:

   - :literal:`getTudatRootPath` (Tudat)

      This function returns the root-path corresponding with the root-directory of the Tudat library as a string with trailing slash included.

   - :literal:`getSpiceKernelPath` (Tudat)

      This function returns the path where the Spice Kernels are located as a string with trailing slash included. This is very useful when writing code that will work cross-platform. For instance, if you have a file inputData.txt that is located at::

         C:/Software/tudatBundle/tudat/Tudat/Astrodynamics/
    
      you will want to avoid writing this full path in your code, as it will not run on any other system where Tudat is installed in a different directory. Instead, what you can do:

      .. code-block:: cpp

         std::string filePath = input_output::getTudatRootPath( ) + "Astrodynamics/";
         std::string fileName = "inputData.txt" 

   - :literal:`getAtmosphereTablesPath` (Tudat)

      This function returns the path where the tabulated atmosphere files are located. You can use this repository to store and easily access your custom atmosphere tables or to access the built-in atmosphere files.

   .. note:: The full path/name of the input data file is now simply :literal:`filePath + fileName`, which will correctly identify the file location on any system.

File-readers
************

If needed, you can use this functionality as a starting point to create your own file-reader for your specific file type. If you run into issues when doing so, please contact the Tudat support team. However, there is a dedicated file-reader available in the Tudat library:

   - :literal:`readMatrixFromFile` (Tudat)
      
      This function can be used to read a simple text file with separated numbers. The documentation includes the options for a separator and the character used to skip a line (e.g. a comment-line). Note that the new-line character (:literal:`\\n`) is reserved to split the lines of the matrix.

   - :literal:`MultiArrayFileReader< NumberOfDimensions >::readMultiArray` (Tudat)

      You can use this function to read and access the information stored in a file with multiple dimensional data. The output is a variable of type :literal:`boost::multi_array< double, NumberOfDimensions >`, whose number of dimensions is specified by the template parameter :literal:`NumberOfDimensions`.

      .. note:: It is important to keep in mind that to use this function, you also need to use the class specifier :literal:`MultiArrayFileReader< NumberOfDimensions >`, where you replace the template argument with the actual number of dimensions yoiu expect from the file.

      You can also use the function :literal:`MultiArrayFileReader< NumberOfDimensions >::readMultiArrayAndIndependentVariables` to access both the multi-array data and a set of independent variables. This is the format used, e.g., to store aerodynamic coefficients and the atmosphere tables.

      .. tip:: You can find a description of how the data is expected to be stored in the file, by looking at the end of :ref:`tudatTabulatedAtmosphere`.

File-writers
************
The :literal:`InputOutput` directory in the Tudat library also contains functionality to write data (e.g. the propagation history of a satellite) to a file. The following function is available:

   - :literal:`writeDataMapToTextFile` (Tudat)
      This function writes data stored in a map to a text file. A number of overloads exists for this function based on the input given to the function. Furthermore, the data-map can store different types of data (e.g. doubles and Eigen vectors, which are typical types for the propagation history). The following overload is most relevant:

      .. code-block:: cpp

         writeDataMapToTextFile( dataMap, 
      	                         outputFileName,
      	                         outputDirectory,
      	                         fileHeader,
      	                         precisionOfKeyType,
      	                         precisionOfValueType,
      	                         delimiter )

      If only the :literal:`dataMap` and file-name are provided, the default :literal:`KeyType`-precision and :literal:`ValueType`-precision (:literal:`digits10` from :literal:`limits` standard library), output directory (Tudat root-path), and delimiter (space) are used. If also the :literal:`outputDirectory` is given an empty file header, a precision of 16 significant digits and a tab as delimiter are used, but these can be user-specified.

   - :literal:`writeMatrixToTextFile` (Tudat)
      This function writes data stored in a :literal:`Eigen::Matrix` to a text file. The input required are the matrix itself and the file-name. Note that any :literal:`scalarType` and number of rows and collumns can be used.

   - :literal:`MultiArrayFileWriter< NumberOfDimensions, NumberOfCoefficients >::writeMultiArrayToFile` (Tudat)

      .. warning:: This function is not yet available, but you can use :literal:`MultiArrayFileWriter< NumberOfDimensions, NumberOfCoefficients >::writeMultiArrayAndIndependentVariablesToFiles` to write both a multi-array and a set of independent variables to a file.

Examples
~~~~~~~~
Text File to MatrixXd
*********************

An example of reading data from a text file to a :literal:`Eigen::MatrixXd` is shown in detail in :ref:`walkthroughsUseOfThrustUserDefinedThrustVector`. A small overview is presented here:

For example a file named :literal:`.txt` contains data structured as follows::

    0       0 0 5
    6068    0 1 5
    6097    1.0 0 5
    6097.5  0.8 0 5
    6098    0.6 0.1 5
    6099    0.1 0.5 5
    12192   0.2 1.0 4.5
    18288   0.3 1.5 4.0
    243575  0.4 2.0 3.0
    3.999e6 1.0 1.0 2.0
    4e6     1.1 5.0 1.0

thus 4 columns spaced with tabs. This file can be read with the following code::

   Eigen::MatrixXd thrustForceMatrix =
            tudat::input_output::readMatrixFromFile( cppFolder + "nameOfFile.txt" , " \t", "#" );

where the first argument is the relative path to the :literal:`.txt` file, the second argument indicates the type(s) of separator(s) used (multiple seperators possible). The last argument indicates the character used for lines to be skipped. 

Data-map (double,double) to Text File
*************************************
A data map is a template class that is defined by its key-type and value-type:

.. code-block:: cpp

    std::map< key-type, value-type >

Using this, a data map, where the type of the key is a :literal:`double`, and the type of the value is also a :literal:`double`, can be defined as:

.. code-block:: cpp

    std::map< double, double > keyDoubleValueDoubleMap;

Each entry in the data map consists of a key and a value and is entered using:

.. code-block:: cpp

    keyDoubleValueDoubleMap[ key ] = value;

As an example, three entries are stored in this data map:

.. code-block:: cpp

    keyDoubleValueDoubleMap[ std::sqrt( 3.0 ) ] = 1.0 / std::sqrt( 2.0 );
    keyDoubleValueDoubleMap[ 4.5 ] = 56.89;
    keyDoubleValueDoubleMap[ 12.65 ] = 1.0 / 3.0;

Now, this data-map can be stored to a file using:

.. code-block:: cpp

    tudat::input_output::writeDataMapToTextFile(
                keyDoubleValueDoubleMap, "keyDoubleValueDoubleMapDataFileWithDefaults" );

Data-map (double,Vector3d) to Text File
***************************************
An example of a data map, where the type of the key is a :literal:`double`, and the type of the value is an :literal:`Eigen::Vector3d`:

.. code-block:: cpp

    std::map< double, Eigen::Vector3d > keyDoubleValueVector3dMap;
    keyDoubleValueVector3dMap[ 1.1 ] = Eigen::Vector3d( 0.0, 1.3, -6.54 );
    keyDoubleValueVector3dMap[ 6.5 ] = Eigen::Vector3d( -4.56, 1.23, -9.98 );
    keyDoubleValueVector3dMap[ 10.9 ] = Eigen::Vector3d( -46.13, 1.0 / 3.0, std::sqrt( 2.0 ) );

This data-map can be stored to a file using:

.. code-block:: cpp

    tudat::input_output::writeDataMapToTextFile(
                keyDoubleValueVector3dMap, "keyDoubleValueVector3dMapDataFile" );

Data-map (int, Matrix3d) to Text File
***************************************
An example of a data map, where the type of the key is a :literal:`int`, and the type of the value is an :literal:`Eigen::Matrix3d`:

.. code-block:: cpp

   Eigen::Matrix3d threeDimensionalMatrix;
   threeDimensionalMatrix << 0, 1, 2, 3, 4, 5, 6, 7, 8;

   std::map<int, Eigen::Matrix3d> matrixMap;
   matrixMap[14] = threeDimensionalMatrix;

This data-map can be stored to a file using:

.. code-block:: cpp

    tudat::input_output::writeDataMapToTextFile(
                matrixMap, "matrixMap" );

This results in::

   14          0                 1                 2                 3                 4                 5                 6                 7                 8

Note that all matrix entries are put on one line when writing a map to file. This makes the file easy to read by other programs. 

Eigen::Matrix3d to Text File
****************************
An example of a matrix to save:

.. code-block:: cpp

   Eigen::Matrix3d squareDoubleValueMatrix;
   squareDoubleValueMatrix << 0, 1, 2, 3, 4, 5, 6, 7, 8;

This can be saved to a text file using:

.. code-block:: cpp

   tudat::input_output::writeMatrixToTextFile(
   				squareDoubleValueMatrix, "squareDoubleValueMatrixFile")     

This results in::

   0                	 1                	 2                
   3                	 4                	 5                
   6                	 7                	 8                

Note the difference with saving matrices inside a map, which will put all matrix entries on one line.

Storing Propagation History
***************************
A good example on how to store the propagation history in a data map can be found in the example applications in the Tudat Bundle. If you have downloaded the bundle, these examples can be found in::

    tudatBundle/tudatApplications/satellitePropagatorExamples/SatellitePropagatorExamples

On the other hand, an example of saving the propagation history to a data file is described at the end of the :ref:`walkthroughsUnperturbedEarthOrbitingSatellite` walkthrough. 

.. tip:: You can also scroll to the end of :ref:`tudatFeaturesSimulatorCreation`, for an overview on how to access and save the propagation history, as well as the **dependent variables**, from a :class:`DynamicsSimulator` object.
