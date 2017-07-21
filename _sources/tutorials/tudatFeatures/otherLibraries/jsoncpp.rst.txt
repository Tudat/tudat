.. _tudatFeaturesJsoncpp:

External Libraries: JSONCPP
===========================
JSONCPP is a library that allows reading, writing and manipulation of JSON data. JSON is a lightweight data-interchange format. It can represent numbers, strings, ordered sequences of values, and collections of name/value pairs. JSONCPP can be used to read in configuration files or write out data.
 
JSONCPP is available as a submodule for Tudat Bundle. This means that building and downloading gets taken care of by the top-level :literal:`CMakeLists.txt`. We even created an example application demonstrating the use of :literal:`JSON` and :literal:`CSPICE` for your convience. You can follow the instructions below if you do not wish to use Tudat Bundle.

Download and build JSONCPP library
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Clone JSONCPP from the GitHub project page::

    git clone https://github.com/open-source-parsers/jsoncpp jsoncpp

.. warning:: Make sure the folder is named jsoncpp (and not for example jsoncpp-master) as the FindJSONCPP.cmake file searches for the following file: jsoncpp/include/json/json.h.

You can now build the project with the flags :literal:`-DJSONCPP_WITH_PKGCONFIG_SUPPORT=OFF -DJSONCPP_WITH_TESTS=OFF`. Copy :literal:`libjsoncpp.a` to the :literal:`jsoncpp/lib` folder. After building this file can be found by default in the build directory :literal:`build[-qtcreator]/src/lib_json/libjsoncpp.a`. CMake can be configured to build the static library in the proper directory.
Your final file structure should look as follows::

    jsoncpp
    |-- include
    |   |-- json
    |       |-- assertions.h
    |       |-- autolink.h
    |       |-- config.h
    |       |-- features.h
    |       |-- forwards.h
    |       |-- json.h
    |       |-- reader.h
    |       |-- value.h
    |       |-- version.h
    |       `-- writer.h
    |-- lib
    |   |-- libjsoncpp.a

Integrating JSONCPP in your application
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
1. Copy :literal:`FindJSONCPP.cmake` into your :literal:`CMAKE_MODULE_PATH`.
2. Change your project's :literal:`CMakeLists.txt` to find and include the package:

   .. code-block:: cmake

       # Find JSONCPP library on local system.
       find_package(JSONCPP)

       # Include JSONCPP directories.
       if(NOT APPLE)
          include_directories(SYSTEM AFTER "${JSONCPP_INCLUDE_DIR}")
       else( )
          set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${JSONCPP_INCLUDE_DIR}\"")
       endif( )

3. Add the jsoncpp library to your applications :literal:`target_link_libraries` command:

    .. code-block:: cmake

        target_link_libraries(my_project tudat_[libs] jsoncpp ${Boost_LIBRARIES} )

4. Include the :literal:`json` value header in your code and use the jsoncpp library:

    .. code-block:: cmake

        #include <json/json.h>
        #include <json/value.h>

Usage example: setting up a json configuration file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Consider the following configuration file called config.json:

.. code-block:: javascript

    // Configuration options
    {
        // Default encoding for text
        "encoding" : "UTF-8",

        // Plug-ins loaded at start-up
        "plug-ins" : [
            "python",
            "c++",  // trailing comment
            "ruby" 
            ],

        // Tab indent size
        // (multi-line comment)
        "indent" : { /*embedded comment*/ "length" : 3, 

    "use_space": true }
    }

In your application you should include:

.. code-block:: cpp

    #include <json/json.h>
    #include <json/value.h>
    #include <fstream>
    [...]
    int main(){
    [...]
    // Construct Json Database
    Json::Value jsonRoot;
    jsonRoot["test"] = "hello"; // write value
    std::cout << jsonRoot["test"].asString() << std::endl;

    // Read file and overwrite Json database
    std::string filename ("/full/path/to/file.json");
    std::ifstream config_doc(filename.c_str(), std::ifstream::binary);
    config_doc >> jsonRoot; // write file to jsonRoot
    std::cout << jsonRoot["encoding"].asString() << std::endl;
    std::cout << jsonRoot["test"].asString() << std::endl; // this instance is overwritten.

which will return::

    hello
    UTF-8

Usage example: read an Eigen matrix or an Eigen vector
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The following functions can be used to read a matrix or a vector from a Json file.

.. code-block:: cpp

    //// Read Json Array and return as Eigen Matrix
    Eigen::MatrixXd ReadJsonMatrix( Json::Value Array )
    {
        int rows = Array.size( );
        int cols = Array[0].size( );

        Eigen::MatrixXd Matrix( rows, cols );

        for(int i = 0 ; i<rows ; i++)
        {
            for(int j=0 ; j<cols ; j++)
            {
                Matrix(i,j) = Array[i][j].asDouble( );
            }
        }

    //    std::cout << "Matrix size: rows = " << rows << " cols = " << cols << std::endl;
    return Matrix;
    }

    //// Read Json Array and return as Eigen Vector
    Eigen::VectorXd ReadJsonVector( Json::Value Array )
    {
        int size = Array.size( );

        Eigen::VectorXd Vector( size );

        for(int i = 0 ; i < size ; i++)
        {
            Vector(i) = Array[i].asDouble( );
        }

        return Vector;
    }

    Json::Value settings = ReadConfigFile(folder,filename);
    Eigen::MatrixXd A = ReadJsonMatrix(settings["A"]);
    Eigen::VectorXd X0 = ReadJsonVector(settings["X0"]);
    std::cout << A << std::endl;
    std::cout << X0 << std::endl;

where the matrix A is defined in the config file as an array:

.. code-block:: javascript

    "A" : [[1,2,3],[3,4,5],[1,1]],
    "X0" : [0.5,-0.25],

which will return::

    1 2 3
    3 4 5
    1 1 0
      0.5
    -0.25

The function uses the first row to define the number of columns of the matrix A. Note that the matrix A is defined using an Array with varying size. A 0 is inserted at A(3,3), because the array is not defined for this location.
