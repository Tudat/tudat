.. _externalEigen:

Eigen: Basic Concepts
=====================
This page provides an overview of the available information on the Eigen external library. Eigen is a C++ library specifically for linear algebra. This library is extensively used throughout Tudat for all linear algebra computations. Eigen is a pure template library defined in header files only. This means that Eigen does not need to be compiled and its header files can be used right away.

Eigen's documentation is considered to be very good and extensive, so it is recommended that you check out the documentation on the `Eigen website <http://eigen.tuxfamily.org/index.php?title=Main_Page>`_. Not only do they provide documentation on the code, the website also includes a long tutorial.

Why Use Eigen?
~~~~~~~~~~~~~~
Eigen is free software that can handle many different linear algebra operations and also has a geometry framework. Furthermore, the code is mature, well maintained, well tested, and has good documentation.

Vector storage containers from the Eigen library should only be used when linear algebra operations are required on it; for other purposes (e.g. storing data) the vector storage container ``std`` of the standard C++ library should be used.

Common Eigen Types
~~~~~~~~~~~~~~~~~~

A couple of Eigen types that are frequently used in Tudat are presented in this table:

+-------------+-----------------+----------------+---------------------------------------------------------------------------------------------+
|**Type Name**|**Fixed/Dynamic**|**Eigen/Tudat** | **Example Link**                                                                            |
+-------------+-----------------+----------------+---------------------------------------------------------------------------------------------+
| Vector3d    | Fixed           | Eigen          | :ref:`Tutorial Eigen Vectors <externalEigenExamples>`                                       |
+-------------+-----------------+----------------+---------------------------------------------------------------------------------------------+
| Vector6d    | Fixed           | Tudat          | :ref:`Tutorial Eigen Vectors <externalEigenExamples>`                                       |
+-------------+-----------------+----------------+---------------------------------------------------------------------------------------------+
| VectorXd    | Dynamic         | Eigen          | :ref:`Tutorial Eigen Vectors <externalEigenExamples>`                                       |
+-------------+-----------------+----------------+---------------------------------------------------------------------------------------------+
| Matrix3d    | Fixed           | Eigen          | `Eigen Matrix Example <https://eigen.tuxfamily.org/dox/group__TutorialMatrixClass.html>`_   |
+-------------+-----------------+----------------+---------------------------------------------------------------------------------------------+
| MatrixXd    | Dynamic         | Eigen          | `Eigen Matrix Example <https://eigen.tuxfamily.org/dox/group__TutorialMatrixClass.html>`_   |
+-------------+-----------------+----------------+---------------------------------------------------------------------------------------------+

The table indicates whether it is a fixed size type or a dynamically sized type (second column). This is also indicated by the type name: the X indicates that the number of rows (and columns, for a matrix) are dynamically allocated at compile time. There is a difference between the fixed size types and the dynamically sized types. Information on when to use which type can be found below.

The table also indicates whether the type is actually part of the Eigen library itself or if you should look for it in Tudat (third column). This is indicated, because the 6-element vector typedefs, such as Vector6d, are not provided by the Eigen library. However, for Astrodynamics calculations it is very common to use state vectors, and fixed-size type memory allocation is typically faster that dynamic allocation (see below), therefore, typedefs for 6-element vectors have been added to Tudat.

An example on how to use some of these Eigen types can be found under the links in the last column of the table. Where the Tutorial Eigen Vectors is a specific Tudat example and the link to the Eigen Matrix Example takes you to the Eigen website.

Fixed Size vs. Dynamic Size
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The section on frequently used Eigen types introduced three different vector types and two different matrix types. As you may have noticed, the difference between these types is the size, and whether this size is fixed or dynamic.

This difference is related to the way that the memory is allocated. The fixed-size variables are all declared in the source code, their size is known at compile time, and the data is allocated to a reserved region of the memory referred to as the stack (or in some cases to the main memory). The dynamically-sized variables are also declared in the source code, however, their size is unknown at compile time, the data needs to be allocated at runtime to another large region of the memory referred to as the heap.

Details of this memory management (the heap vs. the stack) are beyond the scope of this discussion (an example of a comprehensive explanation can be found here).
In general, you should consider the following guidelines:

    - Use fixed-size types, when you need relatively small vectors (or matrices) of 2,3,4 or 6 elements. This stack-based memory allocation is very simple and typically faster than heap-based memory allocation.
    - Use dynamically-sized types to allocate large vectors or matrices. The heap is typically a large pool of memory, in contrast to the stack, which has a limited size. This to prevent that it eats up most of the stack memory and to prevent stack overflow, which may result from this. Also, in any cases where the size of the vector or matrix is not known at compile time, use dynamically-sized types.

Frame Transformations
~~~~~~~~~~~~~~~~~~~~~
In the Tudat library, Eigen is also used for frame transformations. However, note that the definition of the rotations in Eigen differs from the definition taught at the Faculty of Aerospace Engineering by Mooij and Chu. The difference is that in Eigen the vector itself is rotated, whereas the definition of Mooij and Chu rotates the reference frame. The definitions as used in the Eigen library are described in 'Applied Cartesian Tensors for Aerospace Simulations' by D. Henderson (or in the code itself).

Care should therefore be taken when using directional cosine matrices or quaternions from the Eigen library. A simple work-around, when defining the DCM or quaternion with angle and axis, is to multiply the angle with minus one.

Tips
~~~~

    - When declaring an arbitrary vector of type :literal:`VectorXd` in your code, you should always allocate the size.
    - Eigen :literal:`VectorXd` objects allow for dynamic allocation of size. Not allocating the size though and operating on a :literal:`VectorXd` object can result in very cryptic errors.
    - Access the coefficients of a matrix using :literal:`matrix( i, j )`. Notice that it is required to use round brackets to access the coefficients of a matrix.
    - The index of the first element in an Eigen matrix or vector is 0, as opposed to the convention in mathematics that the first index is 1.
    - Eigen supports Matrices and Array objects. There is a difference between the two: the matrix class is used for linear algebraic operations such as matrix multiplication; the array class is used for coefficient-wise operations. In case you need to use both operations on one object, Eigen does have a solution for this: the matrix class has an :literal:`.array( )` method; and the array class has an :literal:`.matrix( )` method, which allows a member of one class to be treated as a member of the other class.
    - You can initialize all coefficients of a vector, matrix or array to zero by using the method :literal:`Zero( )`. For example: :literal:`Eigen::Vector3d::Zero()` creates a 3 element vector of doubles and initializes all elements to zero. Note that initialization using this method should happen immediately when an Eigen object is declared.
    - You can create a vector of linearly spaced elements using the method :literal:`LinSpaced( size, low, high )`.
