.. _socisIndex:

SOCIS 2019
==========

This page serves as an introduction to the Tudat software, for students considering to participate in the SOCIS 2019 project. First, a brief overview of the history and functionality of Tudat is given, followed by a brief description of the SOCIS 2019 project ideas.

Tudat - Introduction
~~~~~~~~~~~~~~~~~~~~

Development of the Tudat software suite was started in 2010, after it became clear that there was a pressing need for a modular, verified and efficient astrodynamics and optimization toolbox, for use by M.Sc. and Ph.D. students of the Astrodynamics and Space (AS) Missions (then Astrodynamics and Satellite Systems) section at the faculty of Aerospace Engineering of Delft University of Technology. The main problem with previous efforts was the lack of a clear user interface, limited modularity, and no long-term vision in the development of the software. With the development of Tudat, our goal has been to overcome these issues, and develop an open-source software package that addresses the needs of both hands-on education and cutting-edge research. 

Over the 2010-2014 period, the Tudat software went through several iterations, as the development team learned valuable lessons on software design, testing, code dissemination and project management. As of early 2015, the software has been hosted on github (http://www.github.com/tudat), and is freely available. A complete overhaul of the interfaces of Tudat was performed in 2015, and the software has been quickly developing since that time, in terms of both functionality and interface quality. An important contribution has been the development of a JSON interface during SOCIS 2017, which allows selected Tudat functionality (numerical state propagation) to be performed using JSON input files, as opposed to a dedicated custom-built C++ program. 

The developments made since 2015 has lead to Tudat being key in a number of courses in the TU Delft Space Exploration M.Sc. curricullum (Space Project; Numerical Astrodynamics; Propagation and Optimization in Astrodynamics). Moreover, Tudat has been used in dozens of M.Sc. projects, several Ph.D. projects, and has contributed to >10 peer-reviewed publications.

Tudat - Design and Functionality
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Tudat is written entirely in C++, using the CMake build environment. Key Tudat functionality falls into the following categories:

* Numerical propagation of dynamics. Tudat contains a variety of models for the space environment, accelerations, propagation schemes, integration schemes, etc. Although the main application of the numerical part of Tudat is orbit propagation, it is capable of combining the numerical propagation of various types of dynamics (translational, rotational, mass variations, etc.), over a single arc, multiple arcs, or a combination. The design of the software ensure that the numerical propagation of the governing differential equations is done in a fully consistent manner. Moreover, the modular setup of Tudat provides the capabilities to easily switch between different models, both to assess the impact of model selection, and to tune the fidelity of the results to the requirements of a given project.
* Precise orbit determination. Tudat is set up for the modern requirements of tracking data simulation, with a focus for planetary missions. This includes the definition of numerous observation models and their associated partial derivatives (n-way range, n-way Doppler, angular position, etc.). In addition, Tudat provides the automated computation of the state transition and sensitivity matrices, for use in differential correction. The functionality has mostly been used for simulations, but it has also been applied for the orbit determination of the Lunar Reconaissance Orbiter using laser ranging data, in a joint DLR/TUD project, as well aspects of Doppler/VLBI data analysis of Mars and Venus Express in a joint JIVE/TUD project. It is currently being used in preparations for the PRIDE experiment on ESA’s JUICE missions.
* Mission design. Tudat provides various capabilities to design spacecraft orbits, and includes (among others) Lambert targetera, a patched conic (with DSMs) module, a CR3BP module, as well as the capability to make the link between these simplified model, and a full numerical simulation with user-defined propagation settings. Variational equations can be propagated for the full numerical propagation, allowing a differential correction between the preliminary design and full numerical problem to be easily achieved.

To support this functionality, Tudat includes a variety of mathematical tools (integrators, interpolators, root-finders, quadrature methods, etc.), as well as a library of basic astrodynamics functions, such as conversions between various element sets (Cartesian, Keplerian, Modified Equinoctial, Unified State Model, etc.).

Key aspects of the design of the Tudat software are:

* Modular software design: Models and methods can be easily exchanged, and the software is set up with the bare minimum of a priori assumptions on use cases or model values. This has allowed it to be applied to an extremely broad variety of topics, from space-plane ascent optimization, to orbit determination of Lunar spacecraft and the investigation of the dynamics of the Uranian satellite system. The broad applicability is both motivated and strengthened by the fact that the activities in the Tudat team cover many areas of both space science and engineering. This ensures that the code is critically reviewed in a cross-disciplinary manner, providing developers with unique opportunities for input from many (typically independent) communities.
* Rigorous unit testing. For each new piece of functionality that is added to the software, a unit test is added. This not only ensures that new pieces of functionality do not compromise old functionality, but it is also invaluable in testing the software on a broad range of operating systems/compilers. In case of unit test failures, clear output is provided, which is used to identify and solve the errors. 
* Use of numerous external libraries. Instead of trying to reinvent the wheel, Tudat makes use of existing software wherever possible. This includes the C++ boost libraries, which are nearly ubiquitous in modern computing. Also, we include the packages such as SPICE, SOFA and Pagmo2, developed by NASA, IAU and ESA, respectively. Our setup requires that only a single repository be pulled from github, where all libraries are automatically configured when first building the software. This prevents users from having to manually retrieve, build and link many packages, and allows direct access to their full functionality.

Tudat - SOCIS 2019 Projects
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Below, we describe the three project ideas for SOCIS 2019 in some detail. Interested students are encouraged to contact us with their questions and ideas for these projects.

SOCIS 2019 - Mission Design Module
**********************************

As noted above, Tudat includes a number of mission design features, including Lambert targeters, patched conic trajectory design, CR3BP orbit design, as well as full numerical propagation models. Separate modules currenty exist to link these simplified mission design features to the full numerical propagation. However, several key points are presently missing from the toolbox (in order of importance):

* A unified interface that allows users to combine the different trajectory models into a full mission scenario. To design a full mission, users must currently write a substantial amount of code. A key part of this project is to create a single interface that allows this process to be much more automated, including the capabilities to use the simple model as an initial guess for the full numerical propagation, and then performing a differential correction and providing the resulting Delta V.
* The design of low-thrust trajectories. Tudat includes numerous thrust models in the full numerical propagation. However, the mission design capabilities for low-thrust trajectories currently require a substantial amount of code from the user. We foresee that this project will develop interfaces for one or more low-thrust mission design approaches, including the link to the full numerical propagation.
* Several efforts have been made to create a periodic orbit module for Tudat. However, none of these has found its way into the repository, due to a number issues in this code. For this project, either the existing periodic orbit code could be updated for inclusion into Tudat, or a new module could be developed.

The scope of this project is quite large, and we encourage interested students to identify specific aspects that they would like to focus on.

SOCIS 2019 - Optimization Module
********************************

Tudat currently includes an automatic link with ESA's Pagmo2 toolbox. Numerous users have used Pagmo2 (and in the past Pagmo) to optimize space mission design problems using Tudat to compute the performance metrics. Although our repository includes a number of example applications on how to link Pagmo with Tudat, setting up an optimization problem can still be a time-consuming task. For this project, a new interface layer will be created to automate the link between Tudat and optimization toolboxes (primarily Pagmo). Based on settings for their mission scenario, such as:

* Types of legs in the mission (numerical, lambert, low-thrust, ...). The level of detail of this component will depend on the previous project. If this project is done in SOCIS2019, but not the previous one, the optimization layer can still be fully implemented for use with numerical propagation, with general interfaces to allow future extensions to be incorporated into the optimization layer.
* Detailed settings for any numerical propagations (environment/acceleration/propagator/integrator settings).
* Definitions of constraints and objectives, either custom user-defined, or derived from mission performance. For example, the maximum stagnation point heat flux over a re-entry, the closest approach to a body for an interplanetary transfer, or a maximum dynamic pressure for an ascent vehicle, can all be computed using Tudat. 
* Definition of decision variables. A custom option should be available here, but should also include options for e.g. initial conditions, thrust profiles, vehicle orientation profiles, etc.

The scope of this project is quite large, and we encourage interested students to identify specific aspects that they would like to focus on.

SOCIS 2019 - Interplanetary Doppler Data Analysis
*************************************************

Tudat contains a wealth of functionality in terms of orbit determination and parameter estimation. What it does not include, however, is a number of detailed models that are required for the processing of Doppler data (such as that obtained by the Deep Space Network or ESTRACK) for interplanetary missions. For this SOCIS project, the models required for this functionality will be developed further. This can include:

* File reading of DSN tracking data binary files, and conversion of file contents to Tudat-compatible data formats
* Models for ionospheric/tropospheric/solar corona influence on the Doppler shift of a radio signal
* Selected models for ground station properties and position variations (Earth deformation, etc.). Models for this, compatible with an older version of Tudat, are available, but not yet implemented into our repository.
* Definition of a time-varying reference point on a spacecraft (e.g. pointable antenna), including the reading of such data from Spice Kernels.











