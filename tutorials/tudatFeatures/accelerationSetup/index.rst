.. _tudatFeaturesAccelerationIndex:

Acceleration Set-Up
===================

These pages of the wiki will help you build a strong knowledge basis to get started with Tudat. It is mandatory to understand the concepts taught here before proceding any further. Below a graphical representation of the acceleration setup is shown after which a list of content is shown. 

.. graphviz::

   digraph
   {
      # General diagram settings
      rankdir = "TB";
      splines = ortho;    
      compound = true;   


      # general node settings 
      node [shape = box, style = filled, width = 1.5, fixedsize = true, color = lightgrey, fontname = FontAwesome, fontsize = 9];


      # specific node color settings
      NamedBodyMap, bodiesToPropagate, centralBodies, AccelerationSettings [color = lightblue];
      AccelerationModel [color= lightgreen];


      subgraph clusterAccelerationType
      {
         label = AccelerationType;
         fontsize = 9;
         style = dashed;
      
         {rank = same; "Undefined acceleration", "Central gravity", Aerodynamic, "Third body \ncentral gravity"};
         {rank = same; "Third body spherical \nharmonic gravity", "Thrust", "Cannon ball \nradiation pressure", "Spherical harmonic \ngravity" };
         {rank = same; "Third body mutual \nspherical harmonic gravity", "Mutual spherical \nharmonic gravity", "Relativistic correction", Empirical};
         
         "Third body \ncentral gravity" -> "Third body spherical \nharmonic gravity" [style = invis];
         "Third body spherical \nharmonic gravity" -> "Third body mutual \nspherical harmonic gravity" [style = invis];
         "Central gravity" -> "Third body \ncentral gravity" [style = invis];
         "Spherical harmonic \ngravity" -> "Third body spherical \nharmonic gravity" [style = invis];
         "Undefined acceleration" -> "Central gravity" [style = invis];
         "Thrust" -> "Spherical harmonic \ngravity" [style = invis];
      }
          

      # AccelerationSettings input
      "Additional \ninformation" -> AccelerationSettings;
      "Third body mutual \nspherical harmonic gravity" -> AccelerationSettings [ltail = clusterAccelerationType];    

      
      # AccelerationModel input
      AccelerationSettings -> AccelerationModel;
      bodiesToPropagate -> AccelerationModel;
      centralBodies -> AccelerationModel;
      NamedBodyMap -> AccelerationModel;
     

      # Structure the layout
      {rank = same; NamedBodyMap, AccelerationModel, centralBodies, bodiesToPropagate};
      {rank = same; AccelerationSettings, "Additional \ninformation"};

      # Hyperlinks (Sphinx auto referencing not working here, need to link to exact web adres)
      "AccelerationSettings" [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/accelerationSetup/frameworkAcceleration.html#AccelerationSettings", target = "_top"];
      "Central gravity" [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/accelerationSetup/frameworkAcceleration.html#AccelerationSettings", target = "_top"];
      "AccelerationModel" [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/accelerationSetup/frameworkAcceleration.html##SelectedAccelerationMap", target = "_top"];
      NamedBodyMap [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/environmentSetup/index.html#NamedBodyMap", target = "_top"];
      "Undefined acceleration" [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/accelerationSetup/frameworkAcceleration.html#AccelerationSettings", target = "_top"];
      "Spherical harmonic \ngravity" [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/accelerationSetup/frameworkAcceleration.html#SphericalHarmonicAccelerationSettings", target = "_top"];
      "Mutual spherical \nharmonic gravity" [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/accelerationSetup/frameworkAcceleration.html#MutualSphericalHarmonicAccelerationSettings", target = "_top"];
      "Aerodynamic" [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/accelerationSetup/frameworkAcceleration.html#MutualSphericalHarmonicAccelerationSettings", target = "_top"];
      "Cannon ball \nradiation pressure" [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/accelerationSetup/frameworkAcceleration.html#MutualSphericalHarmonicAccelerationSettings", target = "_top"];
      "Thrust" [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/accelerationSetup/frameworkAcceleration.html#ThrustAccelerationSettings", target = "_top"];
      "Relativistic correction" [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/accelerationSetup/frameworkAcceleration.html#RelativisticAccelerationCorrectionSettings", target = "_top"];
      "Empirical" [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/accelerationSetup/frameworkAcceleration.html#EmpiricalAccelerationSettings", target = "_top"];

   }

.. toctree::

   frameworkAcceleration
   aerodynamicGuidance
   thrustModels
