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
      Empirical, "Relativistic correction", "Mutual spherical \nharmonic gravity", "Tidal Dissipation", "Thrust" [color =  darkturquoise]; 
      "Additional \ninformation" [style = dotted, fillcolor = lightgrey, color = black];

      subgraph clusterAccelerationType
      {
         label = AccelerationType;
         fontsize = 9;
         style = dashed;
      
         {rank = same; "Cannon ball \nradiation pressure", "Central gravity", Aerodynamic};
         {rank = same; "Thrust", "Tidal Dissipation", "Spherical harmonic \ngravity" };
         {rank = same; "Mutual spherical \nharmonic gravity", "Relativistic correction", Empirical};
         
         "Mutual spherical \nharmonic gravity" -> "Spherical harmonic \ngravity" -> "Central gravity" [style = invis];
         Empirical -> Thrust -> Aerodynamic [style = invis];
         "Cannon ball \nradiation pressure" -> "Central gravity" [style = invis];

      }
      
      Aerodynamic -> "Additional \ninformation"-> NamedBodyMap [style = invis];
      "Mutual spherical \nharmonic gravity" -> bodiesToPropagate [style = invis];
      AccelerationSettings -> bodiesToPropagate [style = invis];

      # AccelerationSettings input
      "Additional \ninformation" -> AccelerationSettings;
      "Central gravity" -> AccelerationSettings [ltail = clusterAccelerationType];    

      
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
      "Spherical harmonic \ngravity" [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/accelerationSetup/frameworkAcceleration.html#SphericalHarmonicAccelerationSettings", target = "_top"];
      "Mutual spherical \nharmonic gravity" [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/accelerationSetup/frameworkAcceleration.html#MutualSphericalHarmonicAccelerationSettings", target = "_top"];
      "Aerodynamic" [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/accelerationSetup/frameworkAcceleration.html#MutualSphericalHarmonicAccelerationSettings", target = "_top"];
      "Cannon ball \nradiation pressure" [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/accelerationSetup/frameworkAcceleration.html#MutualSphericalHarmonicAccelerationSettings", target = "_top"];
      "Thrust" [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/accelerationSetup/frameworkAcceleration.html#ThrustAccelerationSettings", target = "_top"];
      "Tidal Dissipation" [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/accelerationSetup/frameworkAcceleration.html#DissipationAccelerationSettings", target = "_top"];
      "Relativistic correction" [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/accelerationSetup/frameworkAcceleration.html#RelativisticAccelerationCorrectionSettings", target = "_top"];
      "Empirical" [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/accelerationSetup/frameworkAcceleration.html#EmpiricalAccelerationSettings", target = "_top"];

   }

.. graphviz::

   digraph
   {
      # General diagram settings
      rankdir = "LR";
      splines = ortho;    
      compound = true;  

      subgraph clusterLegend
      {
      rank = min;
      style = dashed;


     	# general node settings 
     	node [shape = box, style = filled, width = 1.25, fixedsize = true, color = lightgrey, fontname = FontAwesome, fontsize = 9];


   	"Object requiring \nadditional information" [ fillcolor = darkturquoise];
     	"Main block" [fillcolor = lightgreen];
     	"Optional input" [style = dotted, fillcolor = lightgrey, color = black];
     	"Input for \nmain block" [fillcolor = lightblue];
     	"Optional input"-> "Object requiring \nadditional information" -> "Input for \nmain block" -> "Main block" [style = invis];
      }
   }
.. toctree::

   frameworkAcceleration
   rotModel
   massModel
   aerodynamicGuidance
   thrustModels
   
