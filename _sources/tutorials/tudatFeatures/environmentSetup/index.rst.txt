.. _tudatFeaturesEnvironmentIndex:

Environment Set-up
==================
Models for the physical environment are one of the cornerstones of a numerical astrodynamics toolbox. Here, we define the environment in the broadest sense, including all physical properties of the solar system, such as atmospheres and gravity fields, but also any models for the orbits and rotations of these bodies.

On this page, we will give an overview of how the environment is represented in Tudat, which models have been implemented, and how to create the environment that is tailored to your needs. A graphical representation of the basic structure for implementing the environment is shown below. For the options within these settings object, please click on the object in the figure.

.. graphviz::

   digraph
   {
      # General diagram settings
      rankdir = "LR";
      splines = ortho;    
      compound = true;  


      # general node settings 
      node [shape = box, style = filled, width = 1.25, fixedsize = true, color = lightgrey, fontname = FontAwesome, fontsize = 9];


      # specific node color settings
      NamedBodyMap [color = lightgreen];
      BodySettings [color = lightblue];


      # Hyperlinks (Sphinx auto referencing not working here, need to link to exact web adres)
      "GravityField\nVariationSettings" [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/environmentSetup/index.html#GravityFieldVariationSettings", target = "_top"];
      "RotationModelSettings" [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/environmentSetup/index.html#RotationModelSettings", target = "_top"];
      "Aerodynamic\nCoefficientSettings" [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/environmentSetup/index.html#AerodynamicCoefficientSettings", target = "_top"];
      "BodyShapeSettings" [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/environmentSetup/index.html#BodyShapeSettings", target = "_top"];
      "RadiationPressure\nInterfaceSettings" [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/environmentSetup/index.html#RadiationPressureInterfaceSettings", target = "_top"];
      "AtmosphereSettings" [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/environmentSetup/index.html#AtmosphereSettings", target = "_top"];
      "GravityFieldSettings" [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/environmentSetup/index.html#GravityFieldSettings", target = "_top"];
      "EphemerisSettings" [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/environmentSetup/index.html#EphemerisSettings", target = "_top"];
      BodySettings [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/environmentSetup/index.html#BodySettings", target = "_top"];
      NamedBodyMap [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/environmentSetup/index.html#NamedBodyMap", target = "_top"];


      # NamedBodyMap input
      BodySettings -> NamedBodyMap;
      "User-defined \nbodies" [style = dotted, fillcolor = lightgrey, color = black];
      "User-defined \nbodies" -> NamedBodyMap;
      {rank = same; BodySettings, NamedBodyMap};


      # BodySettings input
      EphemerisSettings -> BodySettings [ltail = clusterEnvironmentSettings];
      getDefaultBodySettings -> BodySettings;
      

      # getDefaultBodySettings input
      bodyNames -> getDefaultBodySettings;
      "(initial/final) time" -> getDefaultBodySettings;
      "(initial/final) time" [style = dotted, fillcolor = lightgrey, color = black];

		#point0 [shape = point, style = vis, width = 0.1];
      #point0 -> "setGlobalFrame\nBodyEphemerides";
      #point0 -> "NamedBodyMap"; 
      #BodySettings -> point0;

      # Cluster all environment settings derived classes
      subgraph clusterEnvironmentSettings
      {
         # cluster settings
         label = "User-defined environment settings";
         fontsize = 9;
         style = dashed;
         rank = min;

         # Make three collumns
         {rank = same; AtmosphereSettings, EphemerisSettings, GravityFieldSettings};
         {rank = same; "GravityField\nVariationSettings", "RadiationPressure\nInterfaceSettings"};
         {rank = same; "Aerodynamic\nCoefficientSettings", BodyShapeSettings, RotationModelSettings};

         BodyShapeSettings -> EphemerisSettings [style = invis];
         "RadiationPressure\nInterfaceSettings" -> BodyShapeSettings [style = invis];
         "GravityField\nVariationSettings" -> RotationModelSettings [style = invis];
      }
      "setGlobalFrame\nBodyEphemerides" -> NamedBodyMap [dir = both];
      "Frame origin\n& orientation" -> "setGlobalFrame\nBodyEphemerides";
		
		
      {"RadiationPressure\nInterfaceSettings", "GravityField\nVariationSettings" [fillcolor = lightcoral]};
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


   	"List of settings" [ fillcolor = lightcoral];
     	"Main block" [fillcolor = lightgreen];
     	"Optional input" [style = dotted, fillcolor = lightgrey, color = black];
     	"Input for \nmain block" [fillcolor = lightblue];
     	"Optional input"-> "List of settings" -> "Input for \nmain block" -> "Main block" [style = invis];
      }
   }

.. toctree::

   setup
   availableSettings
   duringPropagation
   tabulatedAtmosphere
   aerodynamicCoefficients
