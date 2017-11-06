.. _troubleshootingFrequentMistakesJson:

Frequently made  coding mistakes - Json
=======================================
This section sums op some of the most made mistakes leading to failures when loading a .json file.

Incorrect comma placement
~~~~~~~~~~~~~~~~~~~~~~~~~

   Comma placement in .json file is very exacting. The following code, for instance, is incorrect: ::
      
        "bodies": {
        "Sun": {
            "useDefaultSettings": true
        },
        "Earth": {
            "useDefaultSettings": true
        },
        "Moon": {
            "useDefaultSettings": true
        }
        "asterix": {

   This will give the following error::

        Parse error in file "/jsonFileDirectory/jsonFileName.json" at line 19, col 17.
        terminate called after throwing an instance of 'nlohmann::detail::parse_error'
        what():  [json.exception.parse_error.101] parse error at 426: syntax error - unexpected string literal; expected '}'

   Where the final line :literal:`"asterix": {` correspond to line 19 of the input file. 

   The error occurs because, when defining the list of :literal:`bodies`, each entry should be followed by a comma. Here, the comma after the :literal:`"Moon"` entry is missing, and should read: ::

         "Moon": {
            "useDefaultSettings": true
         },

   So, a comma needed to be added to fix the error. In the case where a comma too *many* was added, for instance: ::

        "bodies": {
        "Sun": {
            "useDefaultSettings": true
        },
        "Earth": {
            "useDefaultSettings": true
        },
        "Moon": {
            "useDefaultSettings": true
        },
        }

   A similar error will be given. The problem here is that, following the last body (here :literal:`"Moon"`), no comma should be given, and the entry should read: ::

         "Moon": {
            "useDefaultSettings": true
         },



