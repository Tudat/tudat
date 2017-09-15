
.. role:: jsontype
.. role:: jsonkey
.. role:: arrow

- :jsontype:`string` :jsonkey:`model` (mandatory). Identifier of the built-in spherical harmonics gravity field model to be used. Models currently available: `EGM96 <https://en.wikipedia.org/wiki/EGM96>`_, `GGM02C <http://www2.csr.utexas.edu/grace/gravity/ggm02/>`_, `GGM02S <http://www2.csr.utexas.edu/grace/gravity/ggm02/>`_ (Earth); `GLGM3150 <https://pds.nasa.gov/ds-view/pds/viewProfile.jsp?dsid=LP-L-RSS-5-GLGM3/GRAVITY-V1.0>`_, `LPE200 <https://pds.nasa.gov/ds-view/pds/viewProfile.jsp?dsid=LP-L-RSS-5-GLGM3/GRAVITY-V1.0>`_ (Moon); `JGMRO120D <FIXME>`_ (Mars). Possible values: :literal:`"egm96"`, :literal:`"ggm02c"`, :literal:`"ggm02s"`, :literal:`"glgm3150"`, :literal:`"lpe200"`, :literal:`"jgmro120d"`.
