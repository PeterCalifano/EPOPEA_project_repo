KPL/FK
 
   FILE: kernels/dss.tf
 
   This file was created by PINPOINT.
 
   PINPOINT Version 3.2.0 --- September 6, 2016
   PINPOINT RUN DATE/TIME:    2021-09-28T17:11:55
   PINPOINT DEFINITIONS FILE: pinpoint/dss.def
   PINPOINT PCK FILE:         kernels/pck00010.tpc
   PINPOINT SPK FILE:         kernels/dss.spk
 
   The input definitions file is appended to this
   file as a comment block.
 
 
   Body-name mapping follows:
 
\begindata
 
   NAIF_BODY_NAME                      += 'DSS-12'
   NAIF_BODY_CODE                      += 399001
 
   NAIF_BODY_NAME                      += 'ESA-NNO'
   NAIF_BODY_CODE                      += 399002
 
\begintext
 
 
   Reference frame specifications follow:
 
 
   Topocentric frame DSS-12_TOPO
 
      The Z axis of this frame points toward the zenith.
      The X axis of this frame points North.
 
      Topocentric frame DSS-12_TOPO is centered at the
      site DSS-12, which has Cartesian coordinates
 
         X (km):                 -0.2350443812000E+04
         Y (km):                 -0.4651980837000E+04
         Z (km):                  0.3665630988000E+04
 
      and planetodetic coordinates
 
         Longitude (deg):      -116.8054869538837
         Latitude  (deg):        35.2999377199896
         Altitude   (km):         0.9625753575335E+00
 
      These planetodetic coordinates are expressed relative to
      a reference spheroid having the dimensions
 
         Equatorial radius (km):  6.3781366000000E+03
         Polar radius      (km):  6.3567519000000E+03
 
      All of the above coordinates are relative to the frame IAU_EARTH.
 
 
\begindata
 
   FRAME_DSS-12_TOPO                   =  1399001
   FRAME_1399001_NAME                  =  'DSS-12_TOPO'
   FRAME_1399001_CLASS                 =  4
   FRAME_1399001_CLASS_ID              =  1399001
   FRAME_1399001_CENTER                =  399001
 
   OBJECT_399001_FRAME                 =  'DSS-12_TOPO'
 
   TKFRAME_1399001_RELATIVE            =  'IAU_EARTH'
   TKFRAME_1399001_SPEC                =  'ANGLES'
   TKFRAME_1399001_UNITS               =  'DEGREES'
   TKFRAME_1399001_AXES                =  ( 3, 2, 3 )
   TKFRAME_1399001_ANGLES              =  ( -243.1945130461163,
                                             -54.7000622800104,
                                             180.0000000000000 )
 
 
\begintext
 
   Topocentric frame ESA-NNO_TOPO
 
      The Z axis of this frame points toward the zenith.
      The X axis of this frame points North.
 
      Topocentric frame ESA-NNO_TOPO is centered at the
      site ESA-NNO, which has Cartesian coordinates
 
         X (km):                 -0.2414024374456E+04
         Y (km):                  0.4907891255477E+04
         Z (km):                 -0.3270602751126E+04
 
      and planetodetic coordinates
 
         Longitude (deg):       116.1910000000000
         Latitude  (deg):       -31.0482000000000
         Altitude   (km):         0.2520000000003E+00
 
      These planetodetic coordinates are expressed relative to
      a reference spheroid having the dimensions
 
         Equatorial radius (km):  6.3781366000000E+03
         Polar radius      (km):  6.3567519000000E+03
 
      All of the above coordinates are relative to the frame IAU_EARTH.
 
 
\begindata
 
   FRAME_ESA-NNO_TOPO                  =  1399002
   FRAME_1399002_NAME                  =  'ESA-NNO_TOPO'
   FRAME_1399002_CLASS                 =  4
   FRAME_1399002_CLASS_ID              =  1399002
   FRAME_1399002_CENTER                =  399002
 
   OBJECT_399002_FRAME                 =  'ESA-NNO_TOPO'
 
   TKFRAME_1399002_RELATIVE            =  'IAU_EARTH'
   TKFRAME_1399002_SPEC                =  'ANGLES'
   TKFRAME_1399002_UNITS               =  'DEGREES'
   TKFRAME_1399002_AXES                =  ( 3, 2, 3 )
   TKFRAME_1399002_ANGLES              =  ( -116.1910000000000,
                                            -121.0482000000000,
                                             180.0000000000000 )
 
\begintext
 
 
Definitions file pinpoint/dss.def
--------------------------------------------------------------------------------
 
begindata
 
   SITES         = ( 'DSS-12',
                     'ESA-NNO' )
 
   DSS-12_CENTER  = 399
   DSS-12_FRAME   = 'IAU_EARTH'
   DSS-12_IDCODE  = 399001
   DSS-12_XYZ     = ( -2350.443812, -4651.980837, +3665.630988 )
   DSS-12_UP      = 'Z'
   DSS-12_NORTH   = 'X'
 
   ESA-NNO_CENTER = 399
   ESA-NNO_FRAME  = 'IAU_EARTH'
   ESA-NNO_IDCODE = 399002
   ESA-NNO_LATLON = ( -31.0482, 116.191, 0.252 )
   ESA-NNO_UP     = 'Z'
   ESA-NNO_NORTH  = 'X'
 
begintext
 
begintext
 
[End of definitions file]
 
