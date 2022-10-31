KPL/FK
 
   FILE: LATs.tf
 
   This file was created by PINPOINT.
 
   PINPOINT Version 3.3.0 --- December 13, 2021
   PINPOINT RUN DATE/TIME:    2022-10-28T15:59:36
   PINPOINT DEFINITIONS FILE: LATs.def
   PINPOINT PCK FILE:         pck00010_mod.tpc
   PINPOINT SPK FILE:         LATs.bsp
 
   The input definitions file is appended to this
   file as a comment block.
 
 
   Body-name mapping follows:
 
\begindata
 
   NAIF_BODY_NAME                      += 'LAT90'
   NAIF_BODY_CODE                      += 699001
 
   NAIF_BODY_NAME                      += 'LAT80'
   NAIF_BODY_CODE                      += 699002
 
   NAIF_BODY_NAME                      += 'LAT60'
   NAIF_BODY_CODE                      += 699003
 
   NAIF_BODY_NAME                      += 'LAT0'
   NAIF_BODY_CODE                      += 699004
 
\begintext
 
 
   Reference frame specifications follow:
 
 
   Topocentric frame LAT90_TOPO
 
      The Z axis of this frame points toward the zenith.
      The X axis of this frame points North.
 
      Topocentric frame LAT90_TOPO is centered at the
      site LAT90, which has Cartesian coordinates
 
         X (km):                  0.4091127471325E-11
         Y (km):                  0.0000000000000E+00
         Z (km):                 -0.5436400000000E+05
 
      and planetodetic coordinates
 
         Longitude (deg):         0.0000000000000
         Latitude  (deg):       -90.0000000000000
         Altitude   (km):         0.4091127471325E-11
 
      These planetodetic coordinates are expressed relative to
      a reference spheroid having the dimensions
 
         Equatorial radius (km):  6.0268000000000E+04
         Polar radius      (km):  5.4364000000000E+04
 
      All of the above coordinates are relative to the frame IAU_SATURN.
 
 
\begindata
 
   FRAME_LAT90_TOPO                    =  1699001
   FRAME_1699001_NAME                  =  'LAT90_TOPO'
   FRAME_1699001_CLASS                 =  4
   FRAME_1699001_CLASS_ID              =  1699001
   FRAME_1699001_CENTER                =  699001
 
   OBJECT_699001_FRAME                 =  'LAT90_TOPO'
 
   TKFRAME_1699001_RELATIVE            =  'IAU_SATURN'
   TKFRAME_1699001_SPEC                =  'ANGLES'
   TKFRAME_1699001_UNITS               =  'DEGREES'
   TKFRAME_1699001_AXES                =  ( 3, 2, 3 )
   TKFRAME_1699001_ANGLES              =  ( -180.0000000000000,
                                            -180.0000000000000,
                                               0.0000000000000 )
 
 
\begintext
 
   Topocentric frame LAT80_TOPO
 
      The Z axis of this frame points toward the zenith.
      The X axis of this frame points North.
 
      Topocentric frame LAT80_TOPO is centered at the
      site LAT80, which has Cartesian coordinates
 
         X (km):                  0.1156213711807E+05
         Y (km):                  0.0000000000000E+00
         Z (km):                 -0.5335419759332E+05
 
      and planetodetic coordinates
 
         Longitude (deg):         0.0000000000000
         Latitude  (deg):       -80.0000000000000
         Altitude   (km):         0.0000000000000E+00
 
      These planetodetic coordinates are expressed relative to
      a reference spheroid having the dimensions
 
         Equatorial radius (km):  6.0268000000000E+04
         Polar radius      (km):  5.4364000000000E+04
 
      All of the above coordinates are relative to the frame IAU_SATURN.
 
 
\begindata
 
   FRAME_LAT80_TOPO                    =  1699002
   FRAME_1699002_NAME                  =  'LAT80_TOPO'
   FRAME_1699002_CLASS                 =  4
   FRAME_1699002_CLASS_ID              =  1699002
   FRAME_1699002_CENTER                =  699002
 
   OBJECT_699002_FRAME                 =  'LAT80_TOPO'
 
   TKFRAME_1699002_RELATIVE            =  'IAU_SATURN'
   TKFRAME_1699002_SPEC                =  'ANGLES'
   TKFRAME_1699002_UNITS               =  'DEGREES'
   TKFRAME_1699002_AXES                =  ( 3, 2, 3 )
   TKFRAME_1699002_ANGLES              =  (    0.0000000000000,
                                            -170.0000000000000,
                                             180.0000000000000 )
 
 
\begintext
 
   Topocentric frame LAT60_TOPO
 
      The Z axis of this frame points toward the zenith.
      The X axis of this frame points North.
 
      Topocentric frame LAT60_TOPO is centered at the
      site LAT60, which has Cartesian coordinates
 
         X (km):                  0.3248953362127E+05
         Y (km):                  0.0000000000000E+00
         Z (km):                 -0.4578817699992E+05
 
      and planetodetic coordinates
 
         Longitude (deg):         0.0000000000000
         Latitude  (deg):       -60.0000000000000
         Altitude   (km):         0.0000000000000E+00
 
      These planetodetic coordinates are expressed relative to
      a reference spheroid having the dimensions
 
         Equatorial radius (km):  6.0268000000000E+04
         Polar radius      (km):  5.4364000000000E+04
 
      All of the above coordinates are relative to the frame IAU_SATURN.
 
 
\begindata
 
   FRAME_LAT60_TOPO                    =  1699003
   FRAME_1699003_NAME                  =  'LAT60_TOPO'
   FRAME_1699003_CLASS                 =  4
   FRAME_1699003_CLASS_ID              =  1699003
   FRAME_1699003_CENTER                =  699003
 
   OBJECT_699003_FRAME                 =  'LAT60_TOPO'
 
   TKFRAME_1699003_RELATIVE            =  'IAU_SATURN'
   TKFRAME_1699003_SPEC                =  'ANGLES'
   TKFRAME_1699003_UNITS               =  'DEGREES'
   TKFRAME_1699003_AXES                =  ( 3, 2, 3 )
   TKFRAME_1699003_ANGLES              =  (    0.0000000000000,
                                            -150.0000000000000,
                                             180.0000000000000 )
 
 
\begintext
 
   Topocentric frame LAT0_TOPO
 
      The Z axis of this frame points toward the zenith.
      The X axis of this frame points North.
 
      Topocentric frame LAT0_TOPO is centered at the
      site LAT0, which has Cartesian coordinates
 
         X (km):                  0.6026800000000E+05
         Y (km):                  0.0000000000000E+00
         Z (km):                  0.0000000000000E+00
 
      and planetodetic coordinates
 
         Longitude (deg):         0.0000000000000
         Latitude  (deg):         0.0000000000000
         Altitude   (km):         0.0000000000000E+00
 
      These planetodetic coordinates are expressed relative to
      a reference spheroid having the dimensions
 
         Equatorial radius (km):  6.0268000000000E+04
         Polar radius      (km):  5.4364000000000E+04
 
      All of the above coordinates are relative to the frame IAU_SATURN.
 
 
\begindata
 
   FRAME_LAT0_TOPO                     =  1699004
   FRAME_1699004_NAME                  =  'LAT0_TOPO'
   FRAME_1699004_CLASS                 =  4
   FRAME_1699004_CLASS_ID              =  1699004
   FRAME_1699004_CENTER                =  699004
 
   OBJECT_699004_FRAME                 =  'LAT0_TOPO'
 
   TKFRAME_1699004_RELATIVE            =  'IAU_SATURN'
   TKFRAME_1699004_SPEC                =  'ANGLES'
   TKFRAME_1699004_UNITS               =  'DEGREES'
   TKFRAME_1699004_AXES                =  ( 3, 2, 3 )
   TKFRAME_1699004_ANGLES              =  (    0.0000000000000,
                                             -90.0000000000000,
                                             180.0000000000000 )
 
\begintext
 
 
Definitions file LATs.def
--------------------------------------------------------------------------------
 
begindata
 
   SITES         = ('LAT90'
                    'LAT80'
                    'LAT60'
                    'LAT0')
 
   LAT90_CENTER = 699
   LAT90_FRAME  = 'IAU_SATURN'
   LAT90_IDCODE = 699001
   LAT90_LATLON = ( -90, 0, 0)
   LAT90_BOUNDS = ( @200-JAN-1, @2100-JAN-1 )
   LAT90_UP     = 'Z'
   LAT90_NORTH  = 'X'
 
   LAT80_CENTER = 699
   LAT80_FRAME  = 'IAU_SATURN'
   LAT80_IDCODE = 699002
   LAT80_LATLON = ( -80, 0, 0)
   LAT80_BOUNDS = ( @200-JAN-1, @2100-JAN-1 )
   LAT80_UP     = 'Z'
   LAT80_NORTH  = 'X'
 
   LAT60_CENTER = 699
   LAT60_FRAME  = 'IAU_SATURN'
   LAT60_IDCODE = 699003
   LAT60_LATLON = ( -60, 0, 0)
   LAT60_BOUNDS = ( @200-JAN-1, @2100-JAN-1 )
   LAT60_UP     = 'Z'
   LAT60_NORTH  = 'X'
 
   LAT0_CENTER = 699
   LAT0_FRAME  = 'IAU_SATURN'
   LAT0_IDCODE = 699004
   LAT0_LATLON = ( 0, 0, 0)
   LAT0_BOUNDS = ( @200-JAN-1, @2100-JAN-1 )
   LAT0_UP     = 'Z'
   LAT0_NORTH  = 'X'
 
begintext
 
begintext
 
[End of definitions file]
 
