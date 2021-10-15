#!/usr/bin/python

import os

lat_max = 90.
lon_min = -180.
lat_min = -90.
lon_max = 180.
dll = 5.

lat = lat_max
while lat > -90:
    print('processing LAT:{:10.4f}...'.format(lat))
    lon = lon_min
    while lon < 180:
        print('./makereg 12 '
                  +str(lat)+' '
                  +str(lon)+' '
                  +str(lat-dll)+' '
                  +str(lon+dll) )
        os.system('./makereg 12 '
                  +str(lat)+' '
                  +str(lon)+' '
                  +str(lat-dll)+' '
                  +str(lon+dll) )
        lon = lon + dll
    lat = lat - dll    
