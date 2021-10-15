#!/usr/bin/python

import numpy
import os

reg = numpy.loadtxt('reg_5x5_', delimiter='_')
cnt = 1

for i in reg:
    print cnt, i[0], i[1], i[2], i[3]
    
    print('./bin/aggregate-map.sh -PoutputLCCSClasses=false -PnumMajorityClasses=2 -PuserPFTConversionTable="./resources/Default_LCCS2PFT_LUT.csv" -PgridName=GEOGRAPHIC_LAT_LON -PnumRows=43200 -Pnorth=' + str(i[0]-0.00001) + ' -Peast=' + str(i[3]-0.00001) + ' -Psouth=' + str(i[2]) + ' -Pwest=' + str(i[1]) + ' ESACCI-LC-L4-LCCS-Map-300m-P1Y-2015-v2.0.7b.nc')
    os.system('./bin/aggregate-map.sh -PoutputLCCSClasses=false -PnumMajorityClasses=2 -PuserPFTConversionTable="./resources/Default_LCCS2PFT_LUT.csv" -PgridName=GEOGRAPHIC_LAT_LON -PnumRows=43200 -Pnorth=' + str(i[0]-0.00001) + ' -Peast=' + str(i[3]-0.00001) + ' -Psouth=' + str(i[2]) + ' -Pwest=' + str(i[1]) + ' ESACCI-LC-L4-LCCS-Map-300m-P1Y-2015-v2.0.7b.nc')

    print('mv ESACCI-LC-L4-LCCS-Map-300m-P1Y-aggregated-0.004167Deg-USER_REGION-2015-v2.0.7.nc '
          + 'RG_' + str(int(i[0])) + '_' + str(int(i[1])) + '_' + str(int(i[2])) + '_' + str(int(i[3])) + '.CCI2015.nc')
    os.system('mv ESACCI-LC-L4-LCCS-Map-300m-P1Y-aggregated-0.004167Deg-USER_REGION-2015-v2.0.7.nc '
          + 'RG_' + str(int(i[0])) + '_' + str(int(i[1])) + '_' + str(int(i[2])) + '_' + str(int(i[3])) + '.CCI2015.nc')
    cnt = cnt + 1
