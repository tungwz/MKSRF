#!/usr/bin/python

import numpy
import os

reg = numpy.loadtxt('reg_5x5', delimiter='_')
cnt = 1

for i in reg:
    print cnt, i[0], i[1], i[2], i[3]
    
    print('ncks'
          + ' -d lat,' + str(int(i[2])) + '.,' + str(int(i[0])) + '.'
          + ' -d lon,' + str(int(i[1])) + '.,' + str(int(i[3])) + '.'
          + ' Beck_KG_V1_present_0p0083.nc -o ~/hard/mksrf/kg_5x5/RG_'
          + str(int(i[0])) + '_' + str(int(i[1])) + '_' + str(int(i[2])) + '_' + str(int(i[3])) + '.KG.nc')

    os.system('ncks'
          + ' -d lat,' + str(int(i[2])) + '.,' + str(int(i[0])) + '.'
          + ' -d lon,' + str(int(i[1])) + '.,' + str(int(i[3])) + '.'
          + ' Beck_KG_V1_present_0p0083.nc -o ~/hard/mksrf/kg_5x5/RG_'
          + str(int(i[0])) + '_' + str(int(i[1])) + '_' + str(int(i[2])) + '_' + str(int(i[3])) + '.KG.nc')
    cnt = cnt + 1
