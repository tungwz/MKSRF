#!/usr/bin/python

import numpy
import os

reg = numpy.loadtxt('reg_5x5')
cnt = 1

month = ["01","02","03","04","05","06","07","08","09","10","11","12"]

for i in reg:
    print(cnt, i[0], i[1], i[2], i[3])

    for mm in month:
        print('ncks'
              + ' -d lat,' + str(int(i[2])) + '.,' + str(int(i[0])) + '.'
              + ' -d lon,' + str(int(i[1])) + '.,' + str(int(i[3])) + '.'
              + ' wc2_30s_prec_' + mm + '.nc -o ~/tera02/mksrf/wc_5x5/RG_'
              + str(int(i[0])) + '_' + str(int(i[1])) + '_' + str(int(i[2])) + '_' + str(int(i[3]))
              + '.PREC'+mm+'.nc')
        print('ncks'
              + ' -d lat,' + str(int(i[2])) + '.,' + str(int(i[0])) + '.'
              + ' -d lon,' + str(int(i[1])) + '.,' + str(int(i[3])) + '.'
              + ' wc2_30s_tavg_' + mm + '.nc -o ~/tera02/mksrf/wc_5x5/RG_'
              + str(int(i[0])) + '_' + str(int(i[1])) + '_' + str(int(i[2])) + '_' + str(int(i[3]))
              + '.TAVG'+mm+'.nc')
        print('ncks'
              + ' -d lat,' + str(int(i[2])) + '.,' + str(int(i[0])) + '.'
              + ' -d lon,' + str(int(i[1])) + '.,' + str(int(i[3])) + '.'
              + ' wc2_30s_tmax_' + mm + '.nc -o ~/tera02/mksrf/wc_5x5/RG_'
              + str(int(i[0])) + '_' + str(int(i[1])) + '_' + str(int(i[2])) + '_' + str(int(i[3]))
              + '.TMAX'+mm+'.nc')
        print('ncks'
              + ' -d lat,' + str(int(i[2])) + '.,' + str(int(i[0])) + '.'
              + ' -d lon,' + str(int(i[1])) + '.,' + str(int(i[3])) + '.'
              + ' wc2_30s_tmin_' + mm + '.nc -o ~/tera02/mksrf/wc_5x5/RG_'
              + str(int(i[0])) + '_' + str(int(i[1])) + '_' + str(int(i[2])) + '_' + str(int(i[3]))
              + '.TMIN'+mm+'.nc')
        os.system('ncks'
              + ' -d lat,' + str(int(i[2])) + '.,' + str(int(i[0])) + '.'
              + ' -d lon,' + str(int(i[1])) + '.,' + str(int(i[3])) + '.'
              + ' wc2_30s_prec_' + mm + '.nc -o ~/tera02/mksrf/wc_5x5/RG_'
              + str(int(i[0])) + '_' + str(int(i[1])) + '_' + str(int(i[2])) + '_' + str(int(i[3]))
              + '.PREC'+mm+'.nc')
        os.system('ncks'
              + ' -d lat,' + str(int(i[2])) + '.,' + str(int(i[0])) + '.'
              + ' -d lon,' + str(int(i[1])) + '.,' + str(int(i[3])) + '.'
              + ' wc2_30s_tavg_' + mm + '.nc -o ~/tera02/mksrf/wc_5x5/RG_'
              + str(int(i[0])) + '_' + str(int(i[1])) + '_' + str(int(i[2])) + '_' + str(int(i[3]))
              + '.TAVG'+mm+'.nc')
        os.system('ncks'
              + ' -d lat,' + str(int(i[2])) + '.,' + str(int(i[0])) + '.'
              + ' -d lon,' + str(int(i[1])) + '.,' + str(int(i[3])) + '.'
              + ' wc2_30s_tmax_' + mm + '.nc -o ~/tera02/mksrf/wc_5x5/RG_'
              + str(int(i[0])) + '_' + str(int(i[1])) + '_' + str(int(i[2])) + '_' + str(int(i[3]))
              + '.TMAX'+mm+'.nc')
        os.system('ncks'
              + ' -d lat,' + str(int(i[2])) + '.,' + str(int(i[0])) + '.'
              + ' -d lon,' + str(int(i[1])) + '.,' + str(int(i[3])) + '.'
              + ' wc2_30s_tmin_' + mm + '.nc -o ~/tera02/mksrf/wc_5x5/RG_'
              + str(int(i[0])) + '_' + str(int(i[1])) + '_' + str(int(i[2])) + '_' + str(int(i[3]))
              + '.TMIN'+mm+'.nc')

    cnt = cnt + 1
