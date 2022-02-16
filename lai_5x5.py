#!/usr/bin/python

import numpy
import os

reg = numpy.loadtxt('reg_5x5', delimiter='_')
cnt = 1

for i in reg:
    print(cnt, i[0], i[1], i[2], i[3])

    print('./makereglai '
          + '25 ' + str(int(i[0])) + ' ' + str(int(i[1])) + ' ' + str(int(i[2])) + ' ' + str(int(i[3])) )
    os.system('./makereglai '
          + '25 ' + str(int(i[0])) + ' ' + str(int(i[1])) + ' ' + str(int(i[2])) + ' ' + str(int(i[3])) )
    cnt = cnt + 1
