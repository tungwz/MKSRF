#!/bin/bash

#- compile:
gfortran -c precision.F90
gfortran -I /usr/include/ -c ncio.F90
gfortran -g -I /usr/include/ -fopenmp -c MakeRegionSurface.F90
gfortran -g -I /usr/include/ -fopenmp -c MakeGlobalSurface.F90
gfortran -o makeregionsurface -fopenmp precision.o ncio.o MakeRegionSurface.o -lnetcdf -lnetcdff
gfortran -o makeglobalsurface -fopenmp MakeGlobalSurface.o -lnetcdf -lnetcdff

#- run:
#./makeregionsurface
#./makeglobalsurface 2005
