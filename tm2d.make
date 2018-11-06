#!/bin/bash


includedir=`nf-config --fflags`
libdir=`nc-config --flibs`
gfortran -O3 -o tm2d tm2d_v0.9.0.f90 ${includedir} ${libdir}

# for debugging:
# gfortran -fbounds-check -g -Wextra -Wall -pedantic -o tm2d085 tm2d_v0.8.5.f90 -I/home/marco/prog/include -L/home/marco/prog/lib -lnetcdff

