#!/bin/bash

gfortran -O3 -o ../tm2d tm2d_v0.8.6.f90 -I/storage/home/rohrer/prog/include -L/storage/home/rohrer/prog/lib -lnetcdff

# for debugging:
# gfortran -fbounds-check -g -Wextra -Wall -pedantic -o tm2d085 tm2d_v0.8.5.f90 -I/home/marco/prog/include -L/home/marco/prog/lib -lnetcdff

