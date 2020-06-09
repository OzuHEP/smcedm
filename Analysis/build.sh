#!usr/bin/bash

COMPILER=$(root-config --cxx)
FLAGS=$(root-config --cflags --libs)

$COMPILER -g -O3 -Wall -Wextra -Wpedantic -o SMCEDM_Analysis SMCEDM_Analysis.cxx $FLAGS -lTreePlayer