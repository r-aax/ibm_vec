#!/bin/sh

set -x

cd src

rm -f ../ibm
rm -f *.o

icc -O3 *.cpp -fopenmp -g -o ../ibm_0
icc -O3 -xmic-avx512 *.cpp -fopenmp -g -o ../ibm_1

cd ..

objdump -D ibm_0 > ibm_0.diss
objdump -D ibm_1 > ibm_1.diss
