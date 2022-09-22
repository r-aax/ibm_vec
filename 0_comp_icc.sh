#!/bin/sh

set -x

cd src

rm -f ../ibm
rm -f *.o

icc -xmic-avx512 *.cpp -fopenmp -g -o ../ibm

cd ..
