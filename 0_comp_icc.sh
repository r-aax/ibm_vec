#!/bin/sh

set -x

cd src

rm -f ../ibm
rm -f *.o

icc *.cpp -fopenmp -g -o ../ibm

cd ..
