#!/bin/sh

set -x

cd src

rm -f ../ibm
rm -f *.o

icc *.cpp -g -o ../ibm

cd ..
