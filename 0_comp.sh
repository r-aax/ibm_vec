#!/bin/sh

set -x

cd src

rm -f ../ibm
rm -f *.o

g++ *.cpp -g -o ../ibm

cd ..
