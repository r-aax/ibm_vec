#!/bin/sh

cd src

rm -f ../ibm
rm -f *.o

g++ *.cpp -o ../ibm

cd ..
