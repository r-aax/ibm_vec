#!/bin/sh

set -x

cd src

rm -f ../ibm
rm -f *.o
rm -f ../*.gcno ../*.gcda

g++ *.cpp -g -pg -fprofile-arcs -ftest-coverage -o ../ibm

cd ..

./ibm

gprof ibm gmon.out -p > prof_p.txt
gprof ibm gmon.out -q > prof_q.txt
gprof ibm gmon.out -A > prof_A.txt
