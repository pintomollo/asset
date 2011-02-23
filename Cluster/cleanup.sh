#!/bin/bash

rm -f bugs.txt

cat *.err > bugs.txt
cat ML-*-evol.dat > ML-results.txt
