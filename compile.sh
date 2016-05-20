#!/bin/bash

echo "initpackage ismdust lmodel_ismdust.dat "`pwd`"\nquit\ny" | xspec

rm *~ *.o
rm *FunctionMap.* lpack_*
rm -f *.mod Makefile
