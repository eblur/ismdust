#!/bin/bash

ln -sf lmodel_ismdust.dat lmodel.dat

echo "initpackage ismdust lmodel.dat "`pwd`"\nquit\ny" | xspec

rm *~ *.o
rm *FunctionMap.* lpack_*
rm -f *.mod Makefile
