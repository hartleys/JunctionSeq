#!/bin/bash

VERSIONNUM=$1

echo "----- R CMD build ($(date))          -----"
R --version                                      
echo "----- R CMD build STARTING ($(date)) -----"
R CMD check ./JctSeqData_$VERSIONNUM.tar.gz --no-build-vignettes
echo "----- R CMD build COMPLETE ($(date)) -----"

