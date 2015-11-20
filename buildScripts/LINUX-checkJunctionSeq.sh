#!/bin/bash

VERSIONNUM=$1

echo "----- R CMD CHECK ($(date)) ----- "        
R --version                                      
echo "----- R CMD CHECK STARTING ($(date)) -----"
R CMD check ./JunctionSeq_$VERSIONNUM.tar.gz --no-build-vignettes
echo "----- R CMD CHECK COMPLETE ($(date)) -----"

