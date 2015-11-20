#!/bin/bash

INFILE=$1

echo "----- R CMD buildBinary ($(date))          -----"      
R --version                                                  
echo "----- R CMD buildBinary STARTING ($(date)) -----"      
R CMD INSTALL $INFILE --build 
echo "----- R CMD buildBinary COMPLETE ($(date)) -----"      

