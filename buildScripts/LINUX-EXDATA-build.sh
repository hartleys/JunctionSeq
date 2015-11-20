#!/bin/bash

R_EX_DIR=$1

echo "----- R CMD build ($(date))          -----"
R --version                                      
echo "----- R CMD build STARTING ($(date)) -----"
R CMD build --resave-data $R_EX_DIR                            
echo "----- R CMD build COMPLETE ($(date)) -----"

