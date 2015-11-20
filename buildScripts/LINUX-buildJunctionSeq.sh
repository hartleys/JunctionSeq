#!/bin/bash

R_JS_DIR=$1

echo "----- R CMD build ($(date))          -----"
R --version                                      
echo "----- R CMD build STARTING ($(date)) -----"
R CMD build $R_JS_DIR                            
echo "----- R CMD build COMPLETE ($(date)) -----"

