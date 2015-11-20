#!/bin/bash

VERSIONNUM=$1

echo "----- R CMD BiocCheck ($(date))          -----"
R --version                                      
echo "----- R CMD BiocCheck STARTING ($(date)) -----"
R CMD BiocCheck ./JctSeqExData2_$VERSIONNUM.tar.gz
echo "----- R CMD BiocCheck COMPLETE ($(date)) -----"

