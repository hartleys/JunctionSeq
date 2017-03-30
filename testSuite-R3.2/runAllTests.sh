#!/bin/bash

set +o posix

#This script is designed for use specifically on a SGE-style environment.
# It is not portable, and not intended for general use.

rm -rf out
mkdir out

rm -rf logs
mkdir logs

rm -rf status
mkdir status

while read TID warnExpected
do
  #Run on an SL6 node:
  qsub -q low6.q -P SL6 -l mem_free=1G,h_vmem=1G -j yes -N JST-$TID qsub.wrapper.Rscript.sh runTest.R test-$TID $warnExpected
  
  #Run on a older node:
  #qsub -l mem_free=1G,h_vmem=1G -j yes -N JST-$TID qsub.wrapper.Rscript.sh runTest.R test-$TID $warnExpected
  
done < <(tail -n+2 test.list.txt)
#(Note: the 'tail' statement skips the first line, which holds the column titles)

#If you don't use SGE, you could execute the tests like this:
#while read line
#do
#  qsub -l mem_free=1G,h_vmem=1G -j yes -N JST-$line qsub.wrapper.Rscript.sh runTest.R test-$line
#done < $(tail -n+2 test.list.txt)
