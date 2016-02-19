#!/bin/sh
#
# Note: This example wrapper-script is designed for use on a linux machine
# runnung SGE. If you are not running your scripts on such an environment,
# or if you don't know what any of that means, then don't use this script.
#
###############################
# SGE settings Here
# Basically, if a line starts with "#$", then you can enter any
# qsub command line flags .  See the qsub(1) man page.
# Redirect the STDOUT and STDERR files:
#$ -o logs/$JOB_NAME_$JOB_ID.log
#$ -e logs/$JOB_NAME_$JOB_ID.log
# Do some validation checking, and bail on errors
#$ -w e
# Operate in the current directory
#$ -cwd
# End SGE Settings
###############################

# Usable variables: MYJOBTITLE, MYTEMPDIR, MYTEMPISILON, HOSTNAME, MYLOGFILE

echo "#############################################"
echo "### "
echo "### USER:     $USER"
echo "### JOB_ID:   $JOB_ID"
echo "### JOB_NAME: $JOB_NAME"
echo "### HOSTNAME: $HOSTNAME"
echo "### "
echo "### STARTING @ $(date)"
echo "#############################################"

#Selects my environment vars and stuff:
source ~/setallsl6.sh

Rscript --no-save --no-restore $1 $2 $3 $4 $5 $6 $7 $8 $9

echo "#############################################"
echo "### DONE     @ $(date)"
echo "#############################################"


