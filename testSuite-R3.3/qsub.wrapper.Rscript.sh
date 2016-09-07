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

export LANG=C

source /etc/profile.d/modules.sh
module use /home/hartleys/modules
module use /opt/sw/modules

module load all/Java/1.7.0_60
module load sl6/libffi
module load sl6/glib
module load sl6/pixman
module load sl6/cairo
module load sl6/fontconfig
module load sl6/harfbuzz
module load sl6/freetype
module load sl6/pango
module load libs/zlib-cpath
module load sl6/bzip2
module load sl6/liblzma
module load sl6/pcre
module load sl6/libcurl
module load sl6/texinfo
module load sl6/libtiff
module load compiler/gcc/4.8.5
module load sl6/R/3.3.0

Rscript --no-save --no-restore $1 $2 $3 $4 $5 $6 $7 $8 $9

echo "#############################################"
echo "### DONE     @ $(date)"
echo "#############################################"


