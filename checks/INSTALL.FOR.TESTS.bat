:: This ".bat" file is a DOS script that builds JunctionSeq in windows.
::      It also compiles a windows binary file.

@ECHO OFF
:: go to the directory containing this .bat file:
cd %~dp0

echo "INSTALLING EXDATA TO R3.2.2:"
C:/PROGRA~1/R/R-3.2.2/bin/x64/R CMD INSTALL ../helpDocs/install/JctSeqExData2_LATEST.tar.gz

echo "INSTALLING EXDATA TO R-DEVEL:"
C:/PROGRA~1/R/R-devel/bin/x64/R CMD INSTALL ../helpDocs/install/JctSeqExData2_LATEST.tar.gz
echo "DONE!"

SET /P VER=Script complete. Hit enter to end.


