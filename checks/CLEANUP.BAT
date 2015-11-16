@ECHO OFF

SET /P VER=Please enter the version number: 

if NOT exist JunctionSeq_%VER%.zip (
  GOTO Error
)


echo "PERFORMING CLEANUP AND COPYING FILES..."
:: go to the directory containing this .bat file:
cd %~dp0

copy /Y JunctionSeq_%VER%.zip "../helpDocs/install/JunctionSeq_LATEST_WIN.zip"

rmdir /s /q JunctionSeq.Rcheck
rmdir /s /q altBuilds
rmdir /s /q InstBinary

GOTO End
:Error
ECHO Invalid version! Bye bye!!
:End

SET /P VER=Script complete. Hit enter to end.
