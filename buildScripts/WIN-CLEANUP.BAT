@ECHO OFF

echo "PERFORMING CLEANUP AND COPYING FILES..."
:: go to the directory containing this .bat file:
cd %~dp0

copy /Y JunctionSeq_0.6.1.zip "../helpDocs/install/JunctionSeq_LATEST_WIN.zip"

rmdir /s /q JunctionSeq.Rcheck
rmdir /s /q altBuilds
rmdir /s /q InstBinary


SET /P VER=Script complete. Hit enter to end.
