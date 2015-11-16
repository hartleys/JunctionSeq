@ECHO OFF

echo "PERFORMING CLEANUP AND COPYING FILES..."

copy /Y altBuilds/JunctionSeq_0.6.1.zip "../helpDocs/install/JunctionSeq_LATEST_WIN.zip"

rmdir /s /q JunctionSeq.Rcheck
rmdir /s /q altBuilds

SET /P VER=Script complete. Hit enter to end.
