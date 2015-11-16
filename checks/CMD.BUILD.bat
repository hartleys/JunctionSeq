:: This ".bat" file is a DOS script that builds JunctionSeq in windows.
::      It also compiles a windows binary file.

@ECHO OFF
SET /P VER=Please enter the version number: 
IF "%VER%"=="" GOTO Error
ECHO Running tests for version= %VER%

:: cd C:\Users\hartleys\work\nihwork\home_copy\projects\ZZZ-JunctionSeq\releases\v%VER%\checks
cd %~dp0
mkdir altBuilds
cd altBuilds

if NOT exist ../../JunctionSeq_%VER%.tar.gz (
  GOTO Error
)

xcopy /E /I ..\..\JunctionSeq JunctionSeq

echo ----- R CMD build STARTING (%TIME%, %DATE%) -----
echo ----- R CMD build STARTING (%TIME%, %DATE%) -----         >  ../R-CMD-build-WIN64-3.2.2.log
C:/PROGRA~1/R/R-3.2.2/bin/x64/R --version                      >> ../R-CMD-build-WIN64-3.2.2.log 2>&1
echo --------------------------------------------------------------         >> ../R-CMD-build-WIN64-3.2.2.log
C:/PROGRA~1/R/R-3.2.2/bin/x64/R CMD build "JunctionSeq"        >> ../R-CMD-build-WIN64-3.2.2.log
echo ----- R CMD build COMPLETE (%TIME%, %DATE%) -----         >> ../R-CMD-build-WIN64-3.2.2.log
echo ----- R CMD build COMPLETE (%TIME%, %DATE%) -----

mkdir InstBinary

echo ----- R CMD buildBinary STARTING (%TIME%, %DATE%) ----- 
echo ----- R CMD buildBinary STARTING (%TIME%, %DATE%) -----  >  ../R-CMD-buildBinary-WIN64-3.2.2.log
C:/PROGRA~1/R/R-3.2.2/bin/x64/R --version                     >> ../R-CMD-buildBinary-WIN64-3.2.2.log 2>&1
echo --------------------------------------------------------------        >> ../R-CMD-buildBinary-WIN64-3.2.2.log
C:/PROGRA~1/R/R-3.2.2/bin/x64/R CMD INSTALL --build -l "InstBinary" "JunctionSeq_%VER%.tar.gz" >> ../R-CMD-buildBinary-WIN64-3.2.2.log
echo ----- R CMD buildBinary COMPLETE (%TIME%, %DATE%)        >> ../R-CMD-buildBinary-WIN64-3.2.2.log
echo ----- R CMD buildBinary COMPLETE (%TIME%, %DATE%) -----

copy /Y JunctionSeq_%VER%.zip "../JunctionSeq_%VER%_WIN64-R3.2.2-BINARY.zip"

echo ----- R CMD CHECK STARTING (%TIME%, %DATE%) -----
echo ----- R CMD CHECK STARTING (%TIME%, %DATE%) -----         >  ../R-CMD-CHECK-WIN64-built-3.2.2.log
C:/PROGRA~1/R/R-3.2.2/bin/x64/R --version                      >> ../R-CMD-CHECK-WIN64-built-3.2.2.log 2>&1
echo --------------------------------------------------------------      >> ../R-CMD-build-WIN64-built-3.2.2.log
C:/PROGRA~1/R/R-3.2.2/bin/x64/R CMD build "JunctionSeq_%VER%.zip"        >> ../R-CMD-CHECK-WIN64-built-3.2.2.log
echo ----- R CMD CHECK COMPLETE (%TIME%, %DATE%) -----         >> ../R-CMD-CHECK-WIN64-built-3.2.2.log
echo ----- R CMD CHECK COMPLETE (%TIME%, %DATE%) -----


GOTO End
:Error
ECHO Invalid version! Bye bye!!
:End

SET /P VER=Script complete. Hit enter to end.


