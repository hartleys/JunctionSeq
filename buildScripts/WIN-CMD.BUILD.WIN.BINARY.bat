:: This ".bat" file is a DOS script that builds JunctionSeq in windows.
::      It also compiles a windows binary file.

@ECHO OFF
SET /P VER=Please enter the version number: 
IF "%VER%"=="" GOTO Error
ECHO Running tests for version= %VER%

:: cd C:\Users\hartleys\work\nihwork\home_copy\projects\ZZZ-JunctionSeq\releases\v%VER%\checks
:: go to the directory containing this .bat file:
cd %~dp0

if NOT exist ../JunctionSeq_%VER%.tar.gz (
  GOTO Error
)

mkdir InstBinary

:: --------------------------------------------------------------------------------------------------------------
:: Build a windows binary:
echo ----- R CMD buildBinary STARTING (%TIME%, %DATE%) ----- 
echo ----- R CMD buildBinary STARTING (%TIME%, %DATE%) -----  >  ../buildLogs/R-CMD-buildBinary-WIN64-R-v3.2.2.log
C:/PROGRA~1/R/R-3.2.2/bin/x64/R --version                     >> ../buildLogs/R-CMD-buildBinary-WIN64-R-v3.2.2.log 2>&1
echo --------------------------------------------------------------        >> ../buildLogs/R-CMD-buildBinary-WIN64-R-v3.2.2.log
C:/PROGRA~1/R/R-3.2.2/bin/x64/R CMD INSTALL --build -l "InstBinary" "../JunctionSeq_%VER%.tar.gz" >> ../buildLogs/R-CMD-buildBinary-WIN64-R-v3.2.2.log
echo ----- R CMD buildBinary COMPLETE (%TIME%, %DATE%)        >> ../buildLogs/R-CMD-buildBinary-WIN64-R-v3.2.2.log
echo ----- R CMD buildBinary COMPLETE (%TIME%, %DATE%) -----

GOTO End
:Error
ECHO Invalid version! Bye bye!!
:End

SET /P VER=Script complete. Hit enter to end.


