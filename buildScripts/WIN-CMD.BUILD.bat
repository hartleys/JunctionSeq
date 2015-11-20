:: This ".bat" file is a DOS script that builds JunctionSeq in windows.
::      It also compiles a windows binary file.

@ECHO OFF
SET /P VER=Please enter the version number: 
IF "%VER%"=="" GOTO Error
ECHO Running tests for version= %VER%

:: cd C:\Users\hartleys\work\nihwork\home_copy\projects\ZZZ-JunctionSeq\releases\v%VER%\checks
:: go to the directory containing this .bat file:
cd %~dp0
mkdir altBuilds
cd altBuilds

:: if NOT exist ../../JunctionSeq_%VER%.tar.gz (
::   GOTO Error
:: )

xcopy /E /I /Q ..\..\JunctionSeq JunctionSeq

:: --------------------------------------------------------------------------------------------------------------
:: Build a package archive:
echo ----- R CMD build STARTING (%TIME%, %DATE%) -----
echo ----- R CMD build STARTING (%TIME%, %DATE%) -----         >  ../../buildLogs/R-CMD-build-WIN64-R-v3.2.2.log
C:/PROGRA~1/R/R-3.2.2/bin/x64/R --version                      >> ../../buildLogs/R-CMD-build-WIN64-R-v3.2.2.log 2>&1
echo --------------------------------------------------------------         >> ../../buildLogs/R-CMD-build-WIN64-R-v3.2.2.log
C:/PROGRA~1/R/R-3.2.2/bin/x64/R CMD build "JunctionSeq"        >> ../../buildLogs/R-CMD-build-WIN64-R-v3.2.2.log
echo ----- R CMD build COMPLETE (%TIME%, %DATE%) -----         >> ../../buildLogs/R-CMD-build-WIN64-R-v3.2.2.log
echo ----- R CMD build COMPLETE (%TIME%, %DATE%) -----

:: --------------------------------------------------------------------------------------------------------------
:: R CMD Check the windows-built version:
echo ----- R CMD CHECK STARTING (%TIME%, %DATE%) -----
echo ----- R CMD CHECK STARTING (%TIME%, %DATE%) -----         >  ../../checks/R-CMD-CHECK-WIN64_built-R-v3.2.2.log
C:/PROGRA~1/R/R-3.2.2/bin/x64/R --version                      >> ../../checks/R-CMD-CHECK-WIN64_built-R-v3.2.2.log 2>&1
echo --------------------------------------------------------------      >> ../../checks/R-CMD-CHECK-WIN64_built-R-v3.2.2.log
C:/PROGRA~1/R/R-3.2.2/bin/x64/R CMD check --no-build-vignettes "JunctionSeq_%VER%.tar.gz"        >> ../../checks/R-CMD-CHECK-WIN64_built-R-v3.2.2.log
echo ----- R CMD CHECK COMPLETE (%TIME%, %DATE%) -----         >> ../../checks/R-CMD-CHECK-WIN64_built-R-v3.2.2.log
echo ----- R CMD CHECK COMPLETE (%TIME%, %DATE%) -----

:: --------------------------------------------------------------------------------------------------------------
:: Perform BiocChecks (3.2.2):
echo Starting BiocCheck (DEVEL):
echo ----- R CMD BiocCheck STARTING (%TIME%, %DATE%) ----- 
echo ----- R CMD BiocCheck STARTING (%TIME%, %DATE%) -----  >   ../../checks/R-CMD-BiocCheck-WIN64-RDEVEL-v2015.11.12.log
C:/PROGRA~1/R/R-3.2.2/bin/x64/R --version                     >>  ../../checks/R-CMD-BiocCheck-WIN64-RDEVEL-v2015.11.12.log 2>&1
echo --------------------------------------------------------------------        >>  ../../checks/R-CMD-BiocCheck-WIN64-RDEVEL-v2015.11.12.log
C:/PROGRA~1/R/R-3.2.2/bin/x64/Rscript ../../buildScripts/WIN-BiocCheck.R altBuilds/JunctionSeq_%VER%.tar.gz %~dp0 ../checks/R-CMD-BiocCheck-WIN64_built-R-v3.2.2.log
echo ----- R CMD BiocCheck COMPLETE (%TIME%, %DATE%) -----  >>  ../../checks/R-CMD-BiocCheck-WIN64-RDEVEL-v2015.11.12.log
echo ----- R CMD BiocCheck COMPLETE (%TIME%, %DATE%) ----- 

GOTO End
:Error
ECHO Invalid version! Bye bye!!
:End

SET /P VER=Script complete. Hit enter to end.


