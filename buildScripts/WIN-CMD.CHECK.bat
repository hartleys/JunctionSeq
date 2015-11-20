:: This ".bat" file is a DOS script that runs the built JunctionSeq archive package
::      through a series of tests and checks across multiple R distributions.
::      in a windows environment.
::      First through R v3.2.2, then through a recent R development version.


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

:: --------------------------------------------------------------------------------------------------------------
:: Check the build package archive (using the R development version):
echo ----- R-DEVEL CMD CHECK STARTING (%TIME%, %DATE%) ----- 
echo ----- R-DEVEL CMD CHECK STARTING (%TIME%, %DATE%) -----  >   ../checks/R-CMD-CHECK-WIN64-RDEVEL-v2015.11.12.log
C:/PROGRA~1/R/R-devel/bin/x64/R --version                     >>  ../checks/R-CMD-CHECK-WIN64-RDEVEL-v2015.11.12.log 2>&1
echo --------------------------------------------------------------------        >>  ../checks/R-CMD-CHECK-WIN64-RDEVEL-v2015.11.12.log
C:/PROGRA~1/R/R-devel/bin/x64/R CMD check --no-build-vignettes "../JunctionSeq_%VER%.tar.gz" >> ../checks/R-CMD-CHECK-WIN64-RDEVEL-v2015.11.12.log
echo ----- R-DEVEL CMD CHECK COMPLETE (%TIME%, %DATE%) -----  >>  ../checks/R-CMD-CHECK-WIN64-RDEVEL-v2015.11.12.log
echo ----- R-DEVEL CMD CHECK COMPLETE (%TIME%, %DATE%) ----- 

:: --------------------------------------------------------------------------------------------------------------
:: Check the build package archive (using R 3.2.2):
echo Starting R3.2.2 check:
echo ----- R CMD CHECK STARTING (%TIME%, %DATE%) ----- 
echo ----- R CMD CHECK STARTING (%TIME%, %DATE%) ----- >  ../checks/R-CMD-CHECK-WIN64-R-v3.2.2.log
C:/PROGRA~1/R/R-3.2.2/bin/x64/R --version              >> ../checks/R-CMD-CHECK-WIN64-R-v3.2.2.log 2>&1
echo -------------------------------------------------------------------- >>  ../checks/R-CMD-CHECK-WIN64-R-v3.2.2.log
C:/PROGRA~1/R/R-3.2.2/bin/x64/R CMD check --no-build-vignettes "../JunctionSeq_%VER%.tar.gz" >> ../checks/R-CMD-CHECK-WIN64-R-v3.2.2.log
echo ----- R CMD CHECK COMPLETE (%TIME%, %DATE%) -----  >> ../checks/R-CMD-CHECK-WIN64-R-v3.2.2.log
echo ----- R CMD CHECK COMPLETE (%TIME%, %DATE%) ----- 

:: --------------------------------------------------------------------------------------------------------------
:: Perform BiocChecks (RDEVEL):
echo Starting BiocCheck (DEVEL):
echo ----- R-DEVEL CMD BiocCheck STARTING (%TIME%, %DATE%) ----- 
echo ----- R-DEVEL CMD BiocCheck STARTING (%TIME%, %DATE%) -----  >   ../checks/R-CMD-BiocCheck-WIN64-RDEVEL-v2015.11.12.log
C:/PROGRA~1/R/R-devel/bin/x64/R --version                     >>  ../checks/R-CMD-BiocCheck-WIN64-RDEVEL-v2015.11.12.log 2>&1
echo --------------------------------------------------------------------        >>  ../checks/R-CMD-BiocCheck-WIN64-RDEVEL-v2015.11.12.log
C:/PROGRA~1/R/R-devel/bin/x64/Rscript WIN-BiocCheck.R ../JunctionSeq_%VER%.tar.gz %~dp0 ../checks/R-CMD-BiocCheck-WIN64-RDEVEL-v2015.11.12.log
echo ----- R-DEVEL CMD BiocCheck COMPLETE (%TIME%, %DATE%) -----  >>  ../checks/R-CMD-BiocCheck-WIN64-RDEVEL-v2015.11.12.log
echo ----- R-DEVEL CMD BiocCheck COMPLETE (%TIME%, %DATE%) ----- 

:: --------------------------------------------------------------------------------------------------------------
:: Perform BiocChecks (R v3.2.2)
echo Starting BiocCheck (3.2.2):
echo ----- R CMD BiocCheck STARTING (%TIME%, %DATE%) ----- 
echo ----- R CMD BiocCheck STARTING (%TIME%, %DATE%) ----- >  ../checks/R-CMD-BiocCheck-WIN64-R-v3.2.2.log
C:/PROGRA~1/R/R-3.2.2/bin/x64/R --version              >> ../checks/R-CMD-BiocCheck-WIN64-R-v3.2.2.log 2>&1
echo -------------------------------------------------------------------- >>  ../checks/R-CMD-BiocCheck-WIN64-R-v3.2.2.log
C:/PROGRA~1/R/R-3.2.2/bin/x64/Rscript WIN-BiocCheck.R ../JunctionSeq_%VER%.tar.gz %~dp0 ../checks/R-CMD-BiocCheck-WIN64-R-v3.2.2.log
echo ----- R CMD BiocCheck COMPLETE (%TIME%, %DATE%) -----  >> ../checks/R-CMD-BiocCheck-WIN64-R-v3.2.2.log
echo ----- R CMD BiocCheck COMPLETE (%TIME%, %DATE%) ----- 

GOTO End
:Error
ECHO Invalid version! Bye bye!!
:End

SET /P VER=Script complete. Hit enter to end.





