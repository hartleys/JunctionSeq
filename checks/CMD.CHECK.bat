:: This ".bat" file is a DOS script that runs the built JunctionSeq archive package
::      through a series of tests and checks across multiple R distributions.
::      in a windows environment.
::      First through R v3.2.2, then through a recent R development version.


@ECHO OFF
SET /P VER=Please enter the version number: 
IF "%VER%"=="" GOTO Error
ECHO Running tests for version= %VER%

:: cd C:\Users\hartleys\work\nihwork\home_copy\projects\ZZZ-JunctionSeq\releases\v%VER%\checks
cd %~dp0

if NOT exist ../JunctionSeq_%VER%.tar.gz (
  GOTO Error
)

echo ----- R-DEVEL CMD CHECK STARTING (%TIME%, %DATE%) ----- 
echo ----- R-DEVEL CMD CHECK STARTING (%TIME%, %DATE%) -----  >   R-CMD-CHECK-WIN64-RDEVEL-v2015.11.12.log
C:/PROGRA~1/R/R-devel/bin/x64/R --version                     >>  R-CMD-CHECK-WIN64-RDEVEL-v2015.11.12.log 2>&1
echo --------------------------------------------------------------------        >>  R-CMD-CHECK-WIN64-RDEVEL-v2015.11.12.log
C:/PROGRA~1/R/R-devel/bin/x64/R CMD check --no-build-vignettes "../JunctionSeq_%VER%.tar.gz" >> R-CMD-CHECK-WIN64-RDEVEL-v2015.11.12.log
echo ----- R-DEVEL CMD CHECK COMPLETE (%TIME%, %DATE%) -----  >>  R-CMD-CHECK-WIN64-RDEVEL-v2015.11.12.log
echo ----- R-DEVEL CMD CHECK COMPLETE (%TIME%, %DATE%) ----- 

echo Starting R3.2.2 check:
echo ----- R CMD CHECK STARTING (%TIME%, %DATE%) ----- 
echo ----- R CMD CHECK STARTING (%TIME%, %DATE%) ----- >  R-CMD-CHECK-WIN64-R-v3.2.2.log
C:/PROGRA~1/R/R-3.2.2/bin/x64/R --version              >> R-CMD-CHECK-WIN64-R-v3.2.2.log 2>&1
echo -------------------------------------------------------------------- >>  R-CMD-CHECK-WIN64-R-v3.2.2.log
C:/PROGRA~1/R/R-3.2.2/bin/x64/R CMD check --no-build-vignettes "../JunctionSeq_%VER%.tar.gz" >> R-CMD-CHECK-WIN64-R-v3.2.2.log
echo ----- R CMD CHECK COMPLETE (%TIME%, %DATE%) -----  >> R-CMD-CHECK-WIN64-R-v3.2.2.log
echo ----- R CMD CHECK COMPLETE (%TIME%, %DATE%) ----- 

echo Starting BiocCheck:
C:/PROGRA~1/R/R-devel/bin/x64/Rscript checker.R %VER% RDEVEL-v2015.11.12 %~dp0 WIN64
C:/PROGRA~1/R/R-3.2.2/bin/x64/Rscript checker.R %VER% R-3.2.2 %~dp0 WIN64

:: echo ----------------------------------------------------------------------------------------------
:: 
:: echo ----- R-DEVEL CMD CHECK STARTING (%TIME%, %DATE%) ----- 
:: echo ----- R-DEVEL CMD CHECK STARTING (%TIME%, %DATE%) -----  >   R-CMD-CHECK-WIN32-RDEVEL-v2015.11.12.log
:: C:/PROGRA~1/R/R-devel/bin/i386/R --version                     >>  R-CMD-CHECK-WIN32-RDEVEL-v2015.11.12.log 2>&1
:: echo -------------------------------------------------        >>  R-CMD-CHECK-WIN32-RDEVEL-v2015.11.12.log
:: C:/PROGRA~1/R/R-devel/bin/i386/R CMD check --no-build-vignettes "../JunctionSeq_%VER%.tar.gz" >> R-CMD-CHECK-WIN32-RDEVEL-v2015.11.12.log
:: echo ----- R-DEVEL CMD CHECK COMPLETE (%TIME%, %DATE%) -----  >>  R-CMD-CHECK-WIN32-RDEVEL-v2015.11.12.log
:: echo ----- R-DEVEL CMD CHECK COMPLETE (%TIME%, %DATE%) ----- 
:: 
:: echo Starting R3.2.2 check:
:: echo ----- R CMD CHECK STARTING (%TIME%, %DATE%) ----- 
:: echo ----- R CMD CHECK STARTING (%TIME%, %DATE%) ----- >  R-CMD-CHECK-WIN32-R-v3.2.2.log
:: C:/PROGRA~1/R/R-3.2.2/bin/i386/R --version              >> R-CMD-CHECK-WIN32-R-v3.2.2.log 2>&1
:: echo ------------------------------------------------- >>  R-CMD-CHECK-WIN32-R-v3.2.2.log
:: C:/PROGRA~1/R/R-3.2.2/bin/i386/R CMD check --no-build-vignettes "../JunctionSeq_%VER%.tar.gz" >> R-CMD-CHECK-WIN32-R-v3.2.2.log
:: echo ----- R CMD CHECK COMPLETE (%TIME%, %DATE%) -----  >> R-CMD-CHECK-WIN32-R-v3.2.2.log
:: echo ----- R CMD CHECK COMPLETE (%TIME%, %DATE%) ----- 
:: 
:: echo Starting BiocCheck:
:: C:/PROGRA~1/R/R-devel/bin/i386/Rscript checker.R %VER% RDEVEL-v2015.11.12 %~dp0 WIN32
:: C:/PROGRA~1/R/R-3.2.2/bin/i386/Rscript checker.R %VER% R-3.2.2 %~dp0 WIN32


GOTO End
:Error
ECHO Invalid version! Bye bye!!
:End

SET /P VER=Script complete. Hit enter to end.





