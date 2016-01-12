@ECHO OFF

SET /P VER=Please enter the OLD version number: 

echo "COPYING OLD FILES..."
:: go to the directory containing this .bat file:
cd %~dp0

mkdir archive\v%VER%
xcopy /E /I /Q doc archive\v%VER%\doc
xcopy /E /I /Q install archive\v%VER%\install
xcopy /E /I /Q javascripts archive\v%VER%\javascripts
xcopy /E /I /Q Rhtml archive\v%VER%\Rhtml
xcopy * archive\v%VER%\

xcopy /E /I /Q archive\v0.5.1\stylesheets archive\v%VER%\stylesheets

SET /P VER=Script complete. Hit enter to end.
