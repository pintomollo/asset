@echo off
rem MEX.BAT
rem
rem    Compile and link script used for building MEX-files
rem
rem ###################################################################
rem # Call Perl with the script name as an argument and pass inputs
rem ###################################################################
rem
SETLOCAL
set PERLLIB=
set PERL5LIB=
set PERL5OPT=
set PERL5SHELL=
"%~dp0\..\sys\perl\win32\bin\perl.exe" ".\mymex.pl" %*
ENDLOCAL
cmd /c exit %errorlevel%


