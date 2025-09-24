@echo off
gfortran -O2 -o binary BBM_modules.f90 BBM_subroutines.f90 BBM_main.f90
if errorlevel 1 (
    echo Compilation failed. Exiting...
    pause
    exit /b
)

echo Solving FLOW
.\binary.exe
if errorlevel 1 (
    echo Execution failed.
    pause
    exit /b
)

echo Deleting Files
del *.mod *.exe
pause
