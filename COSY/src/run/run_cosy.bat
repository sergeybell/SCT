@echo off
echo ========================================
echo RUNNING: electrostatic.fox
echo ========================================
cd ..
cosy.exe electrostatic.fox
cd run
echo.
echo ========================================
echo RUNNING: mapping.fox
echo ========================================
cd ..
cosy.exe mapping.fox
cd run
echo.
