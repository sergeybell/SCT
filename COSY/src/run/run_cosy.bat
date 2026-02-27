@echo off
echo ========================================
echo RUNNING: lattice_sext.fox
echo ========================================
cd ..
cosy.exe lattice_sext.fox
cd run
echo.
echo ========================================
echo RUNNING: mapping.fox
echo ========================================
cd ..
cosy.exe mapping.fox
cd run
echo.
