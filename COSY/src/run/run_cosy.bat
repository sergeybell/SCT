@echo off
echo ========================================
echo RUNNING: magnetic_sextupoles.fox
echo ========================================
cd ..
cosy.exe magnetic_sextupoles.fox
cd run
echo.
echo ========================================
echo RUNNING: mapping.fox
echo ========================================
cd ..
cosy.exe mapping.fox
cd run
echo.
