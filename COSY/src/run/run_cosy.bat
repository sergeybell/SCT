@echo off
echo ========================================
echo RUNNING: magnetic_mod_sext.fox
echo ========================================
cd ..
cosy.exe magnetic_mod_sext.fox
cd run
echo.
echo ========================================
echo RUNNING: mapping.fox
echo ========================================
cd ..
cosy.exe mapping.fox
cd run
echo.
