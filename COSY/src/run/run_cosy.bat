@echo off
echo ========================================
echo RUNNING: magnetic_sextupoles_maps.fox
echo ========================================
cd ..
cosy.exe magnetic_sextupoles_maps.fox
cd run
echo.
echo ========================================
echo RUNNING: Twiss.fox
echo ========================================
cd ..
cosy.exe Twiss.fox
cd run
echo.
