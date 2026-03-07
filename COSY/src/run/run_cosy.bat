@echo off
echo ========================================
echo RUNNING: QFS_deuteron_maps.fox
echo ========================================
cd ..
cosy.exe QFS_deuteron_maps.fox
cd run
echo.
echo ========================================
echo RUNNING: Twiss.fox
echo ========================================
cd ..
cosy.exe Twiss.fox
cd run
echo.
