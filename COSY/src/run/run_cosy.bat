@echo off
echo ========================================
echo RUNNING: seq_Nuclotron_16.fox
echo ========================================
cd ..
cosy.exe seq_Nuclotron_16.fox
cd run
echo.
echo ========================================
echo RUNNING: mapping.fox
echo ========================================
cd ..
cosy.exe mapping.fox
cd run
echo.
