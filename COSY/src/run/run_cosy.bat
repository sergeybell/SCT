@echo off
echo ========================================
echo RUNNING: lattice.fox
echo ========================================
cd ..
cosy.exe lattice.fox
cd run
echo.
echo ========================================
echo RUNNING: Chrom.fox
echo ========================================
cd ..
cosy.exe Chrom.fox
cd run
echo.
