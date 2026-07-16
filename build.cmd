@echo off
rem Build helper for Windows: runs make with the w64devkit toolchain, no
rem PATH changes needed. From the repository root:
rem
rem   build            engine binaries, fsmp.exe + pack.exe (make windows)
rem   build test       regression test suite (uses the GUI venv python)
rem   build bundle     release-layout bundle in dist\
rem
rem The toolchain is looked up in %W64DEVKIT%, then C:\w64devkit.
setlocal
set "W64=%W64DEVKIT%"
if "%W64%"=="" set "W64=C:\w64devkit"
if not exist "%W64%\bin\make.exe" (
  echo w64devkit not found at %W64% - set W64DEVKIT or install to C:\w64devkit
  echo   https://github.com/skeeto/w64devkit
  exit /b 1
)
set "PATH=%W64%\bin;%PATH%"
set "TARGET=%~1"
if "%TARGET%"=="" set "TARGET=windows"
set "PY=%~dp0gui\.venv\Scripts\python.exe"
rem make runs recipes through sh, which eats backslashes: pass the path
rem with forward slashes
set "PY=%PY:\=/%"
if exist "%PY%" (
  make -C "%~dp0." %TARGET% "PYTHON=%PY%"
) else (
  make -C "%~dp0." %TARGET%
)
