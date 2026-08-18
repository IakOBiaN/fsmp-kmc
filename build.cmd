@echo off
rem Build helper for Windows: runs make with the w64devkit toolchain, no
rem PATH changes needed. From the repository root:
rem
rem   build            engine binaries, fsmp.exe + pack.exe (make windows)
rem   build test       regression test suite (needs a Python 3)
rem   build bundle     release archives in dist\ (Studio bundle and -cli)
rem
rem The toolchain is looked up in %W64DEVKIT%, then C:\w64devkit. Python is
rem looked up in the GUI venv, then on PATH, then through the py launcher.
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

set "PY="
if exist "%~dp0gui\.venv\Scripts\python.exe" set "PY=%~dp0gui\.venv\Scripts\python.exe"
if not defined PY for %%P in (python.exe) do if not defined PY set "PY=%%~$PATH:P"
if not defined PY for %%P in (py.exe) do if not defined PY set "PY=%%~$PATH:P"

if not defined PY (
  if /i "%TARGET%"=="windows" goto :nopython
  echo No Python 3 found, and "%TARGET%" needs one.
  echo   install it from https://www.python.org/downloads/ ^(tick "Add to PATH"^),
  echo   or create the Studio environment: py -3 -m venv gui\.venv
  exit /b 1
)
:nopython

rem make runs recipes through sh, which eats backslashes: pass the path
rem with forward slashes
if defined PY set "PY=%PY:\=/%"

if defined PY (
  make -C "%~dp0." %TARGET% "PYTHON=%PY%"
) else (
  make -C "%~dp0." %TARGET%
)
