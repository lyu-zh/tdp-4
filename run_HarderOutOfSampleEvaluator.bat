@echo off
REM Run HarderOutOfSampleEvaluator.java without IDE.

setlocal
set "SCRIPT_DIR=%~dp0"
cd /d "%SCRIPT_DIR%"

where javac >nul 2>&1
if errorlevel 1 (
    echo [ERROR] javac not found. Install JDK.
    pause
    exit /b 1
)
where java >nul 2>&1
if errorlevel 1 (
    echo [ERROR] java not found. Install JDK.
    pause
    exit /b 1
)

set "SRC=%SCRIPT_DIR%src"
set "OUT=%SCRIPT_DIR%out\production\TDP"

if not exist "%OUT%" mkdir "%OUT%"

echo Compiling...
javac -encoding UTF-8 -sourcepath "%SRC%" -d "%OUT%" "%SRC%\HarderOutOfSampleEvaluator.java"
if errorlevel 1 (
    echo [ERROR] Compilation failed.
    pause
    exit /b 1
)
echo Compile OK.
echo.

echo Starting HarderOutOfSampleEvaluator in a new window...
set "RUN_SCRIPT=%TEMP%\run_harder_oos_%RANDOM%.bat"
(
echo @echo off
echo cd /d "%SCRIPT_DIR%"
echo java -cp "%OUT%" HarderOutOfSampleEvaluator
echo.
echo echo Done.
echo pause
) > "%RUN_SCRIPT%"

start "HarderOutOfSampleEvaluator" cmd /k call "%RUN_SCRIPT%"
echo New window opened. Output will appear there.
exit /b 0
