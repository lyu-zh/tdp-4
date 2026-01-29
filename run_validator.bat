@echo off
REM TSP Dual Extractor Validator Runner Script
REM Usage: .\run_validator.bat

echo === TSP Dual Extractor Validator Runner ===
echo.

REM Get script directory
set SCRIPT_DIR=%~dp0
cd /d "%SCRIPT_DIR%"

REM Check if Java is installed
where javac >nul 2>&1
if %errorlevel% neq 0 (
    echo [ERROR] javac not found. Please install Java JDK.
    pause
    exit /b 1
)

REM Set classpath
set GUROBI_JAR=D:\gurobi1300\win64\lib\gurobi.jar
set JAMA_JAR=%SCRIPT_DIR%Jama-1.0.3.jar
set SRC_DIR=%SCRIPT_DIR%src
set OUTPUT_DIR=%SCRIPT_DIR%output

REM Check if required files exist
if not exist "%GUROBI_JAR%" (
    echo [ERROR] gurobi.jar not found at: %GUROBI_JAR%
    pause
    exit /b 1
)

REM Create output directory if it doesn't exist
if not exist "%OUTPUT_DIR%" mkdir "%OUTPUT_DIR%"

REM Set classpath for compilation and execution
set CLASSPATH=%GUROBI_JAR%;%JAMA_JAR%;%SRC_DIR%

REM Compile Java file
echo Compiling validator...
javac -encoding UTF-8 -cp "%CLASSPATH%" -d "%SRC_DIR%" "%SRC_DIR%\TSPDualExtractorValidator.java"

if %errorlevel% neq 0 (
    echo [ERROR] Compilation failed!
    pause
    exit /b 1
)

echo Compilation successful!
echo.

REM Run the validator
echo Running validator...
echo.

java -Xmx4g -cp "%CLASSPATH%" TSPDualExtractorValidator

if %errorlevel% equ 0 (
    echo.
    echo Validator completed successfully!
    echo Check output file: %OUTPUT_DIR%\validation_result.txt
) else (
    echo.
    echo [ERROR] Validator execution failed!
    pause
    exit /b 1
)

pause
