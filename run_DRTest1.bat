@echo off
REM Run DRTest1.java without IDE. Usage: run_DRTest1.bat from project root.

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

if defined GUROBI_JAR (set "GUROBI_JAR=%GUROBI_JAR%") else set "GUROBI_JAR=D:\gurobi1300\win64\lib\gurobi.jar"
if defined GUROBI_HOME (
    if exist "%GUROBI_HOME%\lib\gurobi.jar" set "GUROBI_JAR=%GUROBI_HOME%\lib\gurobi.jar"
)
if defined COPT_JAR (set "COPT_JAR=%COPT_JAR%") else set "COPT_JAR=D:\Program Files\copt80\lib\copt.jar"
if defined COPT_HOME (
    if exist "%COPT_HOME%\lib\copt.jar" set "COPT_JAR=%COPT_HOME%\lib\copt.jar"
)
set "JAMA_JAR=%SCRIPT_DIR%Jama-1.0.3.jar"
set "SRC=%SCRIPT_DIR%src"
set "OUT=%SCRIPT_DIR%out\production\TDP"

if not exist "%GUROBI_JAR%" (
    echo [ERROR] gurobi.jar not found: %GUROBI_JAR%
    echo Set env GUROBI_JAR or GUROBI_HOME, or edit this script.
    pause
    exit /b 1
)
if not exist "%JAMA_JAR%" (
    echo [ERROR] Jama-1.0.3.jar not found: %JAMA_JAR%
    pause
    exit /b 1
)
if not exist "%COPT_JAR%" (
    echo [ERROR] copt.jar not found: %COPT_JAR%
    echo Set env COPT_JAR or COPT_HOME, or edit this script. D1-SDP needs COPT.
    pause
    exit /b 1
)

if not exist "%OUT%" mkdir "%OUT%"

REM Strip UTF-8 BOM from sources if present (javac reports illegal character \ufeff)
where python >nul 2>&1
if not errorlevel 1 (
  for %%F in ("%SRC%\DistributionallyRobustAlgo.java" "%SRC%\DRTest1.java") do (
    python -c "import pathlib; p=pathlib.Path(r'%%~fF'); d=p.read_bytes(); p.write_bytes(d[3:] if d[:3]==bytes([0xEF,0xBB,0xBF]) else d)" 2>nul
  )
)

echo Compiling...
set "CP=%OUT%;%GUROBI_JAR%;%JAMA_JAR%;%COPT_JAR%"
javac -encoding UTF-8 -cp "%CP%" -sourcepath "%SRC%" -d "%OUT%" "%SRC%\DRTest1.java"
if errorlevel 1 (
    echo [ERROR] Compilation failed.
    pause
    exit /b 1
)
echo Compile OK.
echo.

echo Starting DRTest1 in a new window...
set "RUN_SCRIPT=%TEMP%\run_drtest1_%RANDOM%.bat"
(
echo @echo off
echo cd /d "%SCRIPT_DIR%"
echo java -cp "%OUT%;%GUROBI_JAR%;%JAMA_JAR%;%COPT_JAR%" DRTest1
echo.
echo echo Done.
echo pause
) > "%RUN_SCRIPT%"
start "DRTest1" cmd /k call "%RUN_SCRIPT%"
echo New window opened. Output will appear there.
exit /b 0