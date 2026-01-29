@echo off
REM DRTest1 Runner Script (Windows Batch)
REM Usage: Double-click this file or run: .\run_DRTest1.bat
REM This script runs DRTest1.java in a separate window with HIGH priority
REM CPU usage will not be affected by minimizing the IDE

echo === DRTest1 Runner Script ===
echo This script will compile and run DRTest1.java
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

where java >nul 2>&1
if %errorlevel% neq 0 (
    echo [ERROR] java not found. Please install Java JDK.
    pause
    exit /b 1
)

REM Display Java version
echo Java Version:
java -version
echo.

REM Set classpath - try multiple possible locations for JAR files
set GUROBI_JAR=
set JAMA_JAR=

REM Try to find gurobi.jar
if exist "%SCRIPT_DIR%gurobi.jar" (
    set GUROBI_JAR=%SCRIPT_DIR%gurobi.jar
) else if exist "D:\gurobi1300\win64\lib\gurobi.jar" (
    set GUROBI_JAR=D:\gurobi1300\win64\lib\gurobi.jar
) else if exist "C:\gurobi\win64\lib\gurobi.jar" (
    set GUROBI_JAR=C:\gurobi\win64\lib\gurobi.jar
) else (
    echo [WARNING] gurobi.jar not found. Trying to continue...
    echo Please ensure gurobi.jar is in the classpath.
)

REM Try to find Jama-1.0.3.jar
if exist "%SCRIPT_DIR%Jama-1.0.3.jar" (
    set JAMA_JAR=%SCRIPT_DIR%Jama-1.0.3.jar
) else if exist "%SCRIPT_DIR%..\..\Jama-1.0.3.jar" (
    set JAMA_JAR=%SCRIPT_DIR%..\..\Jama-1.0.3.jar
) else (
    echo [WARNING] Jama-1.0.3.jar not found. Trying to continue...
    echo Please ensure Jama-1.0.3.jar is in the classpath.
)

REM Set directories
set SRC_DIR=%SCRIPT_DIR%src
set OUTPUT_DIR=%SCRIPT_DIR%output

REM Create output directory if it doesn't exist
if not exist "%OUTPUT_DIR%" mkdir "%OUTPUT_DIR%"
if not exist "%OUTPUT_DIR%\csv" mkdir "%OUTPUT_DIR%\csv"
if not exist "%OUTPUT_DIR%\csv\improved" mkdir "%OUTPUT_DIR%\csv\improved"
if not exist "%OUTPUT_DIR%\csv\assignment_dependent" mkdir "%OUTPUT_DIR%\csv\assignment_dependent"

REM Set classpath for compilation and execution
set CLASSPATH=%SRC_DIR%
if not "%GUROBI_JAR%"=="" set CLASSPATH=%CLASSPATH%;%GUROBI_JAR%
if not "%JAMA_JAR%"=="" set CLASSPATH=%CLASSPATH%;%JAMA_JAR%

echo Classpath: %CLASSPATH%
echo.

REM Compile Java files
echo Compiling Java files...
echo Compiling Instance.java...
javac -encoding UTF-8 -cp "%CLASSPATH%" -d "%SRC_DIR%" "%SRC_DIR%\Instance.java"
if %errorlevel% neq 0 (
    echo [ERROR] Failed to compile Instance.java
    pause
    exit /b 1
)

echo Compiling DistributionallyRobustAlgo.java...
javac -encoding UTF-8 -cp "%CLASSPATH%" -d "%SRC_DIR%" "%SRC_DIR%\DistributionallyRobustAlgo.java"
if %errorlevel% neq 0 (
    echo [ERROR] Failed to compile DistributionallyRobustAlgo.java
    pause
    exit /b 1
)

echo Compiling DRTest1.java...
javac -encoding UTF-8 -cp "%CLASSPATH%" -d "%SRC_DIR%" "%SRC_DIR%\DRTest1.java"
if %errorlevel% neq 0 (
    echo [ERROR] Failed to compile DRTest1.java
    pause
    exit /b 1
)

echo Compilation successful!
echo.

REM Check if required data directories exist
if not exist "%SCRIPT_DIR%Instances_new" (
    echo [WARNING] Instances_new directory not found
    echo Program may fail if required instance files are missing
)

if not exist "%SCRIPT_DIR%data\test" (
    echo [WARNING] data\test directory not found
    echo Program may fail if assignment-dependent model data is required
)

REM Create a simple batch file to run Java
set RUN_SCRIPT=%TEMP%\run_java_drtest1_%RANDOM%.bat
(
echo @echo off
echo title DRTest1 - High Priority
echo echo ========================================
echo echo DRTest1 - High Priority Mode
echo echo ========================================
echo echo.
echo echo Process Priority: HIGH
echo echo CPU Usage: Will maintain full CPU usage
echo echo Independent: Yes ^(continues even if IDE is minimized^)
echo echo Output Directory: %OUTPUT_DIR%\csv
echo echo.
echo echo You can minimize this window - process continues at full speed.
echo echo.
echo REM Run Java with high priority and increased memory
echo java -Xmx8g -cp "%CLASSPATH%" DRTest1
echo.
echo echo.
echo echo ========================================
echo echo DRTest1 completed!
echo echo Check output files in: %OUTPUT_DIR%\csv
echo echo ========================================
echo pause
) > "%RUN_SCRIPT%"

REM Start the Java program in a new window with HIGH priority
echo Starting DRTest1 in a separate window with HIGH priority...
echo.
echo IMPORTANT INFORMATION:
echo   - The Java process will run in a NEW independent window
echo   - Process priority is set to HIGH for maximum CPU usage
echo   - CPU usage will NOT be affected by minimizing the IDE
echo   - You can minimize or close Cursor - the process will continue running
echo   - You can minimize the Java window - it will keep running at full speed
echo   - The process maintains full CPU usage regardless of window state
echo   - Output files will be saved to: %OUTPUT_DIR%\csv
echo.

REM Start with HIGH priority
start "DRTest1 - High Priority" /HIGH cmd /k "%RUN_SCRIPT%"

REM Wait a moment, then try to set Java process priority to high using wmic
timeout /t 3 /nobreak >nul
for /f "tokens=2" %%a in ('tasklist /FI "IMAGENAME eq java.exe" /FO LIST ^| findstr "PID:"') do (
    wmic process where "ProcessId=%%a" call setpriority "high priority" >nul 2>&1
)

echo.
echo Java process started in a separate window!
echo.
echo Process Details:
echo   - Window Title: "DRTest1 - High Priority"
echo   - Priority: HIGH
echo   - Independent: Yes ^(continues even if Cursor/IDE is closed^)
echo   - Memory: 8GB maximum
echo   - Output: %OUTPUT_DIR%\csv
echo.
echo Monitoring:
echo   - Task Manager: Look for java.exe with high CPU usage
echo   - Output directory: %OUTPUT_DIR%\csv
echo   - Output files will be named: distributionally_robust_results_YYYYMMDD_HHMMSS.csv
echo.
echo To stop the process:
echo   - Task Manager ^(Ctrl+Shift+Esc^) ^> Find java.exe ^> End Task
echo   - Or close the Java window
echo.
echo You can now minimize or close this window.
echo The Java process will continue running independently.
timeout /t 5 >nul
