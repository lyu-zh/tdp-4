@echo off
REM TSP Dual Extractor Batch Runner Script (Windows Batch)
REM Usage: Double-click this file or run: .\run_tsp_dual_extractor_batch.bat

echo === TSP Dual Extractor Batch Runner Script ===
echo This script will process ALL dates in demand_matrix.csv
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

REM Set classpath
set GUROBI_JAR=D:\gurobi1300\win64\lib\gurobi.jar
set JAMA_JAR=%SCRIPT_DIR%Jama-1.0.3.jar
set SRC_DIR=%SCRIPT_DIR%src
set OUTPUT_DIR=%SCRIPT_DIR%output
set DATA_DIR=%SCRIPT_DIR%data

REM Check if required files exist
if not exist "%GUROBI_JAR%" (
    echo [ERROR] gurobi.jar not found at: %GUROBI_JAR%
    pause
    exit /b 1
)

if not exist "%JAMA_JAR%" (
    echo [WARNING] Jama-1.0.3.jar not found at: %JAMA_JAR%
    echo Continuing, but some dependencies may be missing...
)

REM Create output directory if it doesn't exist
if not exist "%OUTPUT_DIR%" mkdir "%OUTPUT_DIR%"

REM Set classpath for compilation and execution
set CLASSPATH=%GUROBI_JAR%;%JAMA_JAR%;%SRC_DIR%

REM Compile Java file
echo Compiling Java files...
javac -encoding UTF-8 -cp "%CLASSPATH%" -d "%SRC_DIR%" "%SRC_DIR%\TSPDualExtractorBatch.java"

if %errorlevel% neq 0 (
    echo [ERROR] Compilation failed! Please check error messages.
    pause
    exit /b 1
)

echo Compilation successful!
echo.

REM Check if data files exist
set POINTS_FILE=%DATA_DIR%\unique_coordinates_list.csv
set CENTERS_FILE=%DATA_DIR%\selected_centers_p200.csv
set DEMAND_MATRIX_FILE=%DATA_DIR%\demand_matrix.csv

if not exist "%POINTS_FILE%" (
    echo [WARNING] Data file not found: %POINTS_FILE%
    echo Program will use default path, may fail if file doesn't exist
)

if not exist "%CENTERS_FILE%" (
    echo [WARNING] Centers file not found: %CENTERS_FILE%
    echo Program will use default path, may fail if file doesn't exist
)

if not exist "%DEMAND_MATRIX_FILE%" (
    echo [WARNING] Demand matrix file not found: %DEMAND_MATRIX_FILE%
    echo Program will use default path, may fail if file doesn't exist
)

REM Run the program in a separate window with high priority
echo Starting TSP Dual Extractor Batch Processor in a separate window with HIGH priority...
echo.
echo IMPORTANT INFORMATION:
echo   - The Java process will run in a NEW independent window
echo   - Process priority is set to HIGH for maximum CPU usage
echo   - This will process ALL dates in demand_matrix.csv
echo   - Each date will generate a separate output file
echo   - You can minimize or close Cursor - the process will continue running
echo   - You can minimize the Java window - it will keep running at full speed
echo   - The process maintains 100%% CPU usage regardless of window state
echo   - Processing time depends on the number of dates in demand_matrix.csv
echo.

REM Create a simple batch file to run Java
set RUN_SCRIPT=%TEMP%\run_java_tsp_batch_%RANDOM%.bat
(
echo @echo off
echo title TSP Dual Extractor Batch - High Priority
echo echo ========================================
echo echo TSP Dual Extractor Batch Processor - High Priority Mode
echo echo ========================================
echo echo.
echo echo Process Priority: HIGH
echo echo CPU Usage: Will maintain 100%% usage
echo echo Processing: ALL dates in demand_matrix.csv
echo echo Output Directory: %OUTPUT_DIR%
echo echo.
echo echo You can minimize this window - process continues at full speed.
echo echo.
echo REM Run Java with high priority
echo java -Xmx4g -cp "%CLASSPATH%" TSPDualExtractorBatch
echo.
echo echo.
echo echo ========================================
echo echo Batch processing completed!
echo echo Output files: %OUTPUT_DIR%\travel_dist_dual_values_by_date\travel_dist_dual_values_p*_*.csv
echo echo ========================================
echo pause
) > "%RUN_SCRIPT%"

REM Start the Java program in a new window with HIGH priority
REM /HIGH sets the window priority class to high
REM This ensures the Java process gets maximum CPU time
start "TSP Dual Extractor Batch - High Priority" /HIGH cmd /k "%RUN_SCRIPT%"

REM Wait a moment, then try to set Java process priority to high using wmic
timeout /t 3 /nobreak >nul
for /f "tokens=2" %%a in ('tasklist /FI "IMAGENAME eq java.exe" /FO LIST ^| findstr "PID:"') do (
    wmic process where "ProcessId=%%a" call setpriority "high priority" >nul 2>&1
)

echo.
echo Java process started in a separate window!
echo.
echo Process Details:
echo   - Window Title: "TSP Dual Extractor Batch - High Priority"
echo   - Priority: HIGH
echo   - Independent: Yes (continues even if Cursor is closed)
echo   - Processing: All dates in demand_matrix.csv
echo.
echo Monitoring:
echo   - Task Manager: Look for java.exe with high CPU usage
echo   - Output directory: %OUTPUT_DIR%\travel_dist_dual_values_by_date
echo   - Output files will be named: travel_dist_dual_values_p200_YYYYMMDD.csv
echo.
echo To stop the process:
echo   - Task Manager (Ctrl+Shift+Esc) ^> Find java.exe ^> End Task
echo   - Or wait for the program to complete all dates
echo.
echo You can now minimize or close this window.
timeout /t 5 >nul
