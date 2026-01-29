
% DistributionallyRobustMain.m
function DistributionallyRobustMain()
    clc;clear;
    % Start timer for overall execution
    totalStartTime = tic;
    timerActive = false;
    timeoutTimer = [];
    
    try
        % Set parameters
        instanceFile = './Instances/2DU60-05-1.dat';  % Input file path
        outputFileName = 'test_drcc';           % Output file name
        outputDir = './output/';
        
        % Set distributionally robust optimization parameters
        gamma = 0.2;      % Risk parameter
        numScenarios = 1000; % Number of scenarios
        seed = 12345678;   % Random seed
        
        % Demand parameters
        E = 50.0;          % Expected value
        RSDValues = [0.125, 0.25, 0.5]; % Relative standard deviation array
        
        % Select specific RSD value
        RSD = RSDValues(1);
        r = 0.2;
        % Other DRCC parameters
        useD1 = true;      % Whether to use D1 fuzzy set
        delta1 = 1;        % D2 fuzzy set parameter
        delta2 = 2;        % D2 fuzzy set parameter
        useJointChance = false; % Whether to use joint chance constraint
        
        % Algorithm timeout (in seconds)
        timeout = 3600; % 1 hour timeout
        
        % Display configuration
        % displayConfig(instanceFile, outputFileName, gamma, numScenarios, E, RSD, useD1, delta1, delta2, useJointChance);
        
        % Validate parameters
        % validateParams(gamma, numScenarios, E, RSD, delta1, delta2);
        
        % Ensure output directory exists
        if ~exist(outputDir, 'dir')
            [success, msg] = mkdir(outputDir);
            if ~success
                error('Failed to create directory %s: %s', outputDir, msg);
            end
            fprintf('Created output directory: %s\n', outputDir);
        end
        
        % Check if instance file exists
        if ~exist(instanceFile, 'file')
            error('Instance file not found: %s', instanceFile);
        end
        
        % Load instance
        fprintf('Loading instance from %s...\n', instanceFile);
        instance = Instance(instanceFile);
        
        if isempty(instance) || instance.getN() == 0
            error('Failed to load instance or instance has no areas');
        end
        
        % fprintf('Instance contains %d areas\n', instance.getN());
        
        % Generate random scenarios using vectorized operations
        % fprintf('Generating %d scenarios...\n', numScenarios);
        scenarios = generateScenarios(instance.getN(), numScenarios, E, RSD, seed);
        
        % Create and run the distributionally robust algorithm
        % fprintf('Initializing algorithm...\n');
        algo = DistributionallyRobustAlgo(instance, scenarios, gamma, seed, useD1, delta1, delta2, useJointChance,r);
        
        if isempty(algo)
            error('Failed to create algorithm object');
        end
        
        % Try to set algorithm timeout
        try
            if isfield(algo, 'timeLimit') || isprop(algo, 'timeLimit')
                algo.timeLimit = timeout;
                % fprintf('Algorithm timeout set to %d seconds\n', timeout);
            end
        catch
            fprintf('Warning: Could not set algorithm timeout\n');
        end
        
        % Run algorithm with timeout monitoring
        % fprintf('\nRunning algorithm...\n');
        algoRunStartTime = tic;
        
        % Try to create a timer for timeout monitoring
        try
            if exist('timer', 'file') && exist('timerfind', 'file')
                % Delete any existing timers with the same name
                existingTimers = timerfind('Name', 'TimeoutTimer');
                if ~isempty(existingTimers)
                    delete(existingTimers);
                end
                
                % Create a new timer
                timeoutTimer = timer('TimerFcn', @(~,~) timeoutCheck(algoRunStartTime, timeout), ...
                       'ExecutionMode', 'fixedRate', ...
                       'Period', 10, ... % Check every 10 seconds
                       'Name', 'TimeoutTimer');
                start(timeoutTimer);
                timerActive = true;
            end
        catch timerError
            fprintf('Warning: Could not create timeout timer: %s\n', timerError.message);
        end
        
        % Run algorithm
        algo.run(outputFileName);
        
        % Stop the timer if it's active
        cleanupTimer();
        
        fprintf('Distributionally robust chance-constrained district design problem solved.\n');
        
        % Get results from the algorithm
        zones = algo.getZones();
        centers = algo.getCenters();
        % Perform out-of-sample performance evaluation
        fprintf('\nPerforming out-of-sample performance evaluation...\n');
        evaluateOutOfSamplePerformance(instance, zones, centers, E, RSD,r);

% Store performance results in output if needed
% You could write these results to a file if desired
        % Visualize results
        
        % Display solution statistics
        % Total runtime
        % totalRunTime = toc(totalStartTime);
        % fprintf('\nTotal run time: %.2f seconds (%.2f minutes)\n', ...
        %     totalRunTime, totalRunTime/60);
        
    catch e
        % Stop the timer if it's active
        cleanupTimer();
        
        fprintf('\n\n');
        fprintf('Error occurred: %s\n', e.message);
        fprintf('Stack trace:\n');
        for i = 1:length(e.stack)
            fprintf('  File: %s, Line: %d, Function: %s\n', ...
                e.stack(i).file, e.stack(i).line, e.stack(i).name);
        end
    end
    
    % Nested function to clean up timer
    function cleanupTimer()
        if timerActive && ~isempty(timeoutTimer) && isvalid(timeoutTimer)
            try
                stop(timeoutTimer);
                delete(timeoutTimer);
            catch
                % Ignore errors in timer cleanup
            end
            timerActive = false;
        end
    end
end

% Function to display configuration
function displayConfig(instanceFile, outputFileName, gamma, numScenarios, E, RSD, useD1, delta1, delta2, useJointChance)
    fprintf('----------------------------------------\n');
    fprintf('Distributionally Robust Optimization Configuration:\n');
    fprintf('----------------------------------------\n');
    fprintf('Instance file: %s\n', instanceFile);
    fprintf('Output file: %s\n', outputFileName);
    fprintf('Risk parameter (gamma): %.4f\n', gamma);
    fprintf('Number of scenarios: %d\n', numScenarios);
    fprintf('Expected demand value (E): %.2f\n', E);
    fprintf('Relative standard deviation (RSD): %.4f\n', RSD);
    
    if useD1
        fprintf('Fuzzy set type: D1\n');
    else
        fprintf('Fuzzy set type: D2\n');
        fprintf('Delta1 parameter: %.2f\n', delta1);
        fprintf('Delta2 parameter: %.2f\n', delta2);
    end
    
    if useJointChance
        fprintf('Constraint type: Joint chance constraint\n');
    else
        fprintf('Constraint type: Individual chance constraint\n');
    end
    fprintf('----------------------------------------\n\n');
end

% Function to validate parameters
function validateParams(gamma, numScenarios, E, RSD, delta1, delta2)
    % Check gamma (risk parameter)
    if ~isnumeric(gamma) || ~isscalar(gamma) || gamma <= 0 || gamma >= 1
        error('Parameter "gamma" must be a scalar between 0 and 1');
    end
    
    % Check numScenarios
    if ~isnumeric(numScenarios) || ~isscalar(numScenarios) || numScenarios <= 0 || mod(numScenarios, 1) ~= 0
        error('Parameter "numScenarios" must be a positive integer');
    end
    
    % Check E (expected value)
    if ~isnumeric(E) || ~isscalar(E) || E <= 0
        error('Parameter "E" must be a positive scalar');
    end
    
    % Check RSD (relative standard deviation)
    if ~isnumeric(RSD) || ~isscalar(RSD) || RSD <= 0
        error('Parameter "RSD" must be a positive scalar');
    end
    
    % Check delta1 and delta2 (D2 fuzzy set parameters)
    if ~isnumeric(delta1) || ~isscalar(delta1) || delta1 <= 0
        error('Parameter "delta1" must be a positive scalar');
    end
    
    if ~isnumeric(delta2) || ~isscalar(delta2) || delta2 <= 0
        error('Parameter "delta2" must be a positive scalar');
    end
    
    if delta2 <= delta1
        error('Parameter "delta2" must be greater than "delta1"');
    end
    
    % Check if lower bound would be negative
    lowerBound = E * (1 - sqrt(3) * RSD);
    if lowerBound <= 0
        warning(['With E=%.2f and RSD=%.4f, the calculated lower bound (%.2f) is non-positive. ' ...
                'This will be adjusted to 1 for demand generation.'], E, RSD, lowerBound);
    end
end

% Function to generate scenarios - vectorized version
function scenarios = generateScenarios(n, numScenarios, E, RSD, seed)
    % Set random seed for reproducibility
    rng(seed);
    
    % Calculate uniform distribution endpoints
    lowerBound = E * (1 - sqrt(3) * RSD);
    upperBound = E * (1 + sqrt(3) * RSD);
    
    % Handle potential negative lower bound
    if lowerBound <= 0
        warning('Lower bound for demand is non-positive. Setting to 1.');
        lowerBound = 1;
    end
    
    % Initialize scenarios matrix
    scenarios = zeros(numScenarios, n);
    
    % For very large instances, generate in chunks to save memory
    maxChunkSize = 1000; % Adjust based on available memory
    
    % If matrix is very large, generate in chunks
    if numScenarios * n > 10^7 % Threshold for chunk processing
        fprintf('Large instance detected. Generating scenarios in chunks...\n');
        
        for chunkStart = 1:maxChunkSize:numScenarios
            % Determine chunk end
            chunkEnd = min(chunkStart + maxChunkSize - 1, numScenarios);
            chunkSize = chunkEnd - chunkStart + 1;
            
            % Generate random values for this chunk
            randomValues = rand(chunkSize, n);
            
            % Scale to the desired range
            scenarios(chunkStart:chunkEnd, :) = lowerBound + randomValues * (upperBound - lowerBound);
            
            % Report progress
            fprintf('Generated scenarios %d to %d of %d\n', chunkStart, chunkEnd, numScenarios);
        end
    else
        % For smaller instances, generate all at once
        randomValues = rand(numScenarios, n);
        scenarios = lowerBound + randomValues * (upperBound - lowerBound);
    end
    
    % Ensure all demands are positive (at least 1)
    scenarios = max(1, scenarios);
end

% Function to check if algorithm execution has exceeded timeout
function timeoutCheck(startTime, timeout)
    elapsedTime = toc(startTime);
    if elapsedTime > timeout
        % This will be caught by the try-catch block in the main function
        error('Algorithm execution exceeded timeout of %d seconds', timeout);
    end
end


% Function to evaluate out-of-sample performance
function performanceResults = evaluateOutOfSamplePerformance(instance, zones, ~, E, RSD, r)
    % If r is not provided, use default value 0.1
    if nargin < 6
        r = 0.1;
    end

    % fprintf('\n----------------------------------------\n');
    % fprintf('Evaluating Out-of-Sample Performance:\n');
    % fprintf('----------------------------------------\n');
    % 
    % Parameters for out-of-sample testing
    numTestScenarios = 1000;
    testSeed = 54321; % Different seed for test data
    
    % Generate test scenarios
    % fprintf('Generating %d test scenarios...\n', numTestScenarios);
    testScenarios = generateScenarios(instance.getN(), numTestScenarios, E, RSD, testSeed);
    
    % Initialize counters for each district and overall performance
    districtViolationCounts = zeros(1, length(zones));
    scenariosSatisfied = 0;
    
    % Evaluate each test scenario
    % fprintf('Evaluating scenarios...\n');
    
    % Pre-compute zone assignment matrix (each column represents a zone)
    n = instance.getN();
    assignmentMatrix = zeros(n, length(zones));
    for i = 1:length(zones)
        if ~isempty(zones{i})
            assignmentMatrix(zones{i}, i) = 1;
        end
    end
    
    % For each scenario
    for s = 1:numTestScenarios
        % Calculate total demand for this specific scenario
        scenarioTotalDemand = sum(testScenarios(s, :));
        
        % Calculate demand upper bound for this specific scenario
        scenarioDemandUpperBound = (1 + r) * (scenarioTotalDemand / instance.k);
        
        % Calculate all zone demands in one matrix operation
        scenarioDemandVector = testScenarios(s, :);
        zoneDemands = scenarioDemandVector * assignmentMatrix;
        
        % Check violations for all zones at once
        violationMask = zoneDemands > scenarioDemandUpperBound;
        districtViolationCounts = districtViolationCounts + violationMask;
        
        % Scenario is satisfied if no violations
        scenarioSatisfied = ~any(violationMask);
        if scenarioSatisfied
            scenariosSatisfied = scenariosSatisfied + 1;
        end
    end
    
    % Calculate satisfaction rates
    districtSatisfactionRates = 1 - (districtViolationCounts / numTestScenarios);
    overallSatisfactionRate = scenariosSatisfied / numTestScenarios;
    
    % Output results
    % fprintf('\nOut-of-Sample Performance Results:\n');
    fprintf('  Overall satisfaction rate: %.4f\n', overallSatisfactionRate);
    performanceResults = overallSatisfactionRate;
    % fprintf('\nIndividual district satisfaction rates:\n');
    
    % for i = 1:length(zones)
    %     fprintf('  District %d: %.4f\n', i, districtSatisfactionRates(i));
    % end
    
    % fprintf('\nMin district satisfaction rate: %.4f\n', min(districtSatisfactionRates));
    % fprintf('----------------------------------------\n');
    % 
    % Return results as a structure
    % performanceResults = struct(...
    %     'overallSatisfactionRate', overallSatisfactionRate, ...
    %     'districtSatisfactionRates', districtSatisfactionRates, ...
    %     'numTestScenarios', numTestScenarios);
end