classdef DistributionallyRobustAlgo < handle
    % This class implements a distributionally robust optimization approach
    % for solving the area division problem with uncertain demands
    properties
        inst                    % Instance object
        centers                 % List of center areas
        zones                   % Array of zone assignments
        r = 0.2                % Tolerance parameter
        gamma                   % Risk parameter
        scenarios               % Scenario data
        numScenarios            % Number of scenarios
        scenarioDemands         % Scenario demands
        meanVector              % Mean vector
        covarianceMatrix        % Covariance matrix
        delta1                  % D2 fuzzy set parameter
        delta2                  % D2 fuzzy set parameter
        useD1                   % Whether to use D1 fuzzy set
        useJointChance          % Whether to use joint chance constraint
        individualGammas        % Individual risk allocations for Bonferroni approx
        demandUpperBound        % Upper bound on demand
        timeLimit = 1200        % Time limit in seconds
        maxIterations = 1000     % Maximum iterations
    end
    
    methods
        % Constructor
        function obj = DistributionallyRobustAlgo(instance, scenarios, gamma, seed, useD1, delta1, delta2, useJointChance,r)
            obj.inst = instance;
            obj.scenarios = scenarios;
            obj.numScenarios = size(scenarios, 1);
            
            % Vectorized demand initialization
            obj.scenarioDemands = round(scenarios);
            
            obj.gamma = gamma;
            obj.r = r;
            rng(seed); % Set random seed
            obj.zones = cell(1, instance.k);
            obj.useD1 = useD1;
            obj.delta1 = delta1;
            obj.delta2 = delta2;
            obj.useJointChance = useJointChance;
            
            % Calculate moment information (mean and covariance)
            obj.calculateMomentInformation();
            
            % Calculate demand upper bound
            totalMeanDemand = sum(obj.meanVector);
            obj.demandUpperBound = (1 + obj.r) * (totalMeanDemand / instance.k);
            % Initialize individual gammas for Bonferroni approximation
            if useJointChance
                obj.individualGammas = repmat(gamma / instance.k, 1, instance.k);
            end
        end
        
        % Calculate moment information (mean and covariance) - Vectorized
        function calculateMomentInformation(obj)
            n = obj.inst.getN();
            
            % Calculate mean vector (vectorized)
            obj.meanVector = mean(obj.scenarios, 1);
            
            % Calculate covariance matrix (vectorized)
            centered = obj.scenarios - obj.meanVector;
            obj.covarianceMatrix = (centered' * centered) / obj.numScenarios;
            
            
            % Check if covariance matrix is positive semidefinite
            [~, p] = chol(obj.covarianceMatrix);
            isPSD = (p == 0);
            
            if ~isPSD
                fprintf('Warning: Covariance matrix is not positive semidefinite, attempting to fix...\n');
                obj.ensurePSDMatrix();
                
                % Check again
                [~, p] = chol(obj.covarianceMatrix);
                isPSD = (p == 0);
                fprintf('Fixed matrix is positive semidefinite: %d\n', isPSD);
            end
        end
        
        % Ensure the covariance matrix is positive semidefinite
        function ensurePSDMatrix(obj)
            % Add a small value to the diagonal
            epsilon = 1e-5;
            obj.covarianceMatrix = obj.covarianceMatrix + epsilon * eye(size(obj.covarianceMatrix));
        end
        
        % Run the algorithm
        function run(obj, filename)
            startTime = tic;
            Best = inf;
            BestZones = cell(1, obj.inst.k);
            
            % Step 1: Construct initial center set
            initialCenters = obj.selectInitialCenters();
            
            % Convert center IDs to Area objects
            obj.centers = [];
            areas = obj.inst.getAreas();
            for i = 1:length(initialCenters)
                centerId = initialCenters(i);
                area = areas{centerId};
                obj.centers = [obj.centers, area];
                area.setCenter(true);
            end
            % Step 2: Generate initial feasible solution
            feasible = obj.generateInitialSolution();
            if ~feasible
                fprintf('Cannot find feasible solution, check model parameters\n');
                return;
            end
            
            % Step 3: Improve initial solution
            change = true;
            cur_value = obj.evaluateObjective();
            iteration = 0;
            
            while change && iteration < obj.maxIterations
                iteration = iteration + 1;
                change = false;
                
                % Check each zone's true center
                newCenters = obj.findTrueCenters();
                
                % If centers changed, update and resolve
                if ~obj.compareCenters(obj.centers, newCenters)
                    obj.centers = newCenters;
                    change = true;
                    feasible = obj.generateInitialSolution();
                    if feasible
                        cur_value = obj.evaluateObjective();
                        % fprintf('Iteration %d: Objective value = %.2f\n', iteration, cur_value);
                    else
                        % fprintf('Iteration %d: No feasible solution found\n', iteration);
                        break;
                    end
                end
            end
            
            % Ensure connectivity
            if feasible
                obj.ensureConnectivity();
                cur_value = obj.evaluateObjective();
            end
            
            % Evaluate final result
            if feasible && cur_value < Best
                Best = cur_value;
                BestZones = obj.zones;
            end
            
            endTime = toc(startTime);
            timeSpentInSeconds = endTime;
            
            % Print results instead of writing to file
            fprintf('Best objective: %.2f\n', Best);
            fprintf('Runtime: %.2f s\n', timeSpentInSeconds);
            fprintf('Fuzzy set type: %s\n', conditional(obj.useD1, 'D_1', 'D_2'));
            
            if ~obj.useD1
                fprintf('delta1: %.2f, delta2: %.2f\n', obj.delta1, obj.delta2);
            end
            
            % fprintf('Constraint type: %s\n', conditional(obj.useJointChance, 'Joint constraint (DRJCC)', 'Individual constraint (DRICC)'));
            % fprintf('Risk parameter: %.2f\n', obj.gamma);
        end
        
        % Select initial centers - Optimized to avoid multiple scenario solving
        function centerIds = selectInitialCenters(obj)
            InitialNum = min(5, obj.numScenarios); % Adjustable parameter
            centerFrequency = zeros(1, obj.inst.getN());
            
            % Process a fixed number of scenarios without repeated sampling
            scenarioIndices = randperm(obj.numScenarios, InitialNum);
            
            for i = 1:InitialNum
                scenarioIndex = scenarioIndices(i);
                % fprintf('Processing scenario %d/%d, scenario index: %d\n', i, InitialNum, scenarioIndex);
                
                % Solve for the specific scenario
                scenarioCenters = obj.solveForScenario(scenarioIndex);
                
                % Update center frequency using vector operations
                validCenters = scenarioCenters(scenarioCenters > 0 & scenarioCenters <= obj.inst.getN());
                for center = validCenters
                    centerFrequency(center) = centerFrequency(center) + 1;
                end
            end
            
            % Sort centers by frequency
            [~, sortedIndices] = sort(centerFrequency, 'descend');
            
            % Select top k centers
            candidateCenters = sortedIndices(1:min(obj.inst.k, length(sortedIndices)));
            
            % Fill with random centers if needed
            while length(candidateCenters) < obj.inst.k
                remainingNodes = setdiff(1:obj.inst.getN(), candidateCenters);
                if isempty(remainingNodes)
                    break;
                end
                randomIdx = randi(length(remainingNodes));
                candidateCenters = [candidateCenters, remainingNodes(randomIdx)];
            end
            
            centerIds = candidateCenters;
        end
        
        % Solve deterministic model for a specific scenario - Improved error handling
        function scenarioCenters = solveForScenario(obj, scenarioIndex)
            localTimeLimit = 60; % seconds
            
            try
                % Create scenario instance
                scenarioInstance = obj.createScenarioInstance(scenarioIndex);
                
                % Create Algo object and set time limit
                algo = Algo(scenarioInstance);
                algo.setTimeLimit(localTimeLimit);
                
                % Get centers - this might fail if the scenario is infeasible
                scenarioCenters = algo.run3();
                
                % If not enough centers, supplement with random ones
                if length(scenarioCenters) < obj.inst.k
                    remainingNodes = setdiff(1:obj.inst.getN(), scenarioCenters);
                    numNeeded = obj.inst.k - length(scenarioCenters);
                    
                    if numNeeded <= length(remainingNodes)
                        randIndices = randperm(length(remainingNodes), numNeeded);
                        additionalCenters = remainingNodes(randIndices);
                        scenarioCenters = [scenarioCenters, additionalCenters];
                    end
                    
                    fprintf('Scenario %d: Not enough centers, randomly supplemented to %d centers\n', ...
                           scenarioIndex, length(scenarioCenters));
                end
                
            catch e
                fprintf('Error solving scenario %d: %s\n', scenarioIndex, e.message);
                
                % Use random selection as fallback
                scenarioCenters = randsample(obj.inst.getN(), obj.inst.k);
                fprintf('Scenario %d: Solution failed, randomly selected %d centers\n', ...
                       scenarioIndex, length(scenarioCenters));
            end
        end
        
        % Create scenario instance
        function scenarioInstance = createScenarioInstance(obj, scenarioIndex)
            scenarioInstance = obj.inst.createScenarioInstance(obj.scenarioDemands(scenarioIndex, :));
        end
        
        % Generate initial feasible solution using SOCP approach - Improved YALMIP usage
        function feasible = generateInitialSolution(obj)
            try
                % Create YALMIP optimization variables
                n = obj.inst.getN();
                p = length(obj.centers);
                
                % Decision variables - binary assignment variables
                x = binvar(n, p, 'full');
                
                % Constraints as arrays (more efficient)
                constraints = [];
                
                % Each area must be assigned to exactly one district (vector operation)
                constraints = [constraints, sum(x, 2) == ones(n, 1)];
                
                % Centers must be assigned to their own districts
                centerIds = zeros(1, p);
                for j = 1:p
                    centerIds(j) = obj.centers(j).getId();
                    constraints = [constraints, x(centerIds(j), j) == 1]; 
                end
                
                % Add distributionally robust chance constraints
                if obj.useD1
                    % D1 fuzzy set constraints - vectorized SOC constraint creation
                    for j = 1:p
                        xj = x(:, j);
                        meanTerm = obj.meanVector * xj;
                        
                        % Factor for SOC constraint
                        if obj.useJointChance
                            factor = sqrt((1 - obj.individualGammas(j)) / obj.individualGammas(j));
                        else
                            factor = sqrt((1 - obj.gamma) / obj.gamma);
                        end
                        
                        varTerm = sqrt(xj' * obj.covarianceMatrix * xj);
                        constraints = [constraints, meanTerm + factor * varTerm <= obj.demandUpperBound];
                    end
                else
                    % D2 fuzzy set constraints
                    for j = 1:p
                        xj = x(:, j);
                        meanTerm = obj.meanVector * xj;
                        
                        if obj.useJointChance
                            gamma_j = obj.individualGammas(j);
                            if obj.delta1 / obj.delta2 <= gamma_j
                                factor = sqrt(obj.delta1) + sqrt((1 - gamma_j) / gamma_j * (obj.delta2 - obj.delta1));
                            else
                                factor = sqrt(obj.delta2 / gamma_j);
                            end
                        else
                            if obj.delta1 / obj.delta2 <= obj.gamma
                                factor = sqrt(obj.delta1) + sqrt((1 - obj.gamma) / obj.gamma * (obj.delta2 - obj.delta1));
                            else
                                factor = sqrt(obj.delta2 / obj.gamma);
                            end
                        end
                        
                        varTerm = sqrt(xj' * obj.covarianceMatrix * xj);
                        constraints = [constraints, meanTerm + factor * varTerm <= obj.demandUpperBound];
                    end
                end
                
                % Objective function - vectorized distance calculation
                distMatrix = zeros(n, p);
                for j = 1:p
                    distMatrix(:, j) = obj.inst.dist(:, centerIds(j));
                end
                objective = sum(sum(distMatrix .* x));
                
                % Set solver options
                options = sdpsettings('solver', 'gurobi', 'verbose', 0, 'cachesolvers', 1, 'usex0', 0, 'gurobi.MIPGap', 0.01);
                
                
                % Solve the problem
                result = optimize(constraints, objective, options);
                if result.problem == 0
                    % Extract solution using logical indexing
                    x_val = value(x);
                    
                    % Populate zones using find function
                    for j = 1:p
                        obj.zones{j} = find(x_val(:, j) > 0.5)';
                    end
                    
                    feasible = true;
                else
                    fprintf('Optimization problem could not be solved: %s\n', yalmiperror(result.problem));
                    feasible = false;
                end
                
                return;
            catch e
                fprintf('Error in generateInitialSolution: %s\n', e.message);
                feasible = false;
                return;
            end
        end
        
        % Find true centers for each zone - Vectorized distance calculations
        function newCenters = findTrueCenters(obj)
            newCenters = [];
            areas = obj.inst.getAreas();
            
            for j = 1:length(obj.centers)
                if isempty(obj.zones{j})
                    newCenters = [newCenters, obj.centers(j)];
                    continue;
                end
                
                zone = obj.zones{j};
                
                % Extract distance submatrix for this zone
                distMatrix = obj.inst.dist(zone, zone);
                
                % Sum distances for each node (vectorized)
                totalDists = sum(distMatrix, 2);
                
                % Find node with minimum total distance
                [~, minIdx] = min(totalDists);
                bestCenter = zone(minIdx);
                
                newCenters = [newCenters, areas{bestCenter}];
            end
        end
        
        % Compare two sets of centers
        function isSame = compareCenters(obj, centers1, centers2)
            if length(centers1) ~= length(centers2)
                isSame = false;
                return;
            end
            
            % Extract IDs and compare as vectors
            ids1 = zeros(1, length(centers1));
            ids2 = zeros(1, length(centers2));
            
            for i = 1:length(centers1)
                ids1(i) = centers1(i).getId();
                ids2(i) = centers2(i).getId();
            end
            
            isSame = all(ids1 == ids2);
        end
        
        % Ensure connectivity for all zones - Improved constraint generation
        function ensureConnectivity(obj)
            allConnected = false;
            iteration = 0;
            obj.maxIterations = 1000;
            
            % Create YALMIP model with base constraints
            n = obj.inst.getN();
            p = length(obj.centers);
            
            % Decision variables
            x = binvar(n, p, 'full');
            
            % Basic constraints arrays
            constraints = [];
            
            % Each area must be assigned to exactly one district (vectorized)
            constraints = [constraints, sum(x, 2) == ones(n, 1)];
            
            % Centers must be assigned to their own districts
            centerIds = zeros(1, p);
            for j = 1:p
                centerIds(j) = obj.centers(j).getId();
                constraints = [constraints, x(centerIds(j), j) == 1];
            end
            
            % Add DRCC constraints
            if obj.useD1
                for j = 1:p
                    xj = x(:, j);
                    meanTerm = obj.meanVector * xj;
                    
                    if obj.useJointChance
                        factor = sqrt((1 - obj.individualGammas(j)) / obj.individualGammas(j));
                    else
                        factor = sqrt((1 - obj.gamma) / obj.gamma);
                    end
                    
                    varTerm = sqrt(xj' * obj.covarianceMatrix * xj);
                    constraints = [constraints, meanTerm + factor * varTerm <= obj.demandUpperBound];
                end
            else
                for j = 1:p
                    xj = x(:, j);
                    meanTerm = obj.meanVector * xj;
                    
                    if obj.useJointChance
                        gamma_j = obj.individualGammas(j);
                        if obj.delta1 / obj.delta2 <= gamma_j
                            factor = sqrt(obj.delta1) + sqrt((1 - gamma_j) / gamma_j * (obj.delta2 - obj.delta1));
                        else
                            factor = sqrt(obj.delta2 / gamma_j);
                        end
                    else
                        if obj.delta1 / obj.delta2 <= obj.gamma
                            factor = sqrt(obj.delta1) + sqrt((1 - obj.gamma) / obj.gamma * (obj.delta2 - obj.delta1));
                        else
                            factor = sqrt(obj.delta2 / obj.gamma);
                        end
                    end
                    
                    varTerm = sqrt(xj' * obj.covarianceMatrix * xj);
                    constraints = [constraints, meanTerm + factor * varTerm <= obj.demandUpperBound];
                end
            end
            
            % Prepare objective function - vectorized distance calculation
            distMatrix = zeros(n, p);
            for j = 1:p
                distMatrix(:, j) = obj.inst.dist(:, centerIds(j));
            end
            objective = sum(sum(distMatrix .* x));
            
            % Set solver options
            options = sdpsettings('solver', 'gurobi', 'verbose', 0, 'cachesolvers', 1, 'usex0', 0, 'gurobi.MIPGap', 0.01);
            
            
            % Initialize connectivity check
            areas = obj.inst.getAreas();
            
            while ~allConnected && iteration < obj.maxIterations
                iteration = iteration + 1;
                
                % Solve the model with current constraints
                result = optimize(constraints, objective, options);
                
                if result.problem ~= 0
                    fprintf('Connectivity iteration %d failed: %s\n', iteration, yalmiperror(result.problem));
                    break;
                end
                
                % Extract solution
                x_val = value(x);
                
                % Populate zones using find (faster than loops)
                for j = 1:p
                    obj.zones{j} = find(x_val(:, j) > 0.5)';
                end
                
                % Check connectivity of all zones
                disconnectedComponents = {};
                allConnected = true;
                
                % Check each zone for connectivity
                for j = 1:p
                    zone = obj.zones{j};
                    components = obj.findConnectedComponents(zone);
                    
                    if length(components) > 1
                        allConnected = false;
                        
                        % Find component containing center
                        centerId = centerIds(j);
                        centerComponentIdx = find(cellfun(@(comp) ismember(centerId, comp), components));
                        
                        if ~isempty(centerComponentIdx)
                            centerComponentIdx = centerComponentIdx(1);
                            
                            % Add disconnected components to list
                            otherComponentIndices = setdiff(1:length(components), centerComponentIdx);
                            for c = otherComponentIndices
                                disconnectedComponents{end+1} = struct(...
                                    'component', components{c}, ...
                                    'zoneIndex', j);
                            end
                        end
                    end
                end
                
                if allConnected
                    % fprintf('All zones are connected after %d iterations\n', iteration);
                    break;
                end
                
                % Add connectivity constraints in batches
                constraintCounter = 0;
                for i = 1:length(disconnectedComponents)
                    component = disconnectedComponents{i}.component;
                    zoneIdx = disconnectedComponents{i}.zoneIndex;
                    
                    % Find neighbors efficiently
                    neighbors = [];
                    for nodeIdx = 1:length(component)
                        nodeId = component(nodeIdx);
                        nodeNeighbors = areas{nodeId}.getNeighbors();
                        neighbors = [neighbors, nodeNeighbors];
                    end
                    
                    % Remove duplicates and nodes in the component
                    neighbors = setdiff(unique(neighbors), component);
                    
                    % Create constraint efficiently
                    if ~isempty(neighbors)
                        constr_expr = 0;
                        
                        % Sum for neighbors
                        for neighborId = neighbors
                            constr_expr = constr_expr + x(neighborId, zoneIdx);
                        end
                        
                        % Subtract for component nodes
                        for nodeId = component
                            constr_expr = constr_expr - x(nodeId, zoneIdx);
                        end
                        
                        constraints = [constraints, constr_expr >= 1 - length(component)];
                        constraintCounter = constraintCounter + 1;
                    end
                end
                
                % fprintf('Connectivity iteration %d: Added %d constraints\n', iteration, constraintCounter);
                
                if constraintCounter == 0
                    fprintf('Warning: No new constraints added, but connectivity issues persist\n');
                    break;
                end
            end
            
            if ~allConnected && iteration >= obj.maxIterations
                fprintf('Warning: Could not ensure connectivity within %d iterations\n', obj.maxIterations);
            end
        end
        
        % Find connected components in a zone - Optimized BFS algorithm
        function components = findConnectedComponents(obj, zone)
            components = {};
            
            if isempty(zone)
                return;
            end
            
            % Use a more efficient approach with logical indexing
            visited = false(1, obj.inst.getN());
            areas = obj.inst.getAreas();
            
            for i = 1:length(zone)
                startNode = zone(i);
                if ~visited(startNode)
                    % Start new component
                    component = [];
                    queue = startNode;
                    visited(startNode) = true;
                    
                    % BFS to find connected areas
                    while ~isempty(queue)
                        current = queue(1);
                        queue(1) = [];
                        component = [component, current];
                        
                        % Get neighbors
                        neighbors = areas{current}.getNeighbors();
                        
                        % Filter valid neighbors and mark as visited
                        validNeighbors = neighbors(ismember(neighbors, zone) & ~visited(neighbors));
                        visited(validNeighbors) = true;
                        queue = [queue, validNeighbors];
                    end
                    
                    components{end+1} = component;
                end
            end
        end
        
        % Calculate the current objective value - Vectorized calculation
        function objVal = evaluateObjective(obj)
            objVal = 0;
            
            for j = 1:length(obj.centers)
                centerId = obj.centers(j).getId();
                zoneNodes = obj.zones{j};
                
                if ~isempty(zoneNodes)
                    objVal = objVal + sum(obj.inst.dist(zoneNodes, centerId));
                end
            end
        end
        
        % Get zones
        function z = getZones(obj)
            z = obj.zones;
        end
        
        % Get centers
        function c = getCenters(obj)
            c = obj.centers;
        end
    end
end

% Helper function for conditional evaluation
function result = conditional(condition, trueValue, falseValue)
    if condition
        result = trueValue;
    else
        result = falseValue;
    end
end