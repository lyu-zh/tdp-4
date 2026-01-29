classdef Algo < handle
    properties
        inst            % Instance object
        centers_ids     % Array of center IDs (instead of Area objects)
        zones           % Array of zone assignments
        r = 0.1         % Tolerance parameter
        timeLimit       % Time limit for solving
        demandUpperBound % Upper bound for demand
        globalConnectivityConstraints = [] % Global set of connectivity constraints
    end
    
    methods
        % Constructor
        function obj = Algo(instance)
            obj.inst = instance;
            obj.centers_ids = [];  % Initialize as empty array instead of empty list
            obj.zones = cell(1, instance.k);
            obj.r = 0.1;
            obj.timeLimit = inf; % Default: no time limit
            globalConnectivityConstraints = [];
            
            % Calculate demand upper bound
            obj.demandUpperBound = (1 + obj.r) * (sum(obj.inst.capacity) / instance.k);
        end
        
        % Set time limit
        function setTimeLimit(obj, seconds)
            obj.timeLimit = seconds;
        end
        
        % Run the algorithm
        function Best=run(obj, filename)
            startTime = tic;
            beta = 0.4;
            Best = inf;
            iter = 0;
            MaxIter = 10;
            alpha = 0;
            delta = 0.01;
            BestZones = cell(1, obj.inst.k);
            
            % Create random number generator
            rng(42); % For reproducibility
            areas = obj.inst.getAreas();
            while iter < MaxIter
                % Reset centers
                obj.centers_ids = [];
                
                % Step 1: Select initial centers using greedy random approach
                startId = randi(obj.inst.getN()); % Random start index (1-based)
                
                obj.centers_ids = [obj.centers_ids, startId]; % Store the ID instead of Area object
                areas{startId}.setCenter(true);
                
                % Continue selecting centers until we have k centers
                while length(obj.centers_ids) < obj.inst.k
                    % Calculate minimum distance from each area to closest center

                    minDistances = inf(1, obj.inst.getN());
                    % Mark which areas are centers in a single operation
                    centerMask = ismember(1:obj.inst.getN(), obj.centers_ids);
                    
                    % For each non-center area, find minimum distance to any center
                    nonCenterMask = ~centerMask;
                    nonCenterIndices = find(nonCenterMask);
                    
                    % Get the distance matrix subset for non-centers to centers
                    distSubMatrix = obj.inst.dist(nonCenterIndices, obj.centers_ids);
                    
                    % Find minimum distance for each non-center
                    [minDistances(nonCenterIndices), ~] = min(distSubMatrix, [], 2);
                    
                    % Find max and min of these minimum distances
                    maxMinDist = max(minDistances(nonCenterIndices));
                    minMinDist = min(minDistances(nonCenterIndices));
                    
                    % Calculate threshold value
                    Thre = maxMinDist - alpha * (maxMinDist - minMinDist);
                    
                    % Build candidate list
                    % Build candidate list using logical indexing
                    candidatesMask = ~centerMask & (minDistances >= Thre);
                    candidates = find(candidatesMask);
                    
                    % Select random candidate as next center
                    nextId = randi(length(candidates));
                    candidateId = candidates(nextId);
                    obj.centers_ids = [obj.centers_ids, candidateId];
                    areas{candidateId}.setCenter(true);
                end
                
                % Step 2: Solve assignment problem with the selected centers
                change = true;
                cur_value = 0.0;
                
                while change
                    change = false;
                    % disp(obj.centers_ids);
                    % Solve assignment problem using centers
                    [obj.zones, objVal] = obj.solveAssignmentProblem();
                    
                    cur_value = objVal;
                    % fprintf('---------------Current objective value: %.2f\n', cur_value);
                    % disp(obj.centers_ids);
                    % Check and adjust centers
                    for z = 1:length(obj.zones)
                        zone = obj.zones{z};
                        oldCenterId = obj.centers_ids(z);

                        % Find the new center that minimizes total distance
                        % Find the new center that minimizes total distance - vectorized
                        distMatrix = obj.inst.dist(zone, zone);  % 提取区域内节点间的距离矩阵
                        sumDists = sum(distMatrix, 2);  % 计算每个节点到其他所有节点的距离和
                        
                        % 找到距离和最小的节点作为新中心
                        [~, minIdx] = min(sumDists);
                        newCenterId = zone(minIdx);
                        
                        % Update center if a better one is found
                        if newCenterId ~= oldCenterId
                            change = true;
                            obj.centers_ids(z) = newCenterId;
                        end
                    end
                    % disp("change");
                    % disp(obj.centers_ids);
                    % disp("==========================================================================================================================");
                end

                
                
                if(~change)
                    iter = iter + 1;
                    if (alpha < beta)
                        alpha = alpha + delta;
                    else
                        alpha = 0;
                    end
                    if(Best>objVal)
                        BestZones = obj.zones;
                        Best = objVal;
                    end
                end
            end
            
            % Calculate total time
            endTime = toc(startTime);
            
            % Save results to file
            % obj.saveResults(filename, BestZones, Best, endTime);
        end
        % Save results to file
        function saveResults(obj, filename, bestZones, bestValue, runTime)
            % Create output file path
            outputFilePath = ['./output/', strrep(filename, '.dat', '.txt')];
            
            % Open file for writing
            fileID = fopen(outputFilePath, 'w');
            
            % Write center and zone information
            for i = 1:length(bestZones)
                if ~isempty(bestZones{i})
                    fprintf(fileID, 'center ID: %d\n', obj.centers_ids(i));
                    for j = 1:length(bestZones{i})
                        fprintf(fileID, '%d ', bestZones{i}(j));
                    end
                    fprintf(fileID, '\n');
                end
            end
            
            % Write objective value and runtime
            fprintf(fileID, 'best objective: %.2f\n', bestValue);
            fprintf(fileID, 'Runtime: %.2f s\n', runTime);
            
            % Close file
            fclose(fileID);
            
            % Also print to console
            fprintf('Results saved to: %s\n', outputFilePath);
            fprintf('Best objective value: %.2f\n', bestValue);
            fprintf('Runtime: %.2f s\n', runTime);
        end
        
        % % Solve the assignment problem with fixed centers
        function [zones, objVal] = solveAssignmentProblem(obj)
            % Initialize zones
            zones = cell(1, length(obj.centers_ids));
            
            % Get dimensions
            n = obj.inst.getN();
            p = length(obj.centers_ids);
            
            % Create YALMIP binary variables for assignment
            x = binvar(n, p, 'full');
            
            % Create constraints
            constraints = [];
            
            % Each area must be assigned to exactly one center
            constraints = [constraints, sum(x, 2) == 1];
            
            % Centers must be assigned to themselves
            for j = 1:p
                centerId = obj.centers_ids(j);
                constraints = [constraints, x(centerId, j) == 1];
            end
            
            % Capacity constraints for each center
            % Get demands vector from instance capacity
            demands = obj.inst.capacity;
            
            % Make sure demands is a column vector
            if size(demands, 1) == 1
                demands = demands'; % Convert row to column if needed
            end
            
            % Add all capacity constraints in one operation
            constraints = [constraints, demands' * x <= obj.demandUpperBound];
            
            % Objective: minimize total distance
            % Fully vectorized objective calculation
            objective = sum(sum(obj.inst.dist(:, obj.centers_ids) .* x));
            options = sdpsettings('solver', 'gurobi', 'verbose', 0);
            
            % combinedConstraints = [constraints, obj.globalConnectivityConstraints];
            combinedConstraints = constraints;
            
            % Solve the problem
            diagnostics = optimize(combinedConstraints, objective, options);
            
            % Initialize connectivity check variables
            allConnected = false;
            maxIterations = 100; % Limit iterations to prevent infinite loop
            iteration = 0;
            
            % Iteratively solve and add connectivity constraints
            while ~allConnected && iteration < maxIterations
                iteration = iteration + 1;
                
                % Solve the problem with current constraints
                diagnostics = optimize(constraints, objective, options);
                
                % Check if solution is feasible
                if diagnostics.problem ~= 0
                    warning('Optimization failed in iteration %d: %s', iteration, yalmiperror(diagnostics.problem));
                    break;
                end
                
                % Extract assignment values
                x_val = value(x);
                
                % Populate zones
                for j = 1:p
                    zones{j} = find(x_val(:, j) > 0.5)';
                end
                
                % Check connectivity of all zones
                disconnectedComponents = {};
                allConnected = true;
                
                for j = 1:p
                    zone = zones{j};
                    centerId = obj.centers_ids(j);
                    
                    % Find connected components in this zone
                    components = obj.findConnectedComponents(zone);
                    
                    if length(components) > 1
                        allConnected = false;
                        
                        % Find which component contains the center
                        centerComponentIdx = -1;
                        for c = 1:length(components)
                            if ismember(centerId, components{c})
                                centerComponentIdx = c;
                                break;
                            end
                        end
                        
                        % Store all disconnected components for this zone
                        if centerComponentIdx ~= -1
                            for c = 1:length(components)
                                if c ~= centerComponentIdx
                                    disconnectedComponents{end+1} = struct('component', components{c}, 'zoneIndex', j);
                                end
                            end
                        end
                    end
                end
                
                % If all zones are connected, exit the loop
                if allConnected
                    % fprintf('All zones are connected after %d iterations\n', iteration);
                    break;
                end
                
                % Add connectivity constraints for disconnected components
                constraintsAdded = 0;
                for i = 1:length(disconnectedComponents)
                    component = disconnectedComponents{i}.component;
                    zoneIdx = disconnectedComponents{i}.zoneIndex;
                    
                    % Find neighbors of this component (areas adjacent to component but not in it)
                    A = obj.inst.getAreas();
                    neighbors = [];
                    for nodeIdx = 1:length(component)
                        nodeId = component(nodeIdx);
                        nodeNeighbors = A{nodeId}.getNeighbors();
                        
                        for n = 1:length(nodeNeighbors)
                            neighborId = nodeNeighbors(n);
                            if ~ismember(neighborId, component)
                                neighbors = [neighbors, neighborId];
                            end
                        end
                    end
                    
                    % Remove duplicates
                    neighbors = unique(neighbors);
                    
                    % Create connectivity constraint: at least one neighbor must be in the same zone,
                    % or not all component nodes are in this zone
                    constr_expr = 0;
                    
                    % Sum assignments of neighbor nodes to this zone
                    for neighborId = neighbors
                        constr_expr = constr_expr + x(neighborId, zoneIdx);
                    end
                    
                    % Subtract assignments of component nodes to this zone
                    for nodeId = component
                        constr_expr = constr_expr - x(nodeId, zoneIdx);
                    end
                    newConstraint = constr_expr >= 1 - length(component);
                    obj.globalConnectivityConstraints = [obj.globalConnectivityConstraints, newConstraint];
                    % Add constraint: component nodes all assigned to zone j => at least one neighbor also assigned
                    constraints = [constraints, constr_expr >= 1 - length(component)];
                    constraintsAdded = constraintsAdded + 1;
                end
                
                % fprintf('Iteration %d: Added %d connectivity constraints\n', iteration, constraintsAdded);
            end
            
            % Final solution evaluation
            if allConnected
                % Calculate objective value
                objVal = value(objective);
            else
                % If we couldn't achieve connectivity, return the current solution
                warning('Could not achieve connectivity after %d iterations. Using best solution found.', maxIterations);
                objVal = value(objective);
            end
        end
      
        % Other methods will be similarly updated
        
        % Modified getCenters method to maintain compatibility with other code
        function areas = getCenters(obj)
            % Convert center IDs back to Area objects if needed by other code
            areas = [];
            A = obj.inst.getAreas();
            for i = 1:length(obj.centers_ids)
                centerId = obj.centers_ids(i); 
                areas = [areas, A{centerId}];
            end
        end
        % Find connected components in a zone
        % Find connected components in a zone
        function components = findConnectedComponents(obj, zone)
            components = {};
            
            if isempty(zone)
                return;
            end
            
            % Since IDs are 1-based, we can use the total number of areas
            visited = false(1, obj.inst.getN());
            
            for i = 1:length(zone)
                nodeId = zone(i);
                if ~visited(nodeId)
                    % Start a new component
                    component = [];
                    queue = nodeId;
                    visited(nodeId) = true;
                    
                    % BFS to find all connected nodes
                    while ~isempty(queue)
                        current = queue(1);
                        queue(1) = [];
                        component = [component, current];
                        
                        % Get neighbors of current node
                        areas = obj.inst.getAreas();
                        nodeNeighbors = areas{current}.getNeighbors();
                        
                        for n = 1:length(nodeNeighbors)
                            neighbor = nodeNeighbors(n);
                            if ismember(neighbor, zone) && ~visited(neighbor)
                                queue = [queue, neighbor];
                                visited(neighbor) = true;
                            end
                        end
                    end
                    
                    components{end+1} = component;
                end
            end
        end
        % Get centers IDs directly 
        function ids = getCenterIds(obj)
            ids = obj.centers_ids;
        end
        
        % Get zones
        function z = getZones(obj)
            z = obj.zones;
        end
        function run2(obj, filename)
            startTime = tic;
            
            % Step 1: Select initial centers using greedy random approach
            obj.selectInitialCenters();
            
            % Step 2: Iteratively solve assignment and update centers until stabilization
            change = true;
            iterationCount = 0;
            maxIterations = 100; % Prevent infinite loops
            
            while change && iterationCount < maxIterations
                iterationCount = iterationCount + 1;
                fprintf('Assignment iteration %d\n', iterationCount);
                
                % Solve assignment problem (without connectivity constraints)
                [zones, objVal] = obj.assignmentOnly();
                obj.zones = zones;
                
                % Check if assignment produced valid result
                if objVal == Inf
                    fprintf('Assignment problem has no feasible solution.\n');
                    break;
                end
                
                % Store old centers for comparison
                oldCenters = obj.centers_ids;
                
                % Update centers based on current zones
                obj.updateCenters();
                
                % Check if centers have changed
                change = ~isequal(oldCenters, obj.centers_ids);
                
                if change
                    fprintf('Centers updated, continuing iterations.\n');
                else
                    fprintf('Centers stabilized after %d iterations.\n', iterationCount);
                end
            end
            
            % Step 3: Ensure connectivity
            fprintf('Ensuring connectivity of all zones...\n');
            [connectedZones, connectedObjVal] = obj.ensureConnectivity();
            
            % Update with final solution
            obj.zones = connectedZones;
            
            % Calculate total time
            endTime = toc(startTime);
            
            % Log results without writing to file
            fprintf('Run2 complete - Runtime: %.2f seconds, Objective: %.2f\n', endTime, connectedObjVal);
        end
        
        function centers = run3(obj)
            startTime = tic;
            
            % Step 1: Select initial centers using greedy random approach
            obj.selectInitialCenters();
            
            % Step 2: Iteratively solve assignment and update centers until stabilization
            change = true;
            iterationCount = 0;
            maxIterations = 100; % Prevent infinite loops
            
            while change && iterationCount < maxIterations
                iterationCount = iterationCount + 1;
                % fprintf('Assignment iteration %d\n', iterationCount);
                
                % Solve assignment problem (without connectivity constraints)
                [zones, objVal] = obj.assignmentOnly();
                obj.zones = zones;
                
                % Check if assignment produced valid result
                if objVal == Inf
                    fprintf('Assignment problem has no feasible solution.\n');
                    break;
                end
                
                % Store old centers for comparison
                oldCenters = obj.centers_ids;
                
                % Update centers based on current zones
                obj.updateCenters();
                
                % Check if centers have changed
                change = ~isequal(oldCenters, obj.centers_ids);
                
                % if change
                %     fprintf('Centers updated, continuing iterations.\n');
                % else
                %     fprintf('Centers stabilized after %d iterations.\n', iterationCount);
                % end
            end
            
            % Step 3: Ensure connectivity
            % fprintf('Ensuring connectivity of all zones...\n');
            [connectedZones, connectedObjVal] = obj.ensureConnectivity();
            
            % Update with final solution
            obj.zones = connectedZones;
            centers = obj.centers_ids;
            % Calculate total time
            endTime = toc(startTime);
            
            % Log results without writing to file
            % fprintf('Run2 complete - Runtime: %.2f seconds, Objective: %.2f\n', endTime, connectedObjVal);
        end
        
        
        
        function [zones, objVal] = assignmentOnly(obj)
            % Initialize zones
            zones = cell(1, length(obj.centers_ids));
            
            % Get dimensions
            n = obj.inst.getN();
            p = length(obj.centers_ids);
            
            % Create YALMIP binary variables for assignment
            x = binvar(n, p, 'full');
            
            % Create constraints
            constraints = [];
            
            % Each area must be assigned to exactly one center
            constraints = [constraints, sum(x, 2) == 1];
            
            % Centers must be assigned to themselves
            for j = 1:p
                centerId = obj.centers_ids(j);
                constraints = [constraints, x(centerId, j) == 1];
            end
            
            % Capacity constraints for each center
            % Get demands vector from instance capacity
            demands = obj.inst.capacity;
            
            % Make sure demands is a column vector
            if size(demands, 1) == 1
                demands = demands'; % Convert row to column if needed
            end
            
            % Add capacity constraints
            constraints = [constraints, demands' * x <= obj.demandUpperBound];
            
            % Objective: minimize total distance
            objective = sum(sum(obj.inst.dist(:, obj.centers_ids) .* x));
            
            % Set options for solver
            options = sdpsettings('solver', 'gurobi', 'verbose', 0);
            
            % Solve the problem
            diagnostics = optimize(constraints, objective, options);
            
            % Check if solution is feasible
            if diagnostics.problem ~= 0
                warning('Assignment optimization failed: %s', yalmiperror(diagnostics.problem));
                objVal = Inf;
                return;
            end
            
            % Extract assignment values
            x_val = value(x);
            
            % Populate zones
            for j = 1:p
                zones{j} = find(x_val(:, j) > 0.5)';
            end
            
            % Calculate objective value
            objVal = value(objective);
        end
        
        function updateCenters(obj)
            % For each zone, find the node that minimizes total distance to all other nodes
            for z = 1:length(obj.zones)
                zone = obj.zones{z};
                
                if isempty(zone)
                    continue;
                end
                
                % Extract distance submatrix for the zone
                distMatrix = obj.inst.dist(zone, zone);
                
                % Calculate total distance from each node to all other nodes in zone
                sumDists = sum(distMatrix, 2);
                
                % Find node with minimum total distance
                [~, minIdx] = min(sumDists);
                newCenterId = zone(minIdx);
                
                % Update center
                obj.centers_ids(z) = newCenterId;
            end
        end
        
        function [zones, objVal] = ensureConnectivity(obj)
            % Start with current zones
            zones = obj.zones;
            
            % Get dimensions
            n = obj.inst.getN();
            p = length(obj.centers_ids);
            
            % Create YALMIP binary variables for assignment
            x = binvar(n, p, 'full');
            
            % Create basic constraints
            constraints = [];
            
            % Each area must be assigned to exactly one center
            constraints = [constraints, sum(x, 2) == 1];
            
            % Centers must be assigned to themselves
            for j = 1:p
                centerId = obj.centers_ids(j);
                constraints = [constraints, x(centerId, j) == 1];
            end
            
            % Capacity constraints for each center
            demands = obj.inst.capacity;
            
            % Make sure demands is a column vector
            if size(demands, 1) == 1
                demands = demands'; % Convert row to column if needed
            end
            
            % Add capacity constraints
            constraints = [constraints, demands' * x <= obj.demandUpperBound];
            
            % Objective: minimize total distance
            objective = sum(sum(obj.inst.dist(:, obj.centers_ids) .* x));
            
            % Set options for solver
            options = sdpsettings('solver', 'gurobi', 'verbose', 0);
            
            % Initialize connectivity check variables
            allConnected = false;
            maxIterations = 100; % Limit iterations to prevent infinite loop
            iteration = 0;
            
            % Iteratively solve and add connectivity constraints
            while ~allConnected && iteration < maxIterations
                iteration = iteration + 1;
                
                % Solve the problem with current constraints
                diagnostics = optimize(constraints, objective, options);
                
                % Check if solution is feasible
                if diagnostics.problem ~= 0
                    warning('Connectivity optimization failed in iteration %d: %s', iteration, yalmiperror(diagnostics.problem));
                    break;
                end
                
                % Extract assignment values
                x_val = value(x);
                
                % Populate zones
                for j = 1:p
                    zones{j} = find(x_val(:, j) > 0.5)';
                end
                
                % Check connectivity of all zones
                disconnectedComponents = {};
                allConnected = true;
                
                for j = 1:p
                    zone = zones{j};
                    centerId = obj.centers_ids(j);
                    
                    % Find connected components in this zone
                    components = obj.findConnectedComponents(zone);
                    
                    if length(components) > 1
                        allConnected = false;
                        
                        % Find which component contains the center
                        centerComponentIdx = -1;
                        for c = 1:length(components)
                            if ismember(centerId, components{c})
                                centerComponentIdx = c;
                                break;
                            end
                        end
                        
                        % Store all disconnected components for this zone
                        if centerComponentIdx ~= -1
                            for c = 1:length(components)
                                if c ~= centerComponentIdx
                                    disconnectedComponents{end+1} = struct('component', components{c}, 'zoneIndex', j);
                                end
                            end
                        end
                    end
                end
                
                % If all zones are connected, exit the loop
                if allConnected
                    % fprintf('All zones are connected after %d iterations\n', iteration);
                    break;
                end
                
                % Add connectivity constraints for disconnected components
                constraintsAdded = 0;
                for i = 1:length(disconnectedComponents)
                    component = disconnectedComponents{i}.component;
                    zoneIdx = disconnectedComponents{i}.zoneIndex;
                    
                    % Find neighbors of this component
                    A = obj.inst.getAreas();
                    neighbors = [];
                    for nodeIdx = 1:length(component)
                        nodeId = component(nodeIdx);
                        nodeNeighbors = A{nodeId}.getNeighbors();
                        
                        for n = 1:length(nodeNeighbors)
                            neighborId = nodeNeighbors(n);
                            if ~ismember(neighborId, component)
                                neighbors = [neighbors, neighborId];
                            end
                        end
                    end
                    
                    % Remove duplicates
                    neighbors = unique(neighbors);
                    
                    % Create connectivity constraint
                    constr_expr = 0;
                    
                    % Sum assignments of neighbor nodes to this zone
                    for neighborId = neighbors
                        constr_expr = constr_expr + x(neighborId, zoneIdx);
                    end
                    
                    % Subtract assignments of component nodes to this zone
                    for nodeId = component
                        constr_expr = constr_expr - x(nodeId, zoneIdx);
                    end
                    
                    % Add constraint: component nodes all assigned to zone j => at least one neighbor also assigned
                    constraints = [constraints, constr_expr >= 1 - length(component)];
                    constraintsAdded = constraintsAdded + 1;
                end
                
                % fprintf('Connectivity iteration %d: Added %d connectivity constraints\n', iteration, constraintsAdded);
            end
            
            % Final solution evaluation
            if allConnected
                objVal = value(objective);
            else
                warning('Could not achieve connectivity after %d iterations. Using best solution found.', maxIterations);
                objVal = value(objective);
            end
        end
        
        function selectInitialCenters(obj)
            rng(42);
            % Reset centers
            obj.centers_ids = [];
            
            % Select first center randomly
            startId = randi(obj.inst.getN());
            obj.centers_ids = [obj.centers_ids, startId];
            
            % Parameter for greedy randomized selection
            alpha = 0.3;
            
            % Continue selecting centers until we have k centers
            while length(obj.centers_ids) < obj.inst.k
                % Calculate minimum distance from each area to closest center
                minDistances = inf(1, obj.inst.getN());
                
                % Mark which areas are centers
                centerMask = ismember(1:obj.inst.getN(), obj.centers_ids);
                
                % For each non-center area, find minimum distance to any center
                nonCenterIndices = find(~centerMask);
                
                % Get the distance matrix subset for non-centers to centers
                distSubMatrix = obj.inst.dist(nonCenterIndices, obj.centers_ids);
                
                % Find minimum distance for each non-center
                [minDistances(nonCenterIndices), ~] = min(distSubMatrix, [], 2);
                
                % Find max and min of these minimum distances
                maxMinDist = max(minDistances(nonCenterIndices));
                minMinDist = min(minDistances(nonCenterIndices));
                
                % Calculate threshold value
                Thre = maxMinDist - alpha * (maxMinDist - minMinDist);
                
                % Build candidate list using logical indexing
                candidatesMask = ~centerMask & (minDistances >= Thre);
                candidates = find(candidatesMask);
                
                % Select random candidate as next center
                if ~isempty(candidates)
                    nextId = randi(length(candidates));
                    candidateId = candidates(nextId);
                    obj.centers_ids = [obj.centers_ids, candidateId];
                else
                    % If no candidates, choose the farthest area
                    [~, farthestId] = max(minDistances);
                    obj.centers_ids = [obj.centers_ids, farthestId];
                end
            end
        end
    end
end