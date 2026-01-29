% GabrielGraphGenerator.m
function GabrielGraphGenerator(numNodes, outputFilename, k, seed)
    % This function generates a Gabriel graph with the specified number of nodes
    % and saves it in the same format as the existing Instance files
    %
    % Inputs:
    %   numNodes: Number of nodes to generate
    %   outputFilename: Name of the output file (e.g., 'myInstance.dat')
    %   k: Number of districts to create
    %   seed: Random seed for reproducibility
    
    % Set random seed for reproducibility
    if nargin >= 4
        rng(seed);
    else
        rng('default');
    end
    
    % Generate random points in [0, 500] × [0, 500] square
    x = rand(numNodes, 1) * 500;
    y = rand(numNodes, 1) * 500;
    
    % Construct the Gabriel graph
    fprintf('Constructing Gabriel graph...\n');
    adjacencyMatrix = constructGabrielGraph(x, y);
    
    % Generate random activeness metrics
    baseDemand = 500; % Base demand value
    randomFactor = 0.2; % Random variation factor
    
    % Generate three activeness metrics for each node
    activeness1 = baseDemand + (rand(numNodes, 1) * 2 * randomFactor - randomFactor) * baseDemand;
    activeness2 = 1500 + rand(numNodes, 1) * 300; % Second metric (e.g., cost)
    activeness3 = 1500 + rand(numNodes, 1) * 300; % Third metric (e.g., time)
    
    % Count number of edges
    numEdges = sum(adjacencyMatrix(:)) / 2; % Divide by 2 because adjacencyMatrix is symmetric
    
    % Create edge list
    edges = zeros(numEdges, 2);
    edgeIdx = 1;
    for i = 1:numNodes
        for j = i+1:numNodes
            if adjacencyMatrix(i, j)
                edges(edgeIdx, :) = [i-1, j-1]; % Convert to 0-based for compatibility
                edgeIdx = edgeIdx + 1;
            end
        end
    end
    
    % Write to file in the required format
    writeInstanceFile(outputFilename, numNodes, x, y, activeness1, activeness2, activeness3, numEdges, edges, k);
    
    % Plot the graph for visualization
    fprintf('Plotting graph...\n');
    plotGraph(x, y, adjacencyMatrix, outputFilename);
    
    fprintf('Instance generation complete. File saved as: %s\n', outputFilename);
end

function adjacencyMatrix = constructGabrielGraph(x, y)
    % Constructs a Gabriel graph from the provided set of points
    % Two nodes are adjacent if and only if no other node lies inside 
    % the circle with the diameter being the line segment between the two nodes
    
    numNodes = length(x);
    adjacencyMatrix = zeros(numNodes, numNodes);
    
    for i = 1:numNodes
        for j = i+1:numNodes
            % Calculate center and radius of the circle with diameter (i,j)
            centerX = (x(i) + x(j)) / 2;
            centerY = (y(i) + y(j)) / 2;
            radius = sqrt((x(i) - x(j))^2 + (y(i) - y(j))^2) / 2;
            
            % Check if any other point is inside this circle
            isGabrielEdge = true;
            for k = 1:numNodes
                if k ~= i && k ~= j
                    % Calculate distance from point k to the center
                    distanceToCenter = sqrt((x(k) - centerX)^2 + (y(k) - centerY)^2);
                    
                    % If point k is inside the circle, (i,j) is not a Gabriel edge
                    if distanceToCenter < radius
                        isGabrielEdge = false;
                        break;
                    end
                end
            end
            
            % If no other point is inside the circle, add the edge
            if isGabrielEdge
                adjacencyMatrix(i, j) = 1;
                adjacencyMatrix(j, i) = 1;
            end
        end
    end
    
    return;
end

function writeInstanceFile(filename, numNodes, x, y, act1, act2, act3, numEdges, edges, k)
    % Writes the instance data to a file in the same format as the existing instances
    
    fileID = fopen(filename, 'w');
    
    % Write number of nodes
    fprintf(fileID, '%d\n', numNodes);
    
    % Write node data (id, x, y, activeness1, activeness2, activeness3)
    for i = 1:numNodes
        fprintf(fileID, '%d %f %f %f %f %f\n', i-1, x(i), y(i), act1(i), act2(i), act3(i));
    end
    
    % Write number of edges
    fprintf(fileID, '%d\n', numEdges);
    
    % Write edge data
    for i = 1:numEdges
        fprintf(fileID, '%d %d\n', edges(i, 1), edges(i, 2));
    end
    
    % Write k (number of districts)
    fprintf(fileID, '%d %f %f %f\n', k, 6, 7, 0.05);
    
    % Write seed
    fprintf(fileID, 'semilla: %d\n', randi(2^31 - 1));
    
    fclose(fileID);
end

function plotGraph(x, y, adjacencyMatrix, outputFilename)
    % Plot the graph for visualization
    
    figure('Position', [100, 100, 800, 800]);
    
    % Plot nodes
    scatter(x, y, 50, 'filled', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k');
    
    % Plot edges
    hold on;
    [i, j] = find(adjacencyMatrix);
    for idx = 1:length(i)
        if i(idx) < j(idx) % To avoid plotting each edge twice
            plot([x(i(idx)), x(j(idx))], [y(i(idx)), y(j(idx))], 'k');
        end
    end
    hold off;
    
    % Add title and labels
    title(['Gabriel Graph with ', num2str(length(x)), ' nodes']);
    xlabel('X Coordinate');
    ylabel('Y Coordinate');
    grid on;
    axis equal;
    
    % Save plot
    [path, name, ~] = fileparts(outputFilename);
    if isempty(path)
        saveas(gcf, [name, '_visualization.png']);
    else
        saveas(gcf, fullfile(path, [name, '_visualization.png']));
    end
    
    close(gcf);
end