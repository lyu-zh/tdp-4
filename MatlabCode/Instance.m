classdef Instance < handle
    properties
        n               % Number of areas
        areas           % Array of Area objects
        edges           % Adjacency matrix
        k               % Number of districts to create
        dist            % Distance matrix
        capacity        % Array of first activeness metric (demand) for all areas
    end
    
    methods
        % Constructor
        function obj = Instance(filepath)
            if nargin > 0
                obj.loadFromFile(filepath);
            end
        end
        
        % Load instance data from file
        function loadFromFile(obj, filepath)
            % Open the file
            fileID = fopen(filepath, 'r');
            if fileID == -1
                error('Cannot open file: %s', filepath);
            end
            
            % Read number of areas
            obj.n = fscanf(fileID, '%d', 1);
            
            % Initialize areas array
            obj.areas = cell(1, obj.n);
            
            % Initialize capacity array
            obj.capacity = zeros(1, obj.n);
            
            % Read area data
            for i = 1:obj.n
                id = fscanf(fileID, '%d', 1);
                x = fscanf(fileID, '%f', 1);
                y = fscanf(fileID, '%f', 1);
                
                activeness = zeros(1, 3);
                for j = 1:3
                    activeness(j) = fscanf(fileID, '%f', 1);
                end
                
                % Adjust ID to be 1-based right from the start
                id = id + 1;
                
                % Create Area object with the adjusted 1-based ID
                % and store it at the corresponding index
                obj.areas{id} = Area(id, x, y, activeness);
                
                % Store the first activeness metric (demand) in capacity array
                obj.capacity(id) = activeness(1);
            end
            
            % Initialize edges matrix
            obj.edges = zeros(obj.n, obj.n);
            
            % Read number of edges
            m = fscanf(fileID, '%d', 1);
            
            % Read edge data
            for i = 1:m
                a = fscanf(fileID, '%d', 1);
                b = fscanf(fileID, '%d', 1);
                
                % Adjust indices to be 1-based right from the start
                a = a + 1;
                b = b + 1;
                
                % Add neighbors using 1-based IDs
                obj.areas{a}.addNeighbor(b);
                obj.areas{b}.addNeighbor(a);
                
                % For the adjacency matrix, we use 1-based indices
                obj.edges(a, b) = 1;
                obj.edges(b, a) = 1;
            end
            
            % Read k (number of districts)
            obj.k = fscanf(fileID, '%d', 1);
            
            % Close the file
            fclose(fileID);
            
            % Generate distance matrix
            obj.generateDistanceMatrix();
        end
        
        % Create a clone of the current instance
        function cloned = clone(obj)
            cloned = Instance();
            cloned.n = obj.n;
            cloned.k = obj.k;
            
            % Copy areas
            cloned.areas = cell(1, obj.n);
            for i = 1:obj.n
                area = obj.areas{i};
                cloned.areas{i} = Area(area.getId(), area.getX(), area.getY(), area.getActiveness());
                cloned.areas{i}.setCenter(area.getIsCenter());
                
                % Copy neighbors
                for neighbor = area.getNeighbors()
                    cloned.areas{i}.addNeighbor(neighbor);
                end
            end
            
            % Copy edges
            cloned.edges = obj.edges;
            
            % Copy distance matrix
            cloned.dist = obj.dist;
        end
        
        % Create instance with specific scenario demands
        function inst = createScenarioInstance(obj, scenarioDemands)
            inst = obj.clone();
            inst.capacity = scenarioDemands;
        end
        
        % Generate distance matrix
        function generateDistanceMatrix(obj)
            obj.dist = zeros(obj.n, obj.n);
            for i = 1:obj.n
                for j = i+1:obj.n
                    area_i = obj.areas{i};
                    area_j = obj.areas{j};
                    distance = sqrt((area_i.getX() - area_j.getX())^2 + ...
                                   (area_i.getY() - area_j.getY())^2);
                    obj.dist(i, j) = distance;
                    obj.dist(j, i) = distance;
                end
            end
        end
        
        % Getters
        function value = getN(obj)
            value = obj.n;
        end
        
        function value = getAreas(obj)
            value = obj.areas;
        end
        
        function value = getEdges(obj)
            value = obj.edges;
        end
        
        function value = getCapacity(obj)
            value = obj.capacity;
        end
        
        % Output instance information
        function output(obj)
            fprintf('Number of areas: %d\n', obj.n);
            fprintf('Area information:\n');
            for i = 1:length(obj.areas)
                area = obj.areas{i};
                activeness = area.getActiveness();
                fprintf('Area ID: %d, x: %.2f, y: %.2f, Activeness: [%.2f, %.2f, %.2f]\n', ...
                    area.getId(), area.getX(), area.getY(), ...
                    activeness(1), activeness(2), activeness(3));
            end
            
            fprintf('Edge information:\n');
            for i = 1:obj.n
                for j = 1:obj.n
                    fprintf('%d ', obj.edges(i, j));
                end
                fprintf('\n');
            end
        end
    end
end