classdef Area < handle
    properties
        id              % Area ID (1-based)
        x               % x-coordinate
        y               % y-coordinate
        activeness      % Activity metrics array
        isCenter = false % Flag to indicate if this is a center
        neighbors = []  % List of neighboring area IDs (using 1-based IDs)
    end
    
    methods
        % Constructor
        function obj = Area(id, x, y, activeness)
            if nargin > 0
                obj.id = id;  % Store the original ID (0-based)
                obj.x = x;
                obj.y = y;
                obj.activeness = activeness;
                obj.neighbors = [];
            end
        end
        
        % Add a neighbor to this area (using original 0-based ID)
        function addNeighbor(obj, neighborId)
            obj.neighbors = [obj.neighbors, neighborId];
        end
        
        % Getter for ID
        function value = getId(obj)
            value = obj.id;  % Return the original ID (0-based)
        end
        
        % Getter for x-coordinate
        function value = getX(obj)
            value = obj.x;
        end
        
        % Getter for y-coordinate
        function value = getY(obj)
            value = obj.y;
        end
        
        % Getter for activeness
        function value = getActiveness(obj)
            value = obj.activeness;
        end
        
        % Setter for ID
        function setId(obj, id)
            obj.id = id;
        end
        
        % Setter for x-coordinate
        function setX(obj, x)
            obj.x = x;
        end
        
        % Setter for y-coordinate
        function setY(obj, y)
            obj.y = y;
        end
        
        % Setter for activeness
        function setActiveness(obj, activeness)
            obj.activeness = activeness;
        end
        
        % Getter for isCenter
        function value = getIsCenter(obj)
            value = obj.isCenter;
        end
        
        % Setter for isCenter
        function setCenter(obj, center)
            obj.isCenter = center;
        end
        
        % Getter for neighbors
        function value = getNeighbors(obj)
            value = obj.neighbors;
        end
    end
end