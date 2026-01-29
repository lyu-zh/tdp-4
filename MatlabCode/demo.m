% demo.m - Simplified script to run Algo on all instance files

clear;
clc;

% Get all .dat files in the Instances directory
instanceFiles = dir('Instances/*.dat');
fprintf('Found %d instance files\n', length(instanceFiles));

% Create output directory if it doesn't exist
if ~exist('output', 'dir')
    mkdir('output');
end

% Print header
fprintf('%-30s %-20s %-15s\n', 'Instance', 'Objective Value', 'Runtime (s)');
fprintf('%-30s %-20s %-15s\n', '---------', '---------------', '-----------');

% Process each instance file
for i = 1:length(instanceFiles)
    instanceFilePath = fullfile(instanceFiles(i).folder, instanceFiles(i).name);
    
    try
        % Load instance and create algorithm object
        instance = Instance(instanceFilePath);
        algo = Algo(instance);
        
        % Set time limit (optional)
        algo.setTimeLimit(300); % 5 minutes
        % disp(instanceFiles(i).name)
        % Run algorithm and measure time
        tic;
        objValue= algo.run(instanceFiles(i).name);
        runtime = toc;
        fprintf('%-30s %-20.2f %-15.2f\n', instanceFiles(i).name, objValue, runtime);
        % Get objective value from saved file
        % outputFile = ['./output/', strrep(instanceFiles(i).name, '.dat', '.txt')];
        % if exist(outputFile, 'file')
        %     fileContent = fileread(outputFile);
        %     objMatch = regexp(fileContent, 'best objective: ([\d\.]+)', 'tokens');
        %     if ~isempty(objMatch)
        %         objValue = str2double(objMatch{1}{1});
        %         fprintf('%-30s %-20.2f %-15.2f\n', instanceFiles(i).name, objValue, runtime);
        %     else
        %         fprintf('%-30s %-20s %-15.2f\n', instanceFiles(i).name, 'N/A', runtime);
        %     end
        % else
        %     fprintf('%-30s %-20s %-15.2f\n', instanceFiles(i).name, 'No output file', runtime);
        % end
        
    catch e
        fprintf('%-30s %-20s %-15s\n', instanceFiles(i).name, 'ERROR', e.message);
    end
end