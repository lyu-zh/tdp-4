% demo_run2.m - Simple demo for testing the run2 function

% Clear workspace
clear all;

% Define instance file path
instanceFile = './Instances/2DU60-05-1.dat';

% Load the instance
instance = Instance(instanceFile);

% Create algorithm object
algo = Algo(instance);

% Run the algorithm using run2
fprintf('Running algorithm with run2 function...\n');
tic;
algo.run2(instanceFile);
elapsedTime = toc;

% Calculate objective value
objVal = 0;
centers = algo.centers_ids;
zones = algo.zones;
for j = 1:length(centers)
    for i = zones{j}
        objVal = objVal + instance.dist(i, centers(j));
    end
end

% Report results
fprintf('\nAlgorithm completed in %.2f seconds\n', elapsedTime);
fprintf('Final objective value: %.2f\n', objVal);