% test_gabriel_graph.m
% Script to test the Gabriel graph instance generator

% Clear workspace
clear all;
close all;
clc;

% Parameters
numNodes = 30;        % Number of nodes to generate
outputFilename = './Instances/GabrielGraph100.dat';  % Output file path
k = 4;                 % Number of districts to create
seed = 12345;          % Random seed for reproducibility

% Generate Gabriel graph instance
fprintf('Generating Gabriel graph instance with %d nodes...\n', numNodes);
GabrielGraphGenerator(numNodes, outputFilename, k, seed);

% Test if we can load the generated instance
fprintf('Testing if the generated instance can be loaded...\n');
try
    instance = Instance(outputFilename);
    fprintf('Instance loaded successfully!\n');
    fprintf('Number of areas: %d\n', instance.getN());
    fprintf('Number of edges: %d\n', sum(sum(instance.getEdges()))/2);
    fprintf('Number of districts to create: %d\n', instance.k);

    
catch e
    fprintf('Error: %s\n', e.message);
end

fprintf('\nProcess completed.\n');