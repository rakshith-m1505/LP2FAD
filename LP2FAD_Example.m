%% Initialisation
clc;
clear all;

%% Generate Target LPs and Define N (INPUT)

% Target In-Plane LPs [V1,V2,V3,V4]
LP_targets = [0.523205; 0; 0.150000; 0]; % Example case for Symmetric-Balanced FAD

% LP_targets= [0.4964; 0.2098; -0.2500; 0.3464]; % Example case for Symmetric-Unbalanced FAD

% Desired angle multiples in Layup <45 / 30 / 15>
angle_multiples = 15;

% Design Guidelines
symmetry = true;
balanced = true;

% Number of Layers to be designed
D_N = 10;

% Maximum no. of closely matching solutions to be saved
num_of_sols = 5;

% Print Results <true/false>
print_opt = true;

%% Generate Output

% Use LP2FAD function to generate outputs
[FADs,mismatch_errs] = LP2FAD(LP_targets,D_N,angle_multiples,symmetry,balanced,num_of_sols,print_opt);