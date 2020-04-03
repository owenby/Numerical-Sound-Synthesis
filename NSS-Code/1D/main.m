clear sound
close all
%%% THINGS THAT NEED RESOLVING %%%
% I am not confident in the boundarySwitch, especially 'xylo'
% Can't get the ideal bar to be stable



% This script uses finite difference schemes to realistically simulate a 
% string, with given properties, with an end goal of synthesising sounds 
% from the string simulation

%%% Script Variables %%%

Boundary = 'fixed';    % 'free' both ends of the string free
                       % 'fixed' both ends fixed
                       % 'LfreeRfixed' left side free, right side fixed
                       % 'RfreeLfixed' right side free, left side fixed
                       
Starting = false;      % 'true' starts the string from a hann function shape
                       % 'false' starts the string from 0 amplitude everywhere 
                        
Force = 'pluck';  % 'off' for nothing
                       % 'pluck'
                       % 'hardpluck'
                       % 'strike'
                       % 'softstrike' for bar/stiff

Effect = 'realistic';      % 'none' for nothing
                       % 'loss' has just one type of loss
                       % 'realistic' has two
                       % 'bar' is a bar
                       % 'stiff' is a soft bar
                       % 'fdl' is frequency dependent loss
                       
Synthesize = 'play';   % Whether you see string ('plot'), hear string ('play'), or see the sound ('plotsound')


%%% Variable Setup %%%

% General variables
Fs = 44100;        % Sample rate
k = 1/Fs;          % Time step
T = 1;             % end time
L = 1;             % Length of string
c = 400;           % Wavespeed
h = c*k;           % Grid spacing

% Loss variables
b0 = 0.4;       % First loss term
b1 = 0.05;       % Second loss term

% Stiffness variables
kappa = 20;     % Stiffness term
theta = 1;      % Free parameter, must be > 0.5

% Stability Conditions
switch Effect
    case 'bar'
        mu = sqrt(2*theta - 1)/2;
        h = sqrt(k*kappa/mu);
    case 'stiff'
        mu = sqrt(2*theta - 1)/2;
        h = sqrt((c^2*k^2 + sqrt(c^4*k^4 + 16*kappa^2*k^2*(2*theta-1)))/(2*(2*theta-1)));
    case 'fdl'
        h = sqrt((c^2*k^2 + sqrt(c^4*k^4 + 16*kappa^2*k^2*(2*theta-1)))/(2*(2*theta-1)));
end

Ns = floor(T/k);   % Number of samples
N = floor(L/h);    % Number of string chunks
h = L/N;           % Redefine h so it matches with N


% Starting position of the string
u = zeros(Ns,N+1);

% Starting the string from a hann function shape
if Starting 
    u(1,(round(N/2-5):round(N/2+5))) = hann(11);
    u(2,(round(N/2-5):round(N/2+5))) = hann(11);
end


%%% Running Functions For Synthesis %%%

% Initialising the matrices for system
[A, B, C] = effectSwitch(Effect, b0, b1, c, h, k, N, kappa, theta);

% Force switch
[f, a] = forceSwitch(Force, h, N, Ns);

% Boundary condition switch
B = boundarySwitch(B, N, Boundary); 

% Switch to synthesize the sound
synthesizeSwitch(Synthesize, A, B, C, u, f, a, h, N, Ns);





 