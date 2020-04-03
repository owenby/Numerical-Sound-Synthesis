clear sound
close all

% This script uses finite difference schemes to realistically simulate a 
% string and hammer interaction

%%% Script Variables %%%

Boundary = 'fixed';    % 'free' both ends of the string free
                       % 'fixed' both ends fixed
                       % 'LfreeRfixed' left side free, right side fixed
                       % 'RfreeLfixed' right side free, left side fixed
                       
Starting = false;      % 'true' starts the string from a hann function shape
                       % 'false' starts the string from 0 amplitude everywhere 
                       
Effect = 'hamlss';      % Sets the effect on the string
                       % 'hamstr' for lossless hammer-string interaction
                       % 'hamlss' for single loss term ham str interaction
                       % 'hamlss2' for 2 loss term ham string interaction
                       % 'hamstf' for stiff ham str interaction
                       
Synthesize = 'plot';   % Whether you see string ('plot'), hear string ('play'), or see the sound ('plotsound')


%%% Variable Setup %%%

% General variables
Fs = 44100;        % Sample rate
k = 1/Fs;          % Time step
T = 1.5;           % end time
L = 1;             % Length of string
c = 2*L*261;       % Wavespeed
h = c*k;           % Grid spacing  
Ns = floor(T/k);   % Number of samples
N = floor(L/h);    % Number of string chunks
h = L/N;           % Redefine h so it matches with N

% Loss variables
b0 = 1;          % Loss term
b1 = 0.1;          % Second Loss term

% Stiffness variables
E = 2;                       % Youngs modulus of material
p = 2.31;                    % Density
A = 0.007303;                % Cross-sectional area
I = sqrt(A/pi)/2;            % Bar moment of inertia (I = R/2 for cylinder)
kappa = sqrt((E*I)/(p*A));   % Stiffness parameter

theta = 1;                   % 'free parameter'

%Hammer-String variables
hampnt = round((N+1)/2);     % Hammer position relative to string 
alpha = 2;                   % exp for phi   
uh = zeros(Ns,1);            % Hammer position vector
g = zeros(N+1,1);            % Vector to apply force at hampnt   
m = 0;                       % initialising value
v0 = 0.5;                     % velocity of hammer   
Tns=15;                     % Tension
ps= p*A;                     % Set linear mass density
c=Tns/ps;
Kham=10^5;
M= 1;                      % Mass of hammer                   
uh(1)= -0.05;
uh(2)= uh(1)+v0*k;           % init conds of hammer  

h=c*k;

switch Effect                % Coefficient of fn in rn equation
    case {'hamstr','hamstf'}                     
        m=(k^2)/M+(k^2)/(ps*h);          
    case 'hamlss'             
        m=(k^2)/(h*ps*(1+(b0*k)/(2*ps)))+k^2/M;     
end

%Hammer distribution
g(hampnt)=1/h;   

% Starting position of the string
u = zeros(Ns,N+1);

% Starting the string from a hann function shape
if Starting 
    u(1,(round(N/2-5):round(N/2+5))) = hann(11);
    u(2,(round(N/2-5):round(N/2+5))) = hann(11);
end


%%% Running Functions For Synthesis %%%

% Initialising the matrices for system
[A, B, C] = effectSwitch(Effect, b0, b1, c, h, k, N, kappa, ps, g);

% Force and where
f = zeros(Ns,1);
a=(k^2/ps)*g;

% Boundary condition switch
B = boundarySwitch(B, N, Boundary); 

% Switch to synthesize the sound
synthesizeSwitch(Synthesize, Effect, A, B, C, u, uh, g, f, a, c, N, Ns, hampnt, m, kappa, b0, b1, h, k, ps,Kham, M);

 