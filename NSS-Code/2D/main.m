%%%%% Plate and Membrane Modelling %%%%%
% To create a plate reverb effect, insert a sound.mp3 file of the audio you
% wish to add reverb to into this directory. Ensure the simulation is long
% enough for the audio input.

close all
clear sound

%%% Script Variables %%%

Simulation = 'complex';  % What kind of simulation takes place 
                         % 'membrane', losselss. Need lambda <= 1/root(2)
                         % 'plate', parameterised plate
                         % 'complex', 2 loss terms and 1 stiffness term

Starting = 'cosine';     % Starting position of the membrane
                         % 'none'     % No starting amplitude on any point
                         % 'cosine'   % Plate starts with cosine shape
                       
Output = 'play';         % How the output is presented
                         % 'plot'     % Animation of the plate moving
                         % 'play'     % Plays sound generated from the plate
                       
Force = 'none';          % Force applied to membrane
                         % 'none'     % No force acting on the plate
                         % 'reverb'   % Creates a reverb effect off an audio file

%%% Set Variables %%%

Fs = 40000;          % Sample rate
c = 200;             % Wave speed
T = 4;               % End time
Lx = 1;              % Length of x axis (less than 5 if you value your PC fans)
Ly = 1;              % Length of y axis
% Loss and Stiffness
b0 = 0.001;
b1 = 0.005;
kappa = 10;          % Greater than 2 for best results
% Excitation parameters
epoint = [0.3,0.5];  % Excitation point
size = 0.15;          % Excitation size
% Readout
rpoint = [0.2,0.5];  % Readout point
% Parameters for the Plate Model
varphi = 0; 
phi = 0;


%%% Determined Variables %%% 
 
k = 1/Fs;           % Time step
lambda = 1/sqrt(2); % Courant number
hx = c*k/lambda;    % x grid spacing
hy = c*k/lambda;    % y grid spacing
h = c*k/lambda;     % Same grid spacing

% Some stability conditions
switch Simulation
    case 'plate'
        mu = (varphi - sqrt(varphi^2 - (4*phi-1)))/(4*(4*phi-1));
        h = sqrt(kappa*k/mu);
        hx = h;
        hy = h;
    case 'complex' 
        varphi = 0;
        phi = 0;
        mu = (varphi - sqrt(varphi^2 - (4*phi-1)))/(4*(4*phi-1));
        h = sqrt(kappa*k/mu);
        hx = h;
        hy = h;
end


Nx = floor(Lx/hx);  % Number of x grid points        
Ny = floor(Ly/hy);  % Number of y grid points       
Ns = floor(T/k);    % Number of samples
N = (Nx+1)*(Ny+1);  % Size of big matrices
rp_index = (Ny-1)*floor(rpoint(1)*Nx)+floor(rpoint(2)*Ny); % Readout index

% Setting initial conditions
switch Starting
    case 'none'
        u1 = zeros(Nx+1, Ny+1);
        u2 = zeros(Nx+1, Ny+1);
        
    case 'cosine'
        [X, Y] = meshgrid([0:Nx]*hx, [0:Ny]*hy);
        dist = sqrt((X-epoint(1)).^2 +(Y-epoint(2)).^2);
        ind = sign(max(-dist+size/2,0)); 
        Raised_Cosine = 0.5*ind'.*(1+cos(2*pi*dist'/size));  
        
        u1 = Raised_Cosine;
        u2 = Raised_Cosine;
        
        u0=0;
        v0=1;
        u1 = u0*Raised_Cosine; 
        u2 = (u0+k*v0)*Raised_Cosine;
        
end

% Establishing force matrices for the plate reverb
switch Force 
    case 'none'
        a = zeros(Nx+1, Ny+1);
        a = a(:);
        [x, xFs] = audioread('sound.mp3');
    case 'reverb'
        [x, xFs] = audioread('sound.mp3');
        x = x(1:Ns-3);         % Sound input
        a = zeros(Nx+1, Ny+1); % Location of excitation
        a(round(length(a)*epoint(1)), round(length(a)*epoint(2))) = 1;
        a = a(:);     
end
        

% Getting the matrices to solve the scheme
D = laplacian(hx, hy, Nx, Ny, N);
[A, B, C] = effect(Simulation, c, k, b0, b1, kappa, varphi, phi, D, N);

% Unrolling initial vectors, sparse optimality
u1 = u1(:);   
u2 = u2(:);

% Running the simulation
synthesize(Output, u1, u2, x, a, A, B, C, Ns, Nx, Ny, Fs, xFs, rp_index);




 
