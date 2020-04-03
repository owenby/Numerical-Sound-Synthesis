%%%%% Repeated Hammer Strikes of a Stiff, Lossy Membrane Model %%%%%

close all
clear all

%%% Script Variables %%%
                       
Output = 'play';         % How the output is presented
                         % 'plot'     % Animation of the plate moving
                         % 'play'     % Plays sound generated from the plate

%%% Set Variables %%%

Fs = 44100;          % Sample rate
c = 400;             % Wave speed
T = 4;               % End time
Lx = 1;              % Length of x axis (less than 5, depending on your PC)
Ly = 1;              % Length of y axis
% Loss and Stiffness
b0 = 0.05;
b1 = 0.05;
kappa = 10;          % Greater than 5
% Excitation
epoint = [0.5,0.5];  % Excitation point
size = 0.5;          % Excitation size
% Readout
rpoint = [0.4,0.1];  % Readout point

%%% Determined Variables %%% 
 
k = 1/Fs;           % Time step
lambda = 1/sqrt(2); % Courant number
hx = c*k/lambda;    % x grid spacing
hy = c*k/lambda;    % y grid spacing
h = c*k/lambda;     % Same grid spacing
varphi = 0;
phi = 0;

% Some stability conditions
mu = (varphi - sqrt(varphi^2 - (4*phi-1)))/(4*(4*phi-1));
h = sqrt(kappa*k/mu);
hx = h;
hy = h;

Nx = floor(Lx/hx);  % Number of x grid points        
Ny = floor(Ly/hy);  % Number of y grid points       
Ns = floor(T/k);    % Number of samples
N = (Nx+1)*(Ny+1);  % Size of big matrices
rp_index = (Ny-1)*floor(rpoint(1)*Nx)+floor(rpoint(2)*Ny); % Readout index


u1 = zeros(Nx+1, Ny+1);
u2 = zeros(Nx+1, Ny+1);   


% Getting the matrices to solve the scheme
D = laplacian(hx, hy, Nx, Ny, N);
I = eye(N);
A = (1+b0*k/2)*I;
B = 2*I - kappa*k^2*D^2 + (c^2*k^2 + b1*k)*D ;
C = ((b0/2)*k-1)*I - b1*k*D;
A = sparse(A); 
B = sparse(B);
C = sparse(C);

% Unrolling initial vectors, sparse optimality
u1 = u1(:);   
u2 = u2(:);

% Hammer Variables
M = 10^(9);           % Mass of hammer
K = 10^(2);        % Stiffness of hammer
alpha = 2;       % Power variable
v = 10^(3);  % Velocity of hammer


% Hammer vectors
g = zeros(Nx+1, Ny+1);  % Distribution of hammer
index1 = round(epoint(1)*(Nx+1));
index2 = round(epoint(2)*(Ny+1));
g(index1, index2) = 1/h;        % Excitation point
g=g(:);                         % Excitation in vector form

uh = zeros(Ns,1);
uh(1) = -1;             % Initial position of hammer
uh(2) = uh(1) + v*k;    % Initial conditions

% To update the scheme (hammer)
m = -(k^2)/M - kappa^2/((1+b0*k/2));
eta = zeros(Ns,1);
eta(1) = uh(1) - u1((index2-1)*(Nx+1) + index1);
f = zeros(Ns,1);      % Hold the force values
G=1;
eps = 10^(-6);
b = 2*uh(2) - 2*uh(1) - ((2+c^2*k^2*(-4)+b1*k*(-4)-kappa*k^2*(16))/(1+b0*k/2))*u2((index2-1)*(Nx+1) + index1) - ((b0*k/2 - 1 - b1*k*(-4))/(1+b0*k/2) - 1)*u1((index2-1)*(Nx+1) + index1);
r = b;

% Running the simulation
out = zeros(Ns,1);

    switch Output
        case 'plot'
            for i=2:Ns-1
                eta(i) = uh(i) - u2((index2-1)*(Nx+1) + index1);
                a = eta(i-1);
                phi_a = K/(alpha+1)*(max(0,a)^(alpha+1));             

                b = 2*uh(i) - 2*uh(i-1) - ((2+c^2*k^2*(-4)+b1*k*(-4)-kappa*k^2*(16))/(1+b0*k/2))*u2((index2-1)*(Nx+1) + index1) - ((b0*k/2 - 1 - b1*k*(-4))/(1+b0*k/2) - 1)*u1((index2-1)*(Nx+1) + index1);
                
                while abs(G)>eps
                    phi_ra = K/(alpha+1)*((max(0,r+a))^(alpha+1));
                    G = r + m*(phi_ra - phi_a)/(r) - b;
                    Gp = 1 + (m*(K*max(0,r+a)^alpha-phi_ra+phi_a))/(r^2);
                    r = r - G/Gp;
                end 
                f(i) = (b-r)/m;
                r=b+f(i)*m;
                
                % Update scheme vectors
                u = A\(B*u2 + C*u1 + k^2*f(i)*g);
                u1 = u2;
                u2 = u;
                uh(i+1)=2*uh(i)-uh(i-1)-(k^2/M)*f(i);
      
                % Visuals
                mesh(reshape(u,[Nx+1,Ny+1]), 'FaceColor', 'interp');
                hold on
                plot3(index1, index2,  -uh(i), 'o')
                hold off
                axis([-10 Nx+10 -10 Ny+10 -1 1])
                drawnow
                
            end
            
            
        case 'play'
            for i=2:Ns-1
                
                eta(i) = uh(i) - u2((index2-1)*(Nx+1) + index1);
                a = eta(i-1);
                phi_a = K/(alpha+1)*(max(0,a)^(alpha+1));             

                b = 2*uh(i) - 2*uh(i-1) - ((2+c^2*k^2*(-4)+b1*k*(-4)-kappa*k^2*(16))/(1+b0*k/2))*u2((index2-1)*(Nx+1) + index1) - ((b0*k/2 - 1 - b1*k*(-4))/(1+b0*k/2) - 1)*u1((index2-1)*(Nx+1) + index1);
                
                while abs(G)>eps
                    phi_ra = K/(alpha+1)*((max(0,r+a))^(alpha+1));
                    G = r + m*(phi_ra - phi_a)/(r) - b;
                    Gp = 1 + (m*(K*max(0,r+a)^alpha-phi_ra+phi_a))/(r^2);
                    r = r - G/Gp;
                end 
                f(i) = (b-r)/m;
                r=b+f(i)*m;
                
                % Update scheme vectors
                u = A\(B*u2 + C*u1 + k^2*f(i)*g);
                u1 = u2;
                u2 = u;
                uh(i+1)=2*uh(i)-uh(i-1)-(k^2/M)*f(i);

                % Sound
                out(i) = u(rp_index); 

            end
            soundsc(out, Fs);
            
            % Uncomment to create the sound as an audiofile:
            %out = out/max(abs(out));
            %audiowrite('Hammer_Strikes.wav', out, Fs);
    end 




 
