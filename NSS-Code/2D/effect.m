
function [A, B, C] = effect(Simulation, c, k, b0, b1, kappa, varphi, phi ,D, N)

    I = eye(N);

    switch Simulation
        case 'membrane'
            A = I;
            B = 2*I + c^2*k^2*D;
            C = -I;
        case 'plate' 
            A = I + varphi*k*kappa*D + phi*k^2*kappa^2*D^2;
            B = 2*A - k^2*kappa^2*D^2;
            C = -A;
        case 'complex'
            A = (1+b0*k/2)*I;
            B = 2*I - kappa*k^2*D^2 + (c^2*k^2 + b1*k)*D ;
            C = (b0*k/2-1)*I - b1*k*D;
    end

    A = sparse(A);
    B = sparse(B);
    C = sparse(C);
    

end

