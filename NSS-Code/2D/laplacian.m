function D = laplacian(hx, hy, Nx, Ny, N)

    % returns 5+ point laplacian
   
    nOnesX = ones(Nx+1, 1) ;
    nOnesY = ones(N,1);
    
    Dx = -2*diag(nOnesX, 0) + diag(nOnesX(1:Nx), -1) + diag(nOnesX(1:Nx), 1);
    
    % Fixed Boundaries
    Dx(1,:) = 0;
    Dx(end,:) = 0;
    
    Dy = -2*eye(N) + diag(nOnesY(1:N-(Nx+1)), Nx+1) + diag(nOnesY(1:N-(Nx+1)), -(Nx+1));
    
    % Fixed Boundaries
    Dy(1:Nx+1,1:end) = 0;
    Dy(N-Nx:end,1:end) = 0;
    
    Dxx = (1/hx^2)*kron(eye(Ny+1), Dx) ;
    Dyy = (1/hy^2)*Dy;
    
    D = Dxx + Dyy;
    
    D = sparse(D);
    
end

