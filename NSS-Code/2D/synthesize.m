function out = synthesize(Output, u1, u2, x, a, A, B, C, Ns, Nx, Ny, Fs, xFs, rp_index)

out = zeros(Ns,1);

    switch Output
        case 'plot'
            for i=3:Ns-1
    
                % Update scheme
                u = A\(B*u2 + C*u1) + a*x(i-2);
                u1 = u2;
                u2 = u;

                % Visuals
                mesh(reshape(u,[Nx+1,Ny+1]), 'FaceColor', 'interp');
                axis([-10 Nx+10 -10 Ny+10 -0.0001 0.0001])
                drawnow
                %axis equal 
                
            end 
            
        case 'play'
            for i=3:Ns-1

                % Update scheme
                u = A\(B*u2 + C*u1) + a*x(i-2);
                u1 = u2;
                u2 = u;

                % Sound
                out(i) = u(rp_index);   

            end
            soundsc(out, Fs);
            
            % Uncomment to create the sound as an audiofile:
            %out = out/max(abs(out));
            %audiowrite('2D.wav', out, xFs);
    end 
end

