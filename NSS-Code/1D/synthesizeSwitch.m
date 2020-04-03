function out = synthesizeSwitch(Synthesize, A, B, C, u, f, a, h, N, Ns)

% This functions runs the finite difference scheme, in 3 cases, 'plot'
% where a plot of the moving string is shown. 'play' where the oscillation
% of a point on the string is measured to get a wave to play. 'plotsound'
% plots the graph of the sound produced by 'play'.


rp = 0.2;                % Position of reading
rp_int = 1+floor(N*rp);  % rounded grid index for readout 
rp_frac = 1+rp/h-rp_int; % fractional part of readout location


switch Synthesize
    case 'plot'
        % Loop to show the string over time
        for i=2:Ns-1
            % Calculates next step in time
            u(i+1,:) = A\(B*u(i,:)'+C*u(i-1,:)' + f(i)*a);
            % Plots string in real time
            plot(u(i,:), 'LineWidth',2);
            axis([0,N+2,-1,1]);
            drawnow
        end
        
    case 'play'
         
        out = zeros(1,Ns-1);  % To hear the sound
        % Loop to show the string over time
        for i=2:Ns-1
            % Calculates next step in time
            u(i+1,:) = A\(B*u(i,:)'+C*u(i-1,:)' + f(i)*a);
            % Vibration of one point on string
            %out(i) = u(i,round(N/7));
            out(i) = (1-rp_frac)*u(i,rp_int)+rp_frac*u(i,rp_int+1);
            
        end
        soundsc(out)    % To hear sound
        
       
    case 'plotsound'
        out = zeros(1,Ns-1);  % To hear the sound
        % Loop to show the string over time
        for i=2:Ns-1
            % Calculates next step in time
            u(i+1,:) = A\(B*u(i,:)'+C*u(i-1,:)' + f(i)*a);
            % Vibration of one point on string
            %out(i) = u(i,round(N/3));
            out(i) = (1-rp_frac)*u(i,rp_int)+rp_frac*u(i,rp_int+1);
        end
        plot(out)       % To plot the sound 
end

end

