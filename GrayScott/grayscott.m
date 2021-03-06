% //////////////////////////////////////////////////
% THE GRAY-SCOTT EQUATIONS
% //////////////////////////////////////////////////

hold off
clear all
clf

% 0 for fast, 1 for accurate.
fastOrAccurate = 0;

% Animate the evolution.
doMovie   = 0;
drawPlots = 1;

% Do a scan over h and dt?
doStepScan = 0;

if (doStepScan>0)
    
    stepsize  = [0.01, 0.005, 0.002, 0.001];

    fprintf('\\begin{tabular}{c c c c} \\hline\\hline \n');
    for iparam = [1 2]
        % Re-initialize.
        u_p25_p25     = -999;
        u_p25_p25_old = -999;
        for tmax = [20]
            fprintf('\\multicolumn{4}{c}{Param set %i, tmax %i} \\\\ \\hline \n',...
                    iparam,tmax);
            for istep = [stepsize]
                u_p25_p25 = RunGrayScott(iparam,...
                                         0,...
                                         0,...
                                         tmax,...
                                         istep,...
                                         istep,...
                                         0);
                if (u_p25_p25_old ~= -999)
                    fprintf('%10.5f & %10.5f & %10.5f & %10.5f \\\\ \n',...
                            istep,istep,u_p25_p25,u_p25_p25_old- ...
                            u_p25_p25);
                else
                    fprintf('%10.5f & %10.5f & %10.5f & - \\\\ \n',...
                            istep,istep,u_p25_p25); 
                end
                
                % Update "old" measurement.
                u_p25_p25_old = u_p25_p25;
                
            end
        end
    end

    
    fprintf('\\hline\\hline \\end{tabular} \n');
        
    % That's all we really want,
    % when the step size scan has run,
    % the code should signal it's done.
    fprintf('Finished step size scan!\n');
    return;
    
end

if (fastOrAccurate>0)
    % Should give ~4 digits of accuracy.
    fprintf('Running the accurate simulation!\n');
    dt   = 0.001;
    h    = 0.001;
else
    fprintf('Running the fast simulation!\n');
    dt   = 0.01;
    h    = 0.01;
end

% Fixing h = dt to guarantee stability.
for ip = [1 2]
    for tmax = [500 3000]
            u_p25_p25 = RunGrayScott(ip,...
                                     doMovie,...
                                     drawPlots,...
                                     tmax,...
                                     dt,...
                                     h,...
                                     1);
            fprintf('dt: %10.5f, h: %10.5f, u(0.25,0.25): %7.5f \n',...
                    dt,h,u_p25_p25);
    end
end

function [u_p25_p25] = RunGrayScott(iparam,doMovie,drawPlots,tMax,dt,h,verbose)

    useMat    = 1;

    % Frame rate of movie.
    frameRate = 500;
    
    % Definte parameters of system.
    eu     = 5E-5;
    ev     = 2E-5;

    if (iparam==1)
        c = 0.065;
        F = 0.06;
    else
        c = 0.065;
        F = 0.03;
    end

    % Define x, y, t vectors.
    x = 0:h:0.5;
    y = 0:h:0.5;
    t = 0:dt:tMax;

    % Define # of entries.
    N = length(x);

    % //////////////////////////////////////////////////
    % Step 1: Define initial conditions.
    % //////////////////////////////////////////////////

    u = zeros(length(x),length(y));
    v = zeros(length(x),length(y));

    for ix=1:length(x)
        for iy=1:length(y)
            u(iy,ix) = min(1, 10*sqrt((x(ix)-0.2)^2 + (y(iy)-0.2)^2));
            v(iy,ix) = max(0, 1 - 10*sqrt((x(ix)-0.3)^2 + 2*(y(iy)-0.3)^2));
        end
    end

    % Define previous time-step matrices.
    u_old = u;
    v_old = v;

    % Define some helpful constants.
    lambda_u = (eu*dt)/h^2;
    lambda_v = (ev*dt)/h^2;

    % Prefix for saving plots.
    prefix = ['./paramset_' num2str(iparam) '_h_' num2str(h) ...
              '_dt_' num2str(dt) '_tmax_' num2str(tMax)];
    if (drawPlots>0)
        clf
        subplot(1,2,1)
        contourf(x,y,u);
        axis equal;
        xlabel('x');
        ylabel('y');
        drawnow
        subplot(1,2,2)
        contourf(x,y,v);
        axis equal;
        xlabel('x');
        ylabel('y');
        drawnow
        SavePlot([prefix '_init_cond_matrix']);
    end

    % //////////////////////////////////////////////////
    % Step 2: Propagate the system using Euler method.
    % //////////////////////////////////////////////////    
    
    if (useMat>0)
        
        % Convert to column vectors.
        u_np1 = u(:);
        v_np1 = v(:);
        
        u_n = u_old(:);
        v_n = v_old(:);
        
        % Make some vector quantities for
        % easier propagation.
        [l,w] = size(u_n);
        FF    = F*ones(l,w);
        cc    = c*ones(l,w);
        
        % Define toeplitz matrices.
        I   = speye(N);
        I   = sparse(I);
        L   = toeplitz([-2 1 zeros(1,N-2)]);
        L   = sparse(L);
        
        % NOTE: Need to add PBC terms.
        L(1,end) = 1;
        L(end,1) = 1;
        
        Lx = kron(L,I);
        Ly = kron(I,L);
        II = kron(I,I);
        
        A_u  = II+lambda_u*(Lx+Ly);
        A_v  = II+lambda_v*(Lx+Ly);

        A_u = sparse(A_u);
        A_v = sparse(A_v);
        
        % Taking a few calculations out of the loop
        % for speed purposes.
        Fdt = FF * dt;
        
        % For a number of time steps, update x and y.
        for it=1:length(t)

            % Start timing, for posterity.
            if (it==1)
                tic
            end
            
            if (verbose>0)
                if (mod(t(it),(tMax/10))==0) && (it>1)
                    fprintf('Time step: %f, elapsed (s): %f \n',...
                            t(it),toc);
                end
            end
            
            if (doMovie>0) && (mod(t(it),frameRate*dt)==0)
                
                % Fold u back into a 2D array for plotting.
                uu = reshape(u_n,N,N);
                vv = reshape(v_n,N,N);
                
                subplot(1,2,1)
                contourf(x,y,uu);
                axis equal;
                xlabel('x');
                ylabel('y');
                drawnow

                subplot(1,2,2)
                contourf(x,y,vv);
                axis equal;
                xlabel('x');
                ylabel('y');
                drawnow
                
            end
            
            % Break up the full operation into three
            % components:
            % 1. Lx, Ly matrices w/ periodic boundary
            % conditions
            % 2. c, F-related constant terms (add array)
            % 3. u_ij, v_ij related source terms (add array)

            % //////////////////////////////////////////////////
            % Propagate u.
            % //////////////////////////////////////////////////
            
            % Step 1: kron-based propagation matrices.
            u_np1 = A_u*u_n;
            
            % Step 2: Source terms.
            u_s = -u_n .* (dt*(v_n .* v_n + FF));
            
            % Sum them all together.
            u_np1 = u_np1 + u_s + Fdt;

            % //////////////////////////////////////////////////
            % Propagate v.
            % //////////////////////////////////////////////////

            % Step 1: kron-based propagation matrices.
            v_np1 = A_v*v_n;

            % Step 2: Source terms.
            v_s = v_n .* (dt*(u_n .* v_n - cc - FF));
            
            % Sum them all together.
            v_np1 = v_np1 + v_s;
            
            % //////////////////////////////////////////////////
            % Iterate u_n and v_n.
            % //////////////////////////////////////////////////
            
            u_n = u_np1;
            v_n = v_np1;

            % Calculate approximate time to complete.
            if (verbose>0)
                if (it==10)
                    tend = toc;
                    for tMaxTest = [10 500 3000]
                        tval = (tMaxTest / (10*dt))*(tend);
                        if (tval < 60)
                            fprintf('Time: %i, expected time to complete (s): %10.5f\n',...
                                    tMaxTest,tval);
                        else
                            fprintf('Time: %i, expected time to complete (m): %10.5f\n',...
                                    tMaxTest,tval/60);
                        end
                    end
                    pause(2);
                end
            end
            
        end
        
    end
    
    % Save final u, v arrays.
    u = u_n;
    v = v_n;

    % Fold u, v back into a 2D array for plotting.
    uu = reshape(u_n,N,N);
    vv = reshape(v_n,N,N);
    
    % Plot the final conditions.
    if (drawPlots>0)
        clf
        subplot(1,2,1)
        contourf(x,y,uu);
        axis equal;
        xlabel('x');
        ylabel('y');
        drawnow

        subplot(1,2,2)
        contourf(x,y,vv);
        axis equal;
        xlabel('x');
        ylabel('y');
        drawnow
        SavePlot([prefix '_final_cond_matrix']);
    end

    x(((length(uu)-1)/2)+1);
    y(((length(uu)-1)/2)+1);
    
    % Return the value of u(0.25,0.25).
    u_p25_p25 = uu(((length(uu)-1)/2)+1,((length(uu)-1)/2)+1);

end

% //////////////////////////////////////////////////
% Helper routines.
% //////////////////////////////////////////////////

function [] = SavePlot(filename)
    fprintf('Saving plot: %s.eps \n',filename)
    set(gcf, 'Renderer', 'painters');
    saveas(gca,[filename '.eps'],'epsc')
end
