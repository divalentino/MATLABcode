% //////////// EigenShootLargerLimit.m /////////////
% Finding eigenfunctions and eigenvalues of an 
% aharmonic oscillator using shooting techniques.
% Uses RK45.m from previous project.
% //////////////////////////////////////////////////

hold off
clear all
clf

% Choose to plot or not.
doPlot = 0;

% Number of eigenvalues to break after.
nbreak = [6 6];

% Scale x0, x1 by some amount.
x_sf = 1.0;

for iparam=1
    
    % Define function to be minimized.
    f = @(E)MinimizeMatchingCond(E,x_sf,iparam);
    
    % Start by doing a coarse step size and looking at
    % d(MC). When we see d(MC) flip sign between
    % iterations, then we pick that narrow range
    % and do the root finding.
    if (iparam==1)
        E_vals = [1.8 2.2 6.7 7.2 12.7 13.2];
    else
        E_vals = [-19.2 -18.8 -9.2 -8.8 -1.2 -0.8];
    end
    
    E_eigenvals     = [];

    for ie=1:(length(E_vals)-1)

        % Check for change in sign of matching condition.
        [dmatch1] = MinimizeMatchingCond(E_vals(ie),x_sf,iparam);
        [dmatch2] = MinimizeMatchingCond(E_vals(ie+1),x_sf,iparam);

        fprintf('E1: %f, E2: %f, DM1: %f, DM2: %f \n',...
                E_vals(ie),E_vals(ie+1),dmatch1,dmatch2);
        
        % We cross zero in dMC, so find the minimum.
        if (dmatch1*dmatch2 < 0)

            fprintf('Trying root finding ...\n');
            
            [E_eigen,n_iter]=bisect(f,E_vals(ie),E_vals(ie+1),1e-5,1);

            fprintf('Got potential eigenvalue, E = %10.5f\n', ...
                    E_eigen);
            
            %Check that we don't have a duplicate.
            skipE = 0;
            for ivec=1:length(E_eigenvals)
                if (abs(E_eigen-E_eigenvals(ivec))<1e-5)
                    fprintf('Found a duplicate - skipping! \n')
                    skipE = 1;
                end
            end
            
            % //////////////////////////////////////////////////
            % Check for smoothness - how? Can we approximate
            % the derivative to check for kinks?
            % //////////////////////////////////////////////////
            
            if (skipE>0)
                continue;
            end
            
            %Save everything for later.
            E_eigenvals     = [E_eigenvals; E_eigen];

            % //////////////////////////////////////////////////
            % Break after finding a number of eigenvalues.
            if (length(E_eigenvals)>nbreak(iparam))
                break;
            end
            % //////////////////////////////////////////////////
            
        end
    end
    
    if (doPlot>0)
        if (iparam==1)
            overall_sf = 10;
        else
            overall_sf = 10;
        end
        PlotWaveFunction(E_eigenvals, iparam, overall_sf);
    end
    
    % Save the plot.
    set(gcf, 'Renderer', 'painters');
    saveas(gcf,[ 'plots/E_psi_vs_x_paramset_' num2str(iparam) ...
                 '_larger_limit.eps' ],'epsc');
        
end

% //////////////////////////////////////////////////
% Get the potential for a given set of points.
% //////////////////////////////////////////////////

function [Vx] = V(alpha,beta,gamma,x)
    Vx = alpha*x.^2 + gamma*x.^3 + beta*x.^4;
end

% //////////////////////////////////////////////////
% Shooting functions.
% //////////////////////////////////////////////////

% Propagate forward, backward solutions to midpoint.
function [dmatch,x_a,u_a,x_b,u_b] = MinimizeMatchingCond(E,x_sf,iparam)

    [x_a,u_a] = PropagatePsi(E,x_sf,iparam,+1);
    [x_b,u_b] = PropagatePsi(E,x_sf,iparam,-1);
    
    % Calculate the difference in matching condition.
    dmatch = (u_a(end,2)/u_a(end,1)) - (u_b(end,2)/u_b(end,1));
    
end
    
% Propagate wave function.
function [x,u] = PropagatePsi(E,x_sf,iparam,dir)

    % Define potential well.
    if (iparam==1)
        parameters = [50, 2500, 0];
    else
        parameters = [50, 2500, 1500];
    end

    alpha = parameters(1);
    beta  = parameters(2);
    gamma = parameters(3);
    
    % Constant factor of 2/h^2 
    % to scale potential, energy.
    thh_fac = 2 / 0.076199682;    
    
    % Initial parameters for root finding.
    u0      = [0; 1E-5];
    if (iparam==1)
        x0      = -0.6;
        x_end   = 0.1;
        if (dir<0)
            x0    = 0.7;
            x_end = 0.1;
        end
    else
        x0      = -1.5;
        x_end   = -0.38;
        if (dir<0)
            x0    =  1.5;
            x_end = -0.38;
        end        
    end
    
    % Redefine V, E as V*(2m/h^2), E*(2m/h^2).
    f_fwd = @(x,u)[u(2); (thh_fac*V(alpha,beta,gamma,x)-thh_fac*E)*u(1)];

    % Do a test to propagate the solution and get
    % an approximate wave function for some choice of 
    tol = 1e-10;
    [x,u] = RK45(f_fwd,[x_sf*x0 x_end],u0,tol,0);
    x = x';
    u = u';
    
end

% //////////////////////////////////////////////////
% Bisection algorithm.
% //////////////////////////////////////////////////

% Using the bisection method to compute the root
% based on some input guess.
function [rn, niter, relerr] = bisect(f, r0, r1, tolerance, useDebug)

    max_iter = 2000; %maximum number of iterations
    
    %Dummy return values.
    rn     = 0.0;
    n_iter = 0;
    rmid   = (r0+r1)/2.0;
    
    % *** Using code from C1_1.m (J. Krich)
    err    = 1; %current estimate of error
    n_iter = 0; %current iteration count
    
    while err>tolerance
        n_iter=n_iter+1;
        rmid=(r0+r1)/2.0;
        if f(rmid)*f(r0)<0
            r1=rmid;
        elseif f(rmid)*f(r1)<0
            r0=rmid;
        else
            error('no sign change')
        end

        %Want an absolute error to five decimal places now.
        err = abs(r1-r0);

        if (useDebug>0)
            fprintf('Guess: %15.10f, error: %15.10f\n',...
                    rmid,err);
        end
        
        if n_iter>max_iter
            warning('max iterations exceeded')
            break
        end
    end
    
    %Return values.
    rn     = rmid;
    niter  = n_iter;
    relerr = abs(r1-r0)/rmid;
    
end
