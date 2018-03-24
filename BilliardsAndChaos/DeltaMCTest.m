% ///////////////// EigenShoot.m ///////////////////
% Finding eigenfunctions and eigenvalues of an 
% aharmonic oscillator using shooting techniques.
% Uses RK45.m from previous project.
% //////////////////////////////////////////////////

hold off
clear all
clf

% Choose to plot or not.
doPlot = 1;

% Number of eigenvalues to break after.
nbreak = [6 6];

for iparam=1
    
    % Define function to be minimized.
    f = @(E)MinimizeMatchingCond(E,iparam);
    
    % Start by doing a coarse step size and looking at
    % d(MC). When we see d(MC) flip sign between
    % iterations, then we pick that narrow range
    % and do the root finding.
    if (iparam==1)
        E_vals = [0:0.5:30];
    else
        E_vals = [-30:1:20];
    end
    
    E_eigenvals     = [];
    x_fwd_eigenvals = [];
    u_fwd_eigenvals = [];
    x_bwd_eigenvals = [];
    u_bwd_eigenvals = [];
    
    dmc1_vals = [];
    dmc2_vals = [];

    for ie=1:(length(E_vals)-1)

        % Save ourselves some computational efficiency.
        dmatch1 = MinimizeMatchingCond(E_vals(ie),iparam);
        dmatch2 = MinimizeMatchingCond(E_vals(ie+1),iparam);
        
        dmc1_vals = [dmc1_vals; dmatch1];
        dmc2_vals = [dmc2_vals; dmatch2];
        
        fprintf('E1: %f, E2: %f, DM1: %f, DM2: %f \n',...
                E_vals(ie),E_vals(ie+1),dmatch1,dmatch2);
        
    end
    
    plot(E_vals(1:(end-1)),dmc1_vals,'red-');
    hold on
    plot(E_vals(1:(end-1)),dmc2_vals,'blue-');
    plot(E_vals(1:(end-1)),dmc1_vals.*dmc2_vals,'green-');
    
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
function [dmatch,x_a,u_a,x_b,u_b] = MinimizeMatchingCond(E,iparam)

    [x_a,u_a] = PropagatePsi(E,iparam,+1);
    [x_b,u_b] = PropagatePsi(E,iparam,-1);
    
    % Calculate the difference in matching condition.
    dmatch = (u_a(end,2)/u_a(end,1)) - (u_b(end,2)/u_b(end,1));
    
end
    
% Propagate wave function.
function [x,u] = PropagatePsi(E,iparam,dir)

    % Define potential well.
    if (iparam==1)
        parameters = [50, 2500, 0];
    else
        parameters = [50, 2500, 1500];
    end

    alpha = parameters(1);
    beta  = parameters(2);
    gamma = parameters(3);

    % Scale x0, x1 by some amount.
    x_sf = 1.0;
    
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
    %[x,u] = ode45(f_fwd,[x0 x_end],u0,odeset('RelTol',1e-10));
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
function [rn, niter, relerr] = bisect(f, r0, r1, tolerance)

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

        %err=abs(r1-r0)/abs(rmid)
        err = abs(r1-r0);
        
        fprintf('Guess: %15.10f, error: %15.10f\n',...
                rmid,err);
        
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