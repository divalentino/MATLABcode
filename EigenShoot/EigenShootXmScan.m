% ///////////////// EigenShoot.m ///////////////////
% Finding eigenfunctions and eigenvalues of an 
% aharmonic oscillator using shooting techniques.
% Uses RK45.m from previous project.
% //////////////////////////////////////////////////

hold off
clear all
clf

% Choose to plot or not.
doPlot = 0;

% Plot the d_MC function vs. E.
doPlotdMCvsE = 1;

% Number of eigenvalues to break after.
nbreak = [2 2];

% Scale x0, x1 by some amount.
x_sf = 1.75;

for iparam=2
    hold off
    clf
            
    E_eigenvals     = [];
    
    %for idx=-0.05:0.05:0.05
    for idx=0.0
    
        if (iparam==1)
            x_a = -0.8;
            x_m =  0.1 + idx;
            x_b =  0.8;
        else
            x_a = -0.8;
            x_m = -0.38 + idx;
            x_b =  0.6;
        end
        
        % Define function to be minimized.
        f = @(E)MinimizeMatchingCond(E,x_sf,x_a,x_m,x_b,iparam);
        
        % Start by doing a coarse step size and looking at
        % d(MC). When we see d(MC) flip sign between
        % iterations, then we pick that narrow range
        % and do the root finding.
        if (iparam==1)
            E_vals = [0:0.1:15];
        else
            E_vals = [-20:0.5:10];
        end

        for ie=1:(length(E_vals)-1)

            % Check for change in sign of matching condition.
            [dmatch1] = MinimizeMatchingCond(E_vals(ie),x_sf,x_a,x_m,x_b,iparam);
            [dmatch2] = MinimizeMatchingCond(E_vals(ie+1),x_sf,x_a,x_m,x_b,iparam);
            
            plot(E_vals(ie), dmatch1, 'black .','MarkerSize',20);
            hold on
            drawnow
            
            % We cross zero in dMC, so find the minimum.
            if (dmatch1 > 0) && (dmatch2 < 0)

                plot(E_vals(ie), dmatch1, 'red x', 'MarkerSize', 30);
                drawnow
                
                %[E_eigen,n_iter]=bisect(f,E_vals(ie),E_vals(ie+1), ...
                %                        1e-5,0);
                [E_eigen,n_iter] = fzero(f,E_vals(ie));
                                
                fprintf('Got potential eigenvalue, E = %10.6f\n', ...
                        E_eigen);

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
    end

    set(gcf, 'Renderer', 'painters');
    if (doPlotdMCvsE>0)
        xlabel('E (eV)')
        ylabel('\delta_{mc}')
        set(gca,'fontsize',20);
        axis([E_vals(1),E_vals(end),-300,300]);
        saveas(gcf,[ 'plots/dMC_vs_E_paramset_' num2str(iparam) ...
                     '.eps' ],'epsc');
    end
    
    % Plot the wave functions.
    if (doPlot>0)
        if (iparam==1)
            overall_sf = 10;
        else
            overall_sf = 10;
        end
        PlotWaveFunction(E_eigenvals, iparam, overall_sf);
    end
    
    % Save the plot.
    saveas(gcf,[ 'plots/E_psi_vs_x_paramset_' num2str(iparam) ...
                 '.eps' ],'epsc');
    
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
function [dmatch,x_a,u_a,x_b,u_b] = MinimizeMatchingCond(E,x_sf,x_a,x_m,x_b,iparam)

    [x_a,u_a] = PropagatePsi(E,x_sf,x_a,x_m,x_b,iparam,+1);
    [x_b,u_b] = PropagatePsi(E,x_sf,x_a,x_m,x_b,iparam,-1);
    
    % Calculate the difference in matching condition.
    dmatch = (u_a(end,2)/u_a(end,1)) - (u_b(end,2)/u_b(end,1));
    
end
    
% Propagate wave function.
function [x,u] = PropagatePsi(E,x_sf,x_a,x_m,x_b,iparam,dir)

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
    x0    = x_a;
    x_end = x_m;
    if (dir<0)
        x0 = x_b;
    end
    
    % Redefine V, E as V*(2m/h^2), E*(2m/h^2).
    f_fwd = @(x,u)[u(2); (thh_fac*V(alpha,beta,gamma,x)-thh_fac*E)*u(1)];

    % Do a test to propagate the solution and get
    % an approximate wave function for some choice of 
    tol = 1e-10;
    %[x,u] = RK45(f_fwd,[x_sf*x0 x_end],u0,tol,0);
    %x = x';
    %u = u';
    
    [x,u] = ode45(f_fwd,[x0 x_end],u0,odeset('RelTol',1e-10));
    
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
