
% //////////////////////////////////////////////////
% Finite difference method for nonlinear ODE
% //////////////////////////////////////////////////

hold off
clf
clear all

% Initial conditions.
n0    = 1;
n_end = 0;

r0    = 0;
r_end = 100;

N = 250;
h = (r_end-r0)/N;
r = r0:h:r_end;

% Initialize 0th order guess.
n = (n_end-n0)*(0:N)/N+n0;
[n(1), n(end)];

n_old  = n;
n_orig = n;

% Solution tolerance.
tol = 1e-5;

% Max # of iterations.
max_iter  = 10000;

alpha_vec = 0.0 : 0.01 : 0.8;
alpha_vec = [alpha_vec 0.81 : 0.005 : 0.97];
num_iter  = [];

for ia = 1:length(alpha_vec)

    % Re-initialize parameters.
    n         = n_orig;
    n_old     = n;

    iteration = 0;
    done      = 0;
    
    % Over-relaxation parameter.
    alpha = alpha_vec(ia);

    % Iterate to the desired tolerance.
    while done==false
        iteration=iteration+1;
        if iteration>max_iter
            display('max iterations exceeded')
            break
        end
        for i=2:N
            
            % Calculate updated estimate using Gauss-Seidel.
            n(i) = (sqrt(r(i))/h^2) * (n_old(i+1) + n(i-1)) / ...
                   (sqrt(n_old(i)) + (2/h^2)*sqrt(r(i)));
            
            % Over-relaxation.
            n(i) = n(i) + alpha * (n(i) - n_old(i));
            
        end
        
        % Tolerance check.
        rel_err=max(abs(n-n_old)./abs(n));
        rel_err(isnan(rel_err))=0; %set rel_err to zero in places where n=0
        if all(rel_err<tol)
            done=true;
        end    
        
        % Update the previous measurement.
        n_old = n;

    end
    

    if (mod(alpha,0.1)==0)

        clf
        hold off
        plot(r,zeros(1,length(r)),'black-');
        hold on
        plot(r,n_orig,'blue-','LineWidth',2.5);
        plot(r,n,'red-','LineWidth',2.5)
        axis([0,100,-1,1.2]);
        xlabel('r');
        ylabel('n(r)');
        aleg = legend('Zero', 'Initial guess', 'Converged solution');
        set(gca,'fontsize',20);
        set(aleg,'fontsize',15);
        hold off
        
        saveas(gca,[ 'plots/Q2/n_vs_r_alpha' num2str(alpha) '.eps' ...
                     ],'epsc');
        
        % Also zoom in on the interesting region.
        axis([0,20,0,1]);
        saveas(gca,[ 'plots/Q2/n_vs_r_alpha' num2str(alpha) ...
                     '_zoomed.eps' ],'epsc');
        
    end

    % Save the # of iterations taken.
    num_iter = [num_iter; iteration];

end

plot(alpha_vec, num_iter, 'black .', 'MarkerSize', 12);
axis([0,1,0,max(num_iter)]);
hold on;
drawnow;

xlabel('Relaxation parameter (\alpha)');
ylabel('# of iterations');
set(gca,'fontsize',20);
saveas(gca,'plots/Q2/niter_vs_alpha.eps','epsc');

hold off
