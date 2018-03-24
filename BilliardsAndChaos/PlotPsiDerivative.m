function [] = PlotPsiDerivative(E_eigen, iparam, overall_sf)

    clf
    hold off
    
    if (iparam==1)
        parameters = [50, 2500, 0];
    else
        parameters = [50, 2500, 1500];
    end

    % Define potential well.
    alpha = parameters(1);
    beta  = parameters(2);
    gamma = parameters(3);

    if (iparam==1)
        xlow  = -0.6;
        xmid  =  -0.05;
        xhigh =  0.7;
    else
        xlow  = -2.1;
        xmid  = -0.2;
        xhigh =  1.5;
    end
    
    % Draw the potential well.            
    if (iparam==1)
        xrange = -0.5:0.01:0.5;
    else
        xrange = -0.85:0.01:0.5;
    end
    plot(xrange, V(alpha,beta,gamma,xrange), 'black');
    hold on;

    indices = 1:length(E_eigen);

    for ie = indices
        
        % Get the forward, backward wave functions.
        [x_fwd, u_fwd, x_bwd, u_bwd] = DoFinalPropagate(E_eigen(ie),xlow,xmid,xhigh,iparam,+1);
        
        % Normalize them individually.
        psi_fwd = u_fwd(:,1) / norm(u_fwd(:,1));
        psi_bwd = u_bwd(:,1) / norm(u_bwd(:,1));
        
        % Normalize psi'.
        psip_fwd = u_fwd(:,2) / norm(u_fwd(:,2));
        psip_bwd = u_bwd(:,2) / norm(u_bwd(:,2));
        
        % Stitch the wave forms together at the matching point.
        scale_sf_psip = psip_fwd(end) / psip_bwd(end);
        scale_sf      = psi_fwd(end) / psi_bwd(end)
        
        % ** Try to test for a kink somehow? i.e. Look for sudden
        % changes in psi(i+1) - psi(i) / delta_x. ***

        plot([x_fwd;x_bwd], E_eigen(ie)*ones(1,length([x_fwd;x_bwd])), ...
             'black-', 'LineWidth',2);
        hold on
        plot(x_fwd, E_eigen(ie) + overall_sf*psip_fwd,'green-', ...
             'LineWidth',2); 
        plot(x_bwd, E_eigen(ie) + overall_sf*scale_sf_psip*psip_bwd, ...
             'green-','LineWidth',2); 

        plot(x_fwd, E_eigen(ie) + overall_sf*psi_fwd,'red-', ...
             'LineWidth',2); 
        plot(x_bwd, E_eigen(ie) + overall_sf*scale_sf*psi_bwd,'blue-','LineWidth',2); 

        if (iparam==1)
            axis([-0.5,0.5,0,25]) 
        else
            axis([-1.2,0.5,-30,10]) 
        end
        
        xlabel('x (nm)');
        ylabel('E (eV)');
        drawnow
        
    end

end

function [x_fwd, u_fwd, x_bwd, u_bwd] = DoFinalPropagate(E,xlow,xmid,xhigh,iparam,dir)

    if (iparam==1)
        parameters = [50, 2500, 0];
    else
        parameters = [50, 2500, 1500];
    end
    
    % Define potential well.
    alpha = parameters(1);
    beta  = parameters(2);
    gamma = parameters(3);
    
    u0    = [0; 1E-5];
    
    % Define functions for root finding.
    thh_fac = 2 / 0.076199682;
    f_fwd = @(x,u)[u(2); (thh_fac*V(alpha,beta,gamma,x)-thh_fac*E)*u(1)];
    
    [x_fwd,u_fwd] = RK45(f_fwd,[xlow xmid],u0,1e-10,0);
    x_fwd = x_fwd';
    u_fwd = u_fwd';
    
    [x_bwd,u_bwd] = RK45(f_fwd,[xhigh xmid],u0,1e-10,0);
    x_bwd = x_bwd';
    u_bwd = u_bwd';
    
end

% Get the potential for a given set of points.
function [Vx] = V(alpha,beta,gamma,x)
    Vx = alpha*x.^2 + gamma*x.^3 + beta*x.^4;
end