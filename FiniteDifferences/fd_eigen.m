
% //////////////////////////////////////////////////
% Finite differences for eigenvalues and 
% eigenfunctions
% //////////////////////////////////////////////////

clear all
hold off
clf

% Initial configuration.
drawPlot = 0;

for iparam=1:2
    hold off
    clf

    % Define potential well.
    if (iparam==1)
        parameters = [50, 2500, 0];
        overall_sf = 15;
    else
        parameters = [50, 2500, 1500];
        overall_sf = 25;
    end
    alpha = parameters(1);
    beta  = parameters(2);
    gamma = parameters(3);

    % Constant to make eq'n nice.
    thh_fac = 2 / 0.076199682;

    % Define boundary conditions.
    x0    = -2.5;
    x_end =  2.5;
    N     =  31250;

    % ******************************************************
    % Note: Drawing plots with large N is somehow consuming
    % all the RAM on my system and basically wiping out 
    % my free space with a gigantic swap file,
    % forcing a reboot. So, this is necessary ...
    % ******************************************************

    if (drawPlot>0)
        N = 2000;
    end

    % Define step conditions.
    h      = (x_end-x0)/N;
    x      = x0:h:x_end;
    x(1)   = [];
    x(end) = [];

    % //////////////////////////////////////////////////
    % Make the finite difference matrix.
    % //////////////////////////////////////////////////

    D = 2/h^2*eye(N-1) + diag(Vm(alpha,beta,gamma,x));  % Diagonals
    D = D - 1/h^2*diag(ones(N-2,1),1);                  % Upper diag
    D = D - 1/h^2*diag(ones(N-2,1),-1);                 % Lower diag
    D = sparse(D);

    [psi,e]=eigs(D,10,'sm');
    eigenvalues=diag(e);

    % Derive the energies.
    % Recall: Energy term is actually 2E / hbar^2.
    energies = sort(eigenvalues) / thh_fac;
    evals = [energies(1); energies(2); energies(3)];

    for ie=1:length(evals)
        fprintf('E(%i) = %10.5f\n',ie,evals(ie));
    end

    hold off
    clf

    % //////////////////////////////////////////////////
    % Plot the potential and eigenfunctions.
    % //////////////////////////////////////////////////

    % Want positive even functions, if possible (just a
    % sign choice, really.)
    psi_1 = psi(:,end) / norm(psi(:,end));
    psi_2 = psi(:,end-1) / norm(psi(:,end-1));
    psi_3 = psi(:,end-2) / norm(psi(:,end-2));

    % Make the odd function positve from the left.
    [max, max_pos] = max(psi_2);
    [min, min_pos] = min(psi_2);
    if (max_pos > min_pos)
        psi_2 = -1*psi_2;
    end
    clear max
    clear min
    clear max_pos
    clear min_pos
    
    psi   = [sign(sum(psi_1))*psi_1 psi_2 sign(sum(psi_3))*psi_3];

    if (drawPlot>0)
        if (iparam==1)
            xrange=-0.5:0.01:0.5;
        else
            xrange=-1.0:0.01:0.4;
        end
        plot(xrange,Vo(alpha,beta,gamma,xrange),'black-');
        hold on;
        for ie=1:3
            plot(x, evals(ie)*ones(1,length(x)), 'black','LineWidth',2);
            plot(x, overall_sf*psi(:,ie) + evals(ie)*ones(1,length(x)), 'red','LineWidth',2);
        end

        if (iparam==1)
            axis([xrange(1),xrange(end),0,20])
        else
            axis([xrange(1),xrange(end),-30,10])
        end
        xlabel('x (nm)');
        ylabel('E (eV)');
        set(gca,'fontsize',20);

        set(gcf, 'Renderer', 'painters');
        saveas(gca,['plots/Q1/psi_vs_x_paramset_' num2str(iparam) ...
                    '.eps'],'epsc');
    end
end

% //////////////////////////////////////////////////
% Potential function.
% //////////////////////////////////////////////////

function [V] = Vm(alpha,beta,gamma,x)
    thh_fac = 2 / 0.076199682;
    V = thh_fac*(alpha*x.^2 + gamma*x.^3 + beta*x.^4);
end

function [Vout] = Vo(alpha,beta,gamma,x)
    Vout = (alpha*x.^2 + gamma*x.^3 + beta*x.^4);
end
