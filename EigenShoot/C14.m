
% //////////////////////////////////////////////////////
% C14_6
% Shooting method for particle in a box eigenfunctions
% and eigenvalues
% Use C14_4 to find psi for particle in a box
% //////////////////////////////////////////////////////

clear all
clf

E=9.87;
%V=@(x)0; %no potential energy

V1=-5;
V2=50;
V=@(x)C14_7(x,V1,V2);

[x,u]=C14_4(V,E);
figure(1),clf
plot(x,u(:,1))
ylabel('\psi')
xlabel('x')
grid on

% Look at boundary condition (psi(1)) as a function of E
f=@(E)C14_5(V,E); %function handle to get psi(a)

EE=linspace(0,50,500);
psi_a=0*EE;
for i=1:numel(EE)
    psi_a(i)=f(EE(i));
end
figure(2),clf
plot(EE,psi_a)
grid on
xlabel('E')
ylabel('\psi(1)')

% Find energy eigenvalues as a function of E using root finding
clear E n_iter
[E(1),n_iter(1)]=fzero(f,0);
[E(2),n_iter(2)]=fzero(f,7);
[E(3),n_iter(3)]=fzero(f,15);
[E(4),n_iter(4)]=fzero(f,30);
E
E/pi^2*4

% Find the eigenfunctions at each of these eigenvalues
figure(3),clf,hold on
for i=1:4
    [x,u]=C14_4(V,E(i));
    psi=u(:,1)/norm(u(:,1)); %normalized wave function
    plot(x,18*psi+E(i));
    line([-1 1],E(i)*[1,1],'color','k','linewidth',0.5)
end
xlabel('x')
ylabel('\psi')
legend('1','2','3','4','location','best')

xx=linspace(-1,1,2000);
plot(xx,V(xx),'.k')
xlabel('x')
ylabel('V_{new}')

% //////////////////////////////////////////////////
% End function.
% //////////////////////////////////////////////////

return

%C14_4
%Particle in a box, assuming a=1
%V is a function handle for a potential defined on [-1,1]
%Solve the IVP with u=[psi;psi'] and u(-1)=[0;1]
function [x,u]=C14_4(V,E)
    
    %NOTE: This function takes hbar^2 / 2 = 1, not necessarily
    %true.
    u0=[0;1]; %initial conditions, at x=-1
    f=@(x,u)[u(2);(V(x)-E)*u(1)];
    [x,u]=ode45(f,[-1 1],u0);
    
end


%C14_5
%Call C14_4 and extract psi(a)
function psi_a=C14_5(V,E)
    [x,u]=C14_4(V,E);
    psi_a=u(end,1);
end

function [] = C14_6()
    
    % Use C14_4 to find psi for particle in a box
    E=9.87;
    V=@(x)0; %no potential energy
    [x,u]=C14_4(V,E);
    figure(1),clf
    plot(x,u(:,1))
    ylabel('\psi')
    xlabel('x')
    grid on

    % Look at boundary condition (psi(1)) as a function of E
    f=@(E)C14_5(V,E); %function handle to get psi(a)

    EE=linspace(0,50,500);
    psi_a=0*EE;
    for i=1:numel(EE)
        psi_a(i)=f(EE(i));
    end
    figure(2),clf
    plot(EE,psi_a)
    grid on
    xlabel('E')
    ylabel('\psi(1)')
    %% Find energy eigenvalues as a function of E using root finding
    clear E n_iter
    [E(1),n_iter(1)]=P1_rootfinder(f,0,7);
    [E(2),n_iter(2)]=P1_rootfinder(f,7,15);
    [E(3),n_iter(3)]=P1_rootfinder(f,15,30);
    [E(4),n_iter(4)]=P1_rootfinder(f,30,45);
    E
    E/pi^2*4
    
    % Find the eigenfunctions at each of these eigenvalues
    figure(3),clf,hold on
    for i=1:4
        [x,u]=C14_4(V,E(i));
        psi=u(:,1)/norm(u(:,1)); %normalized wave function
        plot(x,18*psi+E(i));
        line([-1 1],E(i)*[1,1],'color','k','linewidth',0.5)
    end
    xlabel('x')
    ylabel('\psi')
    legend('1','2','3','4','location','best')

    % Consider a different potential function. Still in a box

    xx=linspace(-1,1,2000);
    figure(4),clf
    V1=-5;
    V2=50;
    Vnew=@(x)C14_7(x,V1,V2);
    plot(xx,Vnew(xx),'.k')
    xlabel('x')
    ylabel('V_{new}')

    % Look at psi(1) boundary condition to get rough roots
    f=@(E)C14_5(Vnew,E); %function handle to get psi(a)
    EE=linspace(min([V1,0]),2*max([V1,V2]),200);
    dE=EE(2)-EE(1);
    psi_a=0*EE;
    for i=1:numel(EE)
        psi_a(i)=f(EE(i));
    end
    figure(5),clf
    plot(EE,psi_a)
    grid on

    ylim([-10 10])
    xlabel('E')
    ylabel('\psi(1)')
    rough_root_inds=find(psi_a(2:end).*psi_a(1:end-1)<0);
    rough_roots=EE(rough_root_inds)

    % Find energy eigenvalues as a function of E using root finding
    clear E n_iter
    f=@(E)C14_5(Vnew,E);
    for i=1:min([numel(rough_roots),8])
        [E(i),n_iter(i)]=P1_rootfinder(f,rough_roots(i)-dE,rough_roots(i)+dE);
    end
    E

    figure(6),clf, hold on
    clear p
    for i=1:numel(E)
        [x,u]=C14_4(Vnew,E(i));
        u(:,1)=u(:,1)/norm(u(:,1)); %normalize wave functions

        p(i)=plot(x,30*u(:,1)+E(i));
        plot(xx,Vnew(xx),'k','linewidth',3)
        l=line([-1 1],E(i)*[1 1],'linewidth',0.5,'color','k');
        xlabel('x')
        ylabel('\psi')
        title(['\psi_' num2str(i) ', E=' num2str(E(i),'%0.3g')])
        pause
        delete(l)
    end
end

%C14_7
%More complicated potential function for particle in a box
%V=0 for -1<x<-2/3
%V=V1 for 1/2<x<1
%between -2/3 and 1/2, V is a quadratic function matching the boundary
%conditions and having maximum height V2
function V=C14_7(x,V1,V2)
%Write the function so that it can take a vector of x values
    V=0*x;
    V(x>=1/2)=V1;
    if V1==V2 && V1==0,    return, end
    inds = x>-2/3 & x<1/2;
    alpha = 36/49*(V1-2*V2-2*sqrt(V2*(V2-V1)));
    beta = 1/2-6*V1/7/alpha;
    V(inds)=alpha*(x(inds)+2/3).*(x(inds)-beta);
    V=real(V);
end