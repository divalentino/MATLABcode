
% Using the bisection method to compute the root
% based on some input guess.
function [rn, niter, relerr] = bisect(r0, r1, tolerance)

    max_iter = 100; %maximum number of iterations
    
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
        if dFx(rmid)*dFx(r0)<0
            r1=rmid;
        elseif dFx(rmid)*dFx(r1)<0
            r0=rmid;
        else
            error('no sign change')
        end
        err=abs(r1-r0)/rmid;
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