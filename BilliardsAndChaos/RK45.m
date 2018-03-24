function [xx,yy,n_failures,xsave,ysave]=RK45(f,xrange,y0,tol,useDebug)

% Define start and end points for RK process.
x0    = xrange(1);
x_end = xrange(end);
xsave = x0;
ysave = y0;

useRange = 0;
if (length(xrange)>2)
    fprintf(['NOTE: Returning RK evaluated for a specific set ' ...
             'of points!\n']);
    useRange = 1;
end

c = [0   1/4    3/8   12/13 1 1/2];
b5 =[16/135 0   6656/12825  28561/56430     -9/50   2/55];
b4 =[25/216 0   1408/2565   2197/4104   -1/5 0];
A=zeros(6);
A(2,1)=1/4;
A(3,1:2)=[3/32 9/32];
A(4,1:3)=[1932/2197   -7200/2197  7296/2197];
A(5,1:4)=[439/216     -8      3680/513    -845/4104];
A(6,1:5)=[-8/27   2    -3544/2565     1859/4104   -11/40];

%ensure that y0 is a column vector
y0=y0(:);

xx=x0; %initialize
yy=y0; %initialize

%Find the direction of propagation
if x_end>=x0
    xdir=1;
else
    xdir=-1;
end

%guess an initial step size
h=min([abs(x_end-x0)/10,0.1]); %h is always positive 

x=x0;y=y0;
done=false; 
n_failures=0;
dimy=numel(y0);

%Main Loop
xiter=2;

while ~done
    minh=16*eps(x);     %minimum acceptable value of h
    failures=false;     %no step-size failures yet
    
    %Make sure to hit last step exactly
    if abs(x_end-x)<=h
        h=abs(x_end-x);
        done=true;
    end
    
    if (useRange>0) && (xiter >= length(xrange))
        break;
    end
    
    %Loop for advancing one step
    while true
        fn=zeros(6,dimy);                
        
        % ////////////////////////////////////////////////////////////
        % See pg. 223 of textbook.
        % Recursively calculates f(n+1) based on the results of f(n).
        % ////////////////////////////////////////////////////////////
        
        % Please study these lines carefully. How do they work? 
        for i=1:6
            fn(i,:)=f(x+xdir*h*c(i),y + xdir*h* (A(i,:)*fn).'  ).';
        end        
        y5=y+xdir*h*(b5*fn).';
        y4=y+xdir*h*(b4*fn).';

        % Relative error between 4th, 5th order evaluations.
        err_rel = norm((y5-y4)./y5);
        
        if err_rel<tol
            
            if (useDebug>0)
                fprintf('Accepting step size!\n')
            end
            
            % New: If we're within a range that one of our
            % desired x values, take h to be abs(xval-x),
            % and save that intermediary (x,y) pair.
            if (useRange>0)
                if (xrange(xiter) > x) && ...
                        (xrange(xiter) < (x+h))
                    if (useDebug>0)
                        fprintf(['Got point: %i, t = %20.15f ' ...
                                 '\n'],xiter,xrange(xiter));
                        pause;
                    end
                    if (mod(xiter,500)==0)
                        fprintf('Processed %i / %i points \n',xiter,length(xrange));
                    end
                    hs    = abs(xrange(xiter)-x);
                    y5s   = y+xdir*hs*(b5*fn).';
                    y4s   = y+xdir*hs*(b4*fn).';
                    xiter = xiter+1;
                    xs    = x+xdir*hs;
                    ys    = y5s;
                    xsave = [xsave,xs];
                    ysave = [ysave,ys];
                end
            end
            
            %accept step and update step size. Do not increase step size by
            %more than a factor of ten nor reduce it by more than a factor
            %of 2.
            x=x+xdir*h;
            y=y5;
            xx=[xx,x];
            yy=[yy,y];
            
            % //////////////////////////////////////////////////
            % Start the next iteration with a coarser 
            % or equivalent stepsize.
            % //////////////////////////////////////////////////
            
            hnew = 0.9 * h * (tol/abs(err_rel))^0.2;
            if (hnew >= h)
                h = hnew;
            end
            break;
            
        else

            % *** Need this line to make sure we don't stop 
            % short of the end point.
            done = false;
            
            if (useDebug>0)
                fprintf('Rejecting step size!\n')
            end
            
            %reject step. Reduce h (by at most a factor of 2) and try again 
            if ~failures
                
                % //////////////////////////////////////////////////
                % Step size is too big, calculate a smaller value.
                % //////////////////////////////////////////////////
                
                if (err_rel>0)
                    hnew = 0.9 * h * (tol/abs(err_rel))^0.2;
                else
                    hnew = h;
                end
                if (hnew < h)
                    h = hnew;
                else
                    h = h/2;
                end
                
                %Make sure we don't go beyond the minimum.
                if (h < minh)
                    h = minh;
                end
                
                failures   = true;
                n_failures = n_failures+1;
            else
                h=h/2;
            end
            
            if (useDebug>0)
                fprintf('New smaller step size: %20.15f \n',h)
            end
        end
    end
end


