function q = flowfun(xv)

rknots = [-2 -0.5 -0.25 0.25 0.5 2]';
theta = [-58950 -8000 -2950 2950 8000 58950]';
Phi = hatfun(xv, rknots);
q = Phi*theta;

    function Phi = hatfun(xv, rknots)
        Phi = zeros(length(xv),length(rknots));
        
        % left extrapolation
        r1 = rknots(1);
        r2 = rknots(2);
        mask = xv <= r1;
        h = r2 - r1;
        Phi(mask,1) = (r2-xv(mask))/h;
        Phi(mask,2) = (xv(mask)-r1)/h;

        % right extrapolation
        r1 = rknots(end-1);
        r2 = rknots(end);
        mask = xv > r2;
        h = r2 - r1;
        Phi(mask,end) = (xv(mask)-r1)/h;
        Phi(mask,end-1) = (r2-xv(mask))/h;

        % interior points
        for n = 1:length(rknots)-1
            r1 = rknots(n);
            r2 = rknots(n+1);
            mask = (xv > r1) & (xv <= r2);
            h = r2 - r1;
            Phi(mask,n) = (r2 - xv(mask))/h;
            Phi(mask,n+1) = (xv(mask) - r1)/h;
        end
    end

end

% k1 = 1.18e4;
% k2 = 3.2e4;
% delt = 0.25;

% q = zeros(size(xv));
% 
% q(xv<=delt) = k1*xv(xv<=delt);
% q(xv>delt) = k1*delt + k2*(xv(xv>delt) - delt);
% q(xv<-delt) = -k1*delt + k2*(xv(xv<-delt) + delt);
