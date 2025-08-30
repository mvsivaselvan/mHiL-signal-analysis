function xdot = nonlinear_ode(t, x, u, dt, m, kact, bet, alph, ...
                              Ac, Bc, Cc, Dc)

lu = size(u,2);
if t<0 || t>(lu-1)*dt
    u_ = zeros(2,1);
else
    idx = floor(t/dt)+1;
    u_ = u(:,idx) + (u(:,idx+1)-u(:,idx))/dt*(t - (idx-1)*dt);
end

xdot = zeros(6,1);
xdot(1) = x(2);
xdot(2) = (x(3) + u_(1) - 40*tanh(100*x(2)))/m;
xdot(3) = -kact*x(2) - bet*x(3) + 20600*x(4); % flowfun(x(4));
xdot(4) = alph*(Dc*x([1 3]) - x(4) + Cc*x(5:6) + u_(2));
xdot(5:6) = Ac*x(5:6) + Bc*x([1 3]);
