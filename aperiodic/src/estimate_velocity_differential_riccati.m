function [x_est, u_est] = estimate_velocity_differential_riccati(t, x1_m, u_m, Q)
% Finite horizon velocity estimation using differential Riccati equation
% 
% Estimate velocity from measured displacement and acceleration
% by solving the time-varying Riccati differential equation (finite horizon).
%
% Inputs:
%   t    - uniform time vector (Nx1)
%   x1_m - measured displacement (Nx1)
%   u_m  - measured acceleration (Nx1)
%   Q    - scalar penalty on displacement tracking error
%
% Outputs:
%   x_est - estimated [displacement, velocity] (Nx2)
%   u_est - estimated acceleration (Nx1)

A = [0 1; 0 0];
B = [0; 1];
C = [1 0];
R = 1;

N = length(t);
dt = t(2) - t(1);

x1_m_vec = x1_m(:);
u_m_vec = u_m(:);

% Fast linear interpolation for uniform grid
function val = fast_interp(vec, time)
    if time <= t(1)
        val = vec(1);
        return;
    elseif time >= t(end)
        val = vec(end);
        return;
    end
    alpha = (time - t(1))/dt;
    idx = floor(alpha) + 1;
    frac = alpha - floor(alpha);
    val = (1 - frac)*vec(idx) + frac*vec(idx+1);
end

x1_m_fun = @(time) fast_interp(x1_m_vec, time);
u_m_fun = @(time) fast_interp(u_m_vec, time);

% Riccati ODE for packed P and s:
% z = [p11; p21; p22; s1; s2]
% P = [p11 p21; p21 p22], s = [s1; s2]

function dz = riccati_ode(t_, z)
    p11 = z(1); p21 = z(2); p22 = z(3);
    s1 = z(4); s2 = z(5);
    
    P = [p11 p21; p21 p22];
    s = [s1; s2];
    
    xm = x1_m_fun(t_);
    um = u_m_fun(t_);
    
    dP = -(A'*P + P*A - P*B*(R\B')*P + C'*Q*C);
    ds = -((A - B*(R\B')*P)'*s - C'*Q*xm + P*B*um);
    
    dz = zeros(5,1);
    dz(1) = dP(1,1);
    dz(2) = dP(2,1);
    dz(3) = dP(2,2);
    dz(4) = ds(1);
    dz(5) = ds(2);
end

% Terminal conditions at t = T
zT = zeros(5,1);

opts = odeset('RelTol',1e-6,'AbsTol',1e-9);
t_rev = flipud(t);

[~, Zrev] = ode45(@riccati_ode, t_rev, zT, opts);
Z = flipud(Zrev); % Nx5 matrix of [p11 p21 p22 s1 s2]

% Compute initial state from P(0) and s(0)
P0 = [Z(1,1) Z(1,2); Z(1,2) Z(1,3)];
s0 = Z(1,4:5)';
x0 = -P0 \ s0;

% Helper to interpolate P,s at arbitrary time t_ during forward integration
function val = fast_interp_vec(Zmat, time)
    if time <= t(1)
        val = Zmat(1, :)';
        return;
    elseif time >= t(end)
        val = Zmat(end, :)';
        return;
    end
    alpha = (time - t(1))/dt;
    idx = floor(alpha) + 1;
    frac = alpha - floor(alpha);
    val = (1 - frac)*Zmat(idx, :)' + frac*Zmat(idx+1, :)';
end

% Forward integration of state x(t)
function dx = state_ode(t_, x)
    z_interp = fast_interp_vec(Z, t_);
    P = [z_interp(1) z_interp(2); z_interp(2) z_interp(3)];
    s = z_interp(4:5);
    lambda = P*x + s;
    um = u_m_fun(t_);
    u = um - R \ (B' * lambda);
    dx = A*x + B*u;
end

[~, x_est] = ode45(@state_ode, t, x0, opts);

% Compute control input u_est at sampled points
u_est = zeros(N,1);
for k = 1:N
    z_sample = Z(k, :)';
    Pk = [z_sample(1) z_sample(2); z_sample(2) z_sample(3)];
    sk = z_sample(4:5);
    xk = x_est(k,:)';
    lambda = Pk*xk + sk;
    um = u_m_fun(t(k));
    u_est(k) = um - R \ (B' * lambda);
end

end
