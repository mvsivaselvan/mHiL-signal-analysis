function [x_est, u_est] = estimate_derivative_differential_riccati_analytic(t, x_m, rho)
% Estimate derivative of x_m using a Riccati approach with analytical P(t)
%   Cost = ∫ (x - x_m)^2 + rho*u^2 dt
%   System: dx/dt = u
%
% Inputs:
%   t    - time vector (Nx1), uniformly spaced
%   x_m  - measured signal (Nx1)
%   Q    - state tracking weight (scalar)
%   rho  - control effort weight (scalar)
%
% Outputs:
%   x_est - estimated state (Nx1)
%   u_est - estimated derivative (Nx1)

dt = t(2) - t(1);
x_m_vec = x_m(:);

% === Step 1: Analytical Riccati solution ===
% Riccati equation: dP/dt = -2 + P^2 / (2*rho)
% With terminal condition P(T) = 0

% Let a = sqrt(2/rho), then:
T = t(end);
tau = T - t;  % time to go

% Analytical solution: P(t) = 2\sqrt{\rho}*tanh((T - t)/\sqrt{rho)
P_vec = 2*sqrt(rho)*tanh(tau/sqrt(rho));

% === Step 2: Interpolate x_m for continuous evaluation ===
function xm_val = xm_interp(t_)
    if t_ <= t(1)
        xm_val = x_m_vec(1);
    elseif t_ >= t(end)
        xm_val = x_m_vec(end);
    else
        idx = floor((t_ - t(1)) / dt) + 1;
        alpha = (t_ - t(idx)) / dt;
        xm_val = (1 - alpha)*x_m_vec(idx) + alpha*x_m_vec(idx+1);
    end
end

% === Step 3: Solve for s(t) backward using ode45 ===
function dsdt = s_ode(tau_, s)
    t_ = T - tau_;  % map to tau for P(t)
    P_now = 2*sqrt(rho)*tanh(t_/sqrt(rho));
    xm_now = xm_interp(t_);
    dsdt = -((1/(2*rho)) * P_now * s + 2 * xm_now);
end

% Integrate backward in time
opts = odeset('RelTol',1e-9,'AbsTol',1e-12);
[~, s_rev] = ode15s(@s_ode, t, 0, opts);
s_vec = flipud(s_rev);  % Match original time order

% === Step 4: Compute x0 = -s(0)/P(0) ===
x0 = -s_vec(1) / P_vec(1);

% === Step 5: Forward integrate x(t) using dx/dt = u = -(Px + s)/(2rho) ===
function dxdt = x_ode(t_, x)
    tau_ = T - t_;
    P_now = 2*sqrt(rho)*tanh(tau_/sqrt(rho));

    % Interpolate s(t) at t_
    if t_ <= t(1)
        s_now = s_vec(1);
    elseif t_ >= t(end)
        s_now = s_vec(end);
    else
        idx = floor((t_ - t(1)) / dt) + 1;
        alpha = (t_ - t(idx)) / dt;
        s_now = (1 - alpha)*s_vec(idx) + alpha*s_vec(idx+1);
    end

    dxdt = - (P_now * x + s_now) / (2*rho);
end

[~, x_est] = ode15s(@x_ode, t, x0, opts);
x_est = x_est(:);  % Ensure column

% === Step 6: Compute u_est ===
u_est = - (P_vec .* x_est + s_vec) / (2*rho);

end

