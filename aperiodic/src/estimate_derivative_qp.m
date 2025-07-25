function [x_est, u_est, QP] = estimate_derivative_qp(t, x_m, Q)
% Finite horizon derivative estimation using QP with trapezoidal integration.
%
% Estimate derivative of a signal using noisy position measurements.
%
% INPUTS:
%   t    - time vector (Nx1)
%   x_m  - measured signal (Nx1)
%   Q    - penalty weight on signal tracking error
%
% OUTPUTS:
%   x_est  - estimated state [x; dx] over time (Nx2)
%   u_est  - optimal control input (derivative estimate) (Nx1)
%   QP     - (optional) struct with QP matrices (H, f, Aeq, beq, exitflag)

N = length(t);
dt = t(2) - t(1);

% System matrices: dx/dt = u
A = 0;   % scalar system: dx/dt = u
B = 1;

% Number of decision variables
nx = N;       % x (signal)
nu = N;       % u (derivative)
nz = nx + nu; % total variables

% -------------------
% Construct Aeq and beq for trapezoidal integration
% -------------------
Aeq = zeros(N-1, nz);
beq = zeros(N-1, 1);

for k = 1:N-1
    xk   = k;
    xkp1 = k+1;
    uk   = nx + k;
    ukp1 = nx + k + 1;

    row = k;
    Aeq(row, xkp1) =  1;
    Aeq(row, xk)   = -1;
    Aeq(row, uk)   = -dt/2;
    Aeq(row, ukp1) = -dt/2;
end

% -------------------
% Construct cost function H and f
% -------------------
% Objective: Q*(x - x_m)^2 + u^2
nnzH = N + N;  % N x terms + N u terms
rowsH = zeros(nnzH,1);
colsH = zeros(nnzH,1);
valsH = zeros(nnzH,1);
f     = zeros(nz,1);

idx = 0;

% Signal tracking term
for k = 1:N
    i = k;  % x index
    idx = idx + 1;
    rowsH(idx) = i;
    colsH(idx) = i;
    valsH(idx) = 2 * Q;

    f(i) = -2 * Q * x_m(k);
end

% Derivative penalty term
for k = 1:N
    i = nx + k;  % u index
    idx = idx + 1;
    rowsH(idx) = i;
    colsH(idx) = i;
    valsH(idx) = 2;

    % No linear term for u
end

H = sparse(rowsH, colsH, valsH, nz, nz);

% -------------------
% Solve the QP
% -------------------
opts = optimoptions('quadprog', 'Display', 'off');
[z, ~, exitflag] = quadprog(H, f, [], [], Aeq, beq, [], [], [], opts);

if exitflag ~= 1
    error('QP failed to solve.');
end

% -------------------
% Extract estimates
% -------------------
x_est = z(1:nx);
u_est = z(nx+1:end);

if nargout > 2
    QP.H = H;
    QP.f = f;
    QP.Aeq = Aeq;
    QP.beq = beq;
    QP.exitflag = exitflag;
end

end
