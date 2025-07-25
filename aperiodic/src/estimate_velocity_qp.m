function [x_est, u_est, QP] = estimate_velocity_qp(t, x1_m, u_m, Q)
% Finite horizon velocity estimation using QP with trapezoidal integration.
% 
% Estimate velocity from measured displacement and acceleration
% by solving a QP (finite horizon).
%
% INPUTS:
%   t      - time vector (Nx1)
%   x1_m   - measured displacement (Nx1)
%   u_m    - measured acceleration (Nx1)
%   Q      - penalty weight on displacement tracking error
%
% OUTPUTS:
%   x_est  - estimated state [x1; x2] over time (Nx2)
%   u_est  - optimal control input (acceleration) (Nx1)
%   QP    : (optional) struct with QP matrices (H, f, Aeq, beq, exitflag)

N = length(t);
dt = t(2) - t(1);

% System matrices (double integrator):
% dx/dt = A*x + B*u, where x = [x1; x2], u = acceleration
A = [0 1; 0 0];
B = [0; 1];

% Total number of optimization variables:
% x = [x1; x2] ∈ ℝ^{2N}, u ∈ ℝ^N
nx = 2 * N;
nu = N;
nz = nx + nu;

% -------------------
% Construct Aeq matrix and beq vector for trapezoidal constraints
% -------------------
Aeq = zeros(2*(N-1), nz);
beq = zeros(2*(N-1), 1);

% Matrices for trapezoidal rule
Aplus  = eye(2) + (dt/2)*A;
Aminus = eye(2) - (dt/2)*A;
Bscaled = (dt/2)*B;

for k = 1:N-1
    % Indices
    xk_idx   = 2*(k-1) + (1:2);
    xkp1_idx = 2*k     + (1:2);
    uk_idx   = nx + k;
    ukp1_idx = nx + k + 1;

    % Fill constraint block row k
    % Aminus*x_{k+1} - Aplus*x_k - B*(u_k + u_{k+1})*(dt/2) = 0
    row_idx = 2*(k-1) + (1:2);
    Aeq(row_idx, xkp1_idx) = Aminus;
    Aeq(row_idx, xk_idx)   = -Aplus;
    Aeq(row_idx, uk_idx)   = -Bscaled;
    Aeq(row_idx, ukp1_idx) = -Bscaled;
end

% -------------------
% Construct cost function H anf f
% -------------------
% z = [x; u] ∈ ℝ^{2N + N}
% Objective: Q*(x1 - x1_m)^2 + (u - u_m)^2
% Preallocate nonzeros for H
nnzH = 2*N + N;  % N x1 terms, N u terms
rowsH = zeros(nnzH,1);
colsH = zeros(nnzH,1);
valsH = zeros(nnzH,1);
f     = zeros(nz, 1);

% Index counter
idx = 0;

% Penalize displacement tracking: Q(x1 - x1_m)^2
for k = 1:N
    i = 2*(k-1) + 1;  % x1 index in state vector

    idx = idx + 1;
    rowsH(idx) = i;
    colsH(idx) = i;
    valsH(idx) = 2 * Q;

    f(i) = -2 * Q * x1_m(k);
end

% Penalize acceleration tracking: (u - u_m)^2
for k = 1:N
    i = nx + k;  % u index

    idx = idx + 1;
    rowsH(idx) = i;
    colsH(idx) = i;
    valsH(idx) = 2;

    f(i) = -2 * u_m(k);
end

rowsH = rowsH(1:idx);
colsH = colsH(1:idx);
valsH = valsH(1:idx);
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
x_est = reshape(z(1:nx), [2, N])';
u_est = z(nx+1:end);

if nargout > 2
    QP.H = H;
    QP.f = f;
    QP.Aeq = Aeq;
    QP.beq = beq;
    QP.exitflag = exitflag;
end
end



