addpath ../time-varying-param-estim/
load flowfn.mat
rmpath ../time-varying-param-estim/

import casadi.*

% Build parametrization of flow function, g, in terms of hat basis
% functions
% ---- build knots and step
M = 11;
rknots = quantile_knots(xv,M);

% ---- Precompute symbolic basis evaluation function
x4_sym = MX.sym('x4');
Phi_sym = basis_eval_symbolic(x4_sym, rknots);  % deg+1 elements
Phi_fun = Function('Phi_fun', {x4_sym}, {Phi_sym});

Ndata = length(xv);
Phi_mat = zeros(Ndata, length(rknots));

for i = 1:Ndata
    Phi_mat(i, :) = full(Phi_fun(xv(i)))';
end

theta_fit = Phi_mat \ (g_flow');

figure(101),
    plot(xv, g_flow, 'x', xv, Phi_mat*theta_fit, 'o')
    grid on

p3_nom = 20600;
figure(102),
    plot(rknots, theta_fit, 'x', rknots, p3_nom*rknots(:), 'o');