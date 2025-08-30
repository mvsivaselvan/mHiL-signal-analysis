import casadi.*

addpath ../../data/
load ControllerTEXVS20020.mat % to get kact, beta, d
load datatab.mat % data for Tower E (5.8Hz, 2%) X, 0.75g CERL
rmpath ../../data/

% Some geometric parameters
larm = 15; % moment arm of simulator (in)
h = 80; %% tower height from cruciform, in
mLHE = 15; % lb load head mass
ILH = 2; % (lb-in-s2/rad) load head I
grav = 386.4; % in/s^2
Itower = 1000; % (lb-s^2/in)-in^2
m = Itower/larm^2;

% Sampling rate for measurement
fs_meas = 1000; % Hz

% construct ym_data and um_data
accel = (datatab.accel-datatab.accel(1));
wbar = ((datatab.Fxcalc+mLHE*accel)*h+...
        (datatab.Mycalc+ILH*accel/h))/larm; % conductor force after 
                                            % removing load head inertia
wbar = wbar - wbar(1);
t_meas = (0:length(accel)-1)'/fs_meas;

addpath ../time-varying-param-estim/ % for velocity_initial and prepdata
% construct initial guess for velocity by estimating from measured
% displacement and acceleration
if (~isfile('velocity_initial.mat'))
    addpath ../vel-est-derivative-comp/
    [x_est_ric, ~] = estimate_velocity_differential_riccati(t_meas, ...
                                                  datatab.displ, ...
                                                  accel*larm/h*grav, ...
                                                  100); %Q=100
    vel_initial = x_est_ric(:,2);
    clear x_est_ric
    rmpath ../vel-est-derivative-comp/
    save velocity_initial vel_initial
else
    load('velocity_initial.mat')
end

% construct initial guess for xv by integrating xvdot
sysxv = tf(alph,[1 alph]);
xv_initial = lsim(sysxv, datatab.uv-datatab.uv(1), t_meas);
clear sysxv

% measured data for interpolation
xm_data = [datatab.displ ... % actuator displacement
           vel_initial ... % estimated actuator velocity
           datatab.F ... % actuator force
           xv_initial ... % estimated xv
           datatab.State1 ... % controller state 1
           datatab.State2]'; % controller state 2
ym_data = [datatab.displ ... % actuator displacement
           accel*larm/h*grav ... % actuator accel
           datatab.F ... % actuator force
           datatab.State1 ... % controller state 1
           datatab.State2]'; % controller state 2
um_data = [wbar datatab.uext]';

% pre-processing low pass filter data and remove offset
% low-pass filtering is so that optimization does not attempt "silly" fast
% control input to track high-frequency content that model has no chance of
% tracking - nonlinearity, higher-mode dynamics and noise
xm_data = prepdata(xm_data, 32, fs_meas);
ym_data = prepdata(ym_data, 32, fs_meas);
um_data = prepdata(um_data, 32, fs_meas);
rmpath ../time-varying-param-estim/

% System matrices
Ac = controller.A;
Bc = controller.B;
Cc = controller.C;
Dc = controller.D;
Amat = zeros(6);
Amat(1,2) = 1;
Amat(2,3) = 1/m;
Amat(3,2:4) = [-kact -bet d];
Amat(4,[1 3:6]) = [alph*Dc -alph alph*Cc];
Amat(5:6,5:6) = Ac;
Amat(5:6,[1 3]) = Bc;
Bmat = zeros(6,2);
Bmat(2,1) = 1/m;
Bmat(4,2) = alph;
Cmat = [1 0 0   0 0 0; ... % displacement
        0 0 1/m 0 0 0; ... % acceleration
        0 0 1   0 0 0; ... % force
        0 0 0   0 1 0; ... % controller state 1
        0 0 0   0 0 1];    % controller state 2
Dmat = [0   0; ...   % displacement
        1/m 0; ... % acceleration
        0   0; ...   % force
        0   0; ...   % controller state 1
        0   0];      % controller state 2

xd = @(t,x)(nonlinear_ode(t,x,um_data,0.001,m,kact,bet,alph,Ac,Bc,Cc,Dc));
[T,X] = ode15s(xd, 0:0.001:11, zeros(6,1));
T = T'; X = X';
y_comp = Cmat*X + Dmat*um_data(:,1:11001);
y_comp(2,:) = y_comp(2,:) - 40*tanh(100*X(2,:))/m;
figure(102), 
    plot(t_meas(1:11001), ym_data(2,1:11001), T, y_comp(2,:))
return

% Create CasADi interpolants
xm_fun = interpolant('xm','linear',{t_meas}, xm_data(:));
ym_fun = interpolant('ym','linear',{t_meas}, ym_data(:));
um_fun = interpolant('um','linear',{t_meas}, um_data(:));

% Degree of interpolating polynomial
deg = 3;

% Get collocation points
tau = collocation_points(deg, 'legendre');

% Collocation linear maps
[C,D,B] = collocation_coeff(tau);

% Get an instance of Opti
opti = Opti();

% Dimensions
nx = 6; ny = 5; nu = 2;

% Build parametrization of flow function, g, in terms of hat basis
% functions
% ---- build knots and step
M = 5;
rknots = quantile_knots(xm_data(4,:),M);

% ---- Precompute symbolic basis evaluation function
x4_sym = MX.sym('x4');
Phi_sym = basis_eval_symbolic(x4_sym, rknots);  % deg+1 elements
Phi_fun = Function('Phi_fun', {x4_sym}, {Phi_sym});

% ---- coefficients, theta, of the hat basis functions 
% g(x) = \sum_i=1^M theta(i) phi_i(x)
theta = opti.variable(M,1);
theta0 = d * rknots(:); % initial values correspond to linear function g
opti.set_initial(theta, theta0);

% ---- monotonicity constraints
opti.subject_to( theta(2:end) >= theta(1:end-1) );

% ---- set simple bounds around init to help solver
% tol = 10; % e.g. allow 10x variation; tune as needed
% opti.subject_to( theta >= theta0 - tol*abs(theta0) );
% opti.subject_to( theta <= theta0 + tol*abs(theta0) );

% Cost function matrices
Q = opti.parameter(5,5);
opti.set_value(Q,diag([1e6 1 1 1 1]));

% Time horizon
T = 45;

% discretization
dt = 0.01;
N = T/dt;

% Decision variables and cost
X = cell(N+1,1);
J = 0;

% Initial condition
Xk = opti.variable(nx); X{1} = Xk; opti.set_initial(Xk, xm_fun(0));

for k = 1:N % over all the sample intervals
    if (mod(k,100)==0)
        fprintf('At time step %d\n',k);
    end
    
    % Collocation states
    Xc = cell(deg,1);
    for j = 1:deg
        Xc{j} = opti.variable(nx);
        tau_j = tau(j);
        t_j = (k-1 + tau_j)*dt;
        opti.set_initial(Xc{j}, xm_fun(t_j));
    end

    % Collocation equations
    % xp_j = Pi'(tau_j) = [Xk Xc_1 Xc_2 ... Xx_deg]*C(:,j)
    % Xk_end = [Xk Xc_1 Xc_2 ... Xx_deg]*D
    Xk_end = Xk*D(1);
    for j = 1:deg
        xp = Xk*C(1,j);
        for r = 1:deg
            xp = xp + Xc{r}*C(r+1,j);
        end
        
        t_j = (k-1 + tau(j))*dt;
        u_mj = um_fun(t_j); % fixed measured input
        
        x4_j = Xc{j}(4);  % 4th state is x4
        Phi_x4_j = Phi_fun(Xc{j}(4));
        g_j = theta' * Phi_x4_j;
        
        % Build dynamics manually with g_j in place of p3*x4
        xdot_j = MX(6,1);
        xdot_j(1) = Xc{j}(2);
        xdot_j(2) = (Xc{j}(3) + u_mj(1))/m;
        xdot_j(3) = -kact*Xc{j}(2) - bet*Xc{j}(3) + g_j; % nonlinear
        xdot_j(4) = alph*(Dc*[Xc{j}(1); Xc{j}(3)] ...
                           - Xc{j}(4) + Cc*Xc{j}(5:6) + u_mj(2));
        xdot_j(5:6) = Ac*Xc{j}(5:6) + Bc*[Xc{j}(1); Xc{j}(3)];        
        
        opti.subject_to(dt * xdot_j == xp);
        Xk_end = Xk_end + Xc{j}*D(j+1);
    end

    % New knot point
    Xk = opti.variable(nx); X{k+1} = Xk; opti.set_initial(Xk,xm_fun(k*dt));
    opti.subject_to(Xk == Xk_end);

    % Quadrature cost
    for j = 1:deg
        tau_j = tau(j);
        t_j = (k-1 + tau_j)*dt;

        y_mj = ym_fun(t_j);  % ny × 1
        u_mj = um_fun(t_j);  % nu × 1
        
        y_pred = Cmat * Xc{j} + Dmat * u_mj;
        
        e_y = y_pred - y_mj;
        
        integrand = e_y' * Q * e_y;
        J = J + B(j) * dt * integrand;
    end
end

% % ---- Add smoothness regularization to objective J
% lambda = opti.parameter(1,1);
% opti.set_value(lambda, 1e-2);
% for j=1:M-1
%     J = J + lambda * (theta(j+1) - theta(j))^2;
% end

opti.minimize(J);

opti.solver('ipopt');

sol = opti.solve();

X = [X{:}];
return
%% Plot optimal solution
x_opt = sol.value(X);
theta_opt = sol.value(theta);

% low-pass filter for post-processing measurement
[bLP, aLP] = butter(6, 32/fs_meas, 'low'); % 32 Hz cutoff

% high-pass filter for post-processing optimal sol
[bHP, aHP] = butter(6, 0.5*dt, 'high'); % 0.5 Hz cutoff

y_opt = Cmat*x_opt;

figure(101), 
    plot(t_meas, xm_data(1,:), (0:N)*dt, filtfilt(bHP,aHP,x_opt(1,:)))
    title('Displacement')

figure(102), 
    plot(t_meas, xm_data(2,:), (0:N)*dt, x_opt(2,:))
    title('Velocity')

figure(103), 
    plot(t_meas, ym_data(2,:), (0:N)*dt, y_opt(2,:)), grid on
    title('Acceleration')

figure(104),
    plot(t_meas, ym_data(3,:), (0:N)*dt, y_opt(3,:)), grid on
    title('Actuator force')

figure(105),
    plot(t_meas, ym_data(4,:), (0:N)*dt, y_opt(4,:)), grid on
    title('Controller state 1')

figure(106),
    plot(t_meas, ym_data(5,:), (0:N)*dt, y_opt(5,:)), grid on
    title('Controller state 2')

figure(107),
    plot(t_meas, xm_data(4,:), (0:N)*dt, x_opt(4,:)), grid on
    title('x_v')

Phi_mat = zeros(N+1, length(rknots));
for i = 1:N+1
    Phi_mat(i, :) = full(Phi_fun(x_opt(4,i)))';
end
g_vals = Phi_mat*theta_opt;

figure(301),
    plot(x_opt(4,:), g_vals, 'x');
    xlabel('x_4');
    ylabel('g(x_4)');
    title('Estimated Nonlinear Flow Function');
    grid on;

return
%% Parameter changes
opti.set_value(Q,diag([1e6 100 1 1 1]));
opti.set_value(lambda, 1e-4);
sol = opti.solve();
