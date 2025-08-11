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

% construct initial guess for velocity by estimating from measured
% displcaement and acceleration
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

% Dimensions
nx = 6; ny = 5; nu = 2; np = 1;

% Declare model variables
x1 = MX.sym('x');
x2 = MX.sym('v');
x3 = MX.sym('F');
x4 = MX.sym('xv');
x5 = MX.sym('xc1');
x6 = MX.sym('xc2');
u1 = MX.sym('w');
u2 = MX.sym('uext');
p3 = MX.sym('d');
t = MX.sym('t');
x = [x1; x2; x3; x4; x5; x6];
u = [u1; u2];

% Nominal fixed parameters
p1_nom = kact; 
p2_nom = bet;

% System matrices
Ac = controller.A;
Bc = controller.B;
Cc = controller.C;
Dc = controller.D;
Amat = MX(6,6);
Amat(1,2) = 1;
Amat(2,3) = 1/m;
Amat(3,2:4) = [-p1_nom -p2_nom p3];
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

% Model equation
xdot = Amat*x + Bmat*u;
f = Function('f',{x,u,p3},{xdot});

% Get an instance of Opti
opti = Opti();

% Cost function matrices
Q = opti.parameter(5,5);
opti.set_value(Q,diag([1e6 1 1 1 1]));
rho = opti.parameter(1,1);
opti.set_value(rho, 1e-6);
p3_nom = d;

% Time horizon
T = 45;

% discretization
dt = 0.01;
N = T/dt;

% Decision variables and cost
X = cell(N+1,1); P = cell(N,1);
J = 0;

% Initial condition
Xk = opti.variable(nx); X{1} = Xk; opti.set_initial(Xk, xm_fun(0));

for k = 1:N
    if (mod(k,100)==0)
        fprintf('At time step %d\n',k);
    end

    % Input and parameters for interval (k-1, k)
    Pk = opti.variable(np); P{k} = Pk; opti.set_initial(Pk, p3_nom);
    
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
        fj = f(Xc{j}, u_mj, Pk);        
        opti.subject_to(dt * fj == xp);
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
        e_p = Pk - p3_nom;
        
        integrand = e_y' * Q * e_y + e_p' * rho * e_p;
        J = J + B(j) * dt * integrand;
    end
end

alpha = opti.parameter(1,1);
opti.set_value(alpha, 1e2);
for k = 1:N-1
    J = J + alpha * (P{k+1} - P{k})^2 * dt;
end

opti.minimize(J);

opti.solver('ipopt');

sol = opti.solve();

X = [X{:}];
P = [P{:}];
return
%% Plot optimal solution
x_opt = sol.value(X);
p_opt = sol.value(P);

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

figure(303),
    plot((1:N)*dt, p_opt)
    title('d')
return
%% Parameter changes
opti.set_value(Q,diag([1e6 1 1 1 1]));
% opti.set_value(rho, 1e-6); % This setting does the job, but will silly
                             % fast d variation and apparent poor 
                             % conditioning-many iterations for convergence 
opti.set_value(rho, 1e-6);
opti.set_value(alpha, 0.1);
sol = opti.solve();

%% Continue
opti.set_initial(sol.value_variables());
opti.set_value(Q,diag([1e6 1 1 1 1]));
opti.set_value(R, diag([1 100]));
opti.set_value(rho, diag([1 1 1e-12]));
sol = opti.solve();

%% Forward simulation
fA = casadi.Function('fA',{p},{Amat});
A__ = full(fA(p_nom));
sys = ss(A__,Bmat,[Cmat; 0 0 0 1 0 0],[Dmat; 0 0]);
x_meas = full(xm_fun((0:N)*dt));
u_meas = full(um_fun((1:N)*dt));
y_comp_meas = lsim(sys, u_meas, (1:N)*dt);
y_comp_opt = lsim(sys, u_opt, (1:N)*dt);
y_comp_mix = lsim(sys, [u_meas(1,:); u_opt(2,:)], (1:N)*dt);
y_meas = full(ym_fun((1:N)*dt));
figure(501), 
    plot((1:N)*dt, y_comp_meas(:,2), 'w-', ...
         (1:N)*dt, y_comp_mix(:,2), 'r-', ...
         (1:N)*dt, y_meas(2,:),'b-'), 
    title('Acceleration')
    grid on
uv_opt = Cc*x_opt([5 6],:) + Dc*x_opt([1 3],:) + [u_opt(2,:) 0];
sysxv = tf(alph,[1 alph]);
xv_opt = lsim(sysxv, uv_opt, (0:N)*dt);
figure(502), 
    plot((0:N)*dt, x_opt(4,:), 'b-', ...
         (1:N)*dt, y_comp_opt(:,6), 'r-', ...
         (0:N)*dt, xv_opt, 'k'),
    title('x_v')
    grid on

%% Compare physical valve commands
figure(601),
    % plot((0:N)*dt,x_meas(4,:),(0:N)*dt,x_opt(4,:))
    plot(t_meas,xm_data(4,:),(0:N)*dt,x_opt(4,:))