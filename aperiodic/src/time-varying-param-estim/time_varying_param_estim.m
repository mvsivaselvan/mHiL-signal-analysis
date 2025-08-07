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
t_meas = (0:length(accel)-1)'/fs_meas;

% construct initial guess for velocity by estimating from measured
% displcaement and acceleration
addpath ../vel-est-derivative-comp/
[x_est_ric, ~] = estimate_velocity_differential_riccati(t_meas, ...
                                                  datatab.displ, ...
                                                  accel*larm/h*grav, ...
                                                  100); %Q=100
vel_initial = x_est_ric(:,2);
clear x_est_ric
rmpath ../vel-est-derivative-comp/

% construct initial guess for xv by integrating xvdot
sysxv = tf(alph,[1 alph]);
xv_initial = lsim(sysxv, datatab.uv, t_meas);
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
um_data = [datatab.uext wbar]';

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
nx = 6; ny = 5; nu = 2; np = 3;

% Declare model variables
x1 = MX.sym('x');
x2 = MX.sym('v');
x3 = MX.sym('F');
x4 = MX.sym('xv');
x5 = MX.sym('xc1');
x6 = MX.sym('xc2');
u1 = MX.sym('w');
u2 = MX.sym('uext');
p1 = MX.sym('kact');
p2 = MX.sym('bet');
p3 = MX.sym('d');
t = MX.sym('t');
x = [x1; x2; x3; x4; x5; x6];
u = [u1; u2];
p = [p1; p2; p3];

% System matrices
Ac = controller.A;
Bc = controller.B;
Cc = controller.C;
Dc = controller.D;
Amat = MX(6,6);
Amat(1,2) = 1;
Amat(2,3) = 1/m;
Amat(3,2:4) = [-p1 -p2 p3];
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
f = Function('f',{x,u,p},{xdot});

% Get an instance of Opti
opti = Opti();

% Cost function matrices
Q = opti.parameter(5,5);
opti.set_value(Q,diag([1 1 1 1 1]));
R = opti.parameter(2,2);
opti.set_value(R, diag([1 1]));
rho = opti.parameter(3,3);
opti.set_value(rho, diag([1 1 1]));
p_nom = [kact; bet; d];

% Cost function
e_y = Cmat*x + Dmat*u - ym_fun(t);
e_u = u - um_fun(t);
e_p = p - [kact; bet; d];
L = e_y'*Q*e_y + e_u'*R*e_u + e_p'*rho*e_p;

% Time horizon
T = 45;

% discretization
dt = 0.01;
N = T/dt;

% Decision variables and cost
X = cell(N+1,1); U = cell(N,1); P = cell(N,1);
J = 0;

% Initial condition
Xk = opti.variable(nx); X{1} = Xk; opti.set_initial(Xk, xm_fun(0));

for k = 1:N
    if (mod(k,100)==0)
        fprintf('At time step %d\n',k);
    end

    % Input and parameters for interval (k-1, k)
    Uk = opti.variable(nu); U{k} = Uk; opti.set_initial(Uk, um_fun(k*dt));
    Pk = opti.variable(np); P{k} = Pk; opti.set_initial(Pk, p_nom);
    
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
        fj = f(Xc{j}, Uk, Pk);
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
        
        y_pred = Cmat * Xc{j} + Dmat * Uk;
        
        e_y = y_pred - y_mj;
        e_u = Uk - u_mj;
        e_p = Pk - p_nom;
        
        integrand = e_y' * Q * e_y + ...
                    e_u' * R * e_u + ...
                    e_p' * rho * e_p;
        J = J + B(j) * dt * integrand;
    end
end

opti.minimize(J);

opti.solver('ipopt');

sol = opti.solve();

X = [X{:}];
U = [U{:}];
P = [P{:}];

%% Plot optimal solution
x_opt = sol.value(X);
u_opt = sol.value(U);
p_opt = sol.value(P);

y_opt = Cmat*x_opt;

figure(101), 
    plot(t_meas, datatab.displ, (0:N)*dt, x_opt(1,:))

figure(102), 
    plot(t_meas, vel_initial, (0:N)*dt, x_opt(2,:))

figure(103), 
    plot(t_meas, accel/h*larm*grav, (0:N)*dt, y_opt(2,:)), grid on

figure(201),
    plot(t_meas, um_data(1,:), (1:N)*dt, u_opt(1,:))

figure(202),
    plot(t_meas, um_data(2,:), (1:N)*dt, u_opt(2,:))