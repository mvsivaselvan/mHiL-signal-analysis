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

%% Compare load cell and pressure difference
F_p = (datatab.P1 - datatab.P2)*4.2; % force from pressure measurements
F_p = F_p - mean(F_p(1:100));

F = datatab.F;

% figure(101),
%     plot(t_meas, [F F_p])

r = 1; % increase sample rate by factor r

F_fine = interp(F, r);
F_p_fine = interp(F_p, r);

t_fine = (0:length(F_fine)-1)'/(r*fs_meas);

% figure(102),
%     plot(t_fine, [F_fine F_p_fine])

[c,lags] = xcorr(F_fine, F_p_fine);
[~,idx] = max(c);
shft = lags(idx);

figure(103),
    plot(t_fine(1:end+shft),[F_fine(1:end+shft) F_p_fine(1-shft:end)])
    
vel_fine = interp(vel_initial, r);

F_diff = F_p_fine-F_fine;

%% Beginning and end of the time series, friction is clear here
figure(104), 
    plot(vel_fine(7950:11000), F_diff(7950:11000), ...
         vel_fine(36000:38000), F_diff(36000:38000));
    
figure(105),
    plot(vel_fine(7950:11000), [F_diff(7950:11000) 100*tanh(100*vel_fine(7950:11000))])