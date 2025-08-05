addpath ../vel-est-derivative-comp/
addpath ../../data/

load ControllerTEXVS20020.mat % to get kact, beta, d
load datatab.mat %%% data for Tower E (5.8Hz, 2%) X, 0.75g CERL

%%% tower data
mbar = 0.4; %% (lb-s2/in) 
hCM = 55;%%% in, from the Center of Rotation
Itower = 1000; % (lb-s^2/in)-in^2
IVS = Ibushing;
larm = 15; % moment arm of simulator (in)
kVS = kVSrot/larm^2; % lb/in
mVS = IVS/larm^2; % lb-s^2/in
w = sqrt(kVS/mVS); % rad/s
mtower = Itower/larm^2;
m = mtower;
h = 80; %% tower height from cruciform, in
mLHE = 15; % lb load head mass
ILH = 2; % (lb-in-s2/rad) load head I
%% target VS
AVS = [0 1; -w^2 -2*zet*w];
BVS = [0 0 ;1/mVS -mbar*hCM/mVS/larm]; %% (mbar*hCM/mVS/larm)
CVS = [1 0; 0 1; -w^2 -2*zet*w];
DVS = [0 0; 0 0; 1/mVS -mbar*hCM/mVS/larm];
ssVS = ss(AVS, BVS, CVS, DVS);
ssVS.InputName = {'w','weq'};
ssVS.OutputName = {'xVS', 'vVS','aVS'};
%% AVS model with no fit parameters
Ac = controller.A;
Bc = controller.B;
Cc = controller.C;
Dc = controller.D;
A = zeros(6);
A(1,2) = 1;
A(2,3) = 1/m;
A(3,2:4) = [-kact -bet d];
A(4,[1 3:6]) = [alph*Dc -alph alph*Cc];
A(5:6,5:6) = Ac;
A(5:6,[1 3]) = Bc;
B = zeros(6,2);
B(2,1) = 1/m;
B(4,2) = alph;
C = [eye(6); A(2,:)];
D = [zeros(6,2); B(2,:)];
ssActCont = ss(A, B, C, D);
ssActCont.InputName = {'w','uext'};
ssActCont.OutputName = {'x','v','F','xv','xC1','xC2','a'};

S1 = sumblk('e = xVS - x');

ssAVS=connect(ssVS, ssActCont,S1,{'w','uext'},{'x','v','F','xv','xC1','xC2','e','a'});

%% equation by equation check - estimating velocity from displacement and acceleration
%%%%% not using qp - matrix too large
x_m = datatab.displ;
u_m = (datatab.accel-datatab.accel(1))*larm/h*386.4;
Q = 100;
t = datatab.t;
% [x_est_qp, a_est_qp] = estimate_velocity_qp(t, x_m, u_m, Q); %%% using qp
[x_est_ric, a_est_ric] = estimate_velocity_differential_riccati(t, x_m, u_m, Q); %%% using riccati

figure(201),
    plot(t, [x_m, x_est_ric(:,1)]),
    grid on
    legend('measured', 'ric'),
    title('displacement')
figure(202),
    plot(t, x_est_ric(:,2)),
    grid on
    legend('ric')
    title('velocity')
figure(203),
    plot(t,u_m,'r',t, a_est_ric,'b--'),
    grid on
    legend('measured', 'ric')
    title('acceleration')
%% estimate Fdot
Q = 1e5;
x_m = datatab.F;
% [F_est_qp, Fdot_est_qp] = estimate_derivative_qp(t, x_m, Q); %%% not
% using qp
[F_est_ric, Fdot_est_ric] = estimate_derivative_differential_riccati_analytic(t, x_m, 1/Q);
figure(101),
    plot(t, Fdot_est_ric),
    grid on
    legend('exact', 'ric'),
    title('derivative')
figure(102),
    plot(t, [x_m F_est_ric]),
    grid on
    legend('exact', 'ric')
    title('signal')
%% not calculating the states right now (S2dot is an issue)
Q = 1e5;
% [S1_est_qp, S1dot_est_qp] = estimate_derivative_qp(t, State1, Q);
[S1_est_ric, S1dot_est_ric] = estimate_derivative_differential_riccati_analytic(t, datatab.State1, 1/Q);
figure(401),
    plot(t, S1dot_est_ric),
    grid on
    legend('qp', 'ric'),
    title('derivative')
figure(402),
    plot(t, [State1 S1_est_ric]),
    grid on
    legend('exact', 'ric')
    title('signal')

% [S2_est_qp, S2dot_est_qp] = estimate_derivative_qp(t, State2, Q);
[S2_est_ric, S2dot_est_ric] = estimate_derivative_differential_riccati_analytic(t, datatab.State2, 1/Q);
figure(501),
    plot(t, S2dot_est_ric),
    grid on
    legend('ric'),
    title('derivative')
figure(502),
    plot(t, [State2 S2_est_ric]),
    grid on
    legend('exact','ric')
    title('signal')
%% crude way of checking equation 3
sys = tf(alph,[1 alph]);
xv = lsim(sys, datatab.uv, t);
lhs3 = Fdot_est_ric;
rhs3 = -kact*x_est_ric(:,2)-bet*datatab.F+d*xv;

figure,
plot(t,[lhs3, rhs3-mean(rhs3)]), grid on,
%% equation by equation check for only first 3 equations
sys = tf(alph, [1 alph]);
xv = lsim(sys, datatab.uv, t);
thetaddot = datatab.accel/h;
wbar = ((datatab.Fxcalc+mLHE*datatab.accel)*h+(datatab.Mycalc+ILH*thetaddot))/larm; % conductor force after removing load head inertia
xmeas = [datatab.displ x_est_ric(:,2) datatab.F]; % bvpsol.y(:,1:N)';

xmeasdot = [x_est_ric(:,2) u_m Fdot_est_ric];

usys_ = [wbar xv];
rhs = xmeas(:,1:3)*A(1:3,1:3)' + usys_*[0 1/m 0; 0 0 d];

figure(21), plot(t,xmeasdot(:,1),'r',t, rhs(:,1),'b--'), grid on, title('equation 1')
figure(22), plot(t, xmeasdot(:,2)/386.4, 'r', t, (rhs(:,2)-mean(rhs(:,2)))/386.4,'b--'), grid on, title('equation 2')
figure(23), plot(t, xmeasdot(:,3), 'r',t, rhs(:,3)-mean(rhs(:,3)),'b--'), grid on, title('equation 3')
%% Run this only if state2dot is available
sys = tf(alph, [1 alph]);
xv = lsim(sys, datatab.uv, t);
thetaddot = datatab.accel/h;
xvdot = -alph*(xv-datatab.uv);
wbar = ((datatab.Fxcalc+mLHE*datatab.accel)*h+(datatab.Mycalc+ILH*thetaddot))/larm; % conductor force after removing load head inertia
xmeas = [datatab.displ x_est_ric(:,2) datatab.F xv datatab.State1 datatab.State2]; 

xmeasdot = [x_est_ric(:,2) x_est_ric(:,2) Fdot_est_ric xvdot S1dot_est_ric S2dot_est_ric];

usys_ = [wbar datatab.uext];
rhs = xmeas*A' + usys_*B';

figure(21), plot(t,xmeasdot(:,1),'r',t, rhs(:,1),'b--'), grid on, title('equation 1')
figure(22), plot(t, xmeasdot(:,2)/386.4, 'r', t, (rhs(:,2)-mean(rhs(:,2)))/386.4,'b--'), grid on, title('equation 2')
figure(29), plot([dddash/386.4 (rhs(:,2)-mean(rhs(:,2)))/386.4]), grid on, title('equation 2')
figure(23), plot(t, xmeasdot(:,3), 'r',t, rhs(:,3),'b--'), grid on, title('equation 3')
figure(24), plot(t, xmeasdot(:,4), 'r', t, rhs(:,4), 'b--'), grid on, title('equation 4')
figure(25), plot([xmeasdot(:,5) rhs(:,5)]), grid on, title('equation 5')
figure(26), plot([xmeasdot(:,6) rhs(:,6)]), grid on, title('equation 6')


%% fit parameters
% Afit = [x_est_ric(:,2) datatab.F xv]\(rhs-mean(rhs)); %% Fitting Fdot = -bet*F -kact*xdot + d*xv
% rhscalc = [x_est_ric(:,2) datatab.F xv]*Afit;

Afit = [x_est_ric(:,2) xv]\(rhs3-mean(rhs3)+bet*datatab.F); %%% Fitting Fdot + bet*F = -kact*xdot + d*xv
rhscalc = [x_est_ric(:,2) xv]*Afit;

figure(303), plot([Fdot_est_ric rhscalc-mean(rhscalc)-bet*datatab.F]), grid on,
AA = A;
AA(3,2:4) = [Afit(1) -bet Afit(2)];
fprintf('Frequencies and damping ratios of fit model ...\n')
damp(ss(AA,B,eye(1,6),[]))
[kact bet d]


%% model with fit parameters
ssActCont = ss(AA, B, C, D);
ssActCont.InputName = {'w','uext'};
ssActCont.OutputName = {'x','v','F','xv','xC1','xC2','a'};

S1 = sumblk('e = xVS - x');

ssAVS1=connect(ssVS, ssActCont,S1,{'w','uext'},{'x','v','F','xv','xC1','xC2','e','a'});

%% VeriStand adds arbitrary number of zeros before the actual uext.. 
%%% line up earthquake data to VS to match with uext to AVS or measured
load cerl1gyu.mat %%% this is already 1g
dir = "X";
tow = "1";
tcerl = (0:length(cerl)-1)/512;
fs = 1000;

ttt = (0:1/fs:tcerl(end));

if dir == "X"
    uu = interp1(tcerl,cerl(:,1)-cerl(1,1),ttt);
    wbarbar = ((datatab.Fxcalc+mLHE*datatab.accel)*h+(datatab.Mycalc+ILH*thetaddot))/larm; 
else
    uu = interp1(tcerl,cerl(:,2)-cerl(2,1),ttt);
    wbarbar = ((datatab.Fycalc+mLHE*datatab.accel)*h-(datatab.Mxcalc+ILH*thetaddot))/larm;
end

tshift = 15.225-10.025; %%% eq 0.75g, TEWVS12

load TEXEQControllerVS20020.mat uext %%% external valve command computed with EQ controller
uext1 = 0.75*uext;

tt = (0:length(uext1)-1)/fs;
figure, plot(tt+tshift,uext1,'r',t,datatab.uext,'b--'), grid on, %% check time shift is correct

w12 = [zeros(round(tshift*fs),1); uu'];
w2 = [zeros(round(tshift*fs),1); uu'; zeros(length(datatab.uext)-length(w12),1)];
ttt = (0:length(w2)-1)/fs;

%%
weq = w2(1:length(wbarbar))*0.75*386.4;
yVS = lsim(ssVS, [wbarbar weq], t);

ywVS = lsim(ssVS, [zeros(size(wbarbar)) weq], t);
yAVS = lsim(ssAVS, [wbarbar datatab.uext], t);
yAVS1 = lsim(ssAVS1, [wbarbar datatab.uext], t);
%%%%%%%%%% Displacement
%%% compare VS and AVS (no fit) displacement
figure(404), %% very good comparison, showing that EQ controller is okay
plot(ttt,yVS(:,1),'r',t,yAVS(:,1),'b--'), grid on,
xlabel('Time (s)')
ylabel('Displacement (in)')
legend('Target VS','AVS - no fit','EdgeColor','None','Location','northwest','Color','None')
pbaspect([4 1 1])
set(gcf, 'Position',[100 100 950 250])
xlim([5 46])
title("No Fit")
%%% compare VS and AVS (fit) displacement
figure(405), %% not very good comparison, because d value has changed significantly
plot(ttt,yVS(:,1),'r',t,yAVS1(:,1),'b--'), grid on,
xlabel('Time (s)')
ylabel('Displacement (in)')
legend('Target VS','Achievable VS','EdgeColor','None','Location','northwest','Color','None')
pbaspect([4 1 1])
set(gcf, 'Position',[100 100 950 250])
xlim([5 46])
title("Fit")
%%% compare AVS (no fit) and AVS (fit) displacement
figure(406), %% parameter change do make a difference at the start of the signal
plot(ttt,yAVS(:,1),'r',t,yAVS1(:,1),'b--'), grid on,
xlabel('Time (s)')
ylabel('Displacement (in)')
legend('No Fit','Fit','EdgeColor','None','Location','northwest','Color','None')
pbaspect([4 1 1])
set(gcf, 'Position',[100 100 950 250])
xlim([5 46])
title("AVS")

%%%%%%%%%% Acceleration
%%% compare VS and AVS (no fit) acceleration
figure(505),% very good comparison, showing that EQ controller is okay
plot(ttt,yVS(:,3)/386.4,'r',t,yAVS(:,8)/386.4,'b--'), grid on,
xlabel('Time (s)')
ylabel('Acceleration (g)')
legend('Target VS','Achievable VS','EdgeColor','None','Location','northwest','Color','None')
pbaspect([4 1 1])
set(gcf, 'Position',[100 100 950 250])
xlim([5 46])
%%% compare VS and AVS (fit) acceleration
figure(506), %%% like displacement, comparison not that good
plot(ttt,yVS(:,3)/386.4,'r',t,yAVS1(:,8)/386.4,'b--'), grid on,
xlabel('Time (s)')
ylabel('Acceleration (g)')
legend('Target VS','Achievable VS','EdgeColor','None','Location','northwest','Color','None')
pbaspect([4 1 1])
set(gcf, 'Position',[100 100 950 250])
xlim([5 46])
%%% compare AVS (no fit) and AVS (fit) acceleration
figure(507), %% difference in the start of the signal
plot(ttt,yAVS(:,8)/386.4,'r',t,yAVS1(:,8)/386.4,'b--'), grid on,
xlabel('Time (s)')
ylabel('Acceleration (g)')
legend('No Fit','Fit','EdgeColor','None','Location','northwest','Color','None')
pbaspect([4 1 1])
set(gcf, 'Position',[100 100 950 250])
xlim([5 46])
title("AVS")

%%%%%% Displacement with measurement
%%% comparison of measured with AVS displacement
%%% Not much difference between the fit and no fit displacement
figure(608),% not fit AVS with measurement
plot(ttt,yAVS(:,1)-mean(yAVS(:,1)),'r',t,datatab.displ-datatab.displ(1),'k--'), grid on,
xlabel('Time (s)')
ylabel('Displacement (in)')
legend('AVS - no fit','Measured','EdgeColor','None','Location','northwest','Color','None')
pbaspect([4 1 1])
set(gcf, 'Position',[100 100 950 250])
xlim([5 46])

figure(609),% fit AVS with measurement
plot(ttt,yAVS1(:,1)-mean(yAVS1(:,1)),'r',t,datatab.displ-datatab.displ(1),'k--'), grid on,
xlabel('Time (s)')
ylabel('Displacement (in)')
legend('AVS - fit','Measured','EdgeColor','None','Location','northwest','Color','None')
pbaspect([4 1 1])
set(gcf, 'Position',[100 100 950 250])
xlim([5 46])
%%%%%% Acceleration with measurement
%%% comparison of measured with AVS acceleration
%%% Fit parameters AVS compare well with measurement, espeically at the
%%% start and end signals.. these are the plots I was talking about
figure(605),% not fit AVS with measurement
plot(t,yAVS(:,8)/386.4,'r',t,(datatab.accel-mean(datatab.accel))*larm/h,'b--'), grid on,
xlabel('Time (s)')
ylabel('Acceleration (g)')
legend('AVS - no fit','Measured','EdgeColor','None','Location','northwest','Color','None')
pbaspect([4 1 1])
set(gcf, 'Position',[100 100 950 250])
xlim([5 46])

figure(607), % fit AVS with measurement
plot(t,yAVS1(:,8)/386.4,'r',t,(datatab.accel-mean(datatab.accel))*larm/h,'b--'), grid on,
xlabel('Time (s)')
ylabel('Acceleration (g)')
legend('AVS - fit','Measured','EdgeColor','None','Location','northwest','Color','None')
pbaspect([4 1 1])
set(gcf, 'Position',[100 100 950 250])
xlim([5 46])

%%
rmpath ../vel-est-derivative-comp/
rmpath ../../data/
