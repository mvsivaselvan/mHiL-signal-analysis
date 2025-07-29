
load ControllerTEXVS20020.mat % to get kact, beta, d
load datatab.mat


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
%% target VS
AVS = [0 1; -w^2 -2*zet*w];
BVS = [0 0 ;1/mVS -mbar*hCM/mVS/larm]; %% (mbar*hCM/mVS/larm)
CVS = [1 0; 0 1; -w^2 -2*zet*w];
DVS = [0 0; 0 0; 1/mVS -mbar*hCM/mVS/larm];
ssVS = ss(AVS, BVS, CVS, DVS);
ssVS.InputName = {'w','weq'};
ssVS.OutputName = {'xVS', 'vVS','aVS'};
%% AVS
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

%% equation by equation check
x_m = datatab.displ;
u_m = datatab.accel*larm/h*386.4;
Q = 0.001;
[x_est_qp, a_est_qp] = estimate_velocity_qp(t, x_m, u_m, Q); %%% using qp
[x_est_ric, a_est_ric] = estimate_velocity_differential_riccati(t, x_m, u_m, Q); %%% using riccati

% figure(201),
%     plot(t, [x_m, x_est_qp(:,1), x_est_ric(:,1)]),
%     grid on
%     legend('measured', 'qp', 'ric'),
%     title('displacement')
% figure(202),
%     plot(t, [x_est_qp(:,2), x_est_ric(:,2)]),
%     grid on
%     legend('qp', 'ric')
%     title('velocity')
% figure(203),
%     plot(t, [u_m, a_est_qp, a_est_ric]),
%     grid on
%     legend('measured', 'qp', 'ric')
%     title('acceleration')
figure(201),
    plot(t, [x_m, x_est_ric(:,1)]),
    grid on
    legend('measured', 'ric'),
    title('displacement')
figure(202),
    plot(t, x_est_ric(:,2)),
    grid on
    legend('qp', 'ric')
    title('velocity')
figure(203),
    plot(t, [u_m, a_est_ric]),
    grid on
    legend('measured', 'qp', 'ric')
    title('acceleration')

Q = 1e7;
% [x_est_ric, u_est_ric] = estimate_derivative_differential_riccati(t, x_m, Q);
x_m = F;
[F_est_qp, Fdot_est_qp] = estimate_derivative_qp(t, x_m, Q);
[F_est_ric, Fdot_est_ric] = estimate_derivative_differential_riccati_analytic(t, x_m, 1/Q);
figure(101),
    plot(t, [u_exact Fdot_est_qp Fdot_est_ric]),
    grid on
    legend('exact', 'qp', 'ric'),
    title('derivative')
figure(102),
    plot(t, [x_m F_est_qp F_est_ric]),
    grid on
    legend('exact', 'qp', 'ric')
    title('signal')

[S1_est_qp, S1dot_est_qp] = estimate_derivative_qp(t, State1, Q);
[S1_est_ric, S1dot_est_ric] = estimate_derivative_differential_riccati_analytic(t, State1, 1/Q);
figure(401),
    plot(t, [S1dot_est_qp S1dot_est_ric]),
    grid on
    legend('qp', 'ric'),
    title('derivative')
figure(402),
    plot(t, [State1 S1_est_qp S1_est_ric]),
    grid on
    legend('exact', 'qp', 'ric')
    title('signal')

[S2_est_qp, S2dot_est_qp] = estimate_derivative_qp(t, State2, Q);
[S2_est_ric, S2dot_est_ric] = estimate_derivative_differential_riccati_analytic(t, State2, 1/Q);
figure(501),
    plot(t, [S2dot_est_qp S2dot_est_ric]),
    grid on
    legend('qp', 'ric'),
    title('derivative')
figure(502),
    plot(t, [State2 S2_est_qp S2_est_ric]),
    grid on
    legend('exact', 'qp', 'ric')
    title('signal')

[S3_est_qp, S1dot_est_qp] = estimate_derivative_qp(t, State1, Q);
[S3_est_ric, S1dot_est_ric] = estimate_derivative_differential_riccati_analytic(t, State1, 1/Q);
figure(401),
    plot(t, [S3dot_est_qp S3dot_est_ric]),
    grid on
    legend('qp', 'ric'),
    title('derivative')
figure(402),
    plot(t, [State3 S3_est_qp S3_est_ric]),
    grid on
    legend('exact', 'qp', 'ric')
    title('signal')

sys = tf(alph, [1 alph]);
xv = lsim(sys, u, t);
thetaddot = a/h;
xvdot = -alph*(xv-u);
wbar = ((Fcalc+mLHE*a)*h+(Mcalc+ILH*thetaddot))/larm; % conductor force after removing load head inertia
xmeas = [d v F xv State1 State2 State3]; % bvpsol.y(:,1:N)';

xmeasdot = [x_est_qp(:,2) x_est_qp(:,2) Fdot_est_qp xvdot ];
usys_ = [uextf wbarf];
rhs = xmeas*AA' + usys_*BB';

figure(21), plot([xmeasdot(:,1) rhs(:,1)]), grid on, title('equation 1')
figure(22), plot([xmeasdot(:,2)/386.4 (rhs(:,2)-mean(rhs(:,2)))/386.4]), grid on, title('equation 2')
figure(29), plot([dddash/386.4 (rhs(:,2)-mean(rhs(:,2)))/386.4]), grid on, title('equation 2')

figure(23), plot([xmeasdot(:,3) rhs(:,3)]), grid on, title('equation 3')
figure(24), plot([xmeasdot(:,4) rhs(:,4)]), grid on, title('equation 4')
figure(25), plot([xmeasdot(:,5) rhs(:,5)]), grid on, title('equation 5')
figure(26), plot([xmeasdot(:,6) rhs(:,6)]), grid on, title('equation 6')
figure(27), plot([xmeasdot(:,7) rhs(:,7)]), grid on, title('equation 7')

