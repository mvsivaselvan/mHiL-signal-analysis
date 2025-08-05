%% Driver code to test
% 1. estimate_derivative_differential_riccati
% 2. estimate_derivative_qp
% 3. estimate_velocity_differential_riccati
% 4. estimate_velocity_qp

%% test estimate derivative
t = (0:0.001:1)';
x0 = 5;
f = 1;
x_m = x0*sin(2*pi*f*t); % signal to differentiate
u_exact = 2*pi*f*x0*cos(2*pi*f*t); % exact derivative
Q = 1e7;
% [x_est_ric, u_est_ric] = estimate_derivative_differential_riccati(t, x_m, Q);
[x_est_qp, u_est_qp] = estimate_derivative_qp(t, x_m, Q);
[x_est_ric, u_est_ric] = estimate_derivative_differential_riccati_analytic(t, x_m, 1/Q);
figure(101),
    plot(t, [u_exact u_est_qp u_est_ric]),
    grid on
    legend('exact', 'qp', 'ric'),
    title('derivative')
figure(102),
    plot(t, [x_m x_est_qp x_est_ric]),
    grid on
    legend('exact', 'qp', 'ric')
    title('signal')

%% test estimate velocity
t = (0:0.001:1)';
x0 = 5;
f = 2;
x_m = x0*sin(2*pi*f*t); % measured displacement
u_m = -(2*pi*f)^2*x0*sin(2*pi*f*t); % measured acceleration
v_exact = (2*pi*f)*x0*cos(2*pi*f*t); % exact velocity
Q = 0.001;
[x_est_qp, u_est_qp] = estimate_velocity_qp(t, x_m, u_m, Q);
[x_est_ric, u_est_ric] = estimate_velocity_differential_riccati(t, x_m, u_m, Q);
figure(201),
    plot(t, [x_m, x_est_qp(:,1), x_est_ric(:,1)]),
    grid on
    legend('measured', 'qp', 'ric'),
    title('displacement')
figure(202),
    plot(t, [v_exact, x_est_qp(:,2), x_est_ric(:,2)]),
    grid on
    legend('exact', 'qp', 'ric')
    title('velocity')
figure(203),
    plot(t, [u_m, u_est_qp, u_est_ric]),
    grid on
    legend('measured', 'qp', 'ric')
    title('acceleration')