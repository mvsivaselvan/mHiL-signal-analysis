import casadi.*

fs_meas = 10;
N = 100;
t_meas = (0:N-1)' / fs_meas;   % [100 × 1]

ny = 2;
ym_data = [sin(t_meas)'; cos(t_meas)'];  % [2 × 100]

% Ensure dimensions
disp(size(t_meas))    % Should be [100 1]
disp(size(ym_data))   % Should be [2 100]

ym_fun = interpolant('ym','linear',{t_meas}, ym_data(:));