function Phi = basis_eval_symbolic(x, r)
    import casadi.*
    M = numel(r);
    Phi = MX.zeros(M,1);

    % Determine interval index idx, 0 if x < r(1), M if x >= r(M)
    ge = MX.zeros(M,1);
    for j = 1:M
        ge(j) = if_else(x >= r(j), 1, 0);
    end
    idx = sum(ge);

    Phi_left = MX.zeros(M,1);
    h_left = r(2) - r(1);
    t_left = (x - r(1)) / h_left;
    Phi_left(1) = 1 - t_left;
    Phi_left(2) = t_left;

    Phi_right = MX.zeros(M,1);
    h_right = r(M) - r(M-1);
    t_right = (x - r(M)) / h_right;
    Phi_right(M-1) = -t_right;
    Phi_right(M) = 1 + t_right;

    % Interior - compute idx symbolically
    Phi_interior = MX.zeros(M,1);
    for j = 1:(M-1)
        h = r(j+1) - r(j);
        t = (x - r(j)) / h;
        cond = (idx == j);
        Phi_interior(j)   = Phi_interior(j)   + if_else(cond, 1 - t, 0);
        Phi_interior(j+1) = Phi_interior(j+1) + if_else(cond, t, 0);
    end

    % Now combine all three cases with if_else:
    Phi = if_else(idx == 0, Phi_left, ...
        if_else(idx == M, Phi_right, ...
        Phi_interior));

end
