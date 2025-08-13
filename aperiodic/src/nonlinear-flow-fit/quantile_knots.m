function r = quantile_knots(x, M)
    % x: vector of sample points (e.g., x4 data)
    % M: number of knots
    % Output: r: knot positions

    % Ensure x is column vector
    x = x(:);

    % Compute quantile positions from 0 to 1
    q = linspace(0, 1, M);

    % Use MATLAB's quantile function to get knot positions
    r = quantile(x, q);

    % Ensure knots are strictly increasing (remove duplicates if any)
    r = unique(r, 'stable');

    % If duplicates occur due to repeated values in x, pad with small eps
    if numel(r) < M
        warning('Some knots collapsed due to repeated values in data.');
        r = r + (0:numel(r)-1)' * eps;
    end
end
