% compare CasADi's collocation_points function to spherepack's gaqdm that
% I have used a lot, to understand

Ng = 3; % number of collocation points

% Casadi
xgcas = casadi.collocation_points(Ng, 'legendre');
[C,D,B] = casadi.collocation_coeff(xgcas);

% spherepack
[xg,wg] = gaqdm(Ng); % Gauss points \in (0,pi) and weights
xg = flipud(cos(xg)); % Gauss points \in (-1,1)
xg = (1+xg)/2; % Gauss points \in (0,1)
wg = wg*(1/2); % weights multiplied by appropriate jacobian for integration
               % over (0,1)

fprintf('Quadrature points ...\n')
disp([xgcas' xg]);
fprintf('Quadrature weights ...\n')
disp([B wg])
