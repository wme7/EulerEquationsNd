function [rho, vx, vy, p] = exactSolution_ringlebFlow(psi, V)

% Compute exact solution
mask = abs(psi.*V-1) < 1.0e-15;
theta = asin( psi .* V );
theta(mask) = 0.5*pi; % do not rely on `asin` here

b = sqrt(1-0.2*V.*V);
rho = b.^5;
vx = -V.*cos(theta);
vy = -V.*sin(theta);
p = b.^7;

end % function