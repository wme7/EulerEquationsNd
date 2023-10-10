function [x, y, psi, V, th] = mesh_ringlebFlow(psi_min, psi_max, V_min, nx, ny)
% Generate boundary points

% Build uniform Psi-coordinate
dpsi = (psi_max-psi_min)./nx;
psi = psi_min + dpsi * (0:nx);
psi = repmat(transpose(psi),[1,ny+1]);

% Build uniform Theta-coordinate
th_min = asin(psi(:,1) .* V_min);
th_max = pi/2;
dth = (th_max-th_min)./ny;
th = th_min + dth * (0:ny);

% V-coordiante
V = sin(th)./ psi;

% Geometric variables
b = sqrt(1-0.2*V.^2);
L = 1./b + 1./(3*b.^3) + 1./(5.*b.^5) - 1/2*log((1+b)./(1-b));
rho = b.^5;

% Physical coordinates
x = (1./rho).*(0.5./V./V - psi.*psi) + 0.5*L;
y = psi./(rho.*V) .*sqrt(1-(V.*psi.*V.*psi));

end % function