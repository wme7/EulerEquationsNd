% Test Exact solutions
addpath('../common/');

% Aceptable numerical error to stablish solution validity:
TOL = 1E-15;

% This defines the domain.
  psi_min = 0.69; % suggested value = 0.69_dp
  psi_max = 1.20; % suggested value = 1.20_dp
    V_min = 0.50; % suggested value = 0.50_dp
% This defines the grid size.
       nx = 30;   % the number of nodes in psi-direction
       ny = 60;   % the number of nodes in V-direction

% Build mesh
[x,y,psi,V,theta] = mesh_ringlebFlow(psi_min, psi_max, V_min, nx, ny);

% Compute mesh and Ringleb Flow solution
[rho,u,v,p] = exactSolution_ringlebFlow(psi, V);

% Initialize figure
figure(6);
subplot(221); contourf(x,y,rho); title(sprintf('$\\rho(x,y)$'), Interpreter='latex');
subplot(222); contourf(x,y,u); title(sprintf('$u(x,y)$'), Interpreter='latex');
subplot(223); contourf(x,y,v); title(sprintf('$v(x,y)$'), Interpreter='latex');
subplot(224); contourf(x,y,p); title(sprintf('$p(x,y)$'), Interpreter='latex');
