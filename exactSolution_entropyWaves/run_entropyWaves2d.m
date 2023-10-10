% Test Exact solutions
addpath('../common/');

% Aceptable numerical error to stablish solution validity:
TOL = 1E-12;

% 2D - mesh
[x,y] = meshgrid(linspace(0,1,25));

% Initialize figure
figure(1)

% Compute the exact solution and derivatives
t0=0; dt=0.1; tEnd=1;
for t = t0:dt:tEnd
    
    % Compute exact solution
    [q,qx,qy,qt] = exactSolution2D_entropyWaves(x(:),y(:),t);
    
    % Stop if solution does not meet standard
    if not(evalEulerEquations2d(q,qx,qy,qt,TOL)), return; end

    % Update figure
    plotEulerEquations2d(x,y,q,t);
    drawnow;
end
