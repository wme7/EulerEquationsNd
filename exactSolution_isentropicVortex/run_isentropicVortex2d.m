% Test Exact solutions
addpath('../common/');

% Aceptable numerical error to stablish solution validity:
TOL = 1E-12;

% 2D - mesh
[x,y] = meshgrid(linspace(0,20,50));

% Initialize figure
figure(3)

% Vortex initial position
[x0,y0] = deal(5,5);

% Compute the exact solution and derivatives
t0=0; dt=0.5; tEnd=10;
for t = t0:dt:tEnd

    % Compute exact solution
    [q,qx,qy,qt] = exactSolution_isentropicVortex2d(x(:),y(:),x0,y0,t);
    
    % Stop if solution does not meet standard
    if not(evalEulerEquations2d(q,qx,qy,qt,TOL)), return; end

    % Update figure
    plotContoursEulerEquations2d(x,y,q,t);
    drawnow;
end
