% Test Exact solutions
addpath('../common/');

% Aceptable numerical error to stablish solution validity:
TOL = 9E-11;

% 2D - mesh
[x,y,z] = meshgrid(linspace(0,10,50));

% Initialize figure
figure(4)

% Initial vortex center
[x0,y0,z0] = deal(2.5,2.5,5);

% Compute the exact solution and derivatives
t0=0; dt=0.5; tEnd=5;
for t = t0:dt:tEnd

    % Compute exact solution
    [q,qx,qy,qz,qt] = exactSolution_isentropicVortex3d(x(:),y(:),z(:),x0,y0,z0,t);
    
    % Stop if solution does not meet standard
    if not(evalEulerEquations3d(q,qx,qy,qz,qt,TOL)), return; end

    % Update figure
    plotEulerEquations3d(x,y,z,q,t);
    drawnow;
end

% Iso-surface at final time
figure(5)
plotContoursEulerEquations3d(x,y,z,q,t);
