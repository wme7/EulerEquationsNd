% Test Exact solutions
addpath('../common/');

% Aceptable numerical error to stablish solution validity:
TOL = 9E-11;

% 3D -mesh
[x,y,z] = meshgrid(linspace(0,1,25));

% Initialize figure
figure(2)

% Compute the exact solution and derivatives
t0=0; dt=0.1; tEnd=1;
for t = t0:dt:tEnd

    % Compute exact solution
    [q,qx,qy,qz,qt] = exactSolution3D_entropyWaves(x(:),y(:),z(:),t);

    % Stop if solution does not meet standard
    if not(evalEulerEquations3d(q,qx,qy,qz,qt,TOL)), return; end

    % Update figure
    plotEulerEquations3d(x,y,z,q,t);
    drawnow;
end
