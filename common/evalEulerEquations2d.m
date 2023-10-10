function flag = evalEulerEquations2d(q, qx, qy, qt, TOL)
% Parameters
gamma = 1.4;

% Solutions and derivatives
r = q(:,1);
u = q(:,2);
v = q(:,3);
p = q(:,4);

rx = qx(:,1);
ux = qx(:,2);
vx = qx(:,3);
px = qx(:,4);

ry = qy(:,1);
uy = qy(:,2);
vy = qy(:,3);
py = qy(:,4);

rt = qt(:,1);
ut = qt(:,2);
vt = qt(:,3);
pt = qt(:,4);

% Compute some quantities needed to evaluate the 2D unsteady Euler system
rhout = rt .* u + r .* ut;
rhovt = rt .* v + r .* vt;

rhoux = rx .* u + r .* ux;
rhovy = ry .* v + r .* vy;

rhouux = rx .* u .* u + r .* ux .* u + r .* u .* ux;
rhouvx = rx .* u .* v + r .* ux .* v + r .* u .* vx;

rhouvy = ry .* u .* v + r .* uy .* v + r .* u .* vy;
rhovvy = ry .* v .* v + r .* vy .* v + r .* v .* vy;

rH  = gamma * p  ./ (gamma - 1) + 0.5 *  r .* (u .* u + v .* v);
rHx = gamma * px ./ (gamma - 1) + 0.5 * (rx .* (u .* u + v .* v) + 2 * r .* (u .* ux + v .* vx));
rHy = gamma * py ./ (gamma - 1) + 0.5 * (ry .* (u .* u + v .* v) + 2 * r .* (u .* uy + v .* vy));

rhouHx = ux .* rH + u .* rHx;
rhovHy = vy .* rH + v .* rHy;

rEt = pt ./ (gamma - 1) + 0.5 * (rt .* (u .* u + v .* v) + 2 * r .* (u .* ut + v .* vt));

% Form the 2D unsteady Euler equations
equation(1) = sum(    rt +  rhoux + rhovy );
equation(2) = sum( rhout + rhouux + rhouvy + px );
equation(3) = sum( rhovt + rhouvx + rhovvy + py );
equation(4) = sum(   rEt + rhouHx + rhovHy );

% Display the results. All must be zero.
fprintf('\nSubstitution yields:\n\n');
fprintf('Continuity = %1.12f\n', equation(1));
fprintf('X-momentum = %1.12f\n', equation(2));
fprintf('Y-momentum = %1.12f\n', equation(3));
fprintf('Energy     = %1.12f\n', equation(4));

% Output
flag = all(equation < TOL);
if (flag)
    fprintf('\n Is an exact solution!\n\n'); 
else
    fprintf('\n Not an exact solution!\n\n'); 
end

end % function