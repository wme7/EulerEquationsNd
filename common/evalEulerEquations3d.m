function flag = evalEulerEquations3d(q, qx, qy, qz, qt, TOL)
% Parameters
gamma = 1.4;

% Solutions and derivatives
r = q(:,1);
u = q(:,2);
v = q(:,3);
w = q(:,4);
p = q(:,5);

rx = qx(:,1);
ux = qx(:,2);
vx = qx(:,3);
wx = qx(:,4);
px = qx(:,5);

ry = qy(:,1);
uy = qy(:,2);
vy = qy(:,3);
wy = qy(:,4);
py = qy(:,5);

rz = qz(:,1);
uz = qz(:,2);
vz = qz(:,3);
wz = qz(:,4);
pz = qz(:,5);

rt = qt(:,1);
ut = qt(:,2);
vt = qt(:,3);
wt = qt(:,4);
pt = qt(:,5);

% Compute some quantities needed to evaluate the 2D unsteady Euler system
rhout = rt .* u + r .* ut;
rhovt = rt .* v + r .* vt;
rhowt = rt .* w + r .* wt;

rhoux = rx .* u + r .* ux;
rhovy = ry .* v + r .* vy;
rhowz = rz .* w + r .* wz;

rhouux = rx .* u .* u + r .* ux .* u + r .* u .* ux;
rhouvx = rx .* v .* u + r .* vx .* u + r .* v .* ux;
rhouwx = rx .* w .* u + r .* wx .* u + r .* w .* ux;

rhouvy = ry .* u .* v + r .* uy .* v + r .* u .* vy;
rhovvy = ry .* v .* v + r .* vy .* v + r .* v .* vy;
rhowvy = ry .* w .* v + r .* wy .* v + r .* w .* vy;

rhouwz = rz .* u .* w + r .* uz .* w + r .* u .* wz;
rhovwz = rz .* v .* w + r .* vz .* w + r .* v .* wz;
rhowwz = rz .* w .* w + r .* wz .* w + r .* w .* wz;

rH  = gamma * p  ./ (gamma - 1) + 0.5 *   r .* (u .* u + v .* v + w .* w);
rHx = gamma * px ./ (gamma - 1) + 0.5 * (rx .* (u .* u + v .* v + w .* w) + 2 * r .* (u .* ux + v .* vx + w .* wx));
rHy = gamma * py ./ (gamma - 1) + 0.5 * (ry .* (u .* u + v .* v + w .* w) + 2 * r .* (u .* uy + v .* vy + w .* wy));
rHz = gamma * pz ./ (gamma - 1) + 0.5 * (rz .* (u .* u + v .* v + w .* w) + 2 * r .* (u .* uz + v .* vz + w .* wz));

rhouHx = ux .* rH + u .* rHx;
rhovHy = vy .* rH + v .* rHy;
rhowHz = wz .* rH + w .* rHz;

rEt = pt ./ (gamma - 1) + 0.5 * (rt .* (u .* u + v .* v + w .* w) + 2 * r .* (u .* ut + v .* vt + w .* wt));

% Form the 2D unsteady Euler equations
equation(1) = sum(    rt +  rhoux +  rhovy + rhowz );
equation(2) = sum( rhout + rhouux + rhouvy + rhouwz + px );
equation(3) = sum( rhovt + rhouvx + rhovvy + rhovwz + py );
equation(4) = sum( rhowt + rhouwx + rhowvy + rhowwz + pz );
equation(5) = sum(   rEt + rhouHx + rhovHy + rhowHz );

% Display the results. All must be zero.
fprintf('\nSubstitution yields:\n\n');
fprintf('Continuity = %1.12f\n', equation(1));
fprintf('X-momentum = %1.12f\n', equation(2));
fprintf('Y-momentum = %1.12f\n', equation(3));
fprintf('Z-momentum = %1.12f\n', equation(4));
fprintf('Energy     = %1.12f\n', equation(5));

% Output
flag = all(equation < TOL);
if (flag)
    fprintf('\n Is an exact solution!\n\n'); 
else
    fprintf('\n Not an exact solution!\n\n'); 
end

end % function