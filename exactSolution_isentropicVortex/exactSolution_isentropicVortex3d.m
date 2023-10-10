function [q, qx, qy, qz, qt] = exactSolution_isentropicVortex3d(x, y, z, x0, y0, z0, time)
    % Parameters
    gamma = 1.4;
     
    % Free stream values
    r_inf = 1;
    u_inf = 1;
    v_inf = 1;
    w_inf = 1;
    p_inf = 1 / gamma;
    a_inf = sqrt(gamma * p_inf / r_inf);

    % Some parameters
    Vmax = 5.0;
    alpha = 1.0;

    % Coordinates
    xp = x - x0 - u_inf * time;
    yp = y - y0 - v_inf * time;
    zp = z - z0 - w_inf * time;

    % Delta V
    [a,b,c] = deal(1,0,0);
    rp = [b*zp - c*yp, c*xp - a*zp, a*yp - b*xp]; % [a,b,c] X [xp,yp,zp];
    rp2 = (b*zp - c*yp).^2 + (c*xp - a*zp).^2 + (a*yp - b*xp).^2;
    dV = Vmax/(2 * pi * alpha) .* exp(0.5 * alpha * (1-rp2));

    % Compute the exact solution:
    V = ([u_inf, v_inf, w_inf] + dV .* rp) / a_inf;
    [u,v,w] = deal(V(:,1),V(:,2),V(:,3));
    T = 1 - 0.5*(gamma - 1) * (dV / a_inf).^2;
    rho = T.^(1 / (gamma - 1)) * r_inf;
    p = T.^(gamma / (gamma - 1)) * p_inf / (r_inf * a_inf^2);
    
    % Compute derivatives
    rp2x = -2*b*(a*yp - b*xp) - 2*c*(a*zp - c*xp);
    rp2y =  2*a*(a*yp - b*xp) - 2*c*(b*zp - c*yp);
    rp2z =  2*a*(a*zp - c*xp) + 2*b*(b*zp - c*yp);
    ux = dV/a_inf .* (rp(:,1) .* rp2x);
    uy = dV/a_inf .* (rp(:,1) .* rp2y - c);
    uz = dV/a_inf .* (rp(:,1) .* rp2z + b);
    vx = dV/a_inf .* (rp(:,2) .* rp2x + c);
    vy = dV/a_inf .* (rp(:,2) .* rp2y);
    vz = dV/a_inf .* (rp(:,2) .* rp2z - a);
    wx = dV/a_inf .* (rp(:,3) .* rp2x - b);
    wy = dV/a_inf .* (rp(:,3) .* rp2y + a);
    wz = dV/a_inf .* (rp(:,3) .* rp2z);
    Tx = 0.5 * (gamma - 1) * (dV/a_inf).^2 .* rp2x;
    Ty = 0.5 * (gamma - 1) * (dV/a_inf).^2 .* rp2y;
    Tz = 0.5 * (gamma - 1) * (dV/a_inf).^2 .* rp2z;
    rhox = (1 / (gamma - 1)) * T.^(1 / (gamma - 1) - 1) .* Tx * r_inf;
    rhoy = (1 / (gamma - 1)) * T.^(1 / (gamma - 1) - 1) .* Ty * r_inf;
    rhoz = (1 / (gamma - 1)) * T.^(1 / (gamma - 1) - 1) .* Tz * r_inf;
    px = (gamma / (gamma - 1)) * T.^(gamma / (gamma - 1) - 1) .* Tx * p_inf / (r_inf * a_inf^2);
    py = (gamma / (gamma - 1)) * T.^(gamma / (gamma - 1) - 1) .* Ty * p_inf / (r_inf * a_inf^2);
    pz = (gamma / (gamma - 1)) * T.^(gamma / (gamma - 1) - 1) .* Tz * p_inf / (r_inf * a_inf^2);
    
    % Compute time derivative
    ut = ux .* (-u_inf) + uy .* (-v_inf) + uz .* (-w_inf);
    vt = vx .* (-u_inf) + vy .* (-v_inf) + vz .* (-w_inf);
    wt = wx .* (-u_inf) + wy .* (-v_inf) + wz .* (-w_inf);
    Tt = Tx .* (-u_inf) + Ty .* (-v_inf) + Tz .* (-w_inf);
    rhot = (1 / (gamma - 1)) .* T.^(1 / (gamma - 1) - 1) .* Tt;
    pt = (gamma / (gamma - 1)) .* T.^(gamma / (gamma - 1) - 1) .* Tt * p_inf / (r_inf * a_inf^2);
    
    % Store results in output arrays
    q = [rho, u, v, w, p];
    qx = [rhox, ux, vx, wx, px];
    qy = [rhoy, uy, vy, wy, py];
    qz = [rhoz, uz, vz, wz, pz];
    qt = [rhot, ut, vt, wt, pt];
end
