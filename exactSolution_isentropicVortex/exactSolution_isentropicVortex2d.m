function [w, wx, wy, wt] = exactSolution_isentropicVortex2d(x, y, x0, y0, time)
    % Parameters
    gamma = 1.4;
     
    % Free stream values
    r_inf = 1;
    u_inf = 1;
    v_inf = 1;
    p_inf = 1 / gamma;
    a_inf = sqrt(gamma * p_inf / r_inf);

    % Some parameters
    Vmax = 5.0;
    alpha = 1.0;

    % Coordinates
    xp = x - x0 - u_inf * time;
    yp = y - y0 - v_inf * time;
    zp = zeros(size(x));

    % Delta V
    [a,b,c] = deal(0,0,1);
    rp = [b*zp - c*yp, c*xp - a*zp, a*yp - b*xp]; % [a,b,c] X [xp,yp,zp];
    rp2 = (b*zp - c*yp).^2 + (a*zp - c*xp).^2 + (a*yp - b*xp).^2;
    dV = Vmax/(2 * pi * alpha) .* exp(0.5 * alpha * (1-rp2));
    
    % Compute the exact solution
    V = ([u_inf, v_inf, 0] + dV .* rp) / a_inf;
    [u,v] = deal(V(:,1),V(:,2));
    T = 1 - 0.5*(gamma - 1) * (dV / a_inf).^2;
    rho = T.^(1 / (gamma - 1));
    p = T.^(gamma / (gamma - 1)) * p_inf / (r_inf * a_inf^2);
    
    % Compute derivatives
    rp2x = -2*b*(a*yp - b*xp) - 2*c*(a*zp - c*xp);
    rp2y =  2*a*(a*yp - b*xp) - 2*c*(b*zp - c*yp);
    ux = dV/a_inf .* (rp(:,1) .* rp2x);
    uy = dV/a_inf .* (rp(:,1) .* rp2y - c);
    vx = dV/a_inf .* (rp(:,2) .* rp2x + c);
    vy = dV/a_inf .* (rp(:,2) .* rp2y);
    Tx = 0.5 * (gamma - 1) * (dV/a_inf).^2 .* rp2x;
    Ty = 0.5 * (gamma - 1) * (dV/a_inf).^2 .* rp2y;
    rhox = (1 / (gamma - 1)) * T.^(1 / (gamma - 1) - 1) .* Tx * r_inf;
    rhoy = (1 / (gamma - 1)) * T.^(1 / (gamma - 1) - 1) .* Ty * r_inf;
    px = (gamma / (gamma - 1)) * T.^(gamma / (gamma - 1) - 1) .* Tx * p_inf / (r_inf * a_inf^2);
    py = (gamma / (gamma - 1)) * T.^(gamma / (gamma - 1) - 1) .* Ty * p_inf / (r_inf * a_inf^2);
    
    % Compute time derivative
    ut = ux .* (-u_inf) + uy .* (-v_inf);
    vt = vx .* (-u_inf) + vy .* (-v_inf);
    temperature_t = Tx .* (-u_inf) + Ty .* (-v_inf);
    rhot = (1 / (gamma - 1)) .* T.^(1 / (gamma - 1) - 1) .* temperature_t;
    pt = (gamma / (gamma - 1)) .* T.^(gamma / (gamma - 1) - 1) .* temperature_t * p_inf / (r_inf * a_inf^2);
    
    % Store results in output arrays
    w = [rho, u, v, p];
    wx = [rhox, ux, vx, px];
    wy = [rhoy, uy, vy, py];
    wt = [rhot, ut, vt, pt];
end
