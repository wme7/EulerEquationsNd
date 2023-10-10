function [q, qx, qy, qz, qt] = exactSolution_entropyWaves3d(x, y, z, time)
    % Parameters
    gamma = 1.4;

    % Free stream values
    r_inf = 1.0;
    u_inf = 0.6;
    v_inf = 1.4;
    w_inf = 1.0;
    p_inf = 1.0 / gamma;
    Q_inf = u_inf + v_inf + w_inf;

    % Compute the exact solution
    Amplitude = 0.5;
    coef = 2*pi;

    r = r_inf + Amplitude * sin(coef * (x + y + z - Q_inf * time));
    u = u_inf * ones(size(x));
    v = v_inf * ones(size(x));
    w = w_inf * ones(size(x));
    p = p_inf * ones(size(x));

    % Compute derivatives
    rx = Amplitude * coef * cos(coef * (x + y + z - Q_inf * time));
    ux = zeros(size(x));
    vx = zeros(size(x));
    wx = zeros(size(x));
    px = zeros(size(x));

    ry = Amplitude * coef * cos(coef * (x + y + z - Q_inf * time));
    uy = zeros(size(y));
    vy = zeros(size(y));
    wy = zeros(size(y));
    py = zeros(size(y));

    rz = Amplitude * coef * cos(coef * (x + y + z - Q_inf * time));
    uz = zeros(size(z));
    vz = zeros(size(z));
    wz = zeros(size(z));
    pz = zeros(size(z));

    rt = - Amplitude * coef * Q_inf .* cos(coef * (x + y + z - Q_inf * time));
    ut = zeros(size(y));
    vt = zeros(size(y));
    wt = zeros(size(y));
    pt = zeros(size(y));

    % Store results in output arrays
    q = [r, u, v, w, p];
    qx = [rx, ux, vx, wx, px];
    qy = [ry, uy, vy, wy, py];
    qz = [rz, uz, vz, wz, pz];
    qt = [rt, ut, vt, wt, pt];

end % funtion