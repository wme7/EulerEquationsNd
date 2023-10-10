function [q, qx, qy, qt] = exactSolution_entropyWaves2d(x, y, time)
    % Parameters
    gamma = 1.4;

    % Free stream values
    r_inf = 1.0;
    u_inf = 0.6;
    v_inf = 1.4;
    p_inf = 1.0 / gamma;
    Q_inf = u_inf + v_inf;

    % Compute the exact solution
    Amplitude = 0.5;
    coef = 2*pi;

    r = r_inf + Amplitude * sin(coef * (x + y - Q_inf * time));
    u = u_inf * ones(size(x));
    v = v_inf * ones(size(x));
    p = p_inf * ones(size(x));

    % Compute derivatives
    rx = Amplitude * coef * cos(coef * (x + y - Q_inf * time));
    ux = zeros(size(x));
    vx = zeros(size(x));
    px = zeros(size(x));

    ry = Amplitude * coef * cos(coef* (x + y - Q_inf * time));
    uy = zeros(size(y));
    vy = zeros(size(y));
    py = zeros(size(y));

    rt = - Amplitude * coef * Q_inf .* cos(coef*(x + y - Q_inf * time));
    ut = zeros(size(y));
    vt = zeros(size(y));
    pt = zeros(size(y));

    % Store results in output arrays
    q = [r, u, v, p];
    qx = [rx, ux, vx, px];
    qy = [ry, uy, vy, py];
    qt = [rt, ut, vt, pt];

end % funtion