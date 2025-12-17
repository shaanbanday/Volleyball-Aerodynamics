% Appendix B: MATLAB Code for Volleyball Trajectory
% (Reconstructed from report appendix images)
% Script with local functions at the bottom.

clear; clc; close all;

%% Physical constants and ball properties
rho = 1.225;        % air density [kg/m^3]
mu  = 1.80e-5;      % dynamic viscosity [Pa*s]
g   = 9.81;         % gravity [m/s^2]
m   = 0.266;        % volleyball mass [kg]
r   = 0.105;        % volleyball radius [m]
A   = pi*r^2;       % cross-sectional area [m^2]

%% Court geometry
x0     = 1.5;       % serve contact location [m]
y0     = 3.4;       % serve height [m]
x_net  = 9.0;       % net position [m]
h_net  = 2.43;      % net height (men's) [m]
x_end  = 18.0;      % endline [m]

%% Initial serve parameters (user-editable)
V0    = 38.0;       % initial speed [m/s]
theta = 7.0;        % launch angle [deg]
omega = 90.0;       % spin rate [rad/s] (topspin)

% Initial velocity components
vx0 = V0*cosd(theta);
vy0 = V0*sind(theta);

% Spin vector (topspin into the page, negative z)
omega_vec = [0; 0; -omega];

%% Initial state vector
X0 = [x0; y0; vx0; vy0];

%% ODE integration
opts = odeset( ...
    'RelTol', 1e-8, ...
    'AbsTol', 1e-9, ...
    'Events', @(t,X) stopEvent(t, X, x_end) );

[t, X] = ode45(@(t,X) eom(t, X, rho, mu, A, m, r, omega_vec), [0 5], X0, opts);

x  = X(:,1);
y  = X(:,2);

%% Net clearance check
y_at_net   = NaN;
clears_net = false;

if any(x >= x_net)
    y_at_net   = interp1(x, y, x_net);
    clears_net = (y_at_net > h_net);
end

%% Landing location (interpolate where y crosses 0)
x_land = x(end);

idx = find(y <= 0, 1, 'first');
if ~isempty(idx) && idx > 1
    x_land = interp1(y(idx-1:idx), x(idx-1:idx), 0);
end

%% Print results
fprintf('Inputs: V0 = %.1f m/s, theta = %.1f deg, omega = %.1f rad/s\n', V0, theta, omega);
fprintf('Initial spin ratio Sp = %.3f\n', abs(omega)*r/max(V0, eps));

if ~isnan(y_at_net)
    fprintf('Net height at x = 9m: y = %.2f m (%s)\n', ...
        y_at_net, ternary(clears_net, 'clears net', 'hits net'));
end

fprintf('Landing: x = %.2f m\n', x_land);

if clears_net && x_land >= 9 && x_land <= 18
    fprintf('Verdict: IN\n');
else
    fprintf('Verdict: OUT\n');
end

%% Plot trajectory
figure; hold on; grid on; box on;

plot(x, y, 'LineWidth', 2);

xlabel('x [m]');
ylabel('y [m]');
title('Volleyball Topspin Serve Trajectory');

% Court + net
yline(0, 'k');
xline(x_net, 'k:');
plot([x_net x_net], [0 h_net], 'r', 'LineWidth', 3);

% Landing marker
plot(x_land, max(y(end), 0), 'ms', 'MarkerFaceColor', 'm');

% Net marker
if ~isnan(y_at_net)
    plot(x_net, y_at_net, 'bo', 'MarkerFaceColor', 'b');
end

xlim([0 19]);
ylim([-0.2 4]);

legend('Trajectory', 'Floor', 'Net', 'Location', 'best');

%% ---------------- Local functions ----------------

function dXdt = eom(~, X, rho, mu, A, m, r, omega_vec)
    vx = X(3);
    vy = X(4);

    v = sqrt(vx^2 + vy^2) + 1e-9;
    v_hat = [vx; vy; 0] / v;

    % Reynolds number and drag coefficient
    Re = rho * v * (2*r) / mu;
    Cd = Cd_sphere(Re);

    % Magnus coefficient (as function of spin ratio)
    Sp = norm(omega_vec) * r / v;
    Cm = Cm_quad(Sp);

    % Aerodynamic forces
    Fd = -0.5 * rho * Cd * A * v^2 * v_hat(1:2);

    if norm(omega_vec) < 1e-12
        FM = [0; 0];
    else
        omega_hat = omega_vec / norm(omega_vec);
        FM_3d = 0.5 * rho * Cm * A * v^2 * cross(omega_hat, v_hat);
        FM = FM_3d(1:2);
    end

    % Gravity
    Fg = [0; -m*9.81];

    % Accelerations
    ax = (Fd(1) + FM(1)) / m;
    ay = (Fd(2) + FM(2) + Fg(2)) / m;

    dXdt = [vx; vy; ax; ay];
end

function [value, isterminal, direction] = stopEvent(~, X, x_end)
    % Stop when ball hits ground (y = 0) OR passes endline (x = x_end)
    value      = [X(2); X(1) - x_end];   % ground or past end line
    isterminal = [1; 1];
    direction  = [-1; 1];
end

function Cd = Cd_sphere(Re)
    % Drag correlation used in the report (Equation 6 form). :contentReference[oaicite:2]{index=2}
    Cd = 24./Re ...
        + (2.6*(Re/5)) ./ (1 + (Re/5).^1.52) ...
        + (0.411*(Re/2.63e5).^(-7.94)) ./ (1 + (Re/2.63e5).^(-8.00)) ...
        + (0.25*(Re/1e6)) ./ (1 + (Re/1e6));
    Cd = max(Cd, 0.05);
end

function Cm = Cm_quad(Sp)
    % Quadratic Magnus coefficient fit from the report (Equation 7). :contentReference[oaicite:3]{index=3}
    Cm = -0.0398*Sp.^2 + 0.8260*Sp - 0.0005;
end

function out = ternary(cond, a, b)
    if cond
        out = a;
    else
        out = b;
    end
end
