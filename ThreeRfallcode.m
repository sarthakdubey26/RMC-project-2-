
clear; clc; close all;

% Model
L = [0.7; 0.9; 0.6];   % link lengths [m]
m = [1.0; 0.8; 0.5];   % point masses at tips [kg]
g = 9.81;              % gravity [m/s^2], +y downward

dt     = 0.001;        % timestep [s]
Tfinal = 4.0;          % total simulation time [s]
N      = round(Tfinal/dt);

% Initial state 
theta  = [0; 0; 0];    % initial angles [rad]
dtheta = [0; 0; 0];    % initial angular velocities [rad/s]

thhistory        = zeros(3, N+1);
thhistory(:, 1)  = theta;

% Euler integration loop 
for k = 1:N
    tau = [0; 0; 0];   % no actuation: free fall

    ddtheta = Forwarddynamics(theta, dtheta, tau, L, m, g);

    theta  = theta  + dtheta * dt;
    dtheta = dtheta + ddtheta * dt;

    thhistory(:, k+1) = theta;
end

time = 0:dt:Tfinal;

% joint angles vs time
figure;
plot(time, thhistory(1,:), 'LineWidth', 1.5); hold on;
plot(time, thhistory(2,:), 'LineWidth', 1.5);
plot(time, thhistory(3,:), 'LineWidth', 1.5);
grid on;
xlabel('Time [s]');
ylabel('Joint angle [rad]');
legend('\theta_1','\theta_2','\theta_3','Location','best');
title('Joint angles vs time (3R arm falling under gravity)');

% End-effector trajectory plot 
xe = zeros(1, N+1);
ye = zeros(1, N+1);

for k = 1:N+1
    [~, ~, p3] = ForwardKinematics3R(thhistory(:,k), L);
    xe(k) = p3(1);
    ye(k) = p3(2);
end

figure;
plot(xe, ye, 'LineWidth', 1.5);
axis equal;
grid on;
xlabel('x [m]');
ylabel('y [m]');
title('End-effector trajectory in the plane');

% Figure 3: animation 
figure;
axis equal;
axis([-sum(L) sum(L) -sum(L) sum(L)]);
grid on;
xlabel('x [m]');
ylabel('y [m]');
title('3R planar arm falling under gravity');

for k = 1:20:N+1
    th = thhistory(:,k);
    [p1, p2, p3] = ForwardKinematics3R(th, L);

    x0 = 0; y0 = 0;
    x1 = p1(1); y1 = p1(2);
    x2 = p2(1); y2 = p2(2);
    x3 = p3(1); y3 = p3(2);

    cla;
    hold on;
    plot([x0 x1], [y0 y1], 'r-', 'LineWidth', 2);
    plot([x1 x2], [y1 y2], 'g-', 'LineWidth', 2);
    plot([x2 x3], [y2 y3], 'b-', 'LineWidth', 2);
    plot(x3, y3, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k');

    axis equal;
    axis([-sum(L) sum(L) -sum(L) sum(L)]);
    grid on;
    xlabel('x [m]');
    ylabel('y [m]');
    title(sprintf('Time = %.2f s', time(k)));
    drawnow;
end

% Helper functions

function [p1, p2, p3] = ForwardKinematics3R(theta, L)
    th1 = theta(1); th2 = theta(2); th3 = theta(3);
    L1 = L(1); L2 = L(2); L3 = L(3);

    x1 = L1 * cos(th1);
    y1 = L1 * sin(th1);
    p1 = [x1; y1];

    x2 = x1 + L2 * cos(th1 + th2);
    y2 = y1 + L2 * sin(th1 + th2);
    p2 = [x2; y2];

    x3 = x2 + L3 * cos(th1 + th2 + th3);
    y3 = y2 + L3 * sin(th1 + th2 + th3);
    p3 = [x3; y3];
end

function [p1, p2, p3, Jv1, Jv2, Jv3] = KinematicsJacobians3R(theta, L)
    th1 = theta(1); th2 = theta(2); th3 = theta(3);
    L1 = L(1); L2 = L(2); L3 = L(3);

    th12  = th1 + th2;
    th123 = th1 + th2 + th3;

    c1   = cos(th1);   s1   = sin(th1);
    c12  = cos(th12);  s12  = sin(th12);
    c123 = cos(th123); s123 = sin(th123);

    x1 = L1 * c1; y1 = L1 * s1;
    p1 = [x1; y1];

    x2 = x1 + L2 * c12; y2 = y1 + L2 * s12;
    p2 = [x2; y2];

    x3 = x2 + L3 * c123; y3 = y2 + L3 * s123;
    p3 = [x3; y3];

    Jv1 = [ -L1*s1,          0,           0;
             L1*c1,          0,           0 ];

    Jv2 = [ -L1*s1 - L2*s12, -L2*s12,     0;
             L1*c1 + L2*c12,  L2*c12,     0 ];

    Jv3 = [ -L1*s1 - L2*s12 - L3*s123, -L2*s12 - L3*s123, -L3*s123;
             L1*c1 + L2*c12 + L3*c123,  L2*c12 + L3*c123,  L3*c123 ];
end

function M = Massmatrix(theta, L, m)
    [~, ~, ~, Jv1, Jv2, Jv3] = KinematicsJacobians3R(theta, L);
    m1 = m(1); m2 = m(2); m3 = m(3);
    M = m1*(Jv1')*Jv1 + m2*(Jv2')*Jv2 + m3*(Jv3')*Jv3;
end

function G = Gravitytorque(theta, L, m, g)
    [~, ~, ~, Jv1, Jv2, Jv3] = KinematicsJacobians3R(theta, L);
    m1 = m(1); m2 = m(2); m3 = m(3);

    F1 = [0; m1*g];
    F2 = [0; m2*g];
    F3 = [0; m3*g];

    G = Jv1' * F1 + Jv2' * F2 + Jv3' * F3;
end

function ddtheta = Forwarddynamics(theta, ~, tau, L, m, g)
    M = Massmatrix(theta, L, m);
    G = Gravitytorque(theta, L, m, g);
    ddtheta = M \ (tau - G);
end
