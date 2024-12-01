%%PSC for 6dof rigid body using DQ 
%%Rotation sequence - 3,2,1
%% Using Lagrange ploynomial to approximate states and control are approximated by Lagrange ploynomial
%% Collocation points are generated as the CGL nodes and clenshaw curtis weights are used for approximating the integration.
clc; clear; close all;
set(0, 'defaultFigureWindowState', 'maximized');
set(0, 'defaultAxesFontSize', 20);
set(0, 'defaultAxesLineWidth', 1.5);
set(groot, 'defaultAxesFontName','Century');
set(groot, 'defaultLegendFontName','Century');
set(groot, 'defaultTextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultColorbarTickLabelinterpreter','latex'); 
set(groot, 'defaultLineLineWidth', 1.5);

addpath("RK/","dq_fcn/") 
%% initial conds and final conds and bounds
tau0 = 0;
tau_f = 650; %10 minutes
omega = 0.0012; % orbital velocity in rad/s 
mc = 20;
alpha = pi/3;   %FOV

%initial condition 

 
q_b_i = quaternion(1,[0 0 0]); % initial quaternion of body frame w.r.t inertial frame
q_b_i_conj = conj(q__i);

q_dq_b_i = dualquaternion();
q_dq_b_i.qr = q_b_i; 
r_b = quaternion();
q_dq_b_i.qd =  0.5*q_dq_b_i.qr*r_b; % rotation first followed by translation

xini = 0.0;
xfin = 0.0;
yini = 0.0;
yfin = 0.0;
zini = 0.0;
zfin = 0.0;

xdotini = 0.0;
xdotfin = 0.005;
ydotini = 0.0;
ydotfin = 0.005;
zdotini = 0.0;
zdotfin = 0.005;

% xini = 0.05;
% xfin = 0.02;
% yini = 0.05;
% yfin = 0.05;
% zini = -10.05;
% zfin = -0.1;
% 
% xdotini = 0.0;
% xdotfin = 0.015;
% ydotini = 0.0;
% ydotfin = 0.005;
% zdotini = 0.0;
% zdotfin = 0.05;


bc = [xini; xdotini; yini; ydotini; zini; zdotini; xfin; xdotfin; yfin; ydotfin; zfin; zdotfin];


Txmin = -0.3;
Txmax = 0.3;
Tymin = -0.3;
Tymax = 0.3;
Tzmin = -0.3;
Tzmax = 0.3;

Fxmin = -0.3;
Fxmax = 0.3;
Fymin = -0.3;
Fymax = 0.3;
Fzmin = -0.3;
Fzmax = 0.3;

%% Node distribution, Clenshaw Curtis weights and D matrix
a = -1;
b = 1;
N = 70;
tk = (((b - a) / 2) .* cos(linspace(0, pi, N + 1)) + (b + a) / 2)';
D = -Dmatrix_CGL(tk);
w = flip(cc_quad_weights(N));
tk = flip(tk);
tau = ((tau_f-tau0).*tk+(tau_f+tau0))./2;

%% bounds
lb = [xmin.*ones(N+1, 1); xdotmin.*ones(N+1, 1); ymin.*ones(N+1, 1); ydotmin.*ones(N+1, 1); zmin.*ones(N+1, 1); zdotmin.*ones(N+1, 1); Txmin.*ones(N+1, 1); Tymin.*ones(N+1, 1); Tzmin.*ones(N+1, 1)];
ub = [xmax.*ones(N+1, 1); xdotmax.*ones(N+1, 1); ymax.*ones(N+1, 1); ydotmax.*ones(N+1, 1); zmax.*ones(N+1, 1); zdotmax.*ones(N+1, 1); Txmax.*ones(N+1, 1); Tymax.*ones(N+1, 1); Tzmax.*ones(N+1, 1)];

%% equality constraints and linear constraints
A = []; B = []; Aeq = []; Beq = [];

%% initial guess for decision vector - [X(0)...X(N) Y(0)...Y(N) Z(0)...Z(N) 
% Xdot(0)...X(N) Ydot(0)...Ydot(N) Zdot(0)...Zdot(N) 
% Tx(0)...Tx(N) Ty(0)...Ty(N) Tz(0)...Tz(N)]
load("guess.mat");
%load("position_lvlh.mat");
Pos = resample(Position, tau);
xguess = Pos.Data(:, 1);
yguess = Pos.Data(:, 2);
zguess = Pos.Data(:, 3);

%load("velocity_lvlh.mat");
Vel = resample(Velocity, tau);
xdotguess = Vel.Data(:, 1);
ydotguess = Vel.Data(:, 2);
zdotguess = Vel.Data(:, 3);

%load('thrust_pid_lvlh.mat');
Th = resample(Thrust, tau);
Txguess = Th.Data(:, 1);
Tyguess = Th.Data(:, 2);
Tzguess = Th.Data(:, 3);

DV0 = awgn([xguess; xdotguess; yguess; ydotguess; zguess; zdotguess; Txguess; Tyguess; Tzguess],650);
%DV0 = [xguess; xdotguess; yguess; ydotguess; zguess; zdotguess; Txguess; Tyguess; Tzguess];

%% Optimization options
options =  optimoptions ('fmincon','Display','Iter','OptimalityTolerance',...
1e-4, 'ConstraintTolerance', 1e-1, 'MaxIterations', 2000,'MaxFunctionEvaluations',...
500000,'Algorithm','sqp');

[DV, costval, exitflag, output] = fmincon(@(DV)costfunc(DV, w, tau0, tau_f), DV0, A, B,...
    Aeq, Beq, lb, ub, @(DV)nonlcon(DV, D, bc, w, tau0, tau_f, omega, mc, alpha),options);

exitflag
output