%% Quadrotr Inertial Parameters
P.Jx = 4.43;    % kg.m2
P.Jy = 6.01;
P.Jz = 9.68;
P.m = 24;         % kg
P.g = 9.806650;
P.k = 1;    % k: lift constant
P.l = 1.0105;   % l: distance between rotor and center of mass
P.b = 0.2;  % b: drag constant

%% Initial conditions
P.pn0 = 0.0;        % initial X position, m
P.pe0 = 0.0;        % initial y position, m
P.pd0 = 100;        % initial z position, m
P.phi0 = 10 * pi/180;       % initial roll angle
P.theta0 = 0 * pi/180;     % initial pitch angle
P.epsi0 = 0 * pi/180;      % initial yaw angle
P.ub0 = 0.0;        % initial u velocity
P.vb0 = 0.0;        % initial v velocity
P.wb0 = 0.0;        % initial z velocity
P.p0 = 0.0;         % initial roll rate angle
P.q0 = 0.0;         % initial pitch rate angle
P.r0 = 0.0;         % initial yaw rate angle

%% Simulation Parameters
P.t_start = 0.0;  % Start time of simulation
P.t_end = 20.0;   % End time of simulation
P.Ts = 0.1;      % sample time for simulation

hover = (P.m*P.g)/4;

%% Quad Model
% Create symbolic functions for time-dependent states
syms x(t) y(t) z(t) phi(t) theta(t) psi(t)

% Transformation matrix for angular velocities from inertial frame
% to body frame
W = [ 1,  0,        -sin(theta);
      0,  cos(phi),  cos(theta)*sin(phi);
      0, -sin(phi),  cos(theta)*cos(phi) ];

% Rotation matrix R_ZYX from body frame to inertial frame 
R = rotationMatrixEulerZYX(phi,theta,psi);

% Create symbolic variables for diagonal elements of inertia matrix
syms Ixx Iyy Izz

% Jacobian that relates body frame to inertial frame velocities
I = [Ixx, 0, 0; 0, Iyy, 0; 0, 0, Izz];
J = W.'*I*W;

% Coriolis matrix
dJ_dt = diff(J);
h_dot_J = [diff(phi,t), diff(theta,t), diff(psi,t)]*J;
grad_temp_h = transpose(jacobian(h_dot_J,[phi theta psi]));
C = dJ_dt - 1/2*grad_temp_h;
C = subsStateVars(C,t);

% k: lift constant
% l: distance between rotor and center of mass
% b: drag constant
% ui: squared angular velocity of rotor i as control input
syms k l m b g u1 u2 u3 u4

% Torques in the direction of phi, theta, psi
tau_beta = [l*k*(u2+u3-u1-u4); l*k*(u1+u2-u3-u4); b*(u1-u2+u3-u4)];
% Total thrust
T = k*(u1+u2+u3+u4);

% Create state variables consisting of positions, angles,
% and their derivatives
state = [x; y; z; phi; theta; psi; diff(x,t); diff(y,t); ...
    diff(z,t); diff(phi,t); diff(theta,t); diff(psi,t)];
state = subsStateVars(state,t);

f = [state(7:12); -g*[0;0;1] + R*[0;0;T]/m; inv(J)*(tau_beta - C*state(10:12))];
f = subsStateVars(f,t);

f = subs(f, [Ixx Iyy Izz k l m b g], ...
    [P.Jx P.Jy P.Jz P.k P.l P.m P.b P.g]);
f = simplify(f);

y = [state(1:12)];
y = subsStateVars(y,t);

% Calculate Jacobians for nonlinear prediction model
A = jacobian(f,state);
control = [u1; u2; u3; u4];
B = jacobian(f,control);
C = jacobian(y,state);

%% Save Funtions
% Create mStateFcn.m with current state and control
% vectors as inputs and the state time-derivative as outputs
matlabFunction(f,"File","mStateFcn", ...
    "Vars",{state,control});

% Create mStateJacobian.m with current state and control
% vectors as inputs and the Jacobians of the state time-derivative
% as outputs
matlabFunction(A,B,"File","mStateJacobian", ...
    "Vars",{state,control});

% Create mOutFcn.m with current state and control
% vectors as inputs and the output time-derivative as outputs
matlabFunction(y,"File","mOutFcn", ...
    "Vars",{state,control});

% Create mOutJacobianFcn.m with current state and control
% vectors as inputs and the output Jacobians time-derivative as outputs
matlabFunction(C,"File","mOutJacobianFcn", ...
    "Vars",{state,control});


%% Helppers
function [Rz,Ry,Rx] = rotationMatrixEulerZYX(phi,theta,psi)
% Euler ZYX angles convention
    Rx = [ 1,           0,          0;
           0,           cos(phi),  -sin(phi);
           0,           sin(phi),   cos(phi) ];
    Ry = [ cos(theta),  0,          sin(theta);
           0,           1,          0;
          -sin(theta),  0,          cos(theta) ];
    Rz = [cos(psi),    -sin(psi),   0;
          sin(psi),     cos(psi),   0;
          0,            0,          1 ];
    if nargout == 3
        % Return rotation matrix per axes
        return;
    end
    % Return rotation matrix from body frame to inertial frame
    Rz = Rz*Ry*Rx;
end

function stateExpr = subsStateVars(timeExpr,var)
    if nargin == 1 
        var = sym("t");
    end
    repDiff = @(ex) subsStateVarsDiff(ex,var);
    stateExpr = mapSymType(timeExpr,"diff",repDiff);
    repFun = @(ex) subsStateVarsFun(ex,var);
    stateExpr = mapSymType(stateExpr,"symfunOf",var,repFun);
    stateExpr = formula(stateExpr);
end

function newVar = subsStateVarsFun(funExpr,var)
    name = symFunType(funExpr);
    name = replace(name,"_Var","");
    stateVar = "_" + char(var);
    newVar = sym(name + stateVar);
end

function newVar = subsStateVarsDiff(diffExpr,var)
    if nargin == 1 
      var = sym("t");
    end
    c = children(diffExpr);
    if ~isSymType(c{1},"symfunOf",var)
      % not f(t)
      newVar = diffExpr;
      return;
    end
    if ~any([c{2:end}] == var)
      % not derivative wrt t only
      newVar = diffExpr;
      return;
    end
    name = symFunType(c{1});
    name = replace(name,"_Var","");
    extension = "_" + join(repelem("d",numel(c)-1),"") + "ot";
    stateVar = "_" + char(var);
    newVar = sym(name + extension + stateVar);
end