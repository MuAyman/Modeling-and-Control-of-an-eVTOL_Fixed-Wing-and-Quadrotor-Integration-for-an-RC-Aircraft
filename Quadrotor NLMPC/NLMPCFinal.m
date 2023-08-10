close all;

% Create a nonlinear MPC object with 12 states, 12 outputs, and 4 inputs.
nx = 12;
ny = 12;
nu = 4;
nlmpcobj = nlmpc(nx, ny, nu);


nlmpcobj.Model.StateFcn = @mStateFcn;
nlmpcobj.Jacobian.StateFcn = @mStateJacobian;
...nlmpcobj.Model.OutputFcn = @mOutFcn; % in case of ny < nx
% Test consistency
rng(0)
validateFcns(nlmpcobj,rand(nx,1),rand(nu,1));


% Specify a sample time of 0.1 seconds, a prediction horizon of 18 steps, and control horizon of 2 steps. 
Ts = P.Ts;
p = 18;
m = 2;
nlmpcobj.Ts = Ts;
nlmpcobj.PredictionHorizon = p;
nlmpcobj.ControlHorizon = m;

nlmpcobj.MV = struct( ...
    Min={0;0;0;0}, ...
    Max={70;70;70;70}, ...
    RateMin={-2;-2;-2;-2}, ...
    RateMax={2;2;2;2} ...
    );

nlmpcobj.Weights.OutputVariables = [1 1 1 1 1 1 0 0 0 0 0 0];
nlmpcobj.Weights.ManipulatedVariables = [0.1 0.1 0.1 0.1];
nlmpcobj.Weights.ManipulatedVariablesRate = [0.1 0.1 0.1 0.1];


% Specify the initial conditions
x = [P.pn0; P.pe0; P.pd0; P.phi0; P.theta0; P.epsi0; P.ub0; P.vb0; P.wb0; P.p0; P.q0; P.r0];

% Nominal control target (average to keep quadrotor floating)
nloptions = nlmpcmoveopt;
nloptions.MVTarget = [hover hover hover hover];      % P.m * P.g
mv = nloptions.MVTarget;

% Simulation duration in seconds
Duration = 70;

% Display waitbar to show simulation progress
hbar = waitbar(0,"Simulation Progress");

% MV last value is part of the controller state
lastMV = mv;

% Store states for plotting purposes
xHistory = x';
uHistory = lastMV;

% Simulation loop
for k = 1:(Duration/Ts)

    % Set references for previewing
    t = linspace(k*Ts, (k+p-1)*Ts,p);
    yref = QuadrotorReferenceTrajectory(t);

    % Compute control move with reference previewing
    xk = xHistory(k,:);
...    xk = awgn(xk,10,'measured');     % add noise
    [uk,nloptions,info] = nlmpcmove(nlmpcobj,xk,lastMV,yref',[],nloptions);

    % Store control move
    uHistory(k+1,:) = uk';
    lastMV = uk;

    % Simulate quadrotor for the next control interval (MVs = uk) 
    ODEFUN = @(t,xk) mStateFcn(xk,uk);
    [TOUT,XOUT] = ode45(ODEFUN,[0 Ts], xHistory(k,:)');

    % Update quadrotor state
    xHistory(k+1,:) = XOUT(end,:);

    % Update waitbar
    waitbar(k*Ts/Duration,hbar);
end

% Close waitbar 
close(hbar)

plotTrajectory
.... animateTrajectory




