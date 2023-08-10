function [ xdesired ] = QuadrotorReferenceTrajectory( t )


% the trajectory points

x = zeros(1,length(t)); ...3.5*sin(t); ...3*sin(t/3);
y = zeros(1,length(t)); ... 0.5*sin(t); ... -9*sin(t/3).*cos(t/3);
z = 100*ones(1,length(t)); ...0*cos(t/3);
phi = zeros(1,length(t));
theta = zeros(1,length(t));
psi = zeros(1,length(t));
xdot = zeros(1,length(t));
ydot = zeros(1,length(t));
zdot = zeros(1,length(t));
phidot = zeros(1,length(t));
thetadot = zeros(1,length(t));
psidot = zeros(1,length(t));

xdesired = [x;y;z;phi;theta;psi;xdot;ydot;zdot;phidot;thetadot;psidot];
end

