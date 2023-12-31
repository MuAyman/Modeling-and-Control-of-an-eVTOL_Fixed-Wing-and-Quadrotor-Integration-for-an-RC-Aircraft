
% Plot the closed-loop response.
time = 0:Ts:Duration;
yreftot = QuadrotorReferenceTrajectory(time)';

% Plot the states.
figure('Name','States')

subplot(2,3,1)
hold on
plot(time,xHistory(:,1))
plot(time,yreftot(:,1))
grid on
xlabel('time (s)')
ylabel('x (m)')
legend('actual','reference','Location','southeast')
title('Quadrotor x position')

subplot(2,3,2)
hold on
plot(time,xHistory(:,2))
plot(time,yreftot(:,2))
grid on
xlabel('time (s)')
ylabel('y (m)')
legend('actual','reference','Location','southeast')
title('Quadrotor y position')

subplot(2,3,3)
hold on
plot(time,xHistory(:,3))
plot(time,yreftot(:,3))
grid on
xlabel('time (s)')
ylabel('z (m)')
legend('actual','reference','Location','southeast')
title('Quadrotor z position')

subplot(2,3,4)
hold on
plot(time,xHistory(:,4))
plot(time,yreftot(:,4))
grid on
xlabel('time (s)')
ylabel('phi (rad)')
legend('actual','reference','Location','southeast')
title('Quadrotor phi angle')

subplot(2,3,5)
hold on
plot(time,xHistory(:,5))
plot(time,yreftot(:,5))
grid on
xlabel('time (s)')
ylabel('theta (rad)')
legend('actual','reference','Location','southeast')
title('Quadrotor theta angle')

subplot(2,3,6)
hold on
plot(time,xHistory(:,6))
plot(time,yreftot(:,6))
grid on
xlabel('time (s)')
ylabel('epsi (rad)')
legend('actual','reference','Location','southeast')
title('Quadrotor psi angle')

% Plot the manipulated variables.
figure('Name','Control Inputs')

subplot(2,2,1)
hold on
stairs(time,uHistory(:,1))
...ylim([-0.5,12.5])
plot(time,nloptions.MVTarget(2)*ones(1,length(time)))
grid on
xlabel('time (s)')
ylabel('omega (RPM^2)')
legend('actual','reference')
title('Motor 1')

subplot(2,2,2)
hold on
stairs(time,uHistory(:,2))
...ylim([-0.5,12.5])
plot(time,nloptions.MVTarget(2)*ones(1,length(time)))
grid on
xlabel('time (s)')
ylabel('omega (RPM^2)')
title('Motor 2')
legend('actual','reference')

subplot(2,2,3)
hold on
stairs(time,uHistory(:,3))
...ylim([-0.5,12.5])
plot(time,nloptions.MVTarget(2)*ones(1,length(time)))
grid on
xlabel('time (s)')
ylabel('omega (RPM^2)')
title('Motor 3')
legend('actual','reference')

subplot(2,2,4)
hold on
stairs(time,uHistory(:,4))
...ylim([-0.5,12.5])
plot(time,nloptions.MVTarget(2)*ones(1,length(time)))
grid on
xlabel('time (s)')
ylabel('omega (RPM^2)')
title('Motor 4')
legend('actual','reference')


