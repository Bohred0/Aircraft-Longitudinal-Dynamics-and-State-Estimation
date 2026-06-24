% Longitudinal aircraft state-space model obtained from 
% S. M. Özer, "Delay Effects in Stability Augmentation System of Aircraft Longitudinal Dynamics: A Case Study on Boeing 747,"

A = [ ...
    -0.00276  0.0389   -62.1    -32.1;
    -0.0654  -0.3191   771.51   -2.5994;
    0.0002   -0.001013 -0.4285   0.0003;
    0         0         1         0];

B = [1.44;
    -18.021;
    -1.1579;
    0];
C=eye(4);
D=0;
sys = ss(A, B, C, D);
pole=eig(sys);

figure
plot(real(poles), imag(poles), 'x');
xlabel('Real Part');
ylabel('Imaginary Part');
title('Pole-Zero Plot');
grid on;

t = 0:0.01:500;
x0 = [0;
    0;
    0;
    deg2rad(5)];
q0=[0;
    0;
    deg2rad(5);
    0];

[y1,t,x1] = initial(sys,x0,t);
[y2,t,x2] = initial(sys,q0,t);
figure(1)
subplot(2,2,1)
plot(t,x1(:,1),'LineWidth',1.5)
grid on
title('Forward Velocity Perturbation, u')
xlabel('Time (s)')
ylabel('u (m/s)')

subplot(2,2,2)
plot(t,x1(:,2),'LineWidth',1.5)
grid on
title('Vertical Velocity Perturbation, w')
xlabel('Time (s)')
ylabel('w (m/s)')

subplot(2,2,3)
plot(t,rad2deg(x1(:,3)),'LineWidth',1.5)
grid on
title('Pitch Rate, q')
xlabel('Time (s)')
ylabel('q (deg/s)')

subplot(2,2,4)
plot(t,rad2deg(x1(:,4)),'LineWidth',1.5)
grid on
title('Pitch Angle, \theta')
xlabel('Time (s)')
ylabel('\theta (deg)')

sgtitle('Open-Loop Response to Initial Pitch Disturbance')

figure(2)

subplot(2,2,1)
plot(t,x2(:,1),'LineWidth',1.5)
grid on
title('Forward Velocity Perturbation, u')
xlabel('Time (s)')
ylabel('u (m/s)')

subplot(2,2,2)
plot(t,x2(:,2),'LineWidth',1.5)
grid on
title('Vertical Velocity Perturbation, w')
xlabel('Time (s)')
ylabel('w (m/s)')

subplot(2,2,3)
plot(t,rad2deg(x2(:,3)),'LineWidth',1.5)
grid on
title('Pitch Rate, q')
xlabel('Time (s)')
ylabel('q (deg/s)')

subplot(2,2,4)
plot(t,rad2deg(x2(:,4)),'LineWidth',1.5)
grid on
title('Pitch Angle, \theta')
xlabel('Time (s)')
ylabel('\theta (deg)')

sgtitle('Open-Loop Response to Initial Pitch Rate Disturbance')

rank(ctrb(A,B));
%controllable

%LQR
%To tune our LQR controller, we use Bryson's Rule to normalize states
%Allowable deviations: u:10m/s w:5m/s q:10deg/s theta:5deg controller:20deg
Q = diag([
    1/(10)^2,...
    1/(5)^2,...
    1/(deg2rad(10))^2,...
    1/(deg2rad(5))^2
]);

R = 1/(deg2rad(20))^2;
K = lqr(sys, Q, R);
Acl = A - B*K;
sys_cl = ss(Acl, B, C, D);

[y_cl_1,t,x_cl_1] = initial(sys_cl,x0,t);
[y_cl_2,t,x_cl_2] = initial(sys_cl,q0,t);

figure(1)
subplot(2,2,1)
hold on
plot(t,rad2deg(x_cl_1(:,3)),'LineWidth',1.5)
plot(t,rad2deg(x1(:,3)),'LineWidth',1.5)
grid on
title('Pitch Rate, q')
xlabel('Time (s)')
ylabel('q (deg/s)')
legend("LQR","Natural")

subplot(2,2,2)
hold on
plot(t,rad2deg(x_cl_1(:,4)),'LineWidth',1.5)
plot(t,rad2deg(x1(:,4)),'LineWidth',1.5)
grid on
title('Pitch Angle, \theta')
xlabel('Time (s)')
ylabel('\theta (deg)')
legend("LQR","Natural")


figure(2)
subplot(2,2,1)
hold on
plot(t,rad2deg(x_cl_2(:,3)),'LineWidth',1.5)
plot(t,rad2deg(x2(:,3)),'LineWidth',1.5)
grid on
title('Pitch Rate, q')
xlabel('Time (s)')
ylabel('q (deg/s)')
legend("LQR","Natural")

subplot(2,2,2)
hold on
plot(t,rad2deg(x_cl_2(:,4)),'LineWidth',1.5)
plot(t,rad2deg(x2(:,4)),'LineWidth',1.5)
grid on
title('Pitch Angle, \theta')
xlabel('Time (s)')
ylabel('\theta (deg)')
legend("LQR","Natural")

eig(sys_cl)

u_control_1= -K*x_cl_1';
u_control_2= -K*x_cl_2';

max(abs(rad2deg(u_control_1)))
max(abs(rad2deg(u_control_2)))

%Full State Estimation
C_obs = [0 0 1 0;
     0 0 0 1];
rank(obsv(A,C_obs));

xhat0=[0;0;0;0];

%Measurement Noise:
sigma_q     = deg2rad(0.1);   % rad/s
sigma_theta = deg2rad(0.05);  % rad

rng(1) % for repeat-ability

noise_q     = sigma_q     * randn(length(t),1);
noise_theta = sigma_theta * randn(length(t),1);

%Process Noise:
wind_u     = 0.5*randn(length(t),1);
wind_w     = 0.5*randn(length(t),1);
gust_q     = deg2rad(0.2)*randn(length(t),1);

w_proc = [wind_u wind_w gust_q];

%Simulating process noise in true system:
G=[1 0 0;
    0 1 0;
    0 0 1;
    0 0 0];
A_true= Acl;

sys_true=ss(A_true, G, eye(4), 0*G);

% Simulate the true system with process noise
[~, ~, x_true] = lsim(sys_true, w_proc, t, x0);
y_true=x_true*C_obs';

y_noisy=y_true + [noise_q noise_theta];

%Pole Placement(2x Poles):
L2 = place(Acl',C_obs',(2*eig(Acl)))';

sys_obs_2 = ss(Acl-L2*C_obs,...
    L2,...
    eye(4),...
    zeros(4,2));

xhat_2 = lsim(sys_obs_2,y_noisy,t,xhat0);

e2 = x_true - xhat_2;

%Pole Placement(5x Poles):
L5 = place(Acl',C_obs',(5*eig(Acl)))';

sys_obs_5 = ss(Acl-L5*C_obs,...
    L5,...
    eye(4),...
    zeros(4,2));

xhat_5 = lsim(sys_obs_5,y_noisy,t,xhat0);

e5 = x_true - xhat_5;
%%Pole Placement(10x Poles):
L10 = place(Acl',C_obs',(10*eig(Acl)))';

sys_obs_10 = ss(Acl-L10*C_obs,...
    L10,...
    eye(4),...
    zeros(4,2));

xhat_10 = lsim(sys_obs_10,y_noisy,t,xhat0);

e10 = x_true - xhat_10;

%Kalman Filter
%Noise Covariance Matrices:
Q_kf = diag([0.25, 0.25, deg2rad(0.2)^2]); % Process noise covariance
R_kf = diag([deg2rad(0.1)^2,deg2rad(0.05)^2]); % Measurement noise covariance

Kf = lqe(Acl, G ,C_obs,Q_kf,R_kf);

Aobs_Kf=Acl-Kf*C_obs;

Bobs_Kf=Kf;

% Define the observer system with the Kalman filter
sys_kf = ss(Aobs_Kf, Bobs_Kf, eye(4), 0*Bobs_Kf);

[xhat_kf, t_kf] = lsim(sys_kf, y_noisy, t, xhat0);

% Calculate the estimation error for the Kalman filter
e_kf = x_true - xhat_kf;


figure

subplot(2,2,1)
plot(t,x_true(:,1),'LineWidth',1.5)
hold on
plot(t,xhat_kf(:,1),'--','LineWidth',1.5)
grid on
title('u')
legend('True','Estimated(Kf)')

subplot(2,2,2)
plot(t,x_true(:,2),'LineWidth',1.5)
hold on
plot(t,xhat_kf(:,2),'--','LineWidth',1.5)
grid on
title('w')
legend('True','Estimated(Kf)')

subplot(2,2,3)
plot(t,rad2deg(x_true(:,3)),'LineWidth',1.5)
hold on
plot(t,rad2deg(xhat_kf(:,3)),'--','LineWidth',1.5)
grid on
title('q')
legend('True','Estimated(Kf)')

subplot(2,2,4)
plot(t,rad2deg(x_true(:,4)),'LineWidth',1.5)
hold on
plot(t,rad2deg(xhat_kf(:,4)),'--','LineWidth',1.5)
grid on
title('\theta')
legend('True','Estimated(Kf)')

sgtitle('Observer State Estimation with Kalman Filter')

figure
subplot(2,2,1)
title("2x Closed Loop Poles")
grid on
xlabel('Time (s)')
ylabel('Theta Error (deg)')
plot(t,rad2deg(e2(:,4)),'LineWidth',1.5)

subplot(2,2,2)
title("5x Closed Loop Poles")
grid on
xlabel('Time (s)')
ylabel('Theta Error (deg)')
plot(t,rad2deg(e5(:,4)),'LineWidth',1.5)

subplot(2,2,3)
title("10x Closed Loop Poles")
grid on
xlabel('Time (s)')
ylabel('Theta Error (deg)')
plot(t,rad2deg(e10(:,4)),'LineWidth',1.5)

subplot(2,2,4)
title("Kalman Filter")
grid on
xlabel('Time (s)')
ylabel('Theta Error (deg)')
plot(t,rad2deg(e_kf(:,4)),'LineWidth',2)

sgtitle('Pitch Angle Estimation Error Comparison')

figure
subplot(2,2,1)
title("2x Closed Loop Poles")
grid on
xlabel('Time (s)')
ylabel('q Error (deg/s)')
plot(t,rad2deg(e2(:,3)),'LineWidth',1.5)

subplot(2,2,2)
title("5x Closed Loop Poles")
grid on
xlabel('Time (s)')
ylabel('q Error (deg/s)')
plot(t,rad2deg(e5(:,3)),'LineWidth',1.5)

subplot(2,2,3)
title("10x Closed Loop Poles")
grid on
xlabel('Time (s)')
ylabel('q Error (deg/s)')
plot(t,rad2deg(e10(:,3)),'LineWidth',1.5)

subplot(2,2,4)
title("Kalman Filter")
grid on
xlabel('Time (s)')
ylabel('q Error (deg/s)')
plot(t,rad2deg(e_kf(:,3)),'LineWidth',2)

sgtitle('Pitch Rate Estimation Error Comparison')

figure

subplot(2,2,1)
plot(t,rad2deg(y_true(:,1)))
title('True q')

subplot(2,2,2)
plot(t,rad2deg(y_noisy(:,1)),'--')
title('Measured q (Kalman Filter)')

subplot(2,2,3)
plot(t,rad2deg(y_true(:,2)))
title('True \theta')

subplot(2,2,4)
plot(t,rad2deg(y_noisy(:,2)),'--')
title('Measured \theta (Kalman Filter)')

rms_u = rms(e_kf(:,1));
rms_w = rms(e_kf(:,2));
rms_q = rms(rad2deg(e_kf(:,3)));
rms_theta = rms(rad2deg(e_kf(:,4)));

fprintf('RMS u error = %.3f m/s\n',rms_u)
fprintf('RMS w error = %.3f m/s\n',rms_w)
fprintf('RMS q error = %.3f deg/s\n',rms_q)
fprintf('RMS theta error = %.3f deg\n',rms_theta)

fprintf('\n')
fprintf('Theta RMS Errors\n')
idx = t > 50;

fprintf('KF theta RMS after 50s = %.4f deg\n', ...
    rad2deg(rms(e_kf(idx,4))))
fprintf('2x poles theta RMS after 50s = %.4f deg\n', ...
    rad2deg(rms(e2(idx,4))))
fprintf('5x poles theta RMS after 50s = %.4f deg\n', ...
    rad2deg(rms(e5(idx,4))))
fprintf('10x poles theta RMS after 50s = %.4f deg\n', ...
    rad2deg(rms(e10(idx,4))))


fprintf('\n')
fprintf('Pitch Rate RMS Errors\n')
fprintf('KF q RMS after 50s = %.4f deg\n', ...
    rad2deg(rms(e_kf(idx,3))))
fprintf('2x poles q RMS after 50s = %.4f deg\n', ...
    rad2deg(rms(e2(idx,3))))
fprintf('5x poles q RMS after 50s = %.4f deg\n', ...
    rad2deg(rms(e5(idx,3))))
fprintf('10x poles q RMS after 50s = %.4f deg\n', ...
    rad2deg(rms(e10(idx,3))))

%LQG
%% LQG

% True measured outputs
y_true = x_true*C_obs';

% Add sensor noise
y_noisy = y_true + [noise_q noise_theta];

% Augmented state:
% X = [x;
%      xhat]

A_aug = [A        -B*K;
         Kf*C_obs  A-B*K-Kf*C_obs];

B_aug = [G              zeros(4,2);
         zeros(4,3)     Kf];

C_aug = eye(8);
D_aug = zeros(8,5);

sys_lqg = ss(A_aug,B_aug,C_aug,D_aug);

% Initial condition
X0 = [x0;
      xhat0];

% Inputs:
% [process noise  measurement noise]
U_lqg = [w_proc noise_q noise_theta];

[~,~,x_lqg] = lsim(sys_lqg,U_lqg,t,X0);

% Extract states
x         = x_lqg(:,1:4);
xhat_lqg  = x_lqg(:,5:8);

% Estimation error
e = x - xhat_lqg;

% RMS errors
fprintf('LQG RMS u error = %.3f m/s\n',rms(e(:,1)))
fprintf('LQG RMS w error = %.3f m/s\n',rms(e(:,2)))
fprintf('LQG RMS q error = %.3f deg/s\n',rad2deg(rms(e(:,3))))
fprintf('LQG RMS theta error = %.3f deg\n',rad2deg(rms(e(:,4))))

%figure

subplot(2,2,1)
plot(t,x(:,1),'LineWidth',1.5)
hold on
plot(t,xhat_lqg(:,1),'--','LineWidth',1.5)
grid on
title('Forward Velocity Perturbation, u')
xlabel('Time (s)')
ylabel('u (m/s)')
legend('True','Estimated')

subplot(2,2,2)
plot(t,x(:,2),'LineWidth',1.5)
hold on
plot(t,xhat_lqg(:,2),'--','LineWidth',1.5)
grid on
title('Vertical Velocity Perturbation, w')
xlabel('Time (s)')
ylabel('w (m/s)')
legend('True','Estimated')

subplot(2,2,3)
plot(t,rad2deg(x(:,3)),'LineWidth',1.5)
hold on
plot(t,rad2deg(xhat_lqg(:,3)),'--','LineWidth',1.5)
grid on
title('Pitch Rate, q')
xlabel('Time (s)')
ylabel('q (deg/s)')
legend('True','Estimated')

subplot(2,2,4)
plot(t,rad2deg(x(:,4)),'LineWidth',1.5)
hold on
plot(t,rad2deg(xhat_lqg(:,4)),'--','LineWidth',1.5)
grid on
title('Pitch Angle, \theta')
xlabel('Time (s)')
ylabel('\theta (deg)')
legend('True','Estimated')

sgtitle('LQG State Estimation Performance')

figure

subplot(2,2,1)
plot(t,e(:,1),'LineWidth',1.5)
grid on
title('u Error')
xlabel('Time (s)')
ylabel('Error (m/s)')

subplot(2,2,2)
plot(t,e(:,2),'LineWidth',1.5)
grid on
title('w Error')
xlabel('Time (s)')
ylabel('Error (m/s)')

subplot(2,2,3)
plot(t,rad2deg(e(:,3)),'LineWidth',1.5)
grid on
title('q Error')
xlabel('Time (s)')
ylabel('Error (deg/s)')

subplot(2,2,4)
plot(t,rad2deg(e(:,4)),'LineWidth',1.5)
grid on
title('\theta Error')
xlabel('Time (s)')
ylabel('Error (deg)')

sgtitle('LQG Estimation Errors')

figure

plot(t,rad2deg(x_cl_1(:,3)),'LineWidth',1.5)
hold on
plot(t,rad2deg(x(:,3)),'--','LineWidth',1.5)

grid on

title('LQR vs LQG Pitch Rate Response')
xlabel('Time (s)')
ylabel('q (deg/s)')

legend('Ideal LQR (Full State)','LQG (Estimated State)')

figure

plot(t,rad2deg(x_cl_1(:,4)),'LineWidth',1.5)
hold on
plot(t,rad2deg(x(:,4)),'--','LineWidth',1.5)

grid on

title('LQR vs LQG Pitch Angle Response')
xlabel('Time (s)')
ylabel('\theta (deg)')

legend('Ideal LQR (Full State)','LQG (Estimated State)')