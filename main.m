%% AA 273 Project
% Code to generate different types of noise and use different types of
% filters to estimate true trajectory measurements.

clear
clc

%% Simulation conditions.
t0 = 0; tf = 300; dt = 0.1;
tVec = t0:dt:tf;
x0 = zeros(6,1);
mass = 1;

%% Control inputs 
% At the moment its assumed that u = 0 is hover
numSteps = size(tVec,2);
u = [sin(0.05*tVec);cos(0.05*tVec);0.02*ones(1,numSteps)];

%% Process Noise Generation
Qgauss = diag([0,0.1,0,0.1,0,0.001]);
Qbi1 = diag([0,0.075,0,0.025,0,0.001]);
Qbi2 = diag([0,0.025,0,0.1,0,0.001]);
mubi = [0;-0.5;0;-0.5;0;-0.05];

new_noise = 1; % 1 to generate new noise. 0 to load same noise as prev run
generate_noise(Qgauss,Qbi1,Qbi2,mubi,numSteps,new_noise); 
load('noise')

%% True (noisy) Trajectories 
% Using double integrator \ddot{x} = u/mass
% State X = [x;xdot,y,ydot,z,zdot]
A = zeros(6,6); A(1,2) = 1; A(3,4) = 1; A(5,6) = 1;
B = zeros(6,3); B(2,1) = 1/mass; B(4,2) = 1/mass; B(6,3) = 1/mass;
C = zeros(3,6); C(1,1) = 1; C(2,3) = 1; C(3,5) = 1;

xNoNoise = zeros(6,numSteps); xNoNoise(:,1) = x0;
xGaussian = xNoNoise;
xUniform = xNoNoise;
xExp = xNoNoise;
xBrn = xNoNoise;
xCauchy = xNoNoise;
xBiMod = xNoNoise;
for i = 2:numSteps
    xNoNoise(:,i) = act_dyn(xNoNoise(:,i-1),u(:,i-1), dt, A, B);
    xGaussian(:,i) = act_dyn(xGaussian(:,i-1),u(:,i-1), dt, A, B) + gaussNoise(:,i); 
    xUniform(:,i) = act_dyn(xUniform(:,i-1),u(:,i-1), dt, A, B) + uniformNoise(:,i); 
    xExp(:,i) = act_dyn(xExp(:,i-1),u(:,i-1), dt, A, B) + expNoise(:,i); 
    xBrn(:,i) = act_dyn(xBrn(:,i-1),u(:,i-1), dt, A, B) + brnNoise(:,i); 
    xCauchy(:,i) = act_dyn(xCauchy(:,i-1),u(:,i-1), dt, A, B) + cauchNoise(:,i); 
    xBiMod(:,i) = act_dyn(xBiMod(:,i-1),u(:,i-1), dt, A, B) + biModalNoise(:,i); 
end

%% Measurent Noise Generation
% If y has smaller dimension than x, then need to redo this a little
Rgauss = 0.1*eye(6);
Qbi1 = diag([0,0.075,0,0.025,0,0.001]);
Qbi2 = diag([0,0.025,0,0.1,0,0.001]);
mubi = [0;-0.5;0;-0.5;0;-0.05];
new_noise = 1; % 1 to generate new noise. 0 to load process noise
generate_noise(Qgauss,Qbi1,Qbi2,mu,numSteps,new_noise);

%% Generate measurements. Use measurement model defined in Estimation
%Quadcopter Paper. Assumes measurement come from GPS and Magnetometer
yNoNoise = zeros(6,numSteps);
yGaussian = yNoNoise;
yUniform = yNoNoise;
yExp = yNoNoise;
yBrn = yNoNoise;
yCauchy = yNoNoise;
yBiMod = yNoNoise;

H = eye(6);

for i = 1:numSteps
    yNoNoise(:,i) = meas_model(H,xNoNoise(:,i));
    yGaussian(:,i) = meas_model(H,xGaussian(:,i)) + gaussNoise(:,i);
    yUniform(:,i) = meas_model(H,xUniform(:,i))+ uniformNoise(:,i);
    yExp(:,i) = meas_model(H,xExp(:,i)) + expNoise(:,i); 
    yBrn(:,i) = meas_model(H,xBrn(:,i)) + brnNoise(:,i);
    yCauchy(:,i) = meas_model(H,xCauchy(:,i)) + cauchNoise(:,i); 
    yBiMod(:,i) = meas_model(H,xBiMod(:,i)) + biModalNoise(:,i);
end

%% Add filter stuff
Cmean = [0;0;0;0;0;0]; % Default: Gauss, Uniform, Cauchy, Bimodal
%Cmean = [0.002;0.002;0.002;0.002;0.0004;0.0004]; % Exp
%Cmean = mubi; % Bimodal
Q = Qgauss; % Default: Gauss
R = Rgauss; % Default: Gauss

% [mu,sigma] = kalmanFilter(yGaussian,u,A,B,H,Q,R,Cmean,dt);
% plotting(mu,xGaussian,'Kalman Filter');

[mu,sigma] = ExtendedkalmanFilter(yGaussian,u,A,B,H,Q,R,dt);
plotting(mu,xGaussian,'Extended Kalman Filter');

numParticles = 1e3;
dim = 6;
[muPF] = particleFilter(yGaussian,u,numParticles,dim,Q,R,dt,A,B,C);
plotting(muPF,xGaussian,'Particle Filter');

%% Plot Trajectories
figure
plot3(xNoNoise(1,:),xNoNoise(3,:),xNoNoise(5,:),'DisplayName','No Noise')
hold on
plot3(xGaussian(1,:),xGaussian(3,:),xGaussian(5,:),'DisplayName','Gaussian noise')
plot3(xUniform(1,:),xUniform(3,:),xUniform(5,:),'DisplayName','Uniform noise')
plot3(xExp(1,:),xExp(3,:),xExp(5,:),'DisplayName','Exponetial noise')
plot3(xBrn(1,:),xBrn(3,:),xBrn(5,:),'DisplayName','Brownian noise')
plot3(xCauchy(1,:),xCauchy(3,:),xCauchy(5,:),'DisplayName','Cauchy noise')
plot3(xBiMod(1,:),xBiMod(3,:),xBiMod(5,:),'DisplayName','Bi-Model noise')
grid on
xlabel('X')
ylabel('Y')
zlabel('Z')
legend
title('Trajectories')