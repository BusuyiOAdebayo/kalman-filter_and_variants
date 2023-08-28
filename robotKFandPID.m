clear all; close all; clc
% Initialization
Ts = 0.1; % Sampling time
tinit = 0;
tfin = 100;
time = tinit:Ts:tfin;
N = length(time)-1;

% For plant
x0 = [0;0]; % Position, velocity
x(:,1) = x0;

% For Kalman Filter
sigma_a = 0.1;
sigma_y = 0.01;
mu = 0;
x_hat(:,1) = x0;
P_k{1} = 0*eye(2,2);
P_k_norm(1) = norm(P_k{1});

% Velocity Profile
Velocity(1:20) = 0; % m/s
Velocity(20:100) = 0.5;
Velocity(100:300) = 1.2;
Velocity(300:500) = 1.0;
Velocity(500:700) = 0;
Velocity(700:800) = 0.7;
Velocity(800:1001) = 0;

%% Plant Model
F = [1 Ts; 0 1];
G = [Ts^2/2; Ts];
H = [1 0]; % We have only positions measurements

% Check system controllability
if(rank(ctrb(F,G)) == 2)
    disp('System is controllable!')
end
% Check system observability
if(rank(obsv(F,H)) == 2)
    disp('System is observable!')
end

%% PID Controller
Kp = 0.8;
Ki = 1e-3;

[Kp,Ki] = mydialog;

%% Velocity estimation using Kalman filter
Q = G*G'*sigma_a^2;
R = sigma_y^2;
noise = sigma_y*randn(1,N+1)+mu; % Generate white noise

%% Simulation
% Error = zeros(1,length(N+1));
% Prop = zeros(1,length(N+1));
% Integral = zeros(1,length(N+1));
% Int = zeros(1,length(N+1));
Error = zeros(1,N+1);
Prop = zeros(1,N+1);
Integral = zeros(1,N+1);
%Int = zeros(1,N+1);

for t = 1:N
    % PID controller
    Error(t+1) = Velocity(t)-x_hat(2,t);

    Prop(t+1) = Error(t+1);
    Int(t+1) = Error(t+1)*Ts;
    Integral(t+1) = sum(Int);

    a(t+1) = Kp*Prop(t+1)+Ki*Integral(t+1);

    % Plant open loop simulation
    x(:,t+1) = F*x(:,t)+G*a(t+1);
    y(t+1) = H*x(:,t+1)+noise(t+1); % Add noisse to the measurements

    % Estimate the velocity using Kalman filter
    [x_up,P_k{t+1}] = Kalman_Filter(x_hat(:,t),P_k{t},y(t),F,H,Q,R);
    x_hat(:,t+1) = x_up;
end

function [x_up,P_up] = Kalman_Filter(x_k_1,P_k_1,y,F,H,Q,R)
x_k = F*x_k_1;
P_k = F*P_k_1*F'+Q; % Project the state covariance ahead
Kk = P_k*H'*inv(H*P_k*H'+R); % Kalman gain
x_up = x_k+Kk*(y-H*x_k); % Update measurement with measurement
P_up = P_k-Kk*H*P_k; % Update the error covariance
end

function [Kp,Ki] = mydialog
prompt = {'Proportional gain (Kp):', 'Integral gain (Ki):'};
dlgtitle = 'PID controller gains';
dims = [1 50;1 50];
definput = {'0.8', '0.001'};
values = inputdlg(prompt,dlgtitle,dims,definput);
%values = str2double(answer);

Kp = values(1);
Ki = values(2);
end