%Initialization
close all; clear all; clc;
sympref('AbbreviateOutput',false);
%Variable and parameter definition
syms X G E t % X G E P t
syms Y_gx Y_ge Y_ex mu1 mu2 mue1 K_M_G K_M_E
%Variables/Parameters
%initX = [2.5; 6; 0.2; 0.15; 0.08]; %Initial state (Biomass, Glucose, Ethanol)
initX = [1.75;0.00;1.75;0.00;4.94;2.46;0.00;5.00;0.00;0.00;0.00]; %Initial state (Biomass, Glucose, Ethanol)
initP = diag([0.1,0.02,0.02,0.2,0.02,0.1,0.02,0.02,0.2,0.02,0.02]); % Initial process estimation covariance matrix
% I changed from [initX, initP(:)] to [initX; initP(:)] due to this error: Error using horzcat, Dimensions of arrays being concatenated are not consistent.
init = [initX; initP(:)]; % Combined initial value vector for the odesolver
H = [1,1,1,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]; % Observation matrix
Q = diag([0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001]); % Process noise covariance matrix;
R = 0.05; % Measurement noise covariance matrix
K1 = 0.1; % Monod constant glucose
K2 = 0.1; % Monod constant ethanol

% Estimated parameter values
Ygx = 0.15; % Yield glucose -> biomass
Yge = 0.34; % Yield glucose -> ethanol
Yex = 0.43; % Yield ethanol -> biomass

% Process model
% Monod terms
mue1 = mu1*G/(G+K_M_G);
mue2 = mu2*E/(E+K_M_E)*(1-mue1/mu1);
% mue1_t = eval(subs(mue1, {G, mu1}, {G_t, mu1_t}));
% mue2_t = eval(subs(mue2, {G, mu2}, {G_t, mu2_t}));

% Model OD
dS = sym(X*[ ...
    (mue1+mue2); ...                    % Biomass
    -mue1/Y_gx; ...                     % Glucose
    (mue1/Y_gx*Y_ge-mue2/Y_ex); ...     % Ethanol
    0; ...                              % mue1
    0; ...                              % mue2
    ]);

% Jacobian of Model with respect to state variables
F = jacobian(dS, [X, G, E, mu1, mu2]);

% P matrix
P = sym('P', [5,5]);
dP = F*P+P*F'+Q;
% Simulation/State prediction and filtering
% Replace all symbolic parameters with their respective numeric values
F = subs(F, [Y_gx Y_ge Y_ex K_M_G K_M_E], [Ygx Yge Yex K1 K2]);
dS = subs(dS, [Y_gx Y_ge Y_ex K_M_G K_M_E], [Ygx Yge Yex K1 K2]);
dP = subs(dP, [Y_gx Y_ge Y_ex K_M_G K_M_E], [Ygx Yge Yex K1 K2]);

% Assemble all differential equations into a vector of 12 elements (3x state, 9x P)
OdeSys = matlabFunction([dS(:);dP(:)], 'Vars', {t, [X; G; E; mu1; mu2; P(:)]});

% Load measurement values from file:
%load BC2_eth_pred.mat % Ethanol sensor measurements (ME)
ME = [0.2;0.266666667;0.333333333;0.4;0.466666667;0.533333333;0.6;0.666666667;0.733333333;0.8;0.866666667;0.933333333;1;1.066666667;1.133333333;1.2;1.266666667;1.333333333;1.4;1.466666667;1.533333333;1.6;1.666666667;1.733333333;1.8;1.7625;1.725;1.6875;1.65;1.6125;1.575;1.5375;1.5;1.4625;1.425;1.3875;1.35;1.3125;1.275;1.2375;1.2;1.1625;1.125;1.0875;1.05;1.0125;0.975;0.9375;0.9;0.8625;0.825;0.7875;0.75;0.7125;0.675;0.6375;0.6;0.5625;0.525;0.4875;0.45;0.4125;0.375;0.3375;0.3;0.2625;0.225;0.1875;0.15;0.1125;0.075;0.0375;0];
%load BC2.mat % Offline values for biomass, glucose and ethanol (M)
time = [0.75;0.9;1;1.1;1.3;1.5;1.7;1.9;2.1;2.2;2.3;2.4;2.5;5.5;6];
M = zeros(length(time),3);
M(:,1) = [2.6;2.9;3;3.1;3.15;3.16;3.17;3.4;3.45;3.55;3.65;3.68;3.75;4.25;4.5];
M(:,2) = [3.65;3;2.4;2.1;1.75;1.2;0.85;0.4;0.2;0.1;0.05;0.025;0.005;0;0];
M(:,3) = [0.85;1;1.1;1.3;1.4;1.5;1.65;1.8;1.9;1.85;1.75;1.7;1.65;0.3;0];
% mat = Microsoft Access Table
% Simulate the process from one ethanol gas measurement time to the next:
t0 = 0;
SimTime = []; 
%SimTime = zeros(0,1);
timeE = [0.083333333;0.166666667;0.25;0.333333333;0.416666667;0.5;0.583333333;0.666666667;0.75;0.833333333;0.916666667;1;1.083333333;1.166666667;1.25;1.333333333;1.416666667;1.5;1.583333333;1.666666667;1.75;1.833333333;1.916666667;2;2.083333333;2.166666667;2.25;2.333333333;2.416666667;2.5;2.583333333;2.666666667;2.75;2.833333333;2.916666667;3;3.083333333;3.166666667;3.25;3.333333333;3.416666667;3.5;3.583333333;3.666666667;3.75;3.833333333;3.916666667;4;4.083333333;4.166666667;4.25;4.333333333;4.416666667;4.5;4.583333333;4.666666667;4.75;4.833333333;4.916666667;5;5.083333333;5.166666667;5.25;5.333333333;5.416666667;5.5;5.583333333;5.666666667;5.75;5.833333333;5.916666667;6];
% ME = [0.2;0.266666667;0.333333333;0.4;0.466666667;0.533333333;0.6;0.666666667;0.733333333;0.8;0.866666667;0.933333333;1;1.066666667;1.133333333;1.2;1.266666667;1.333333333;1.4;1.466666667;1.533333333;1.6;1.666666667;1.733333333;1.8;1.7625;1.725;1.6875;1.65;1.6125;1.575;1.5375;1.5;1.4625;1.425;1.3875;1.35;1.3125;1.275;1.2375;1.2;1.1625;1.125;1.0875;1.05;1.0125;0.975;0.9375;0.9;0.8625;0.825;0.7875;0.75;0.7125;0.675;0.6375;0.6;0.5625;0.525;0.4875;0.45;0.4125;0.375;0.3375;0.3;0.2625;0.225;0.1875;0.15;0.1125;0.075;0.0375;0];
MF = zeros(0,5); % Store filtered states in these variables
SimState = zeros(0,5);
for i = 1:numel(timeE)
    tspan = [t0 timeE(i)];
    [Time,State] = ode45(OdeSys, tspan, init); % Simulate/solve model

    PS = State(end,1:5)'; % Predicted state
    MS = ME(i); % Measured state
    P = reshape(State(end,6:end),5,5); % Process covariance matrix
    K = P*H'/(H*P*H'+R); % Kalmann gain matrix
    FS = PS+K*(MS-PS(3)); % Filtered state
    PF = P-K*H*P; % Filtered process covariance matrix
    init = [FS; PF(:)]; % New initial condition
    t0 = timeE(i); % New starting time for next iteration

    % Save intermediate states for plotting
    MF = [MF;FS'];
    %State(end,1:3) = NaN;
    SimState = [SimState; State(:,1:5)];
    SimTime = [SimTime;Time];
end

% mue1_t = eval(subs(mue1, {G, E, mu1}, {MF(:,2), MF(:,3), MF(:,4)}));
% mue2_t = eval(subs(mue2, {G, E, mu2}, {MF(:,2), MF(:,3), MF(:,5)}));

% Results
% Plot the results in a presentable figure and save file to disk
f = figure("Position", [0,0,1600,640]);
subplot(1,2,1);
h = plot([0;time], [initX(1:3)'; M],'.','MarkerSize',20); % Plot measurements
set(h,{'color'},{'r';'g';'b'});
hold on;
h = plot(SimTime,SimState(:,1:3)); % Plot simulated values
set(h,{'color'},{'r';'g';'b'});
plot([0;timeE],ME,'+b','MarkerSize',8); % Plot ethanol gas sensor values
hold off;
ax = gca;
ax.FontSize = 14;
ax.FontName = 'Times';
ax.Position = [0.05 0.1 0.4 0.85];
%ax.Position = [0.05 0.1 0.5 1.0];
ax.ActivePositionProperty = 'outerposition';
ax.GridLineStyle = ':';
ax.GridAlpha = 0.7;
xlabel('time $/h$','interpreter','Latex','FontSize',16);
ylabel('concentration $/\frac{g}{L}$','interpreter','Latex','FontSize',16); ylim([0 8]);
grid on; box off; grid(gca,'minor');
legend('Biomass Offline','Glucose Offline','Ethanol Offline','Biomass Kalman','Glucose Kalman','Ethanol Kalman','Ethanol Online (Gas Sensor)','interpreter','Latex','FontSize',12,'Color',[0.9 0.91 0.9]);
subplot(1,2,2);
h = plot(SimTime,SimState(:,4:5),timeE,mue1_t,timeE,mue2_t); % Plot mu values over time
ax = gca;
ax.FontSize = 14;
ax.FontName = 'Times';
ax.Position = [0.55 0.1 0.4 0.85];
%ax.Position = [0.55 0.1 0.5 1.0];
ax.ActivePositionProperty = 'outerposition';
ax.GridLineStyle = ':';
ax.GridAlpha = 0.7;
ytickformat('%.2f')
set(h,{'color'}, {'r';'k'});
xlabel('time $/h$','interpreter','Latex','FontSize',16);
ylabel('$\mu$ value $/\frac{1}{h}$','interpreter','Latex','FontSize',16);
grid on; box off; grid(gca,'minor');
legend('$\mu_{max,G}$','$\mu_{max,E}$','interpreter','Latex','FontSize',12,'Color',[0.9 0.91 0.9]);
annotation('arrow',[0.55 0.97],[0.1 0.1])
annotation('arrow',[0.05 0.47],[0.1 0.1])
annotation('arrow',[0.05 0.05],[0.1 0.98])
annotation('arrow',[0.55 0.55],[0.1 0.98])
saveas(f,'KalmanPred.svg','svg'); % Save a copy of figure to file

% Calculate Errors
SimTime = SimTime+((1:numel(SimTime))*1e-10)';
SimValues = interp1(SimTime,SimState(:,1:3),time);
SSE = sum((SimValues-M).^2);
RMSE = sqrt(SSE/numel(time)); % Actually RMSEP
SQT = sum((M-mean(M)).^2);
RSq = 1-SSE./SQT(1:3);
Performance_Table = table('Size',[3,3], 'VariableTypes',{'double','double','double'},'VariableNames',{'Biomass','Glucose','Ethanol'},'RowNames',{'RMSEP (g/L)','Error (%)','R^2'});
Performance_Table(1,:) = num2cell(RMSE);
Performance_Table(2,:) = num2cell(SSE./(max(M)-min(M))*100);
Performance_Table(3,:) = num2cell(RSq);