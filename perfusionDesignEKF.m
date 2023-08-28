% Define the state vector and the measurement vector
% x = [c(1) c(2) c(3) c(4) c(5) c(6) c(7) c(8) c(9) c(10) c(11) c(12)];
% y = [c(1) c(2) c(3)];
% x = [c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12];
% y = [c1 c2 c3];

% Define the non-linear state transition function
F = @(x) eye(12); % @(x) eye(9); 
% F = @(x) [1 0 0 0 0 0 0 0 0;
%     0 1 0 0 0 0 0 0 0;
%     0 0 1 0 0 0 0 0 0;
%     0 0 0 1 0 0 0 0 0;
%     0 0 0 0 1 0 0 0 0;
%     0 0 0 0 0 1 0 0 0;
%     0 0 0 0 0 0 1 0 0;
%     0 0 0 0 0 0 0 1 0;
%     0 0 0 0 0 0 0 0 1];

% Define the non-linear measurement function
H = @(x) [1 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0;
    1 1 1 0 0 0 0 0 0];

% Define the Jacobian of the state transition function
F_jac = @(x) eye(12); % @(x) eye(9);

% Define the Jacobian of the measurement function
H_jac = @(x) [1 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0;
    1 1 1 0 0 0 0 0 0];

% Define the process noise covariance matrix
q1 = 0.001;
q2 = 0.001;
q3 = 0.001;
q4 = 0.001;
q5 = 0.001;
q6 = 0.001;
q7 = 0.001;
q8 = 0.001;
q9 = 0.001;
q10 = 0.001;
q11 = 0.001;
q12 = 0.001;

Q = diag([q1 q2 q3 q4 q5 q6 q7 q8 q9 q10 q11 q12]);

% Define the measurement noise covariance matrix
r1 = 0.05;
r2 = 0.05;
r3 = 0.05;

R = diag([r1 r2 r3]);

% Initialize the state estimate and the error covariance matrix
c1initial = 1.75;
c2initial = 1e-2;
c3initial = c1initial+c2initial;
c4initial = 4.94;
c5initial = 0.00;
c6initial = 0.25*c4initial;
c7initial = 0.51;
c8initial = 0.96;
c9initial = 0.00;
c10initial = 2500;
c11initial = 2.46;
c12initial = (c1initial/c3initial)*100;
xhat = [c1initial c2initial c3initial c4initial c5initial c6initial c7initial c8initial c9initial c10initial c11initial c12initial]; % c13initial c14initial c15initial c16initial];
p1 = 0.01;
p2 = 0.01;
p3 = 0.01;
p4 = 0.01;
p5 = 0.01;
p6 = 0.01;
p7 = 0.01;
p8 = 0.01;
p9 = 0.01;
p10 = 0.01;
p11 = 0.01;
p12 = 0.01;
P = diag([p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12]);

% Perform Extended Kalman filtering
t = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30];
for i = 1:length(t)
    % Predict the state estimate and the error covariance matrix
    xhat_ = F(xhat);
    P_ = F_jac(xhat)'*P*F_jac(xhat) + Q;

    % Update the state estimate and the error covariance matrix
    K = P_*H'/(H*P_*H' + R);
    xhat = xhat_ + K*(y(i,:)' - H(xhat_));
    P = P_ - K*H*P_;

    % Store the estimated state variables
    % c1_est(i) = xhat(10);
    % c2_est(i) = xhat(11);
    % c3_est(i) = xhat(12);
    c4_est(i) = xhat(1);
    c5_est(i) = xhat(2);
    c6_est(i) = xhat(3);
    c7_est(i) = xhat(4);
    c8_est(i) = xhat(5);
    c9_est(i) = xhat(6);
    c10_est(i) = xhat(7);
    c11_est(i) = xhat(8);
    c12_est(i) = xhat(9);
end

% Plot the estimated state variables
figure;
subplot(4,3,1);
plot(t,c4_est,'b');
xlabel('t, [days]');
ylabel('[GLC] est, [g/L]');
legend('[GLC] est vs t');

subplot(4,3,2);
plot(t,c5_est,'g');
xlabel('t, [days]');
ylabel('[GLN] est, [mM]');
legend('[GLN] est vs t');

subplot(4,3,3);
plot(t,c6_est,'r');
xlabel('t, [days]');
ylabel('[AMI] est, [mM]');
legend('[AMI] est vs t');

subplot(4,3,4);
plot(t,c7_est,'b');
xlabel('t, [days]');
ylabel('[LAC] est, [g/L]');
legend('[LAC] est vs t');

subplot(4,3,5);
plot(t,c8_est,'g');
xlabel('t, [days]');
ylabel('[AMM] est, [mM]');
legend('[AMM] est vs t');

subplot(4,3,6);
plot(t,c9_est,'r');
xlabel('t, [days]');
ylabel('Titer est, [g/L]');
legend('Titer est vs t');

subplot(4,3,7);
plot(t,c10_est,'b');
xlabel('t, [days]');
ylabel('V est, [mL]');
legend('V est vs t');

subplot(4,3,8);
plot(t,c11_est,'g');
xlabel('t, [days]');
ylabel('[GLU] est, [mM]');
legend('[GLU] est vs t');

subplot(4,3,9);
plot(t,c12_est,'r');
xlabel('t, [days]');
ylabel('Viability est, [%]');
legend('Viability est vs t');

% c1initial = 1.75;
% c2initial = 1e-2;
% c3initial = c1initial+c2initial;
% c4initial = 4.94;
% c5initial = 0.00;
% c6initial = 0.25*c4initial;
% c7initial = 0.51;
% c8initial = 0.96;
% c9initial = 0.00;
% c10initial = 2500;
% c11initial = 2.46;
% c12initial = (c1initial/c3initial)*100;
tspan = [0 30];
cinitial = [c1initial c2initial c3initial c4initial c5initial c6initial c7initial c8initial c9initial c10initial c11initial c12initial];
[t,c] = ode45(@perfusionDesignCase1_1, tspan, cinitial);
function dcdt = perfusionBioreactor(~,c,MV,theta,Xv_target,cfeed,V,CSPR,PR,proportioninBleed,proportioninHarvest)
        dcdt = zeros(8,1);

        Fbleed = MV(1);
        cfeedGlc = MV(2);
        cfeedGlu = MV(3);
        %         mumax  = theta(1);
        %         Ksglc = theta(2);
        %         Ksglu = theta(3);
        %         Kilac  = theta(4);
        %         Kiamm = theta(5);
        %         qXd  = theta(6);
        %         qglc  = theta(7);
        %         qglu = theta(8);
        %         qgln1 = theta(9);
        %         qgln2 = theta(10);
        %         qgln3 = theta(11);
        %         qlac1 = theta(12);
        %         qlac2  = theta(13);
        %         qamm = theta(14);
        %         qmAb  = theta(15);

        mumax  = theta(1);
        qXd  = theta(2);
        qglc  = theta(3);
        qglu = theta(4);
        qgln1 = theta(5);
        qgln2 = theta(6);
        qgln3 = theta(7);
        qlac1 = theta(8);
        qlac2  = theta(9);
        qamm = theta(10);
        qmAb  = theta(11);

        mu = mumax;%(c(3)/(Ksglc+c(3)))*(c(4)/(Ksglu+c(4)))*(Kmlac/(Kmlac+c(5)))*(Kmamm/(Kmamm+c(5)));
        % mu = mumax*(c(4)/(Ksglc+c(4)))*(c(5)/(Ksglu+c(5)));
        % mu = mumax*(c(4)/(Ksglc+c(4)))*(c(5)/(Ksglu+c(5)))*(c(7)/(Ksami+c(7)))*(Kmlac/(Kmlac+c(8)))*(Kmamm/(Kmamm+c(9)));%*(c(6)/(Ksgln+c(6)));
        % mu = mumax*(c(4)/(Ksglc+c(4)))*(c(5)/(Ksglu+c(5)))*(c(6)/(Ksgln+c(6)))*(c(7)/(Ksami+c(7)))*(Kmlac/(Kmlac+c(8)))*(Kmamm/(Kmamm+c(9)));
        % Control option A:
        %         if c(1) <= Xv_target
        %             Fmedia = CSPR*c(1)*V; % [mL/d]
        %         else
        %             Fmedia = PR*V;
        %         end
        %         Fharvest = Fmedia-Fbleed;

        % Control option B:
        if c(1) <= Xv_target
            Fharvest = CSPR*c(1)*V; % [mL/d]
        else
            Fharvest = PR*V;
        end
        Fmedia = Fharvest+Fbleed;

        r = [(mu-qXd)*c(1);qXd*c(1);-qglc*c(1);-qglu*c(1);qgln1*c(4)-qgln2*c(1)-qgln3*c(5);qlac1*c(3)-qlac2*c(1);qamm*c(5);qmAb*c(1)];

        % Environmental disturances are neglected for now!
        dcdt(1) = (Fmedia/V)*cfeed(1)-proportioninBleed(1)*(Fbleed/V)*c(1)-proportioninHarvest(1)*(Fharvest/V)*c(1)+r(1); % Viable cells
        dcdt(2) = (Fmedia/V)*cfeed(2)-proportioninBleed(2)*(Fbleed/V)*c(2)-proportioninHarvest(2)*(Fharvest/V)*c(2)+r(2); % Dead cells
        dcdt(3) = (Fmedia/V)*cfeedGlc-proportioninBleed(3)*(Fbleed/V)*c(3)-proportioninHarvest(3)*(Fharvest/V)*c(3)+r(3); % Glucose
        dcdt(4) = (Fmedia/V)*cfeedGlu-proportioninBleed(4)*(Fbleed/V)*c(4)-proportioninHarvest(4)*(Fharvest/V)*c(4)+r(4); % Glutamate
        dcdt(5) = (Fmedia/V)*cfeed(5)-proportioninBleed(5)*(Fbleed/V)*c(5)-proportioninHarvest(5)*(Fharvest/V)*c(5)+r(5); % Glutamine
        dcdt(6) = (Fmedia/V)*cfeed(6)-proportioninBleed(6)*(Fbleed/V)*c(6)-proportioninHarvest(6)*(Fharvest/V)*c(6)+r(6); % Lactate
        dcdt(7) = (Fmedia/V)*cfeed(7)-proportioninBleed(7)*(Fbleed/V)*c(7)-proportioninHarvest(7)*(Fharvest/V)*c(7)+r(7); % Ammonia
        dcdt(8) = (Fmedia/V)*cfeed(8)-proportioninBleed(8)*(Fbleed/V)*c(8)-proportioninHarvest(8)*(Fharvest/V)*c(8)+r(8); % mAb product
end