clear all; close all; clc

% Define the system model and parameters
A = [0.95 0; 0 -0.15];
B = [1 0; 0 1];
C = [1 0; 0 1];
D = [0 0; 0 0];
x0 = [1; 0.5];
u0 = [0.5; 0.5];
N = 10;
Ts = 1;
Q = diag([1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ]);%diag([1 1]);
R = diag([0.01 0.01]);%diag([0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01]);%diag([0.01 0.01]);

% Define the initial state and input estimate
xpred = [1.1; 0.4]; % I changed from xhat to xpred
upred = [0.4; 0.4]; % I changed from uhat to upred

% Define the measurement vector
y = [1.2 0.35; % 1st sampling instant
    1.1 0.4; % 2nd sampling instant
    1.3 0.3]; % 3rd sampling instant

% Implement the MHE algorithm
for i = 1:length(y) % Every sampling instant
    % Define the prediction matrices
    G = eye(2); % [1,0;0,1], This is for xpred
    F = eye(2); % [1,0;0,1], This is for upred
    for j = 1:N-1
        G = [G; A^(j-1)*B]; % This is for xpred
        F = [F; A^j]; % This is for upred
    end
    H = [C;zeros(2*(N-1),2)]; % I changed from [C zeros(2*(N-1),2)]; to [C; zeros(2*(N-1),2)];

    % Define the cost function (from P and q) and constraints (from lb and ub)
    P = G'*Q*G + R;
    q = 2*G'*Q*(F*xpred - H'*y(i,:)');
    lb = [0; 0];
    ub = [2; 2];

    % Solve the optimization problem using quadprog

    % quadprog: Quadratic programming
    % Syntax:
    % x = quadprog(H,f)
    % x = quadprog(H,f,A,b)
    % x = quadprog(H,f,A,b,Aeq,beq)
    % x = quadprog(H,f,A,b,Aeq,beq,lb,ub)
    % x = quadprog(H,f,A,b,Aeq,beq,lb,ub,x0)
    % x = quadprog(H,f,A,b,Aeq,beq,lb,ub,x0,options)
    % x = quadprog(problem)
    % [x,fval] = quadprog(___)
    % [x,fval,exitflag,output] = quadprog(___)
    % [x,fval,exitflag,output,lambda] = quadprog(___)
    % [wsout,fval,exitflag,output,lambda] = quadprog(H,f,A,b,Aeq,beq,lb,ub,ws)

    options = optimset('Display', 'off');
    [u,~,exitflag] = quadprog(P, q, [], [], [], [], lb, ub, upred, options);

    % Check for successful optimization
    if exitflag == 1
        upred = u;
        xpred = A*xpred + B*upred;
    end
end

% Plot the estimated states and inputs
t = (0:length(y)-1)*Ts;
subplot(2,1,1)
plot(t, xpred(1,:), 'b', t, y(:,1)', 'r--')
ylabel('Cell concentration')
legend('Estimated', 'Measured')
subplot(2,1,2)
plot(t, xpred(2,:), 'b', t, y(:,2)', 'r--')
ylabel('Dissolved oxygen concentration')
legend('Estimated', 'Measured')
xlabel('Time (s)')