% Particle filtering is a statistical method that can be used to estimate the state of a system over time. 
% Here is a simple example of how particle filtering can be implemented in a bioreactor model using MATLAB:
clear all; close all; clc

% Define constants
Fin = 0.5; % feed flow rate (L/hr) % I changed from Q to Fin
Fout = 0.1; % flow rate (L/hr) % I changed from k to Fout
Y = 0.5; % yield coefficient (g substrate/g biomass)
mu = 0.2; % specific growth rate (1/hr)

% Define time points
t = 0:0.1:10; % time points (hr)

% Define initial conditions for variables (Set initial values for variables)
X0 = 0.2; % initial concentration of substrate (g/L)
S0 = 50; % initial concentration of substrate in feed (g/L)
V0 = 1; % initial volume of bioreactor (L)

% Initialize variables
% X = zeros(size(t)); % concentration of substrate (g/L)
% X(1) = X0;
% S = zeros(size(t)); % concentration of substrate in feed (g/L)
% S(1) = S0;
% V = zeros(size(t)); % volume of bioreactor (L)
% V(1) = V0;
% I changed to:
X = zeros(1,length(t)); % concentration of substrate (g/L)
X(1) = X0;
S = zeros(1,length(t)); % concentration of substrate in feed (g/L)
S(1) = S0;
V = zeros(1,length(t)); % volume of bioreactor (L)
V(1) = V0;

% Initialize particle filter variables
n_particles = 1000; % number of particles
% particles = zeros(n_particles, length(t)); % particles for concentration of substrate
% I changed to:
particles = zeros(n_particles, length(t)); % particles for concentration of substrate
weights = ones(n_particles, 1) / n_particles; % weights for particles

% Define initial conditions for particles (Set initial values for particles)
particles(:,1) = X0;

% Loop through time points
for i = 2:length(t)
    % Calculate new concentrations
    X(i) = (Fin*S(i-1) + X(i-1)*V(i-1)) / (Fout+V(i-1));
    S(i) = S0 - (Fin*S(i-1)) / Fout;
    V(i) = V0 + Fin - Fout;

    % Propagate particles (pass/move through particles at every point in time)
    for j = 1:n_particles
        particles(j,i) = (Fin*S(i-1) + particles(j,i-1)*V(i-1)) / (Fout+V(i-1));
    end

    % Calculate likelihood of each particle (put in a column vector for all the particles)
    likelihood = normpdf(particles(:,i), X(i), 0.1*X(i)); % Note here that: X(i) = mean and 0.1*X(i) = standard deviation

    % Calculate new weights for particles
    % weights = weights .* likelihood;
    % I changed to:
    weights = likelihood .* weights; % Kind of posterior = likelihood * prior
    weights = weights / sum(weights); % Normalize here

    % Resample particles
    % particles = randsample(particles, weights, n_particles);
    % I replaced with:
    % particles(:,i) = randsample(particles(:,i), weights, n_particles);
    particles(:,i) = randsample(particles(:,i),n_particles,true,weights); % generally, y = randsample(population,k,true,w)
end

% Plot results
% figure(1)
% plot(t, X, '-');
% figure(2)
% plot(t, S, '--');
% figure(3)
% plot(t, V, '-o');
plot(t, X, '-', t, mean(particles,1), '--');
xlabel('Time (hr)');
ylabel('Concentration of substrate (g/L)');
legend('True concentration of substrate', 'Estimated concentration of substrate');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % I changed to this:
% % Particle filtering is a statistical method that can be used to estimate the state of a system over time. 
% % Here is a simple example of how particle filtering can be implemented in a bioreactor model using MATLAB:
% clear all; close all; clc
% 
% % Define constants
% Fin = 0.5; % feed flow rate (L/hr) % I changed from Q to Fin
% Fout = 0.1; % flow rate (L/hr) % I changed from k to Fout
% Y = 0.5; % yield coefficient (g substrate/g biomass)
% mu = 0.2; % specific growth rate (1/hr)
% 
% % Define time points
% t = 0:0.1:10; % time points (hr)
% 
% % Define initial conditions for variables (Set initial values for variables)
% X0 = 0.2; % initial concentration of substrate (g/L)
% S0 = 50; % initial concentration of substrate in feed (g/L)
% V0 = 1; % initial volume of bioreactor (L)
% 
% % Initialize variables
% % X = zeros(size(t)); % concentration of substrate (g/L)
% % X(1) = X0;
% % S = zeros(size(t)); % concentration of substrate in feed (g/L)
% % S(1) = S0;
% % V = zeros(size(t)); % volume of bioreactor (L)
% % V(1) = V0;
% % I changed to:
% X = zeros(1,length(t)); % concentration of substrate (g/L)
% X(1) = X0;
% S = zeros(1,length(t)); % concentration of substrate in feed (g/L)
% S(1) = S0;
% V = zeros(1,length(t)); % volume of bioreactor (L)
% V(1) = V0;
% 
% % Initialize particle filter variables
% n_particles = 1000; % number of particles
% % particles = zeros(n_particles, length(t)); % particles for concentration of substrate
% % I changed to:
% particles = ones(n_particles, 1); % particles for concentration of substrate
% particles0 = X0*ones(n_particles, 1); % particles for concentration of substrate
% weights = ones(n_particles, 1) / n_particles; % weights for particles
% 
% % Define initial conditions for particles (Set initial values for particles)
% % particles(1) = X0;
% 
% % Loop through time points
% for i = 2:length(t)
%     % Calculate new concentrations
%     X(i) = (Fin*S(i-1) + X(i-1)*V(i-1)) / (Fout+V(i-1));
%     S(i) = S0 - (Fin*S(i-1)) / Fout;
%     V(i) = V0 + Fin - Fout;
% 
%     % Propagate particles (pass/move through particles at every point in time)
%     for j = 1:n_particles
%         particles(j) = (Fin*S(i-1) + particles0(j)*V(i-1)) / (Fout+V(i-1));
%     end
%     particles0 = particles(:);
% 
%     % Calculate likelihood of each particle (put in a column vector for all the particles)
%     likelihood = normpdf(particles(:), X(i), 0.1*X(i)); % Note here that: X(i) = mean and 0.1*X(i) = standard deviation
% 
%     % Calculate new weights for particles
%     % weights = weights .* likelihood;
%     % I changed to:
%     weights = likelihood .* weights; % Kind of posterior = likelihood * prior
%     weights = weights / sum(weights); % Normalize here
% 
%     % Resample particles
%     % particles = randsample(particles, weights, n_particles);
%     % I replaced with:
%     % particles(:,i) = randsample(particles(:,i), weights, n_particles);
%     particles(:) = randsample(particles(:),n_particles,true,weights); % generally, y = randsample(population,k,true,w)
% end
% 
% % Plot results
% % figure(1)
% % plot(t, X, '-');
% % figure(2)
% % plot(t, S, '--');
% % figure(3)
% % plot(t, V, '-o');
% plot(t, X, '-', t, mean(particles,1), '--');
% xlabel('Time (hr)');
% ylabel('Concentration of substrate (g/L)');
% legend('True concentration of substrate', 'Estimated concentration of substrate');