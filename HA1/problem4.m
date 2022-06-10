clear all, close all, clc

% Load datasets
load('stations.mat')        % pos_vec
load('RSSI-measurements.mat')

% Params
mu = 0;
sigma = .5;
dt = .5;
alpha = .6;

% Sizes
m = 500;
N = 10000;

% Number of recievers
s = 6;

% Random error function independent of n
W = @(n) normrnd(mu, sigma,[2,N]);
% Driving commands (states)
Z_state = [0 0;3.5 0;0 3.5; 0 -3.5; -3.5 0]';
% Transition matrix
P = 1/20*[16 1 1 1 1;1 16 1 1 1;1 1 16 1 1;1 1 1 16 1;1 1 1 1 16];
cs_mat = cumsum(P,2);   % Cumulative matrix of each state, rowwize

% Sub matrices
phi_tilde = [1 dt 0.5*dt^2; 0 1 dt;0 0 alpha];
psi_z_tilde = [0.5*dt^2 dt 0]';
psi_w_tilde = [0.5*dt^2 dt 1]';

% Matrices
phi = [phi_tilde zeros(3);zeros(3) phi_tilde];
psi_z = [psi_z_tilde zeros(3,1);zeros(3,1) psi_z_tilde];
psi_w = [psi_w_tilde zeros(3,1);zeros(3,1) psi_w_tilde];

% Defining kernel function
q = @(X,Z,W) phi*X + psi_z*Z + psi_w*W;


%% SISR
% Starting values
% Choose X0 vector according to its distribution for the two components
mu0 = zeros(6,1);
sigma0 = diag([500,5,5,200,5,5]);
X = mvnrnd(mu0,sigma0,N)';      % The col = x1, col 2 = x2
states = randi(5,[N,1]);     % Initially, states uniformly distributed


% Initial weights
w = p(Y(:,1),X([1 4],:)',pos_vec);


% Most common states
mc_states = zeros(1,m+1);
prob = zeros(1,5);
for k = 1:5
    prob(k) = sum(w(states == k));
end
[v,state_idx] = max(prob);
mc_states(1) = state_idx;


% Create tau matrix
tau = zeros(m+1,2);    % col 1 = tau1, col 2 = tau2
tau(1,:) = sum(w.*X([1 4],:)',1)/sum(w);     % initial approximation of tau 
for n = 2:m+1   % Main loop
    idx = randsample(N,N,true,w);
    states = states(idx);
    X = X(:,idx);
    Z = Z_state(:,states);
    X = q(X,Z,W(n));
    Xsub = X([1 4],:)';
    w = p(Y(:,n),Xsub,pos_vec);

    tau(n,:) = sum(w.*Xsub,1)/sum(w);
    
    prob = zeros(1,5);
    next_states = zeros(N,1);
    rv = rand(N,1);
    for k = 1:5
        prob(k) = sum(w(states == k)/sum(w));
        
        % Update driving commands
        vv = cs_mat(states,k);
        ind = rv <= vv;
        rv(ind) = nan;
        next_states(ind) = k;
    end
    [v,state_idx] = max(prob);
    mc_states(n) = state_idx;
    states = next_states;
    
end
figure(1)
subplot(1,2,1)
plot(tau(:,1),tau(:,2))
hold on
plot(pos_vec(1,:),pos_vec(2,:),'hr')
subtitle('SISR, N = 10000')
xlabel('X^1 [m]'); ylabel('X^2 [m]')
legend('Trajectory','Stations')

ml_commands = Z_state(:,mc_states);


% Plot
id1 = mc_states == 1;
id2 = mc_states == 2;
id3 = mc_states == 3;
id4 = mc_states == 4;
id5 = mc_states == 5;

figure(1)
subplot(1,2,2)
plot(tau(id1,1),tau(id1,2),'.b')
hold on
plot(tau(id2,1),tau(id2,2),'.m')
hold on
plot(tau(id3,1),tau(id3,2),'.r')
hold on
plot(tau(id4,1),tau(id4,2),'.g')
hold on
plot(tau(id5,1),tau(id5,2),'.k')
xlabel('X^1 [m]'); ylabel('X^2 [m]')
legend('No added velocity','East','North','South','West')
subtitle('Most common drvining commands')

disp('Effective sample size')
ESS_500 = 1/sum((w/sum(w)).^2)
