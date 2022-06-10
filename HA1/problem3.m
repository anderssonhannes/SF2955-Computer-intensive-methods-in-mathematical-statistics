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


%% SIS
% Starting values
% Choose X0 vector according to its distribution for the two components
mu0 = zeros(6,1);
sigma0 = diag([500,5,5,200,5,5]);
X = mvnrnd(mu0,sigma0,N)';      % The col = x1, col 2 = x2
states = randi(5,[N,1]);     % Initially, states uniformly distributed

% Initial weights
w = p(Y(:,1),X([1 4],:)',pos_vec);
w_plot1 = w;
% Create tau matrix
tau = zeros(m+1,2);    % col 1 = tau1, col 2 = tau2
tau(1,:) = sum(w.*X([1 4],:)',1)/sum(w);     % initial approximation of tau 
for n = 2:m+1   % Main loop
    Z = Z_state(:,states);
    X = q(X,Z,W(n));
    Xsub = X([1 4],:)';
    w = p(Y(:,n),Xsub,pos_vec).*w;

    if n == 15
       w_plot2 = w; 
    elseif n == 70
        w_plot3 = w;
    end
    tau(n,:) = sum(w.*Xsub,1)/sum(w);
    
    next_states = zeros(N,1);
    rv = rand(N,1);
    for k = 1:5
        vv = cs_mat(states,k);
        ind = rv <= vv;
        rv(ind) = nan;
        next_states(ind) = k;
    end
    states = next_states;
end
figure(1)
subplot(1,2,1)
plot(tau(:,1),tau(:,2))
hold on
plot(pos_vec(1,:),pos_vec(2,:),'hr')
title('SIS, N = 10000')
xlabel('X^1 [m]'); ylabel('X^2 [m]')
legend('Trajectory','Stations')

figure(2)
sgtitle('Histogram of SIS')
subplot(1,3,1)
hist(w_plot1)
title('m = 0')
subplot(1,3,2)
hist(w_plot2)
title('m = 15')
subplot(1,3,3)
hist(w_plot3)
title('m = 70')

disp('Effective sample size SIS')
ESS_0  = 1/sum((w_plot1/sum(w_plot1)).^2)
ESS_15 = 1/sum((w_plot2/sum(w_plot2)).^2)
ESS_70 = 1/sum((w_plot3/sum(w_plot3)).^2)


%% SIS with normalized log weights
% Starting values
% Choose X0 vector according to its distribution for the two components
mu0 = zeros(6,1);
sigma0 = diag([500,5,5,200,5,5]);
X = mvnrnd(mu0,sigma0,N)';      % The col = x1, col 2 = x2
states = randi(5,[N,1]);     % Initially, states uniformly distributed



% Initial weights
lw = log(p(Y(:,1),X([1 4],:)',pos_vec));
w_ = exp(lw-max(lw));
w = w_/sum(w_);
w_plot1 = w;
% Create tau matrix
tau = zeros(m+1,2);    % col 1 = tau1, col 2 = tau2
tau(1,:) = sum(w.*X([1 4],:)',1)/sum(w);     % initial approximation of tau 
for n = 2:m+1   % Main loop
    Z = Z_state(:,states);
    X = q(X,Z,W(n));
    Xsub = X([1 4],:)';
    lw = log(p(Y(:,n),Xsub,pos_vec)) + lw;
    L = max(lw);
    w_ = exp(lw-L);
    w = w_/sum(w_);
    if n == 15
       w_plot2 = w; 
    elseif n == 70
        w_plot3 = w;
    end
    tau(n,:) = sum(w.*Xsub,1)/sum(w);
    
    next_states = zeros(N,1);
    rv = rand(N,1);
    for k = 1:5
        vv = cs_mat(states,k);
        ind = rv <= vv;
        rv(ind) = nan;
        next_states(ind) = k;
    end
    states = next_states;
end
figure(1)
subplot(1,2,2)
plot(tau(:,1),tau(:,2))
hold on
plot(pos_vec(1,:),pos_vec(2,:),'hr')
title('SIS with normalized log weights, N=10000')
xlabel('X^1 [m]'); ylabel('X^2 [m]')
legend('Trajectory','Stations')

figure(3)
sgtitle('Histogram of SIS with log weights')
subplot(1,3,1)
hist(w_plot1)
title('m = 0')
subplot(1,3,2)
hist(w_plot2)
title('m = 15')
subplot(1,3,3)
hist(w_plot3)
title('m = 70')

disp('Effective sample size SIS normalized log weights')
ESS_0 = 1/sum((w_plot1/sum(w_plot1)).^2)
ESS_15 = 1/sum((w_plot2/sum(w_plot2)).^2)
ESS_70 = 1/sum((w_plot3/sum(w_plot3)).^2)




