clear all, close all, clc

% Load datasets
load('stations.mat')        % pos_vec
load('RSSI-measurements-unknown-sigma.mat')

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
Np = 20;     % Number of sigmas tested
sig0 = 2.1306;    % Must be larger than 0.3
sig1 = 2.2521;

sig_vec = linspace(sig0,sig1,Np);     % USe this for iteration of sigmas
sig_vec = [2.1872];                   % Use this to plot optimal sigma
ml_vec = zeros(1,length(sig_vec));
for np = 1:length(sig_vec)
    sig = sig_vec(np);
    % Starting values
    % Choose X0 vector according to its distribution for the two components
    mu0 = zeros(6,1);
    sigma0 = diag([500,5,5,200,5,5]);
    X = mvnrnd(mu0,sigma0,N)';      % The col = x1, col 2 = x2
    states = randi(5,[N,1]);     % Initially, states uniformly distributed
    
    w_sum = zeros(m+1,1);
    % Initial weights
    w = p5(Y(:,1),X([1 4],:)',pos_vec,sig);
    w_sum(1) = sum(w)/N;
    
    % Create tau matrix
    tau = zeros(m+1,2);    % col 1 = tau1, col 2 = tau2
    tau(1,:) = sum(w.*X([1 4],:)',1)/sum(w);     % initial approximation of tau
    for n = 2:m+1   % Main loop
        idx = randsample(N,N,true,w);
        states = states(idx);
        Z = Z_state(:,states);
        X = X(:,idx);
        X = q(X,Z,W(n));
        Xsub = X([1 4],:)';
        w = p5(Y(:,n),Xsub,pos_vec,sig);
        w_sum(n) = sum(w)/N;
        tau(n,:) = sum(w.*Xsub,1)/sum(w);
        
        next_states = zeros(N,1);
        rv = rand(N,1);
        for k = 1:5
            
            % Update driving commands
            vv = cs_mat(states,k);
            ind = rv <= vv;
            rv(ind) = nan;
            next_states(ind) = k;
        end
        states = next_states;
    end
    % The log likelihood estimation
    ml_vec(np) = sum(log(w_sum));
end
[ml_est, ml_ind] = max(ml_vec);
ml_sig = sig_vec(ml_ind)

if sig == 2.1872
    figure(1)
    plot(tau(:,1),tau(:,2))
    hold on
    plot(pos_vec(1,:),pos_vec(2,:),'hr')
    title('SISR for \varsigma = 2.1872, N=10000')
    xlabel('X^1 [m]'); ylabel('X^2 [m]')
    legend('Trajectory','Stations')
end