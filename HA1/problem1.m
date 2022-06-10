clear all, close all, clc

% Params
mu = 0;
sigma = .5;
dt = .5;
alpha = .6;

% Random error function independent of n
W = @(n) normrnd(mu, sigma,[2,1]);
% Driving commands (states)
Z = [0 0;3.5 0;0 3.5; 0 -3.5; -3.5 0]';
% Transition matrix
P = 1/20*[16 1 1 1 1;1 16 1 1 1;1 1 16 1 1;1 1 1 16 1;1 1 1 1 16];

% Sub matrices
phi_tilde = [1 dt 0.5*dt^2; 0 1 dt;0 0 alpha];
psi_z_tilde = [0.5*dt^2 dt 0]';
psi_w_tilde = [0.5*dt^2 dt 1]';

% Matrices
phi = [phi_tilde zeros(3);zeros(3) phi_tilde];
psi_z = [psi_z_tilde zeros(3,1);zeros(3,1) psi_z_tilde];
psi_w = [psi_w_tilde zeros(3,1);zeros(3,1) psi_w_tilde];

% Starting values
% Choose X0 vector according to its distribution
mu0 = zeros(6,1);
sigma0 = diag([500,5,5,200,5,5]);
X0 = mvnrnd(mu0,sigma0)';
% Choose random starting position unoformly
state = randi(5);
Z0 = Z(:,state);

m = 3000;
X_res = zeros(6,m+1);
X_res(:,1) = X0;
Z_prev = Z0;
for n = 2:m+1
    
    X_res(:,n) = phi*X_res(:,n-1) + psi_z*Z_prev + psi_w*W(n);
    
    p = rand;
    cs = cumsum(P(state,:));
    for k = 1:length(cs)
        pk = cs(k);
       if p<=pk
           state = k;
           break
       end
    end
    Z_prev = Z(:,state);
end

X1 = X_res(1,:);
X2 = X_res(4,:);
plot(X1,X2)
title('Trajectory')
xlabel('X^1 [m]'); ylabel('X^2 [m]')






