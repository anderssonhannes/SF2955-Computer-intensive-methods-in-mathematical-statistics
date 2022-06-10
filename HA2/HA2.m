clear all, close all, clc

It = load('iran_infected.csv');         % Infected
Rt = load('iran_removed.csv');          % Removed
P = 84*10^6;                            % Population Iran

It = load('germany_infected.csv');      % Infected
Rt = load('germany_removed.csv');       % Removed
P = 83*10^6;                            % Population Germany

St = P-It-Rt;
td = length(It);

n = 10000;          
burn_in = 2000;
N = n + burn_in;

d = 4;

% Initial values
lamb_ = rand(1,d);
t_int = round(linspace(0,td,d+1));
t_ = t_int(2:d);
pir = 0.2;

% Parameters
sigma = 0.02;
M = 3;
a = 2;
b = 3;
b_i = 3*ones(length(lamb_),1);




lamb_res = zeros(length(lamb_),N);
lamb_res(:,1) = lamb_;
t_res = zeros(length(t_),N);
t_res(:,1) = t_;
pir_res = zeros(1,N);
pir_res(1) = full_pir(a,b,St,It);
for k = 1:(N-1)
    
    % Lambdas
    for i = 1:length(lamb_)
        cand = lamb_(i) + sigma*randn;
        if cand <= 0
           lamb_res(:,k+1) = lamb_; 
           continue
        end
        lamb_star = lamb_;
        lamb_star(i) = cand;
        
        alpha = min([0, full_lamb(i,lamb_star,t_,St,It,P,b_i) - full_lamb(i,lamb_,t_,St,It,P,b_i)]); % logarithm of alpha;
        
        if rand <= exp(alpha)
            lamb_res(:,k+1) = lamb_star;
            lamb_ = lamb_star;
        else
            lamb_res(:,k+1) = lamb_;
        end
    end
    
    % Break points 
    for i = 1:length(t_)
        cand = t_(i) + randi([-M M]);
        t_star = t_;
        t_star(i) = cand;
        % Special case when d-1 = 1
        if length(t_star) == 1
            cond = false;
        else
            cond = min(diff(t_star)) <= 0;
        end
        
        if cond || cand >= td || cand <= 0    % Check conditions for t_
            t_res(:,k+1) = t_;
            continue
        end
        
        alpha = min([0, full_t(i,t_star,lamb_,St,It,P) - full_t(i,t_,lamb_,St,It,P)]); % logarithm of alpha
        
        if rand <= exp(alpha)
            t_res(:,k+1) = t_star;
            t_ = t_star;
        else
            t_res(:,k+1) = t_;
        end
    end
    
    pir_res(k+1) = full_pir(a,b,St,It);
    
end 

mean_t = mean(t_res(:,burn_in:N),2)
mean_lamb = mean(lamb_res(:,burn_in:N),2)
mean_pir = mean(pir_res)

% Plots
for k = 1:length(lamb_)
    figure(k)
    plot(lamb_res(k,:))
end



% Example histogram
figure(5)
histogram(lamb_res(1,:))
title('Posterior distribution for \lambda_1')
xlabel('\lambda_1'); ylabel('Frequency')

% Problem 7 plot for germany data
if d == 4
    time = round(mean_t);
    figure(6)
    plot(It)
    hold on
    plot([9+9:9+14],It([9+9:9+14]),'r')
    hold on
    plot([16+9:16+14],It([16+9:16+14]),'g')
    hold on
    plot([23+9:23+14],It([23+9:23+14]),'c')
    hold on
    plot(time,It(time),'kd')
    legend('Infected individuals','Interval 1','Interval 2','Interval 3','Estimated break points')
end


% Defining full conditionals
function prob = full_lamb(j,lamb_,t_,St,It,P,b_i)
phi = 0.995;
td = length(It); 
a_i = 2*ones(length(lamb_),1);      % Recomended to use 2 for all i
lambda = lamb_(j);

    % Not including the last time in each interval to avoid duplicates
    if j == 1                   % If first lambda
        time = 1:t_(j)-1;
    elseif j == length(lamb_)    % If last lambda
        time = t_(end):td-1;
    else                        % Any other lambda
        time = t_(j-1):t_(j)-1;
    end
    dts = St(time,:) - St(time+1,:);
    kappa = (1/phi -1)*St(time).*(1-exp(-lambda*It(time)/P));
    P_dtI = gammaln(dts+kappa) - gammaln(dts+1) - gammaln(kappa) + kappa*log(1-phi) + dts*log(phi);
    
    prob = sum(a_i.*log(b_i) - gammaln(a_i) + (a_i-1)*log(lambda) - b_i*lambda) + sum(P_dtI);
end

function prob = full_t(n,t_,lamb_,St,It,P)
phi = 0.995;
td = length(It); 
t_int = [1 t_ td]; % Add end points and change 0 to 1 to avoid indexing errors. 
prob = 0;
for j=[n,n+1]   % Will always check the two adjacent intervals
    time = t_int(j):(t_int(j+1)-1);
    lambda = lamb_(j);

    dts = St(time,:) - St(time+1,:);
    kappa = (1/phi -1)*St(time).*(1-exp(-lambda*It(time)/P));
    P_dtI = gammaln(dts+kappa) - gammaln(dts+1) - gammaln(kappa) + kappa*log(1-phi) + dts*log(phi);
    
    prob = prob + sum(P_dtI);
    
end
prob = prob + log(sum((min(diff(t_int))>0)));
end



function prob = full_pir(a,b,St,It)
dts = -diff(St);
dti = -diff(It);
pir_alpha = sum(dts+dti) + a;
pir_beta = sum(It(1:end-1)-dts - dti) + b;
prob = betarnd(pir_alpha, pir_beta);
end

