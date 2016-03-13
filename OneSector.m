%     MATLAB program to solve the one-sector model
%     from the article "Solving Nonlinear Dynamic
%     Stochastic Models: An Algorithm Iterating on
%     Value Function by Simulations".
%          
%     1. OneSector.m
%     2. objective.m
%     3. shock.m
%
%     April 5, 2002
%
%     -------------------------------------------------------------------
clear all; 
CPU0 = cputime;

% 1. Initialize the model parameters
% ----------------------------------
alpha   = 0.33;         % Capital share 
delta   = .95;          % Discount factor
d       = .02;          % Depreciation rate
T       = 1000;        % Length of simulation

% 2. The steady state of k and c 
% -----------------------------------------------
ks = ( (1-delta+delta*d) / (alpha*delta) )^(1/(alpha-1) );
cs = ks^alpha - d*ks;												  

% 3. Allocate memory for the simulated series
% -------------------------------------------
k    = zeros(T+1,1)+ks;  % Capital k 
c    = zeros(T,1);       % Consumption
w    = zeros(T-1,1);     % Value function

% 4. Draw the shock
% -------------------------------------------
load shock05;           % loading rho, sigma and shocks
tet = shock(1:T,1); 

% 5. Initialize the Algorithm Parameters
% --------------------------------------
beta    	= [log(cs)/(1-delta)-ks/cs/delta*log(ks); ks/cs/delta; 10^-5; 10^-5];  
criter  	= 1e-5;            	            	 % Convergence criterion
update  	= 0.5;            	                % Updating parameter for homotopy 
maxiter  = 100000;                            % Maximum number of iterations

% 6. The Main Loop 
% -----------------             
ksi      = beta;
iteration  = 1;                               % Initially, iteration is 0
dif	     = 2e-5;					             % Initally, criterion is not satified
kup  = ks*5;                                  % Upper bound on k 
klow = ks/5;                                  % Lower bound on k

while (dif > criter)&(iteration<=maxiter);

% 6.1 Given 'beta', compute the time series
% -----------------------------------------

for t = 1:T
   c(t)  = (k(t)^alpha*tet(t)+k(t)*(1-d))/(1+beta(2)*delta+rho*delta*beta(4)*log(tet(t))); 
   k(t+1)=delta*(beta(2)+beta(4)*log(tet(t))*rho)*c(t);
   k(t+1)=k(t+1)*(k(t+1)>klow)*(k(t+1)<kup)+klow*(k(t+1)<klow)+kup*(k(t+1)>kup);
   c(t)  = k(t)^alpha*tet(t)+k(t)*(1-d)-k(t+1);   
end;
   
% 6.2 Given simulated time series, compute the value function
% -------------------------------------------------------------
for t = 1:T-1
       w(t) = log(c(t))+delta*(beta(1)+beta(2)*log(k(t+1))+beta(3)*rho*log(tet(t))+beta(4)*rho*log(tet(t))*log(k(t+1)));      
end;
    
% 6.3 Recompute 'beta' by using NLLS regression
% ---------------------------------------------
   x = [ones(T-1,1) log( k(1:T-1) ) log( tet(1:T-1) ) log(k(1:T-1) ).*log( tet(1:T-1) )];  % Regressors 
   ksi = nlinfit(x,w,'objective',ksi);                   % NLLS regression
   iteration                                             % Display iteration
   dif = norm(beta-ksi)                                  % Display difference between 
   beta = update*ksi + (1-update)*beta;                  % Update the coefficients (homotopy)
   iteration = iteration+1;			                     % Next iteration
end;

% 7. Statistics  
% ---------------------------------
CPU = cputime-CPU0

% 8. Plot the time series solution  
% ---------------------------------
time=(1:1:T);                         
subplot(2,1,1);
plot (time,k(1:T,1)), xlabel('t'), ylabel('Capital k')
title('Time series solution');
subplot(2,1,2);
plot (time,c(1:T,1)), xlabel('t'), ylabel('Consumption')