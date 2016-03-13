%     MATLAB program to solve the two-sector model
%     from the article "Solving Nonlinear Dynamic
%     Stochastic Models: An Algorithm Iterating on
%     Value Function by Simulations".
%          
%     1. TwoSector.m
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
T       = 1000;         % Length of simulation

% 2. The steady state of k, h and c 
% -----------------------------------------------
ks = ( (1-delta+delta*d) / (alpha*delta) )^(1/(alpha-1) );
hs = ks;
cs = 2*(ks^alpha - d*ks);												  

% 3. Allocate memory for the simulated series
% -------------------------------------------
k    = zeros(T+1,1)+ks;  % Capital k 
h    = zeros(T+1,1)+hs;  % Capital h
c    = zeros(T,1);       % Consumption
w    = zeros(T-1,1);     % Value function

% 4. Draw the shock
% -------------------------------------------
load shock05; 
tetk = shock(1:T,1); 
teth = shock(1:T,2);
prf  = shock(1:T,3); 
dep  = shock(1:T,4);

% 5. Initialize the Algorithm Parameters
% --------------------------------------
beta    	= [log(cs)/(1-delta)-2*ks/cs/delta*log(ks); ks/cs/delta; ks/cs/delta; 10^-5; 10^-5; 10^-5; 10^-5; 10^-5; 10^-5; 10^-5; 10^-5; 10^-5; 10^-5; 10^-5; 10^-5];
criter  	= 1e-5;            	              % Convergence criterion
update  	= .5;            	                 % Updating parameter for homotopy 
maxiter  = 100000;                          % Maximum number of iterations

% 6. The Main Loop 
% -----------------             
ksi      = beta;
iteration  = 1;                               % Initially, iteration is 0
dif	     = 2e-5;					             % Initally, criterion is not satified
kup  = ks*5;  hup  = hs*5;                    % Upper bound on k 
klow = ks/5;  hlow = hs/5;                    % Lower bound on k

while (dif > criter)&(iteration<=maxiter);

% 6.1 Given 'beta', compute the time series
% -----------------------------------------

for t = 1:T
   c(t)  = (k(t)^alpha*tetk(t)+ h(t)^alpha*teth(t)+(1-d*dep(t))*(k(t)+h(t)))/(1+delta*(beta(2)+beta(3)+(beta(8)+beta(12))*log(tetk(t))*rho+(beta(9)+beta(13))*log(teth(t))*rho+(beta(10)+beta(14))*log(prf(t))*rho+(beta(11)+beta(15))*log(dep(t))*rho)/prf(t)); 
   k(t+1)=delta/prf(t)*(beta(2)+beta(8)*log(tetk(t))*rho+beta(9)*log(teth(t))*rho+beta(10)*log(prf(t))*rho+beta(11)*log(dep(t))*rho)*c(t);
   k(t+1)=k(t+1)*(k(t+1)>klow)*(k(t+1)<kup)+klow*(k(t+1)<klow)+kup*(k(t+1)>kup);
   h(t+1)=delta/prf(t)*(beta(3)+beta(12)*log(tetk(t))*rho+beta(13)*log(teth(t))*rho+beta(14)*log(prf(t))*rho+beta(15)*log(dep(t))*rho)*c(t);
   h(t+1)=h(t+1)*(h(t+1)>hlow)*(h(t+1)<hup)+hlow*(h(t+1)<hlow)+hup*(h(t+1)>hup);
   c(t)  = k(t)^alpha*tetk(t)+h(t)^alpha*teth(t) + (1-d*dep(t))*(k(t)+h(t))-k(t+1)-h(t+1);   
end;
   
% 6.2 Given simulated time series, compute the value function
% -------------------------------------------------------------
for t = 1:T-1
       w(t) = log(c(t))+delta*(beta(1)+beta(2)*log(k(t+1))+beta(3)*log(h(t+1))+beta(4)*rho*log(tetk(t))+beta(5)*rho*log(teth(t))+beta(6)*rho*log(prf(t))+beta(7)*rho*log(dep(t))+log(k(t+1))*(beta(8)*rho*log(tetk(t))+beta(9)*rho*log(teth(t))+beta(10)*rho*log(prf(t))+beta(11)*rho*log(dep(t)))+log(h(t+1))*(beta(12)*rho*log(tetk(t))+beta(13)*rho*log(teth(t))+beta(14)*rho*log(prf(t))+beta(15)*rho*log(dep(t))));      
end;
    
% 6.3 Recompute 'beta' by using NLLS regression
% ---------------------------------------------
x = [ones(T-1,1) log( k(1:T-1) ) log( h(1:T-1) ) log( tetk(1:T-1) ) log( teth(1:T-1) ) log( prf(1:T-1) ) log( dep(1:T-1) )...
log(k(1:T-1) ).*log( tetk(1:T-1) ) log(k(1:T-1) ).*log( teth(1:T-1) ) log(k(1:T-1) ).*log( prf(1:T-1) ) log(k(1:T-1) ).*log( dep(1:T-1) ) log(h(1:T-1) ).*log( tetk(1:T-1) ) log(h(1:T-1) ).*log( teth(1:T-1) ) log(h(1:T-1) ).*log( prf(1:T-1) ) log(h(1:T-1) ).*log( dep(1:T-1) )];  % Regressors 
   ksi = nlinfit(x,w,'objective',ksi);                   % NLLS regression
   iteration                                             % Display iteration
   dif = norm(beta-ksi)                                  % Display difference between 
   beta = update*ksi + (1-update)*beta;                  % Update the coefficients (homotopy)
   iteration = iteration+1;			                     % Next iteration
end;

% 7. Statistics  
% ---------------------------------
CPU = cputime-CPU0

% 8. Plot the time series solution y, c and k 
% -------------------------------------------
time=(1:1:T);                         
subplot(3,1,1);
plot (time,k(1:T,1)), xlabel('t'), ylabel('Capital k')
title('Time series solution');
subplot(3,1,2);
plot (time,h(1:T,1)), xlabel('t'), ylabel('Capital h')
subplot(3,1,3);
plot (time,c(1:T,1)), xlabel('t'), ylabel('Consumption')