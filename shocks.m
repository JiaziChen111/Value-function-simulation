%     Constructing shocks for the one- and two-sector 
%     models from the article "Solving Nonlinear Dynamic
%     Stochastic Models: An Algorithm Iterating on
%     Value Function by Simulations".
%
%     April 5, 2002
%     -------------------------------------------------------------------

sigma   = .05;                % Standard deviation for log noise
rho     = 0.95;               % Persistence of log technology shock
T       = 10000;              % Length of simulation
shock = ones(T,4);            % Matrix for shocks
epsi  = randn(T,4)*sigma;     % Innovations
for t = 2:T; 
   shock(t,:) = shock(t-1,:).^rho.*exp(epsi(t,:)); 
end;
save shock05 sigma rho shock;