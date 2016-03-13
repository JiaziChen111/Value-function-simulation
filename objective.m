%     Objective function to solve the one- and two-sector 
%     models from the article "Solving Nonlinear Dynamic
%     Stochastic Models: An Algorithm Iterating on
%     Value Function by Simulations".
%
%     April 5, 2002
%     -------------------------------------------------------------------
% ksiz   initial coefficients 
% x      regressors
% y      explanatory variable 
% ---------------------------

function y = objective(ksiz,x)
y=x*ksiz;