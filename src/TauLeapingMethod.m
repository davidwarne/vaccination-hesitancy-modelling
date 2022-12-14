function [Z,t] = TauLeapingMethod(dsct,T,tau)
%% The Tau Leaping Method
% An approximate stochastic simulation algorithm for a discrete-state 
% continuous-time Markov process.
%
% Inputs:
%    DSCT - model definition structure
%    T    - the end time of the simulation
%    tau   - the timestep
% Outputs:
%    Z    -  time series of copy number vectors
%    t    -  vector of times
%
% Author:
%   David J. Warne (david.warne@qut.edu.au)
%         School of Mathematical Sciences
%         Faculty of Science
%         Queensland University of Technology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialise
Nt = floor(T/tau);
Z = zeros(length(dsct.X0),Nt+1);
t = zeros(1,Nt+1);
Z(:,1) = dsct.X0;
for i=1:Nt
    % compute time-dependent propensities 
    a = dsct.a(Z(:,i),dsct.k,t(i));
    % generate Poisson random variates
    Y = poissrnd(a*tau);
    % update copy numbers
    Z(:,i+1) = Z(:,i) + (dsct.nu') * Y;
    Z(Z(:,i+1) < 0,i+1) = 0;
    t(i+1) = t(i) + tau;
end
