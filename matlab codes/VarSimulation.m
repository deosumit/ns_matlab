% Calculates simulation data
% 
% A function which calculates simulation data 
% for the following system:
%
% y(t) = G1*y(t-1) + impact*e(t)
%
% Where y(t) is a nx1 vector, e(t) is a kx1 vector and 
% G1 and impact are matrices of the appropriate dimensions
%
% Inputs:
%
%    G1 and impact  --- see above
%    shocks         --- Txk vector, shocks
%    relative       --- kx1 vector, relative size of shocks
%    T              --- scalar, number of periods to be simulated
%    Burn           --- scalar, number of periods of burn in
%
% Output:
%
%    X              --- Txn matrix describing the evolution of
%                         y(t)
%
% Created by Jon Steinsson, June 2001 (modified Sept 2005)
%*********************************************************

function X = VarSimulation(G1,impact,shocks,relative,T,Burn)

%*********************************************************
% Calibration of shocks

for j = 1:(T+Burn)
    shocks(j,:) = relative'.*shocks(j,:);
end

%*********************************************************
% Definition of variables

n = size(G1,1);
y = zeros(n,1);
X = zeros(T+Burn,n);

%*********************************************************
% Calculations

for i = 1:(T+Burn)
    y = G1*y + impact*shocks(i,:)';
    X(i,:) = y';
end
X = X((Burn+1):end,:);

    