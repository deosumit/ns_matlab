% Two regions - sticky prices - heterogeneous factor markets 
% home bias
%
% Emi Nakamura and Jon Steinsson -- Feb 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calibration of parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sigma = 1;       % Intertemporal elastisity of substitution
alpha = 0.75;    % Calvo parameter
beta = 0.99;     % Subjective discount factor
nu = 1;          % Frish elasticity of labor suppy
aa = 0.67;        % Curvature of production function
eta = 2;         % Elasticity of substitution between home and foreign goods
theta = 7;      % Elasticity of substitution between varieties

delta = 0.025;   % Rate of depreciation of capital
sigma_a = 0.01;     % Elasticity of capital utilization cost
kappa_i = 2.5;     % investment adjustment cost parameter

nn = 0.1;        % Size of home region
phiH = 0.69;     % Weight of home goods in home consumption basket
phiF = 1-phiH;   % Weight of foreign goods in home consumption basket
phiHstar = phiF*nn/(1-nn);  % Weight of home goods in foreign consumption basket
phiFstar = 1 - phiHstar;    % Weight of foreign goods in foreign consumption basket

if ((min([phiH, phiF, phiHstar, phiFstar]) < 0) || (max([phiH, phiF, phiHstar, phiFstar]) > 1))
    fprintf('Home Bias Error\n')
    fprintf('phiH     = %5.5f\n',phiH)
    fprintf('phiF     = %5.5f\n',phiF)
    fprintf('phiHstar = %5.5f\n',phiHstar)
    fprintf('phiFstar = %5.5f\n',phiFstar)
    error('Home Bias Error')
end

phiPi = 1.5;     % Inflation response in Taylor Rule
phiY  = 0.5;     % Output response in Taylor Rule
rhoii = 0.8;     % Lagged dependence in Taylor Rule
phiG  = 0.0;    % Direct response of monetary policy to fiscal shock
phi_c = 1;       % Consumption elasticity of money demand
phi_i = 0;       % Interest elasticity of money demand

rhoG = 0.933;      % Persistence of government spending shock
corrG = 0.999;       % Correlation of home and foreign gov spending shocks

KoverL = ((1-aa)*(beta/(1-beta*(1-delta)))*(theta-1)/theta)^(1/aa);
Ibar   = delta*KoverL^aa; % Steady state investment-output ratio
Gbar = 0.2;      % Steady state government spending-output ratio
Cbar = 1 - Gbar - Ibar; % Steady state consumption-output ratio

kappa    = (1-alpha)*(1-alpha*beta)/alpha;
psi_v    = aa*(1/nu)/(1+(1/nu)*(1-aa));
zeta     = 1/(1+psi_v*theta);
sigma_c  = sigma^(-1)*(1-aa*((theta-1)/theta)*Cbar^(-1)*(1+1/nu)^(-1))^(-1);
sigma_l  = (Cbar^(-1)*((theta-1)/theta)*aa*sigma_c^(-1))^(-1);
psi_c    = aa/(1+(1/nu)*(1-aa));
psi_k    = ((1/nu)+1)*(1-aa)/(1+(1/nu)*(1-aa));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose the Monetary Policy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TaylorBasic         = 1;
ConstantMoneySupply = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of system
%
% g0*x(t)=g1*x(t-1)+C+Psi*z(t)+Pi*eta(t)
% 
% x_(t)= [
%  [1]    c_t
%  [2]    E_t c_t+1
%  [3]    c*_t
%  [4]    E_t c*_t+1
%  [5]    piH_t
%  [6]    E_t piH_t+1
%  [7]    piF_t
%  [8]    E_t piF_t+1
%  [9]    pi_t
%  [10]    E_t pi_t+1
%  [11]    pi*_t
%  [12]    E_t pi*_t+1
%  [13]    y_t
%  [14]   y*_t
%  [15]   q_t
%  [16]   pH_t
%  [17]   r_nt  (nominal interest rate)
%  [18]   g_Ht
%  [19]   g_Ft
%  [20]   ny_t
%  [21]   ny*_t
%  [22]   L_t
%  [23]   L*_t
%  [24]   w_t - p_Ht (real wage with PPI as deflator)
%  [25]   w*_t - p_Ft + q_t (NB: p_Ft = -(phiH/phiF)p_Ht)
%  [26]   s_t + p_Ht (real MC with PPI as deflator)
%  [27]   s*_t + p_Ft 
%  [28]   kbar_t+1 (home installed capital)
%  [29]   kbar*_t+1 
%  [30]   k_t (home capital input (as oppose to installed capital))
%  [31]   k*_t 
%  [32]   i_t (home investment)
%  [33]   i*_t
%  [34]   u_t (home capital utilization)
%  [35]   u*_t
%  [36]   d_t (home shadow value of installed capital)
%  [37]   d*_t
%  [38]   r_kt (home  capital rental rate)
%  [39]   E_t d_t+1
%  [40]   E_t d*_t+1
%  [41]   E_t r_kt+1
%  [42]   r*_kt
%  [43]   E_t r*_kt+1
%  [44]   E_t u_t+1
%  [45]   E_t u*_t+1
%  [46]   E_t L_t+1
%  [47]   E_t L*_t+1
%  [48]   E_t i_t+1
%  [49]   E_t i*_t+1]

dimg0 = 49;
cc = zeros(dimg0,1);
g0 = zeros(dimg0,dimg0);
g1 = zeros(dimg0,dimg0);
pi = zeros(dimg0,16);  % Expectational Errors
psi = zeros(dimg0,3); % Disturbances

% Home Consumption Euler Equation
g0(1,1) = -1;
g0(1,2) = 1;
g0(1,22) = sigma_c*sigma_l^(-1);
g0(1,46) = -sigma_c*sigma_l^(-1);
g0(1,17) = -sigma_c;
g0(1,10) = sigma_c;

% Backus-Smith Condition
g0(2,1) = 1;
g0(2,3) = -1;
g0(2,22) = -sigma_c*sigma_l^(-1);
g0(2,23) = sigma_c*sigma_l^(-1);
g0(2,15) = -sigma_c;

% Home Phillips Curve
g0(3,5) = -1;
g0(3,6) = beta;
%g0(3,1) = kappa*zeta*sigma^(-1)*psi_c;
g0(3,13) = kappa*zeta*psi_v;
g0(3,16) = -kappa*zeta;
g0(3,38) = kappa*zeta*psi_k;

% Foreign Phillips Curve
g0(4,7) = -1;
g0(4,8) = beta;
%g0(4,3) = kappa*zeta*sigma^(-1)*psi_c;
g0(4,14) = kappa*zeta*psi_v;
g0(4,15) = kappa*zeta*psi_c;
g0(4,16) = phiH/phiF*kappa*zeta;
g0(4,42) = kappa*zeta*psi_k;

% Home Inflation
g0(5,9) = -1;
g0(5,5) = phiH;
g0(5,7) = phiF;

% Foreign Inflation
g0(6,11) = -1;
g0(6,5) = phiHstar;
g0(6,7) = phiFstar;

% Home Resource Constraint
g0(7,13) = -1;
g0(7,1) = Cbar*phiH;
g0(7,3) = Cbar*(1-nn)/nn*phiHstar;
g0(7,15) = eta*(Cbar+Ibar)*(1-nn)/nn*phiHstar;
g0(7,16) = -eta*(Cbar+Ibar)*(phiH+(1-nn)/nn*phiHstar);
g0(7,18) = 1;
g0(7,32) = Ibar*phiH;
g0(7,33) = Ibar*(1-nn)/nn*phiHstar;

% Foreign Resource Constraint
g0(8,14) = -1;
g0(8,1) = Cbar*nn/(1-nn)*phiF;
g0(8,3) = Cbar*phiFstar;
g0(8,15) = eta*(Cbar+Ibar)*phiFstar;
g0(8,16) = eta*(Cbar+Ibar)*(nn/(1-nn)*phiF+phiFstar)*phiH/phiF;
g0(8,19) = 1;
g0(8,32) = Ibar*nn/(1-nn)*phiF;
g0(8,33) = Ibar*phiFstar;

% Home Relative Price
g0(9,16) = 1;
g1(9,16) = 1;
g0(9,5) = -1;
g0(9,9) = 1;

% Real Exchange Rate
g0(10,15) = -1;
g0(10,16) = phiHstar-phiH/phiF*phiFstar;

if TaylorBasic == 1
    g0(11,17) = -1;
    g1(11,17) = -rhoii;
    g0(11,9) = (1-rhoii)*nn*phiPi;
    g0(11,11) = (1-rhoii)*(1-nn)*phiPi;
    g0(11,13) = (1-rhoii)*nn*phiY;
    g0(11,14) = (1-rhoii)*(1-nn)*phiY;
    g0(11,18) = (1-rhoii)*nn*phiG;
    g0(11,19) = (1-rhoii)*(1-nn)*phiG;
    psi(11,3) = 1;
end

if ConstantMoneySupply == 1
    g0(11,:) = zeros(1,dimg0);
    g1(11,:) = zeros(1,dimg0);
    g0(11,17) = phi_i;
    g0(11,9)  = -nn;
    g0(11,11)  = -(1-nn);
    g0(11,1)  = nn*(-phi_c);
    g0(11,3)  = (1-nn)*(-phi_c);
    g1(11,17) = phi_i;
    g1(11,1)  = nn*(-phi_c);
    g1(11,1)  = (1-nn)*(-phi_c);
    psi(11,3) = -1;
end

% Home Government Spending
g0(12,18) = 1;
g1(12,18) = rhoG;
psi(12,1) = 1;

% Foreign Government Spending
g0(13,19) = 1;
g1(13,19) = rhoG;
psi(13,2) = 1;

% Expectation Error Equations
g0(14,1) = 1;
g1(14,2) = 1;
pi(14,1) = 1;

g0(15,3) = 1;
g1(15,4) = 1;
pi(15,2) = 1;

g0(16,5) = 1;
g1(16,6) = 1;
pi(16,3) = 1;

g0(17,7) = 1;
g1(17,8) = 1;
pi(17,4) = 1;

g0(18,9) = 1;
g1(18,10) = 1;
pi(18,5) = 1;

g0(19,11) = 1;
g1(19,12) = 1;
pi(19,6) = 1;

% Home Nominal Output
g0(20,20) = -1;
g0(20,13) = 1;
g0(20,5) = 1;

% Foreign Nominal Output
g0(21,21) = -1;
g0(21,14) = 1;
g0(21,7) = 1;

% Home production function
g0(22,22) = -aa;
g0(22,13) = 1;
g0(22,30) = -(1-aa);

% Foreign production function
g0(23,23) = -aa;
g0(23,14) = 1;
g0(23,31) = -(1-aa);

% Home production real wage
g0(24,24) = 1;
%g0(24,1) = -sigma^(-1);
g0(24,22) = -nu^(-1);
g0(24,16) = 1;

% Foreign production real wage
g0(25,25) = 1;
%g0(25,3) = -sigma^(-1);
g0(25,23) = -nu^(-1);
g0(25,15) = -1;
g0(25,16) = -phiH/phiF;

% Home marginal cost
g0(26,26) = 1;
g0(26,1) = -sigma^(-1)*psi_c;
g0(26,13) = -psi_v;
g0(26,16) = 1;
g0(26,38) = -psi_k;

% Foreign marginal cost
g0(27,27) = 1;
g0(27,3) = -sigma^(-1)*psi_c;
g0(27,14) = -psi_v;
g0(27,15) = -psi_c;
g0(27,16) = -phiH/phiF;
g0(26,42) = -psi_k;

% Home capital accumulation
g0(28,28) = 1;
g1(28,28) = (1-delta);
g0(28,32) = -delta;

% Foreign capital accumulation
g0(29,29) = 1;
g1(29,29) = (1-delta);
g0(29,33) = -delta;

% Home capital input
g0(30,30) = 1;
g1(30,28) = 1;
g0(30,34) = -1;

% Foreign capital input
g0(31,31) = 1;
g1(31,29) = 1;
g0(31,35) = -1;

% Home investment
g0(32,36) = 1;
g0(32,32) = -kappa_i*(1+beta);
g1(32,32) = -kappa_i;
g0(32,48) = beta*kappa_i;
g0(32,1) = sigma_c^(-1);
g0(32,22) = -sigma_l^(-1);

% Foreign investment
g0(33,37) = 1;
g0(33,33) = -kappa_i*(1+beta);
g1(33,33) = -kappa_i;
g0(33,49) = beta*kappa_i;
g0(33,3) = sigma_c^(-1);
g0(33,23) = -sigma_l^(-1);

% Home capital utilization
g0(34,34) = sigma_a;
g0(34,38) = -1;

% Foreign capital utilization
g0(35,35) = sigma_a;
g0(35,42) = -1;
g0(35,15) = 1;

% Home shadow value of installed capital
g0(36,36) = -1;
g0(36,39) = beta*(1-delta);
g0(36,41) = (1-beta*(1-delta));
g0(36,2) = -sigma_c^(-1)*(1-beta*(1-delta));
g0(36,46) = sigma_l^(-1)*(1-beta*(1-delta));

% Foreign shadow value of installed capital
g0(37,37) = -1;
g0(37,40) = beta*(1-delta);
g0(37,43) = (1-beta*(1-delta));
g0(37,10) = (1-beta*(1-delta));
g0(37,15) = -(1-beta*(1-delta));
g0(37,12) = -(1-beta*(1-delta));
g0(37,4) = -sigma_c^(-1)*(1-beta*(1-delta));
g0(37,47) = sigma_l^(-1)*(1-beta*(1-delta));

% Home capital demand
g0(38,38) = -1;
g0(38,22) = nu^(-1)+1;
g0(38,30) = -1;

% Expectations Errors
g0(39,36) = 1;
g1(39,39) = 1;
pi(39,7) = 1;

g0(40,37) = 1;
g1(40,40) = 1;
pi(40,8) = 1;

g0(41,38) = 1;
g1(41,41) = 1;
pi(41,9) = 1;

% Foreign capital demand
g0(42,42) = -1;
g0(42,23) = nu^(-1)+1;
g0(42,31) = -1;
g0(42,15) = 1;

% Expectations Error
g0(43,42) = 1;
g1(43,43) = 1;
pi(43,10) = 1;

g0(44,34) = 1;
g1(44,44) = 1;
pi(44,11) = 1;

g0(45,35) = 1;
g1(45,45) = 1;
pi(45,12) = 1;

g0(46,22) = 1;
g1(46,46) = 1;
pi(46,13) = 1;

g0(47,23) = 1;
g1(47,47) = 1;
pi(47,14) = 1;

g0(48,32) = 1;
g1(48,48) = 1;
pi(48,15) = 1;

g0(49,33) = 1;
g1(49,49) = 1;
pi(49,16) = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run Gensys
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[G1,C,impact,fmat,fwt,ywt,gev,eu] = gensys(g0,g1,cc,psi,pi);

% Shock to Goverment Spending
shock = [1;0;0];

% Impulse Response
period = 40; % specify number of period
irf = VarImpulse(G1,impact,shock,period);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Relative Regressions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation
NSim = 200000;
NBurn = 50;
mu = zeros(NSim+NBurn,3);
sigma_shock = [1 corrG 0; corrG 1 0; 0 0 1];
shocks = mvnrnd(mu,sigma_shock);
relative = [1; 10^(-6); 10^(-6)];
XSim = VarSimulation(G1,impact,shocks,relative,NSim,NBurn);
NSim = size(XSim);
NSim = NSim(1);

% Convert to annual data
NAnnual = floor(NSim/4);
NColumns = size(XSim);
NColumns = NColumns(2);
XAnnual = zeros(NAnnual,NColumns);
for ii = 1:NAnnual 
    XAnnual(ii,:) = sum(XSim((4*(ii-1)+1):(4*ii),:))/4;
end
XSimOrig = XSim;
XSim = XAnnual;
NSim = NAnnual;

% Create annual price level
Nquarters = size(XSimOrig);
Nquarters = Nquarters(1);
Xprices = zeros(Nquarters,4);
Xprices(1,1) = XSimOrig(1,5);
Xprices(1,2) = XSimOrig(1,7);
Xprices(1,3) = XSimOrig(1,9);
Xprices(1,4) = XSimOrig(1,11);
for ii = 2:Nquarters
    Xprices(ii,1) = Xprices(ii-1,1) + XSimOrig(ii,5);
    Xprices(ii,2) = Xprices(ii-1,2) + XSimOrig(ii,7);
    Xprices(ii,3) = Xprices(ii-1,3) + XSimOrig(ii,9);
    Xprices(ii,4) = Xprices(ii-1,4) + XSimOrig(ii,11);
end

XpricesAnnual = zeros(NAnnual,4);
for ii = 2:NAnnual
    XpricesAnnual(ii,:) = sum(Xprices((4*(ii-1)+1):(4*ii),:))/4;
end

% Relative Regressions
dgH = XSim(3:NSim,18) - XSim(1:(NSim-2),18) + (1-Cbar)*(XpricesAnnual(3:NSim,1) - XpricesAnnual(1:(NSim-2),1));
dgF = XSim(3:NSim,19) - XSim(1:(NSim-2),19) + (1-Cbar)*(XpricesAnnual(3:NSim,2) - XpricesAnnual(1:(NSim-2),2));
CC = ones(NSim-2,1);
NXX = [CC, dgH-dgF];

dgH = XSim(3:NSim,18) - XSim(1:(NSim-2),18);
dgF = XSim(3:NSim,19) - XSim(1:(NSim-2),19);
RXX = [CC, dgH-dgF];

fprintf('\n')
fprintf('Relative Regression Results\n')

% Quasi-Real Output Regression
dyH = XSim(3:NSim,13) - XSim(1:(NSim-2),13) + (XpricesAnnual(3:NSim,1) - XpricesAnnual(1:(NSim-2),1)) - (XpricesAnnual(3:NSim,3) - XpricesAnnual(1:(NSim-2),3));
dyF = XSim(3:NSim,14) - XSim(1:(NSim-2),14) + (XpricesAnnual(3:NSim,2) - XpricesAnnual(1:(NSim-2),2)) - (XpricesAnnual(3:NSim,4) - XpricesAnnual(1:(NSim-2),4));
YY = dyH - dyF;
bb = (NXX'*NXX)\NXX'*YY;
fprintf('beta output         = %5.5f\n',bb(2,1))

% CPI Regression
dyH = XpricesAnnual(3:NSim,3) - XpricesAnnual(1:(NSim-2),3);
dyF = XpricesAnnual(3:NSim,4) - XpricesAnnual(1:(NSim-2),4);
YY = dyH - dyF;
bb = (NXX'*NXX)\NXX'*YY;
fprintf('beta cpi            = %5.5f\n',bb(2,1))

fprintf('\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% World Regressions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation
NSim = 200000;
NBurn = 50;
mu = zeros(NSim+NBurn,3);
sigma_shock = [1 corrG 0; corrG 1 0; 0 0 1];
shocks = mvnrnd(mu,sigma_shock);
relative = [1; 1; 10^(-6)];
XSim = VarSimulation(G1,impact,shocks,relative,NSim,NBurn);
NSim = size(XSim);
NSim = NSim(1);

% Convert to annual data
NAnnual = floor(NSim/4);
NColumns = size(XSim);
NColumns = NColumns(2);
XAnnual = zeros(NAnnual,NColumns);
for ii = 1:NAnnual 
    XAnnual(ii,:) = sum(XSim((4*(ii-1)+1):(4*ii),:))/4;
end
XSimOrig = XSim;
XSim = XAnnual;
NSim = NAnnual;

% Create annual price level
Nquarters = size(XSimOrig);
Nquarters = Nquarters(1);
Xprices = zeros(Nquarters,4);
Xprices(1,1) = XSimOrig(1,5);
Xprices(1,2) = XSimOrig(1,7);
Xprices(1,3) = XSimOrig(1,9);
Xprices(1,4) = XSimOrig(1,11);
for ii = 2:Nquarters
    Xprices(ii,1) = Xprices(ii-1,1) + XSimOrig(ii,5);
    Xprices(ii,2) = Xprices(ii-1,2) + XSimOrig(ii,7);
    Xprices(ii,3) = Xprices(ii-1,3) + XSimOrig(ii,9);
    Xprices(ii,4) = Xprices(ii-1,4) + XSimOrig(ii,11);
end

XpricesAnnual = zeros(NAnnual,4);
for ii = 2:NAnnual
    XpricesAnnual(ii,:) = sum(Xprices((4*(ii-1)+1):(4*ii),:))/4;
end

pW = nn*XpricesAnnual(:,3)+(1-nn)*XpricesAnnual(:,4);

% Regressors
dgH = XSim(3:NSim,18) - XSim(1:(NSim-2),18);
dgF = XSim(3:NSim,19) - XSim(1:(NSim-2),19);
XX = [CC, nn*dgH+(1-nn)*dgF];
fprintf('World Regression Results\n')

% Output Regression
dyH = XSim(3:NSim,13) - XSim(1:(NSim-2),13);
dyF = XSim(3:NSim,14) - XSim(1:(NSim-2),14);
YY = nn*dyH + (1-nn)*dyF;
bb = (XX'*XX)\XX'*YY;
fprintf('beta output         = %5.5f\n',bb(2,1))

% CPI Regression
dyH = XpricesAnnual(3:NSim,3) - XpricesAnnual(1:(NSim-2),3);
dyF = XpricesAnnual(3:NSim,4) - XpricesAnnual(1:(NSim-2),4);
YY = nn*dyH + (1-nn)*dyF;
bb = (XX'*XX)\XX'*YY;
fprintf('beta cpi            = %5.5f\n',bb(2,1))

fprintf('\n')
