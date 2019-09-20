% Two regions - sticky prices - heterogeneous factor markets --
% home bias -- GHH preferences
%
% Emi Nakamura and Jon Steinsson -- Oct 2010
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
theta = 7;       % Elasticity of substitution between varieties

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

Gbar = 0.2;     % Steady state government spending - output ratio
Cbar = 1 - Gbar; % Steady state consumption-output ratio

phiPi = 1.5;     % Inflation response in Taylor Rule
phiY  = 0.5;     % Output response in Taylor Rule
rhoii = 0.8;     % Lagged dependence in Taylor Rule
phi_c = 1;       % Consumption elasticity of money demand
phi_i = 1;       % Interest elasticity of money demand
phi_p = 0.025;     % Response to price level in Wickselling rule

rhoG = 0.933; %0.5^(1/8);     % Persistence of government spending shock
corrG = 0.999;       % Correlation of home and foreign gov spending shocks

kappa   = (1-alpha)*(1-alpha*beta)/alpha;
psi_v   = (1+nu*(1-aa))/(nu*aa);
zeta    = 1/(1+psi_v*theta);
sigma_c = sigma^(-1)*(1-aa*((theta-1)/theta)*Cbar^(-1)*(1+1/nu)^(-1))^(-1);
xi_y    = Cbar^(-1)*((theta-1)/theta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose the Monetary Policy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TaylorBasic         = 1;
ConstantRealRate    = 0;
ConstantNominalRate = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose Tax Policy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NoTaxes             = 1;
BalancedBudget      = 0;

% Steady state taxes
if NoTaxes == 1
    taubar = 0;
end
if BalancedBudget == 1
    taubar = theta/(theta-1)*Gbar; 
end

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
%  [17]   i_t
%  [18]   g_Ht
%  [19]   g_Ft
%  [20]   ny_t
%  [21]   ny*_t
%  [22]   L_t
%  [23]   L*_t
%  [24]   E_t y_t+1
%  [25]   p_t (home cpi price level)
%  [26]   p*_t (foreign cpi price level)
%  [27]   w_t - p_Ht (real wage with PPI as deflator)
%  [28]   w*_t - p_Ft + q_t (NB: p_Ft = -(phiH/phiF)p_Ht)
%  [29]   tau_t]

dimg0 = 29;
cc = zeros(dimg0,1);
g0 = zeros(dimg0,dimg0);
g1 = zeros(dimg0,dimg0);
pi = zeros(dimg0,7);  % Expectational Errors
psi = zeros(dimg0,3); % Disturbances

% Home Consumption Euler Equation
g0(1,1) = -1;
g0(1,2) = 1;
g0(1,13) = xi_y;
g0(1,24) = -xi_y;
g0(1,17) = -sigma_c;
g0(1,10) = sigma_c;

% Backus-Smith Condition
g0(2,1) = 1;
g0(2,3) = -1;
g0(2,13) = -xi_y;
g0(2,14) = xi_y;
g0(2,15) = -sigma_c;

% Home Phillips Curve
g0(3,5) = -1;
g0(3,6) = beta;
%g0(3,1) = kappa*sigma_c^(-1)*zeta;
g0(3,13) = kappa*psi_v*zeta;
g0(3,16) = -kappa*zeta;
g0(3,29) = taubar/(1-taubar)*kappa*zeta;

% Foreign Phillips Curve
g0(4,7) = -1;
g0(4,8) = beta;
%g0(4,3) = kappa*sigma_c^(-1)*zeta;
g0(4,14) = kappa*psi_v*zeta;
g0(4,15) = kappa*zeta;
g0(4,16) = phiH/phiF*kappa*zeta;
g0(4,29) = taubar/(1-taubar)*kappa*zeta;

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
g0(7,15) = eta*Cbar*(1-nn)/nn*phiHstar;
g0(7,16) = -eta*Cbar*(phiH+(1-nn)/nn*phiHstar);
g0(7,18) = 1;

% Foreign Resource Constraint
g0(8,14) = -1;
g0(8,1) = Cbar*nn/(1-nn)*phiF;
g0(8,3) = Cbar*phiFstar;
g0(8,15) = eta*Cbar*phiFstar;
g0(8,16) = eta*Cbar*(nn/(1-nn)*phiF+phiFstar)*phiH/phiF;
g0(8,19) = 1;

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
    psi(11,3) = 1;
end

if ConstantRealRate == 1
    a_pi = kappa*psi_v*zeta/(1-beta*rhoG)/(1-Cbar*xi_y);
    g0(11,:) = zeros(1,dimg0);
    g1(11,:) = zeros(1,dimg0);
    g0(11,17) = -1;
    g0(11,9) = nn*phiPi;
    g0(11,11) = (1-nn)*phiPi;
    g0(11,18) = -nn*a_pi*(phiPi-rhoG);
    g0(11,19) = -(1-nn)*a_pi*(phiPi-rhoG);
    psi(11,3) = 1;
end

if ConstantNominalRate == 1
    AA = [-(1-beta*rhoG) kappa*psi_v*zeta;
          Cbar*sigma_c*rhoG -(1-rhoG)*(1-Cbar*xi_y)];
    bb = [0; -(1-rhoG)];
    xx = AA\bb;
    a_pi = xx(1,1);
    g0(11,:) = zeros(1,dimg0);
    g1(11,:) = zeros(1,dimg0);
    g0(11,17) = -1;
    g0(11,9) = nn*phiPi;
    g0(11,11) = (1-nn)*phiPi;
    g0(11,18) = -nn*a_pi*phiPi;
    g0(11,19) = -(1-nn)*a_pi*phiPi;
    psi(11,3) = 1;
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
g0(20,9) = 1;

% Foreign Nominal Output
g0(21,21) = -1;
g0(21,14) = 1;
g0(21,11) = 1;

% Home labor
g0(22,22) = -aa;
g0(22,13) = 1;

% Foreign labor
g0(23,23) = -aa;
g0(23,14) = 1;

% Expectation Error Equations
g0(24,13) = 1;
g1(24,24) = 1;
pi(24,7) = 1;

% Home CPI Price Level
g0(25,25) = 1;
g1(25,25) = 1;
g0(25,9) = -1;

% Foreign CPI Price Level
g0(26,26) = 1;
g1(26,26) = 1;
g0(26,11) = -1;

% Home production real wage
g0(27,27) = 1;
%g0(27,1) = -sigma^(-1);
g0(27,13) = -aa^(-1)*nu^(-1);
g0(27,16) = 1;
g0(27,29) = -taubar/(1-taubar);

% Foreign production real wage
g0(28,28) = 1;
%g0(28,3) = -sigma^(-1);
g0(28,14) = -aa^(-1)*nu^(-1);
g0(28,15) = -1;
g0(28,16) = -phiH/phiF;
g0(28,29) = -taubar/(1-taubar);

% No Taxes
g0(29,29) = 1;

% Balanced Budget
if BalancedBudget == 1
   g0(29,29) = 1;
   g0(29,27) = nn;
   g0(29,28) = 1-nn;
   g0(29,13) = nn*aa^(-1);
   g0(29,14) = (1-nn)*aa^(-1);
   g0(29,18) = -nn;
   g0(29,19) = -(1-nn);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run Gensys
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[G1,C,impact,fmat,fwt,ywt,gev,eu] = gensys(g0,g1,cc,psi,pi);

% Shock to Goverment Spending
shock = [1;0;0];

% Impulse Response
period = 40; % specify number of period
irf = VarImpulse(G1,impact,shock,period);

irfpp = zeros(period+1,1);
for ii = 2:(period+1)
    irfpp(ii,1) = irfpp(ii-1,1) + nn*irf(ii,9) + (1-nn)*irf(ii,11);
end

% Create Price Levels
IrfPrices = zeros(period+1,4);
for ii = 2:(period+1)
    IrfPrices(ii,1) = IrfPrices(ii-1,1) + irf(ii,5);
    IrfPrices(ii,2) = IrfPrices(ii-1,2) + irf(ii,7);
    IrfPrices(ii,3) = IrfPrices(ii-1,3) + irf(ii,9);
    IrfPrices(ii,4) = IrfPrices(ii-1,4) + irf(ii,11);
end

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

% Output Regression
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

