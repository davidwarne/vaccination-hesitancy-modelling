function [Data_s] = simuldata_reg_fA_vax_h(Data,theta)
%% SIMULDATA
% This function simulates stochastic epidemic curve with latent infectous population
% and regulatory mechnisms and vaccine uptake process. 
% Model is a discrete state coninuous time Markov process simulated using a 
% Tau-Leaping method.
%
% Parameters:
%
% Data - observed case data and vaccination data use to initialise the time series.
%
% theta - Vector of model parameters
%
% Returns:
%     Data_s - simulated data using the model with given parameters.
%
% Authors:
%     David J. Warne (david.warne@qut.edu.au)
%           School of Mathematical Sciences
%           Faculty of Science 
%           Queensland University of Technology
%     
%     Abhishek Varghese 
%           Centre for Data Science
%           Queensland University of Technology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract parameters

% for transmission/SEIR + obseration process
k.alpha0 = theta(1)/Data.P;     % residual transmission 
k.alpha = theta(2)/Data.P;      % transmission 
k.beta = theta(3);              % symptom onset rate 
k.gamma = theta(4);             % case detection rate 
k.delta = theta(5);             % case death rate 
k.rho = theta(6);               % case recovery rate
k.eta = theta(7);               % unobserved recovery rate 

% vaccine effect
k.v = theta(8);                 % vaccination rate (1st dose) 
k.alphav1 = theta(9);           % transmission reduction (1st dose) 
k.gammav1 = theta(10);          % detection reduction (1st dose) 
k.rhov1 = theta(11);            % recovery increase (1st dose)
k.deltav1 = theta(12);          % death decrease (1st dose)
k.omega = theta(13);            % rate to 2nd dose 
k.alphav2 = theta(14);          % transmission reduction (2nd dose) 
k.gammav2 = theta(15);          % detection reduction (2nd dose) 
k.rhov2 = theta(16);            % recovery increase (2nd dose) 
k.deltav2 = theta(17);          % death decrease (2nd dose) 

% initial conditions
k.kappa = theta(18);        % initial under-estimate factor 
k.zeta = theta(19);         % R0 = (1-zeta)*(C0-D0)

%Function parameters
k.n = theta(20);            % inibitive strength 
k.w = theta(21);            % active case weight on NPI response
k.tl0 = theta(22);          % day that NPI ceases

k.nv = theta(23);           % promotor strength for vaccination
k.wC = theta(24);           % cases weight on vaccination
k.wD = theta(25);           % deaths weight on vaccination
k.wV = theta(26);           % vaccinations weight on vaccination
k.tv0 = theta(27);          % day the vaccine is available 
start = 1;

%copy initial condition
Data_s.C = Data.C(start:end);   % cumulative confirmed
Data_s.D = Data.D(start:end);   % died
Data_s.V1 = Data.V1(start:end); % Vaccinated (1st dose)
Data_s.V2 = Data.V2(start:end); % Vaccinated (2nd dose)
Data_s.P = Data.P;              % total population

% clear evolution data except initial condition
Data_s.C(2:end) = 0;
Data_s.D(2:end) = 0;
Data_s.V1(2:end) = 0;
Data_s.V2(2:end) = 0;

% response function (inhibitor hill function) set to 1 for t > tl0
g = @(X,k,t) 1/(1+((k.w*(t < k.tl0) )*X(5))^(k.n));
% uptake function (promoter hill function) set to 0 for t < tv0 Sv2 + Av2* + Rv2* + Dv2*
h = @(X,k,t) (t >= k.tv0)*(((k.wC*(X(5)+X(6))+k.wD*X(7)+k.wV*(X(15)+X(19)+X(20)+X(21)))^(k.nv))/(1+((k.wC*(X(5)+X(6))+k.wD*X(7)+k.wV*(X(15)+X(19)+X(20)+X(21)))^(k.nv)))); % currently a simple switch on day tv0 (init of vax rollout.)

% set up model
model = struct();

% dimensions 
model.k = k; % parameters
model.M = 28; % number of interactions
model.N = 21; % number of sub-populations

% state update matrices (M x N)
model.nu_minus = zeros(model.M, model.N);
model.nu_plus = zeros(model.M, model.N);
%    i = 1  2  3  4  5   6   7   8   9   10  11  12   13   14   15  16  17  18  19   20   21
% X(i) = Su Eu Iu Ru Au* Ru* Du* Sv1 Ev1 Iv1 Rv1 Av1* Rv1* Dv1* Sv2 Ev2 Iv2 Rv2 Av2* Rv2* Dv2*
model.nu_minus(1,[1,3]) = 1;    model.nu_plus(1,[2,3]) = 1;    % Su + Iu -> Iu + Eu
model.nu_minus(2,[1,10]) = 1;   model.nu_plus(2,[2,10]) = 1;   % Su + Iv1 -> Iv1 + Eu
model.nu_minus(3,[1,17]) = 1;   model.nu_plus(3,[2,17]) = 1;   % Su + Iv2 -> Iv2 + Eu
model.nu_minus(4,[2]) = 1;      model.nu_plus(4,[3]) = 1;      % Eu -> Iu
model.nu_minus(5,[3]) = 1;      model.nu_plus(5,[4]) = 1;      % Iu -> Ru
model.nu_minus(6,[3]) = 1;      model.nu_plus(6,[5]) = 1;      % Iu -> Au*
model.nu_minus(7,[5]) = 1;      model.nu_plus(7,[6]) = 1;      % Au* -> Ru*
model.nu_minus(8,[5]) = 1;      model.nu_plus(8,[7]) = 1;      % Au* -> Du*
model.nu_minus(9,[1]) = 1;      model.nu_plus(9,[8]) = 1;      % Su -> Sv1
model.nu_minus(10,[4]) = 1;     model.nu_plus(10,[11]) = 1;    % Ru -> Rv1
model.nu_minus(11,[8,3]) = 1;   model.nu_plus(11,[9,3]) = 1;   % Sv1 + Iu -> Iu + Ev1
model.nu_minus(12,[8,10]) = 1;  model.nu_plus(12,[9,10]) = 1;  % Sv1 + Iv1 -> Iv1 + Ev1
model.nu_minus(13,[8,17]) = 1;  model.nu_plus(13,[9,17]) = 1;  % Sv1 + Iv2 -> Iv2 + Ev1
model.nu_minus(14,[9]) = 1;     model.nu_plus(14,[10]) = 1;    % Ev1 -> Iv1
model.nu_minus(15,[10]) = 1;    model.nu_plus(15,[11]) = 1;    % Iv1 -> Rv1
model.nu_minus(16,[10]) = 1;    model.nu_plus(16,[12]) = 1;    % Iv1 -> Av1*
model.nu_minus(17,[12]) = 1;    model.nu_plus(17,[13]) = 1;    % Av1* -> Rv1*
model.nu_minus(18,[12]) = 1;    model.nu_plus(18,[14]) = 1;    % Av1* -> Dv1*
model.nu_minus(19,[8]) = 1;     model.nu_plus(19,[15]) = 1;    % Sv1 -> Sv2
model.nu_minus(20,[11]) = 1;    model.nu_plus(20,[18]) = 1;    % Rv1 -> Rv2
model.nu_minus(21,[15,3]) = 1;  model.nu_plus(21,[16,3]) = 1;  % Sv2 + Iu -> Iu + Ev2
model.nu_minus(22,[15,10]) = 1; model.nu_plus(22,[16,10]) = 1; % Sv2 + Iv1 -> Iv1 + Ev2
model.nu_minus(23,[15,17]) = 1; model.nu_plus(23,[16,17]) = 1; % Sv2 + Iv2 -> Iv2 + Ev2
model.nu_minus(24,[16]) = 1;    model.nu_plus(24,[17]) = 1;    % Ev2 -> Iv2
model.nu_minus(25,[17]) = 1;    model.nu_plus(25,[18]) = 1;    % Iv2 -> Rv2
model.nu_minus(26,[17]) = 1;    model.nu_plus(26,[19]) = 1;    % Iv2 -> Av2*
model.nu_minus(27,[19]) = 1;    model.nu_plus(27,[20]) = 1;    % Av2* -> Rv2*
model.nu_minus(28,[19]) = 1;    model.nu_plus(28,[21]) = 1;    % Av2* -> Dv2*

model.nu = model.nu_plus - model.nu_minus;

% initial state (assume V1 = V2 = 0)
Au0 = k.zeta*(Data_s.C(1)-Data_s.D(1));
Ru0 = (1-k.zeta)*(Data_s.C(1)-Data_s.D(1));
model.X0 = [Data_s.P - Au0 - ceil(3*k.kappa*Au0) - Ru0-Data_s.D(1); % Su0
           ceil(2*k.kappa*Au0); % Eu0
           ceil(k.kappa*Au0); % Iu0
           0; %Ru0
           Au0; %Au0*
           Ru0; %Ru0*
           Data_s.D(1); %Du0*
           zeros(model.N - 7,1)]; % all vax population = 0

% hazard functions
model.a = @(X,k,t) [(k.alpha0 + k.alpha*g(X,k,t))*X(1)*X(3);
                     k.alphav1*(k.alpha0 + k.alpha*g(X,k,t))*X(1)*X(10);
                     k.alphav2*(k.alpha0 + k.alpha*g(X,k,t))*X(1)*X(17);
                     k.beta*X(2);
                     k.eta*X(3);
                     k.gamma*X(3);
                     k.rho*X(5);
                     k.delta*X(5);
                     k.v*h(X,k,t)*X(1);
                     k.v*h(X,k,t)*X(4);
                     k.alphav1*(k.alpha0 + k.alpha*g(X,k,t))*X(8)*X(3);
                     k.alphav1*k.alphav1*(k.alpha0 + k.alpha*g(X,k,t))*X(8)*X(10);
                     k.alphav2*k.alphav1*(k.alpha0 + k.alpha*g(X,k,t))*X(8)*X(17);
                     k.beta*X(9);
                     k.rhov1*k.eta*X(10);
                     k.gammav1*k.gamma*X(10);
                     k.rhov1*k.rho*X(12);
                     k.deltav1*k.delta*X(12);
                     k.omega*X(8);
                     k.omega*X(11);
                     k.alphav2*(k.alpha0 + k.alpha*g(X,k,t))*X(15)*X(3);
                     k.alphav2*k.alphav1*(k.alpha0 + k.alpha*g(X,k,t))*X(15)*X(10);
                     k.alphav2*k.alphav2*(k.alpha0 + k.alpha*g(X,k,t))*X(15)*X(17);
                     k.beta*X(16);
                     k.rhov2*k.eta*X(17);
                     k.gammav2*k.gamma*X(17);
                     k.rhov2*k.rho*X(19);
                     k.deltav2*k.delta*X(19)];

% number of simulated days
T = length(Data_s.C)-1;

% run simulation
tau = 1.0; % timestep (day resolution) 
[Z,t] = TauLeapingMethod(model,T,tau);
% actual observables (init post vax use efficacy) 
for i =1:T
    [J] = find(t <= i);
    Data_s.C(1+i) = Z(5,J(end)) + Z(12,J(end)) + Z(19,J(end)) + Z(6,J(end)) + Z(13,J(end)) + Z(20,J(end)) + Z(7,J(end)) + Z(14,J(end)) + Z(21,J(end)); 
    Data_s.D(1+i) = Z(7,J(end)) + Z(14,J(end)) + Z(21,J(end));
    Data_s.V1(1+i) = Z(8,J(end)) + Z(12,J(end)) + Z(13,J(end)) + Z(14,J(end));
    Data_s.V2(1+i) = Z(15,J(end)) + Z(19,J(end)) + Z(20,J(end)) + Z(21,J(end));
end
%append initial zeros (to correctly match full time-series)
Data_s.C = [zeros(start-1,1); Data_s.C];
Data_s.D = [zeros(start-1,1); Data_s.D];
Data_s.V1 = [zeros(start-1,1); Data_s.V1];
Data_s.V2 = [zeros(start-1,1); Data_s.V2];

Data_s.Z = Z;
Data_s.t = t;
