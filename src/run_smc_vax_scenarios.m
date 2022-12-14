%% Bayesian analysis of vaccine hesitancy scenarios
%  Utilises data provided by from various sources (including Johns-Hopkins University
%   DXY, National Health Commisson (China), Protezione Civile (Italy), 
%   and https://www.corona-data.ch/ (Switzerland)).
%  
% Simulation based model extending a stochastic SIR model with latent infectous
% populatio and regulatory mechnisms.
%
%  Inference of model parameters is perfromed using Approximate Bayesian Computation
%  Within the Sequential Monte Carlo framework of Drovandi and Pettit (2011).
%  
% Authors:
%     David J. Warne (david.warne@qut.edu.au)
%           School of Mathematical Sciences
%           Science and Engineering Faculty 
%           Queensland University of Technology
%     Abhishek Varghese 
%           Centre for Data Science
%           Queensland University of Technology
%
%     Christopher Drovandi (c.drovandi@qut.edu.au)
%           School of Mathematical Sciences
%           Science and Engineering Faculty 
%           Queensland University of Technology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Algorithm Initialisation
%
% Initialise Random number generator for reproducibility
%clear global
rng(1337,'twister')
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex');


%% model fixed parameters
% population
Data.P = 67e6;
M = 12; % Roughly 12 months of simulation time
Data.C = zeros(M*30,1); 
Data.D = zeros(M*30,1); 
Data.V1 = zeros(M*30,1); 
Data.V2 = zeros(M*30,1); 


% Initial data estimated from worldometer.info/coronovirus/uk
% cumulative cases of 337,798 on 1st Sep 2020 in UK
Data.C(1,1) = 337798; 
% cumulative deaths init 41,551 on 1st Sep 2020 in UK
Data.D(1,1) = 41551;

% Transmission and detection model
theta = zeros(21,1);
theta(1) = 0.05;    % residual transmission
theta(2) = 0.4;     % transmission rate 
theta(3) = 0.07;    % symptom onset rate recovery rate (1/2 wks = 0.07)
theta(4) = 0.05;    % case detection rate 
theta(5) = 0.03;    % case death rate
theta(6) = 0.05;    % case recovery rate
theta(7) = 0.1;     % unobserved recovery rate

% vaccination effect model (AstraZeneca)
theta(8) = 0.01;    % max vaccination rate (1st dose) 
theta(9) = 0.693;   % transmission reduction (1st dose) ( AZ vs delta)
theta(10) = 1.0;    % detection reduction (1st dose) 
theta(11) = 1.0;    % recovery increase (1st dose)
theta(12) = 0.08;   % death decrease (1st dose)
theta(13) = 0.048;  % rate to 2nd dose (for AZ 1/21 days = 0.048 )
theta(14) = 0.33;   % transmission reduction (2nd dose) (AZ vs delta)
theta(15) = 1.0;    % detection reduction (2nd dose) 
theta(16) = 1.0;    % recovery increase (2nd dose) 
theta(17) = 0.047;  % death decrease (2nd dose)

% initial constion variables for latent states
theta(18) = 10;     % initial under-estimate factor ( I0 = 10*A0)

% chosen to replicate worldometeer cumulative recoveries of 25,1022 on 1st Sep 2020 in UK
theta(19) = 0.152659773769861; % (proportion of living cases that are active)

% NPI response functions
theta(20) = 10;         % inibitive strength ... like a hill constant
theta(21) = 1/30000;    %  A = 1/w_A -> toggle threshold
theta(22) = 90;         % day that everything opens up (NPI's disabled)

% Hesitancy effect function (example for scenario 1) 
theta(23) = 4;          % promotor strength
theta(24) = 1/1e6;      % w_C
theta(25) = 0;      % w_D
theta(26) = 1/1e5;      % w_V
theta(27) = 90;         % the day vaccine is available

%% Synthetic data generation for identifiability test
Dataset = simuldata_reg_fA_vax_h(Data,theta);

%% set-up ABC-SMC 
% the number of particles
N = 1000; 
% target ABC tolerance. (If zero, use acceptance probability stopping criteria)
epsilon_final = 0; 
% tuning parameters for ABC-SMC -- set to good initial defaults
a = 0.5; 
c = 0.01;
% minimum acceptance probability used as stopping criteria if espilon_final = 0
% if p_acc_min = 0 then use epsilon_final
p_acc_min = 0.05;

% define prior for theta = [nu,nv,-log(wC),-log(wD),-log(wV)]
prior.num_params = 5;
prior.p1 = [0.0, 0,  3, 3,   3];
prior.p2 = [0.02, 10, 8,  8,  8];
prior.sampler = @() [unifrnd(prior.p1,prior.p2)]; 
prior.pdf = @(theta) prod([unifpdf(theta,prior.p1,prior.p2)]);
prior.trans_f = @(theta)[theta];
prior.trans_finv = @(theta) [theta];

% user functions
% stochastic simulation given data template and parameters
sim_func = @(D,params) simuldata_reg_fA_vax_h(D,[theta(1:7)',params(1),theta(9:22)',params(2),10.^(-params(2:4)),theta(27)]);
%discrepancy metric
dist_func = @(S_d,S_s) sqrt(sum((S_d(:) - S_s(:)).^2));
% summary statistic
smry_func = @smry;

% prior labels for plots (set this for plot labels to print correctly
prior_names = {'n_v','-log(w_C)','-log(w_D)','-log(w_V)'};
%% run SMC sampler
j = 1;
start_time = tic;
[part_vals, part_sim, part_s, sims,epsilon_t,p_acc_t] = smc_abc_rw_par(Datasets{j},...
                                       sim_func,dist_func,smry_func,prior,N,...
                                       epsilon_final,a,c,p_acc_min);
experiment_run_time = toc;
%% save results
save(filename,'results-1.mat');

