function [part_vals, part_sim, part_s, sims,dist_t,p_acc_t] = smc_abc_rw(D,...
                                                 sim_func,dist_func,...
                                                 smry_func,prior,N,dist_final,...
                                                 a,c,p_acc_min)
%% Adaptive sequential Monte Carlo for approximate Bayesian computation.
%  
% Parameters: 
%   D         - data struction assumed to contain fields required by user functions.
%   sim_func  - user defined function that simulates the model of interest given
%              a vector of parameter values.
%   dist_func - distance function/discrepancy metric for comaring simulated and 
%               true data.
%   smry_func - summary statistic function for mapping data to a lower dimensional
%               space.
%   prior     - data struction for the prior distribution of the parameters. 
%               Requires the following fields. 
%               
%               num_params - the dimensionality of the parameter space;
%               sampler    - a sampling function that generates random vectors from
%                            the joint prior distribution;
%               pdf        - the joint probability density function;
%               trans_f    - transform of prior paramete space to ensure 
%                            unbounded support for MCMC sampling.
%               trans_finv - inverse of transform.
%  N          - number of particles for SMC sampler.
%  dist_final - target discrepancy threshold. If zero, then p_acc_min is used to
%               determine stopping criteria.
%  a          - tuning parameter for adaptive selection of discrepancy threshold 
%               sequence. 
%  c          - tuning parameter for choosing the number of MCMC iterations in 
%               move step.
%  p_acc_min  - minimum acceptable acceptance rate in the MCMC interations. If
%               zero the dist_final is used to determine stopping criteria. 
%
% Returns:
%    part_val  - parameter values for each particle.
%    part_sim  - summary statistics for each particle.
%    part_c    - discrepancy metric value for each particle.
%    sim       - total number of model simulations performed.
%    dist_t    - smallest discrepacy threshold reached.
%    p_acc_min - smallest MCMC acceptance rate reached.
%
%% Authors:
%     Christopher Drovandi (c.drovandi@qut.edu.au)
%           School of Mathematical Sciences
%           Faculty of Science 
%           Queensland University of Technology
%
%     David J. Warne (david.warne@qut.edu.au)
%           School of Mathematical Sciences
%           Faculty of Science 
%           Queensland University of Technology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute summary statistics for observations
part_obs = smry_func(D); 

% values for adaptive steps
num_drop = floor(N*a);
num_keep = N-num_drop;
mcmc_trials = 5;

% Initialise particle data structures
part_vals = zeros(N,prior.num_params);
part_s = zeros(N,1);
part_sim = zeros(N, length(part_obs));
tic()
% initial prior rejection algorithm
for i = 1:N
    % sample prior 
    part_vals(i,:) = prior.sampler();
    % simulate model
    D_s = sim_func(D,part_vals(i,:));
    % compute summary stats
    part_sim(i,:) = smry_func(D_s);
    % evaluate the discrepancy metric    
    part_s(i) = dist_func(part_obs,part_sim(i,:));
end
toc()
sims = N;

% transform the parameters
for i=1:N
    part_vals(i,:) = prior.trans_f(part_vals(i,:));
end

% sort the particles
[part_s,ix] = sort(part_s); 
part_vals = part_vals(ix,:); 
part_sim = part_sim(ix,:);

% determine next disprepacy threshold
dist_max = part_s(N);
dist_next = part_s(num_keep);
dist_t = dist_next;
p_acc_t = 0;

% interate toward target discrepancy
while (dist_max > dist_final)
    tic()
    % compute the covariance matrix (of particles that remain) required 
    % for the Independent MH move step
    cov_matrix = (2.38^2)*cov(part_vals(1:num_keep,:))/size(part_vals,2);
    
    % resample
    r = randsample(num_keep, N-num_keep, 'true');
    part_vals((num_keep+1):N, :) = part_vals(r,:);
    part_s((num_keep+1):N) = part_s(r);
    part_sim((num_keep+1):N, :) = part_sim(r,:);
    
    i_acc = zeros(N-num_keep,1);
    sims_mcmc = zeros(N-num_keep,1);
    
    % trial MCMC steps 
    for i = (num_keep+1):N
        for r = 1:mcmc_trials
            % Gaussian random walk
            part_vals_prop = mvnrnd(part_vals(i,:),cov_matrix);
            prior_curr = prior.pdf(prior.trans_finv(part_vals(i,:)));
            prior_prop = prior.pdf(prior.trans_finv(part_vals_prop));
            % early rejection (assumes symmetric proposals)
            if (isnan(prior_prop/prior_curr) || rand > prior_prop./prior_curr)  
                continue;
            end
            
            prop = prior.trans_finv(part_vals_prop);
            D_s = sim_func(D,prop);
            part_sim_prop = smry_func(D_s);
            dist_prop = dist_func(part_obs,part_sim_prop);
            
            sims_mcmc(i-num_keep) = sims_mcmc(i-num_keep)+1;
            % ABC part of the acceptance probability
            if (dist_prop <= dist_next) 
                % then the metropolis-hastings ratio is satisfied
                part_vals(i,:) = part_vals_prop; 
                part_s(i) = dist_prop; 
                part_sim(i,:) = part_sim_prop;
                i_acc(i-num_keep) = i_acc(i-num_keep) + 1;
            end
        end
    end
    
    % determine number of MCMC iterations to perfrom 
    acc_rate = sum(i_acc)/(mcmc_trials*(N-num_keep));
    mcmc_iters =   floor(log(c)/log(1-acc_rate)+1);
    fprintf('Total number of mcmc moves for current target is %d, number remaining is %d\n',mcmc_iters,mcmc_iters-mcmc_trials);
    
    % move step
    for i = (num_keep+1):N
        for r = 1:(mcmc_iters-mcmc_trials)
            % Gaussian random walk
            part_vals_prop = mvnrnd(part_vals(i,:),cov_matrix);
            prior_curr = prior.pdf(prior.trans_finv(part_vals(i,:)));
            prior_prop = prior.pdf(prior.trans_finv(part_vals_prop));
            % early rejection (assumes symmetric proposal)
            if (isnan(prior_prop/prior_curr) || rand > prior_prop./prior_curr) 
                continue;
            end
            
            prop = prior.trans_finv(part_vals_prop);
            D_s = sim_func(D,prop);
            part_sim_prop = smry_func(D_s);
            
            dist_prop = dist_func(part_obs,part_sim_prop);
            
            sims_mcmc(i-num_keep) = sims_mcmc(i-num_keep)+1;
            
            % ABC part of the acceptance probability
            if (dist_prop <= dist_next) 
                % then the metropolis-hastings ratio is satisfied
                part_vals(i,:) = part_vals_prop; 
                part_s(i) = dist_prop; 
                part_sim(i,:) = part_sim_prop;
                i_acc(i-num_keep) = i_acc(i-num_keep) + 1;
            end
        end
    end
    
    num_mcmc_iters = max(0, mcmc_iters - mcmc_trials) + mcmc_trials;
    p_acc = sum(i_acc)/(num_mcmc_iters*(N-num_keep));
    
    fprintf('MCMC acceptance probability was %f\n',p_acc);
    
    sims = sims + sum(sims_mcmc);
    mcmc_trials = ceil(mcmc_iters/2);
    
    % compute number of unique particles
    fprintf('The number of unique particles is %d\n',length(unique(part_vals(:,1))));
    
    % compute the next distance and maximum distance
    % sort the particles
    [part_s,ix] = sort(part_s); part_vals = part_vals(ix,:); part_sim = part_sim(ix,:);
    
    % if most of the particles are under the final target then don't
    % drop them
    if (sum((part_s > dist_final)) < num_drop)
        num_drop = sum((part_s > dist_final));
        num_keep = N-num_drop;
    end
   
    % to return information about convergence
    dist_t = dist_next;
    p_acc_t = p_acc;

    dist_max = part_s(N);
    dist_next = part_s(num_keep);
    
    fprintf('The next distance is %f and the maximum distance is %f and the number to drop is %d\n',dist_next,dist_max,num_drop);
    fprintf('The number of sims is %d\n',sims);
   
    if (p_acc < p_acc_min)
        fprintf('Getting out as MCMC acceptance rate is below acceptable threshold\n');
        break;
    end
    toc()
end
for i=1:N
    part_vals(i,:) = prior.trans_finv(part_vals(i,:));
end
            
