% runcovid19vax runs the structural identifiability analysis for the 
% 2nd version of the COVID-19 model described by
% 
%    Warne et al. (2022). Bayesian uncertainty quantification to identify 
%       of population level vaccine hesitancy behaviours, medRxiv.org

% Confirm execution
genssiStartup
%syms nu nv wC wD wV 
syms nu nv wC wD wV n wA alpha0 alpha1 beta1 eta gamma1 rho delta
options.reportCompTime = true;
%genssiMain('covid19vax',10,[nu;nv;wC;wD;wV],options);
genssiMain('covid19regvax',10,[nu;nv;wC;wD;wV;n;wA;alpha0;alpha1;beta1;eta;gamma1;rho;delta],options)