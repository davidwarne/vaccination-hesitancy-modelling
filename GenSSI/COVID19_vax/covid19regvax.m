function model = covid19vax()
  
    
    % Symbolic variables
	syms Su Eu Iu Ru Auo Ruo Duo Sv1 Ev1 Iv1 Rv1 Av1o Rv1o Dv1o Sv2 Ev2 Iv2 Rv2 Av2o Rv2o Dv2o
	syms nu nv wC wD wV n wA alpha0 alpha1 beta1 eta gamma1 rho delta

    % Parameters to estimate
	model.sym.p = [nu;nv;wC;wD;wV];
    
    % known parameters
    %alpha0 = 0.05
    %alpha = 0.4;
    %beta = 0.07;
    %eta = 0.1;
    %gamma = 0.05;
    %rho = 0.05;
    %delta = 0.03;
    omega = 0.048;
    alphav1 = 0.693;
    alphav2 = 0.33;
    deltav1 = 0.08;
    deltav2 = 0.047;
    P = 67e6;
    Su0 = 64529649;
    Eu0 = 121079;
    Iu0 = 40203;
    Ru0 = 1314369;
    Auo0 = 23406;
    Ruo0 = 675428;
    Duo0 = 295867;
    %wA = 1/30000;
    %n = 10;
    % State variables
    model.sym.x = [Su;Eu;Iu;Ru;Auo;Ruo;Duo;Sv1;Ev1;Iv1;Rv1;Av1o;Rv1o;Dv1o;Sv2;Ev2;Iv2;Rv2;Av2o;Rv2o;Dv2o];
    size(model.sym.x)
    % Control vectors (g)
    model.sym.g = zeros(size(model.sym.x));
%%  
    % Autonomous dynamics (f)
    %%
	model.sym.xdot = [-nu*((wC*(Auo+Av1o+Av2o+Ruo+Rv1o+Rv2o)+wD*(Duo+Dv1o+Dv2o)+wV*(Sv2+Ev2+Iv2+Rv2+Av2o+Rv2o+Dv2o)^nv)/(1+(wC*(Auo+Av1o+Av2o+Ruo+Rv1o+Rv2o)+wD*(Duo+Dv1o+Dv2o)+wV*(Sv2+Ev2+Iv2+Rv2+Av2o+Rv2o+Dv2o))^nv))*Su-(alpha0+alpha1*(1/(1+(wA*(Auo+Av1o+Av2o))^n)))*((Iu + alphav1*Iv1+alphav2*Iv2)/P)*Su;
                       (alpha0+alpha1*(1/(1+(wA*(Auo+Av1o+Av2o))^n)))*((Iu + alphav1*Iv1+alphav2*Iv2)/P)*Su-beta1*Eu;
                       beta1*Eu-(eta+gamma1)*Iu;
                       eta*Iu-nu*((wC*(Auo+Av1o+Av2o+Ruo+Rv1o+Rv2o)+wD*(Duo+Dv1o+Dv2o)+wV*(Sv2+Ev2+Iv2+Rv2+Av2o+Rv2o+Dv2o)^nv)/(1+(wC*(Auo+Av1o+Av2o+Ruo+Rv1o+Rv2o)+wD*(Duo+Dv1o+Dv2o)+wV*(Sv2+Ev2+Iv2+Rv2+Av2o+Rv2o+Dv2o))^nv))*Ru;
                       gamma1*Iu-(delta+rho)*Auo;
                       rho*Auo;
                       delta*Auo;
                       nu*((wC*(Auo+Av1o+Av2o+Ruo+Rv1o+Rv2o)+wD*(Duo+Dv1o+Dv2o)+wV*(Sv2+Ev2+Iv2+Rv2+Av2o+Rv2o+Dv2o)^nv)/(1+(wC*(Auo+Av1o+Av2o+Ruo+Rv1o+Rv2o)+wD*(Duo+Dv1o+Dv2o)+wV*(Sv2+Ev2+Iv2+Rv2+Av2o+Rv2o+Dv2o))^nv))*Su-alphav1*(alpha0+alpha1*(1/(1+(wA*(Auo+Av1o+Av2o))^n)))*((Iu+alphav1*Iv1+alphav2*Iv2)/P)*Sv1-omega*Sv1;
                       alphav1*(alpha0+alpha1*(1/(1+(wA*(Auo+Av1o+Av2o))^n)))*((Iu+alphav1*Iv1+alphav2*Iv2)/P)*Sv1-beta1*Ev1;
                       beta1*Ev1-(eta+gamma1)*Iv1;
                       eta*Iv1+nu*((wC*(Auo+Av1o+Av2o+Ruo+Rv1o+Rv2o)+wD*(Duo+Dv1o+Dv2o)+wV*(Sv2+Ev2+Iv2+Rv2+Av2o+Rv2o+Dv2o)^nv)/(1+(wC*(Auo+Av1o+Av2o+Ruo+Rv1o+Rv2o)+wD*(Duo+Dv1o+Dv2o)+wV*(Sv2+Ev2+Iv2+Rv2+Av2o+Rv2o+Dv2o))^nv))*Ru-omega*Rv1;
                       gamma1*Iv1-(deltav1*delta+rho)*Av1o;
                       rho*Av1o;
                       deltav1*delta*Av1o;
                       omega*Sv1-alphav2*(alpha0+alpha1*(1/(1+(wA*(Auo+Av1o+Av2o))^n)))*((Iu+alphav1*Iv1+alphav2*Iv2)/P)*Sv2;
                       alphav2*(alpha0+alpha1*(1/(1+(wA*(Auo+Av1o+Av2o))^n)))*((Iu + alphav1*Iv1+alphav2*Iv2)/P)*Sv2-beta1*Ev2;
                       beta1*Ev2-(eta+gamma1)*Iv2;
                       eta*Iv2+omega*Rv1;
                       gamma1*Iv2-(deltav2*delta+rho)*Av2o;
                       rho*Av2o;
                       deltav2*delta*Av2o];

    % Initial conditions
	model.sym.x0 = [Su0;Eu0;Iu0;Ru0;Auo0;Ruo0;Duo0;zeros(14,1)];

    % Observables
	model.sym.y = [Auo+Av1o+Av2o+Ruo+Rv1o+Rv2o;
                   Duo+Dv1o+Dv2o;
                   Sv1+Ev1+Iv1+Rv1+Av1o+Rv1o+Dv1o;
                   Sv2+Ev2+Iv2+Rv2+Av2o+Rv2o+Dv2o];
end