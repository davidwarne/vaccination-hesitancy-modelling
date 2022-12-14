%% Plot joint posterior distribution estimates
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex');
prior_names =  {'\nu','n_v','-\log(w_C)','-\log(w_D)','-\log(w_V)','\mathcal{D}'};

colourmap = [255,255,204;
            27,158,119;
            217,95,2;
            117,112,179]/255;
theta_true = [0.01,4,1/1e6,0,1/1e5;
        0.01,4,1/1e6,0,0;
        0.01,4,0,1/5e5,1/1e5;
        0.01,4,0,1/5e5,0;
        0.01,4,1/2e6,1/1e6,1/1e5;
        0.01,4,1/2e6,1/1e6,0;
        0.005,4,1,1,1;
        0.002,4,1,1,1];

    theta_true(:,3:5) = -log10(theta_true(:,3:5));

Sigma = cell(8,2);
Spec = cell(8,4);
for j=1:8
    col = lines(4);

    load(['results-',num2str(j),'.mat'])
    [h,ax] = ksdensitymatrix([],[],part_vals,[],[],100,prior_names,'-',col(1,:),theta_true(j,:),0.01)
    print(['post-examples-scen-',num2str(j)],'-dsvg','-painters')
end
