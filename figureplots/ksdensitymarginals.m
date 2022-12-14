function [h,ax] = ksdensitymatrix(h,ax,X,w,lims,res,labels,linspec,colorspec,ref,tol)
%% [h, ax] = ksdensitymatrix(h,ax,X,w,lims,res,labels,linspec,colorspec)
% Generates a kernel smoothing density plot matrix of an M-dimensional 
% multivariate probability density function from N samples.
%
% Inputs:
%    h          - the current figure handle, set h = [] to create a new figure.
%    ax         - axes list handles, ignored if h = []. 
%    X          - N x M array of multivariate samples, X(i,:) is the i-th sample
%                 and X(:,j) is all the samples of the j-th dimension.
%    w          - sample importance weights use empty vector [] to indicate 
%                 uniform weights.
%    lims       - 2 x M list of plot limits for each dimension lims(1,:) is 
%                 lower limits, and lims(2,:) is the upper limits. Defaults 
%                 to min and max of samples if user provides an empty array.
%    labels     - an (M+1) x 1 cell array of labels for each dimension 
%                 and the dataset label 
%    linspec    - line properties for plotting. e.g., '--' or ':'
%    colorspec  - 3 x 1 array pf colour properties form plotting.
%
% Outputs:
%    h          - the figure handle used or created
%    ax         - axes list handles for h. Use h and ax again as inputs 
%                 to overlay results.
%
% Note: Assumes LaTeX Interpreter, use the following before first call,
%       set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
%       set(groot, 'defaultLegendInterpreter','latex');
%       set(groot, 'defaulttextInterpreter','latex');
%
% Author: David J. Warne (david.warne@qut.edu.au)
%         School of Mathematical Sciences
%         Queensland University of Technology
%
% Date: Sep 2018
%
% Example usage:
% X1 = mvnrnd([1,0],[1,0;0,1],1000);
% X2 = mvnrnd([1,2],[1,0.5;0.5,1],1000);
% [h,ax] = ksdensitymatrix([],[],X1,[],[],300,{'X','Y','\mathcal{D}'},'-',[0,0,1]);
% [h,ax] = ksdensitymatrix(h,ax,X2,[],[],300,{'X','Y','\mathcal{D}'},'--',[1,0,0]);
%------------------------------------------------------------------------------

[N,M] = size(X); % N number of samples, M data dimension

[rho,pval] = corr(X,'Type','Spearman');

% set weights to uniform if not provided
if isempty(w)
    w = ones(N,1);
end

% create a new figure if h and ax are not provided
if isempty(h)
    h = figure;
    %for i=1:M
        for j=1:M
            ax(j) = subplot(1,M,j);
        end
    %end
end

% set limits to min and max sample if not provided
if isempty(lims)
    lims = zeros(2,M);
    lims(1,:) = min(X);
    lims(2,:) = max(X);
end

%for i=1:M
    for j=1:M
        axes(ax(j));
        box on
        hold on;
        % draw a univariate marginal
       
            % build grid
            pts = linspace(lims(1,j),lims(2,j),res);
            % generate kernel density estimate
            [p,xi] = ksdensity(X(:,j),pts,'Weights',w);
            % plot the result
            ha = area(xi,p,'LineStyle',linspec,'FaceColor',colorspec,'FaceAlpha',0.2,'LineWidth',1);
            %scatter(xpts,ypts,10,'filled','MarkerFaceAlpha',0.2,'MarkerFaceColor',corcol);
            plot([ref(j),ref(j)],[0,max(p)],'--r');
            xlim([lims(1,j),lims(2,j)]);
            ylabel(['$p(',labels{j},' | ',labels{M+1},')$'])
            xlabel(['$',labels{j},'$']);
      
        hold off;
    end
end
