function [h,ax] = ksdensitymatrix(h,ax,X,w,lims,res,labels,linspec,colorspec,ref,alpha)
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
    h = figure('units','centimeters','Position',[1,1,25,20]);
    for i=1:M
        for j=1:M
            ax(i,j) = subplot(M,M,M*(i-1) + j);
        end
    end
end

% set limits to min and max sample if not provided
if isempty(lims)
    lims = zeros(2,M);
    lims(1,:) = min(X);
    lims(2,:) = max(X);
end

for i=1:M
    for j=i:M
        
        % draw a univariate marginal
        if i==j
            axes(ax(i,i));
            box on
            hold on;
            % build grid
            pts = linspace(lims(1,i),lims(2,i),res);
            % generate kernel density estimate
            [p,xi] = ksdensity(X(:,i),pts,'Weights',w);
            % plot the result
            ha = area(xi,p,'LineStyle',linspec,'FaceColor',colorspec,'FaceAlpha',0.2,'LineWidth',1);
            %scatter(xpts,ypts,10,'filled','MarkerFaceAlpha',0.2,'MarkerFaceColor',corcol);
            plot([ref(j),ref(j)],[0,max(p)],'--r');
            xlim([lims(1,i),lims(2,i)]);
            ylabel(['$p(',labels{i},' | ',labels{M+1},')$'])
            xlabel(['$',labels{i},'$']);
        else % draw a bivariate marginal
            axes(ax(i,j));
            box on
            hold on;
            % build grid
            xpts = X(:,j) % linspace(lims(1,j),lims(2,j),res);
            ypts = X(:,i) %linspace(lims(1,i),lims(2,i),res);
            %[Xp,Yp] = ndgrid(xpts,ypts);
            % generate kernel density estimate
            %[p,xy] = ksdensity(X(:,[j,i]),[Xp(:),Yp(:)],'PlotFcn','contour','Weights',w);
            if pval(i,j) <= alpha
                if rho(i,j) > 0
                    corcol = [150,101,104]./255;
                else 
                    corcol = [101,104,150]./255;
                end
            else
                corcol = [0.5,0.5,0.5];
            end
            scatter(xpts,ypts,10,'filled','MarkerFaceAlpha',0.2,'MarkerFaceColor',corcol);
            % place correlation and pval at mid-point
            %xm = (max(xpts) - min(xpts))/2 + min(xpts);
            %ym = (max(ypts) - min(ypts))/2 + min(ypts);
            %str = ['$\rho = ', num2str(rho(i,j)),'$, $p\rm{-val} = ',num2str(pval(i,j)),'$'];
            %text(xm,ym,str);
            %P = reshape(p,[res,res]);
            %xi = reshape(xy(:,1),[res,res]);
            %yi = reshape(xy(:,2),[res,res]);
            % plot the result as a contour plot
            %[c,h] = contour(xi,yi,P,4,linspec,'color',colorspec,'LineWidth',1);
            %view(2);
            xlim([lims(1,j),lims(2,j)]);
            ylim([lims(1,i),lims(2,i)]);
            ylabel(['$',labels{i},'$']);
            xlabel(['$',labels{j},'$']);
            axes(ax(j,i));
            box off
            hold on;
            if pval(i,j) > 0.05
                text(0.05,0.5,['$\rho = ', num2str(rho(i,j),2),'$ (ns)']);
            elseif pval(i,j) <= 0.001
                text(0.05,0.5,['$\rho = ', num2str(rho(i,j),2),'$ (***)']);
            elseif pval(i,j) <= 0.01
                text(0.05,0.5,['$\rho = ', num2str(rho(i,j),2),'$ (**)']);
            elseif pval(i,j) <= 0.05
                text(0.05,0.5,['$\rho = ', num2str(rho(i,j),2),'$ (*)']);
            end
            set(gca,'xtick',[])
            set(gca,'xticklabel',[])
            set(gca,'ytick',[])
            set(gca,'yticklabel',[])
        end
        hold off;
    end
end
