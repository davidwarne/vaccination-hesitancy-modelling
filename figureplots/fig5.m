%% Demonstration of vaccination effect without any hesitancy

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex');

% initial data
Data.P = 67e6;
M = 12;
Data.C = zeros(M*30,1); 
Data.D = zeros(M*30,1); 
Data.V1 = zeros(M*30,1); 
Data.V2 = zeros(M*30,1); 
Data.C(1,1) = 337798; 
Data.D(1,1) = 41551;

% model parameters
% Transmission and detection model
theta = zeros(21,1);
theta(1) = 0.05; 
theta(2) = 0.4;   
theta(3) = 0.07;
theta(4) = 0.05; 
theta(5) = 0.03;    
theta(6) = 0.05; 
theta(7) = 0.1;  

% vaccination effect model
theta(8) = 0.01;  
theta(9) = 0.693; 
theta(10) = 1.0; 
theta(11) = 1.0;
theta(12) = 0.08; 
theta(13) = 0.048;
theta(14) = 0.33; 
theta(15) = 1.0;   
theta(16) = 1.0; 
theta(17) = 0.047;  

% initial constion variables for latent states
theta(18) = 10;      
theta(19) = 0.152659773769861; 

% NPI response and vaccine uptake functions
theta(20) = 10;        
theta(21) = 1/30000; 
theta(22) = Inf;  

theta(23) = 4; 
theta(24) = 1;  
theta(25) = 1; 
theta(26) = 1;
theta(27) = Inf;
% strategy -- for model demonstration
scen = [10,1/500000,Inf,Inf;10,1/500000,90,Inf;10,1/500000,Inf,90;10,1/500000,90,90];

M = 4;
figure; 
lim_y = [2e7,15e6];
for i=1:4
    theta(20) = scen(i,1);
    theta(21) = scen(i,2);
    theta(22) = scen(i,3);
    theta(27) = scen(i,4);
    sims = cell(M,1);
    for j=1:M
        sims{j} = simuldata_reg_fA_vax_h(Data,theta);
    end
    
    subplot(2,2,i);
    Data_s = sims{1};
    plot([theta(22)/30,theta(22)/30],[0,1e8],'--k','Linewidth',2);
    hold on
    plot([theta(27)/30,theta(27)/30],[0,1e8],':','Color',[67,170,139]/255,'Linewidth',2); 
    for j=1:M
        Data_s = sims{j};
        plot(Data_s.t/30,Data_s.Z(2,:),':b','Color',[85,153,255]/255,'Linewidth',2);
        plot(Data_s.t/30,Data_s.Z(3,:),'--b','Color',[85,153,255]/255,'Linewidth',2);
        plot(Data_s.t/30,Data_s.Z(5,:)+ Data_s.Z(12,:) + Data_s.Z(19,:),'-','Color',[85,153,255]/255,'Linewidth',2);
        plot(Data_s.t/30,Data_s.C,'-','Color',[255,176,52]/255,'Linewidth',2);
        plot(Data_s.t/30,Data_s.D,'-','Color',[255,128,128]/255,'Linewidth',2);
        plot(Data_s.t/30,Data_s.V1,'--','Color',[67,170,139]/255,'Linewidth',2);
        plot(Data_s.t/30,Data_s.V2,'-','Color',[144,190,109]/255,'Linewidth',2);
    end
    ylim([0,lim_y(mod(i,2)+1)])
    legend('$T_d$', '$T_v$','$E$', '$I$','$A^*$','$C^*$','$D^*$','$V_1^*$','$V_2^*$');
    xlabel('Time $t$ (months)');
    ylabel('Counts');
    xlim([0,12]);
end 
