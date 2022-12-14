%% Plot example synthetic datasets 

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

% strategy -- for model demonstration
scen = [0.01,4,1,1,1;0.01,4,1/2e6,1/1e6,0;0.01,4,1/2e6,1/1e6,1/1e5]
    
theta(22) = 90;
theta(27) = 90;
figure;
for i=1:3
    theta(8) = scen(i,1);
    theta(23) = scen(i,2)
    theta(24) = scen(i,3); 
    theta(25) = scen(i,4); 
    theta(26) = scen(i,5);
    Data_s = simuldata_reg_fA_vax_h(Data,theta);
    subplot(1,3,i);
    plot([theta(22)/30,theta(22)/30],[0,1e8],'--k','Linewidth',2);
    hold on;
    plot([theta(27)/30,theta(27)/30],[0,1e8],':','Color',[67,170,139]/255,'Linewidth',2); 
    plot([0,Data_s.t(end)/30],[1/theta(24),1/theta(24)],':','Color',[255,176,52]/255);
    plot([0,Data_s.t(end)/30],[1/theta(25),1/theta(25)],':','Color',[255,128,128]/255,'Linewidth',2);
    plot([0,Data_s.t(end)/30],[1/theta(26),1/theta(26)],':','Color',[144,190,109]/255,'Linewidth',2);
    plot(Data_s.t/30,Data_s.C,'-','Color',[255,176,52]/255,'Linewidth',2);
    plot(Data_s.t/30,Data_s.D,'-','Color',[255,128,128]/255,'Linewidth',2);
    plot(Data_s.t/30,Data_s.V1,'--','Color',[67,170,139]/255,'Linewidth',2);
    plot(Data_s.t/30,Data_s.V2,'-','Color',[144,190,109]/255,'Linewidth',2);
    xlabel('Time $t$ (months)')
    ylabel('Counts');
    legend('$T_d$','$T_v$','$1/w_C$','$1/w_D$','$1/w_V$','$C^*$','$D^*$','$V_1^*$','$V_2^*$');
    ylim([0,3.5e7]);
end    
