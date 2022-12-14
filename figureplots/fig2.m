%% Demonstration of response function
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
scen = [10,1/30000,Inf;3,1/30000,Inf;10,1/500000,Inf;3,1/500000,Inf];

M = 4;
figure;  
for i=1:2
    theta(20) = scen(i,1);
    theta(21) = scen(i,2);
    theta(22) = scen(i,3);
    
    sims = cell(M,1);
    for j=1:M
        sims{j} = simuldata_reg_fA_vax_h(Data,theta);
    end
    
    subplot(2,3,i);
    Data_s = sims{1};
    plot([min(Data_s.t),max(Data_s.t)]/30,[1/theta(21),1/theta(21)],':k','Linewidth',2); 
    hold on;
    for j=1:M
        Data_s = sims{j};
        plot(Data_s.t/30,Data_s.Z(2,:),':b','Color',[85,153,255]/255,'Linewidth',2);
        plot(Data_s.t/30,Data_s.Z(3,:),'--b','Color',[85,153,255]/255,'Linewidth',2);
        plot(Data_s.t/30,Data_s.Z(5,:)+ Data_s.Z(12,:) + Data_s.Z(19,:),'-','Color',[85,153,255]/255,'Linewidth',2);
        plot(Data_s.t/30,Data_s.C,'-','Color',[255,176,52]/255,'Linewidth',2);
        plot(Data_s.t/30,Data_s.D,'-','Color',[255,128,128]/255,'Linewidth',2);
    end
    legend('$1/w_A$', '$E$', '$I$','$A^*$','$C^*$','$D^*$');
    xlabel('Time $t$ (months)');
    ylabel('Counts');
    ylim([0,Inf]);
end 

for i=3:4
    theta(20) = scen(i,1);
    theta(21) = scen(i,2);
    theta(22) = scen(i,3);
    
    sims = cell(M,1);
    for j=1:M
        sims{j} = simuldata_reg_fA_vax_h(Data,theta);
    end
    
    subplot(2,3,i+1);
    Data_s = sims{1};
    plot([min(Data_s.t),max(Data_s.t)]/30,[1/theta(21),1/theta(21)],':k','Linewidth',2); 
    hold on;
    for j=1:M
        Data_s = sims{j};
        plot(Data_s.t/30,Data_s.Z(2,:),':b','Color',[85,153,255]/255,'Linewidth',2);
        plot(Data_s.t/30,Data_s.Z(3,:),'--b','Color',[85,153,255]/255,'Linewidth',2);
        plot(Data_s.t/30,Data_s.Z(5,:)+ Data_s.Z(12,:) + Data_s.Z(19,:),'-','Color',[85,153,255]/255,'Linewidth',2);
        plot(Data_s.t/30,Data_s.C,'-','Color',[255,176,52]/255,'Linewidth',2);
        plot(Data_s.t/30,Data_s.D,'-','Color',[255,128,128]/255,'Linewidth',2);
    end
    legend('$1/w_A$', '$E$', '$I$','$A^*$','$C^*$','$D^*$');
    xlabel('Time $t$ (months)');
    ylabel('Counts');
    ylim([0,Inf]);
end 

subplot(2,3,3);
col = lines(4)
for i=1:4
    g = @(A) (1+(scen(i,2).*A).^scen(i,1)).^-1;
    fplot(g,[0,2/min(scen(:,2))],'-','Color',col(i,:),'Linewidth',2);
    hold on;
end
plot([0,2/min(scen(:,2))],[0.5,0.5],'--k','Linewidth',2);
legend({'$n = 10, w_A = 3\times 10^{-4}$','$n = 3, w_A = 3\times 10^{-4}$','$n = 10, w_A = 5\times 10^{-5}$','$n = 3, w_A = 5\times 10^{-5}$','$g(1/w_A) = 1/2$'})
xlabel('$A^*$');
ylabel('$g(A^*)$')
