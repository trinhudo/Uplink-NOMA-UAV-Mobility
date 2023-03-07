% Paper: Adaptive Decoding Mechanisms for UAV-enabled Double-Uplink Coordinated NOMA
% Authors: Thanh Luan Nguyen, Georges Kaddoum, Tri Nhu Do, Daniel Benevides da Costa, Zygmunt J. Haas
% Journal: IEEE Transactions on Vehicular Technology
% DOI: 10.1109/TVT.2023.3255001

%% Initialize System Parameters

clear all

N = 1e2;

theta1 = 0.8;
theta2 = 0.8;

PmaxdB = 0;
Pmax = db2pow(PmaxdB);

tic;
%
PC1 = theta1*Pmax;
PE1 = (1-theta1)*Pmax;

PC2 = theta2*Pmax;
PU2 = (1-theta2)*Pmax;
%
for iter = 1:N
    system_param_per_iter; %
    %
    SF_C1 = R_F_xC1(PC1) > R_th_C;
    
    EU_CE = hCU2 > hEU2;
        SU_CE= R_ADMUC_U_xC1(PC1,PE1) > R_th_C;
        SU_E = R_ADMUC_U_xE(PC1,PE1) > R_th_E;
    EU_EC = hCU2 < hEU2;
        SU_EC = R_ADMUE_U_xE(PC1,PE1) > R_th_E;

    EF_CU = hCF2 > hUF2;
        SF_CU= R_ADMFC_F_xC2(PC2,PU2) > R_th_C;
        SF_U = R_ADMFC_F_hxE(PC2,PU2) > R_th_E;
    EF_UC = hCF2 < hUF2;
        SF_UC= R_ADMFE_F_hxE(PC2,PU2) > R_th_E;
        SF_C = R_ADMFE_F_xC2(PC2,PU2) > R_th_C;
    
    SF_C2 = R_F_xC2(PC2) > R_th_C;
    %% Simulation ADM
    % Theorem 1:
    OP_U_xE = mean(((~SU_CE)|(~SU_E))&(EU_CE)) + mean((~SU_EC)&(EU_EC));
    OP_F_hXE= mean(((~SF_CU)|(~SF_U))&(EF_CU)) + mean((~SF_UC)&(EF_UC));
    OP_ADM_hxE(iter) = 1 - (1-OP_U_xE) * (1-OP_F_hXE);
    
    % Theorem 2:
    OP_U_xC2 = mean((~SF_CU)&(EF_CU))+mean(((~SF_UC)|(~SF_C))&(EF_UC));
    OP_F_xC2 = mean(~SF_C2);
    OP_ADM_xC2(iter) = (1-OP_U_xE)*OP_U_xC2 + OP_U_xE*OP_F_xC2;
    
    % Theorem 3:
    OP_xC1(iter) = 1 - mean(SF_C1);
    %% NADM
    
    OP_UC_xE = 1-mean(SU_CE & SU_E); % OP at U when i=C
    OP_UE_xE = 1-mean(SU_EC); % OP at U when i=E
    
    % Non-Silent UAV
    OP_FC_xC2_noSilent = 1-mean(SF_CU); % OP of F to decode xC2 when i=C
    OP_FE_xC2_noSilent = 1-mean(SF_UC & SF_C); % OP of F to decode xC2 when i=E
    
    OP_FC_hxE = 1-mean(SF_CU & SF_U); % OP of F to decode hxE when i=C
    OP_FE_hxE = 1-mean(SF_UC); % OP of F to decode hxE when i=E
    
    % Silent UAV
    OP_FC_hxE_Silent = 1;
    OP_FE_xC2_Silent = 1-mean(SF_C2);
    
    % Simulation NADM-CC
    OP_NADM_CC_xE(iter) = 1 - (1-OP_UC_xE)*(1-OP_FC_hxE);
    OP_NADM_CC_xC2(iter)= 1 - (1-OP_UC_xE)*(1-OP_FC_xC2_noSilent) - OP_UC_xE*(1-OP_FE_xC2_Silent);
    
    % Simulation NADM-EC
    OP_NADM_EC_xE(iter) = 1 - (1-OP_UE_xE)*(1-OP_FC_hxE);
    OP_NADM_EC_xC2(iter)= 1 - (1-OP_UE_xE)*(1-OP_FC_xC2_noSilent) - OP_UE_xE*(1-OP_FE_xC2_Silent);
    
    % Simulation nADM3
    OP_NADM_CE_xE(iter) = 1 - (1-OP_UC_xE)*(1-OP_FE_hxE);
    OP_NADM_CE_xC2(iter)= 1 - (1-OP_UC_xE)*(1-OP_FE_xC2_noSilent) - OP_UC_xE*(1-OP_FE_xC2_Silent);
    
    % Simulation nADM4
    OP_NADM_EE_xE(iter) = 1 - (1-OP_UE_xE)*(1-OP_FE_hxE);
    OP_NADM_EE_xC2(iter)= 1 - (1-OP_UE_xE)*(1-OP_FE_xC2_noSilent) - OP_UE_xE*(1-OP_FE_xC2_Silent);
end
toc;

%
Throughput_ADM = R_th_E/2*( 1-OP_ADM_hxE )...
    + R_th_C/2*( 1-OP_ADM_xC2 ) + R_th_C/2*( 1-OP_xC1 );
%
Throughput_NADM_CC = R_th_E/2*( 1-OP_NADM_CC_xE )...
    + R_th_C/2*( 1-OP_NADM_CC_xC2 ) + R_th_C/2*( 1-OP_xC1 );
Throughput_NADM_EC = R_th_E/2*( 1-OP_NADM_EC_xE )...
    + R_th_C/2*( 1-OP_NADM_EC_xC2 ) + R_th_C/2*( 1-OP_xC1 );
Throughput_NADM_CE = R_th_E/2*( 1-OP_NADM_CE_xE )...
    + R_th_C/2*( 1-OP_NADM_CE_xC2 ) + R_th_C/2*( 1-OP_xC1 );
Throughput_NADM_EE = R_th_E/2*( 1-OP_NADM_EE_xE )...
    + R_th_C/2*( 1-OP_NADM_EE_xC2 ) + R_th_C/2*( 1-OP_xC1 );
%
%
blue1=[0, 0.4470, 0.7410];
orange1=[0.8500, 0.3250, 0.0980];
yellow1=[0.9290, 0.6940, 0.1250];
purple1=[0.4940, 0.1840, 0.5560];
green1=[0.4660, 0.6740, 0.1880];
cyan1=[0.3010, 0.7450, 0.9330];
red1=[0.6350, 0.0780, 0.1840];

color_sky = [173 235 255]/255;


figure;
plot(1:N,Throughput_ADM,'-m','markersize',5,'LineWidth',1); hold on;
plot(1:N,Throughput_NADM_CC,'-g','markersize',4,'LineWidth',1); hold on;
plot(1:N,Throughput_NADM_EC,'-b','markersize',4,'LineWidth',1); hold on;
plot(1:N,Throughput_NADM_CE,'-r','markersize',4,'LineWidth',1); hold on;
plot(1:N,Throughput_NADM_EE,'-c','markersize',4,'LineWidth',1); hold on;

xlabel('Location Index');
ylabel('Throughput [bits/s/Hz]');

leg = legend('ADM','NADM, ${d}_{\rm CC}$','NADM, ${d}_{\rm EC}$',...
        'NADM, ${d}_{\rm CE}$','NADM, ${d}_{\rm EE}$',...
        'NumColumns',2,'Interpreter','LaTex');

% leg.ItemTokenSize = [15,18];   
% grid on;
% set(gca,'Fontsize',10);
% set(gcf,'Position',[100 100 300 300]);
% set(gca,'LooseInset',get(gca,'TightInset'));