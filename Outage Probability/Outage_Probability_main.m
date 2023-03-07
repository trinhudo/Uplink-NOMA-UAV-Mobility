% Paper: Adaptive Decoding Mechanisms for UAV-enabled Double-Uplink Coordinated NOMA
% Authors: Thanh Luan Nguyen, Georges Kaddoum, Tri Nhu Do, Daniel Benevides da Costa, Zygmunt J. Haas
% Journal: IEEE Transactions on Vehicular Technology
% DOI: 10.1109/TVT.2023.3255001

%% Initialize System Parameters
close all
clear all

theta1 = 0.8;
theta2 = 0.8;

%
PmaxdB = 0:5:40;
%
for iPmax = 1:length(PmaxdB)
    Pmax = db2pow(PmaxdB(iPmax));

    PC1 = theta1*Pmax;
    PE1 = (1-theta1)*Pmax;

    PC2 = theta2*Pmax;
    PU2 = (1-theta2)*Pmax;

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
    OP_ADM_E_SIM(iPmax) = 1 - (1-OP_U_xE) * (1-OP_F_hXE);

    % Theorem 2:
    OP_U_xC2 = mean((~SF_CU)&(EF_CU))+mean(((~SF_UC)|(~SF_C))&(EF_UC));
    OP_F_xC2 = mean(~SF_C2);
    OP_ADM_C2_SIM(iPmax) = (1-OP_U_xE)*OP_U_xC2 + OP_U_xE*OP_F_xC2;

    % Theorem 3:
    OP_ADM_C1_SIM(iPmax) = 1 - mean(SF_C1);
end

for iPmax = 1:length(PmaxdB)
    Pmax = db2pow(PmaxdB(iPmax));

    PC1 = theta1*Pmax;
    PE1 = (1-theta1)*Pmax;

    PC2 = theta1*Pmax;
    PU2 = (1-theta1)*Pmax;
    %
    %% Analytical
    tau_th_C = 2^(2*R_th_C)-1;
    tau_th_E = 2^(2*R_th_E)-1;

    a1 = tau_th_C/(PC1/noisePow);
    a2 = tau_th_E/(PE1/noisePow);
    alp1= (PE1/noisePow)*a1;
    alp2= (PC1/noisePow)*a2;
    A1 = a1/(1-alp1);
    A2 = a2/(1-alp2);
    A = (A1-a2)/alp2;

    % P1
    K0  = 1;
    chi0= [pCU 1-pCU]; mu0 = [mCU 1];
    Omg0= (waveLen/(4*pi*dCU))^2./[mCU*etaLoS etaNLoS];

    K1  = 1;
    chi1= [pEU 1-pEU]; mu1 = [mEU 1];
    Omg1= (waveLen/(4*pi*dEU))^2./[mEU*etaLoS etaNLoS];

    K2  = 1; V2 = ones(1,K2+1);
    chi2= [1/2 1/2]; mu2 = [1 1];
    Omg2= xi*LCU*V2;

    P1 = EqP1(K0,chi0,Omg0,mu0,...
        K1,chi1,Omg1,mu1,...
        K2,chi2,Omg2,mu2,...
        a1,a2,alp1,alp2,A1,A2,A);

    % P2
    K0  = 1;
    chi0= [pEU 1-pEU]; mu0 = [mEU 1];
    Omg0= (waveLen/(4*pi*dEU))^2./[mEU*etaLoS etaNLoS];

    K1  = 1;
    chi1= [pCU 1-pCU]; mu1 = [mCU 1];
    Omg1= (waveLen/(4*pi*dCU))^2./[mCU*etaLoS etaNLoS];

    P2 = EqP2(K0,chi0,Omg0,mu0,...
        K1,chi1,Omg1,mu1,...
        a2,alp2,A2);

    % P3
    b1 = tau_th_C/(PC2/noisePow);
    b2 = tau_th_E/(PU2/noisePow);
    bet1= (PU2/noisePow)*b1;
    bet2= (PC2/noisePow)*b2;
    B1 = b1/(1-bet1);
    B2 = b2/(1-bet2);
    B = (B1-b2)/bet2;
    hat_B= (B2-b1)/bet1;

    a1= b1; A1= B1; alp1 = bet1;
    a2= b2; A2= B2; alp2 = bet2; A = B;

    alp = 1/(2*bCF)*(2*bCF*mCF/(2*bCF*mCF+OmgCF))^mCF;
    bet = 1/(2*bCF);
    del = OmgCF/(2*bCF*(2*bCF*mCF+OmgCF));
    zet = @(l) (-1).^l .* pochhammer(1-mCF,l) ./ factorial(l) .* del.^l;

    K0  = mCF-1; V0 = ones(1,K0);
    chi0= alp*[zet(1:K0)./(bet-del).^(2:mCF) 1/(bet-del)];
    Omg0= [LCF/(bet-del)*V0 LCF/(bet-del)];
    mu0 = [(2:mCF) 1];

    K1  = 1;
    chi1= [pUF 1-pUF]; mu1 = [mUF 1];
    Omg1= (waveLen/(4*pi*dUF))^2./[mUF*etaLoS etaNLoS];

    K2  = 1;
    chi2= [1/2 1/2];
    Omg2= [xi*LCF xi*LCF];
    mu2 = [1 1];

    P3 = EqP1(K0,chi0,Omg0,mu0,...
        K1,chi1,Omg1,mu1,...
        K2,chi2,Omg2,mu2,...
        a1,a2,alp1,alp2,A1,A2,A);

    % P4 
    K0  = 1;
    chi0= [pUF 1-pUF]; mu0 = [mUF 1];
    Omg0= (waveLen/(4*pi*dUF))^2./[mUF*etaLoS etaNLoS];

    K1  = mCF-1; V1 = ones(1,K1);
    chi1= alp*[zet(1:K1)./(bet-del).^(2:mCF) 1/(bet-del)];
    Omg1= [LCF/(bet-del)*V1 LCF/(bet-del)];
    mu1 = [(2:mCF) 1];

    P4 = EqP2(K0,chi0,Omg0,mu0,...
        K1,chi1,Omg1,mu1,...
        a2,alp2,A2);
    %
    OP_ADM_E_ANA(iPmax) = 1 - (P1+P2)*(P3+P4);

    % P5 
    K0  = 1;
    chi0= [pUF 1-pUF]; mu0 = [mUF 1];
    Omg0= (waveLen/(4*pi*dUF))^2./[mUF*etaLoS etaNLoS];

    K1  = mCF-1; V1 = ones(1,K1);
    chi1= alp*[zet(1:K1)./(bet-del).^(2:mCF) 1/(bet-del)];
    Omg1= [LCF/(bet-del)*V1 LCF/(bet-del)];
    mu1 = [(2:mCF) 1];

    K2  = 1; V2 = ones(1,K2+1);
    chi2= [1/2 1/2]; mu2 = [1 1];
    Omg2= xi*LUF*V2;

    a1 = b2; alp1 = bet2; A1 = B2;
    a2 = b1; alp2 = bet1; A2 = B1; A = hat_B;

    P5 = EqP1(K0,chi0,Omg0,mu0,...
        K1,chi1,Omg1,mu1,...
        K2,chi2,Omg2,mu2,...
        a1,a2,alp1,alp2,A1,A2,A);

    % P6 
    K0  = mCF-1; V0 = ones(1,K0);
    chi0= alp*[zet(1:K0)./(bet-del).^(2:mCF) 1/(bet-del)];
    Omg0= [LCF/(bet-del)*V0 LCF/(bet-del)];
    mu0 = [(2:mCF) 1];

    K1  = 1;
    chi1= [pUF 1-pUF]; mu1 = [mUF 1];
    Omg1= (waveLen/(4*pi*dUF))^2./[mUF*etaLoS etaNLoS];

    P6 = EqP2(K0,chi0,Omg0,mu0,...
        K1,chi1,Omg1,mu1,...
        a2,alp2,A2);

    % P7
    K0  = mCF-1; V0 = ones(1,K0);
    chi0= alp*[zet(1:K0)./(bet-del).^(2:mCF) 1/(bet-del)];
    Omg0= [LCF/(bet-del)*V0 LCF/(bet-del)];
    mu0 = [(2:mCF) 1];

    P7 = 0;
    for k = 1:(K0+1)
        P7 = P7 + chi0(k) * gammainc( b1/Omg0(k),mu0(k),'upper' );
    end

    OP_ADM_C2_ANA(iPmax) = (P1+P2)*(1-P5-P6)+(1-P1-P2)*(1-P7);

    % P8
    a1 = tau_th_C/(PC1/noisePow);
    P8 = 0;
    for k = 1:(K0+1)
        P8 = P8 + chi0(k) * gammainc( a1/Omg0(k),mu0(k),'upper' );
    end
    P8 = 1-P8;

    OP_ADM_C1_ANA(iPmax) = P8;
end

figure;
semilogy(PmaxdB, OP_ADM_E_ANA, 'r-'); hold on;
semilogy(PmaxdB, OP_ADM_E_SIM, 'rs:'); hold on;

semilogy(PmaxdB, OP_ADM_C1_ANA, 'g-'); hold on;
semilogy(PmaxdB, OP_ADM_C1_SIM, 'g^:'); hold on;

semilogy(PmaxdB, OP_ADM_C2_ANA, 'b-'); hold on;
semilogy(PmaxdB, OP_ADM_C2_SIM, 'bd:'); hold on;

xlabel('$P_{\max}$ [dBm]' ,'Interpreter', 'LaTex')
ylabel('OP');
legend('Ana. $\mathrm{OP}_{\mathsf{E}, \mathrm{e2e}}$',...
    'Sim. $\mathrm{OP}_{\mathsf{E}, \mathrm{e2e}}$',...
    'Ana. $\mathrm{OP}_{\mathsf{C}, \mathrm{e2e}}^{[1]}$',...
    'Sim. $\mathrm{OP}_{\mathsf{C}, \mathrm{e2e}}^{[1]}$',...
    'Ana $\mathrm{OP}_{\mathsf{C}, \mathrm{e2e}}^{[2]}$',...
    'Sim. $\mathrm{OP}_{\mathsf{C}, \mathrm{e2e}}^{[2]}$',...
    'Interpreter', 'LaTex', 'Location', 'southwest')

% leg.ItemTokenSize = [10,18];
% grid on;
% set(gca,'Fontsize',10);
% set(gcf,'Position',[100 100 300 300]);
% set(gca,'LooseInset',get(gca,'TightInset'));