
%% General Characteristics
trials = 1e5;


GCdB = 0; % dBi
GFdB = 0; % dBi
fc = 3; % GHz - carrier frequency
waveLen = physconst('LightSpeed')/(fc*10^9);
%
etaLoSdB =1.6; etaLoS = db2pow(etaLoSdB);
etaNLoSdB= 23; etaNLoS= db2pow(etaNLoSdB);
%
NPB = -174; % dBm/Hz
BW = 10*10^6; % Hz
noisePow = 10^(NPB/10)*BW; % mW
%%
% R_th_C = 1.0;
% R_th_E = 0.5;
R_th_C = 0.5;
R_th_E = 0.025;
%% 3D location of nodes
% network area
L = 25;

% Load the mobility model
spacing = 1;
loc_indices = 1:spacing:spacing*N;
pos_C = cell2mat(struct2cell(load('pos_C.mat'))); pos_C = pos_C(loc_indices,:);
pos_E = cell2mat(struct2cell(load('pos_E.mat'))); pos_E = pos_E(loc_indices,:);
pos_U = cell2mat(struct2cell(load('pos_U.mat'))); pos_U = pos_U(loc_indices,:);

% the FC is located as the origin
XF = 0; YF = 0; ZF = 0;

% the UE-C
XC = pos_C(iter,1); YC = pos_C(iter,2); ZC = 0;

X_0 = 10.2947;
Y_0 = 11.7037;
Z_0 = 12.77;
% the UE-E
XE_0 = X_0;
YE_0 = Y_0;
ZE_0 = 0;

XE = XE_0 + pos_E(iter,1);
YE = YE_0 + pos_E(iter,2);
ZE = ZE_0;

% % the UAV 3D Cartesian Coordinates
XU_0 = X_0;
YU_0 = Y_0;
ZU_0 = Z_0;

XU = XU_0 + pos_U(iter,1);
YU = YU_0 + pos_U(iter,2);
ZU = ZU_0;
%%
%% Large-Scale Fading
% ------------- The G2G channel from CU-C to FC ------------- 
dCF = sqrt( (XC-XF)^2+(YC-YF)^2+(ZC-ZF)^2 );
LCF = db2pow( GCdB+GFdB-22.7-26*log10(fc)-36.7*log10(dCF) );

% ------------- The G2A channels from CU-C to UAV ------------- 
dCU = sqrt( (XC-XU)^2+(YC-YU)^2+(ZC-ZU)^2 );
phiCU = 180/pi*asin( abs(ZC-ZU)/dCU );
pCU = 1/(1+12.08*exp(-0.11*(phiCU-12.08)));
LCULoS = (waveLen/(4*pi*dCU))^2/etaLoS;
LCUNLoS= (waveLen/(4*pi*dCU))^2/etaNLoS;
BCU = rand(1,trials) < pCU;
LCU = LCULoS .*pCU + LCUNLoS.*(1-pCU);

% ------------- The G2A channels from CU-E to UAV ------------- 
dEU = sqrt( (XE-XU)^2+(YE-YU)^2+(ZE-ZU)^2 );
phiEU = 180/pi*asin( abs(ZE-ZU)/dEU );
pEU = 1/(1+12.08*exp(-0.11*(phiEU-12.08)));
LEULoS = (waveLen/(4*pi*dEU))^2/etaLoS;
LEUNLoS= (waveLen/(4*pi*dEU))^2/etaNLoS;
BEU = rand(1,trials) < pEU;
LEU = LEULoS .*pEU + LEUNLoS.*(1-pEU);

% ------------- The G2A channels from UAV to FC ------------- 
dUF = sqrt( (XF-XU)^2+(YF-YU)^2+(ZF-ZU)^2 );
phiUF = 180/pi*asin( abs(ZF-ZU)/dUF );
pUF = 1/(1+12.08*exp(-0.11*(phiUF-12.08)));
LUFLoS = (waveLen/(4*pi*dUF))^2/etaLoS;
LUFNLoS= (waveLen/(4*pi*dUF))^2/etaNLoS;
BUF = rand(1,trials) < pUF;
LUF = LUFLoS .*pUF + LUFNLoS.*(1-pUF);

%% joint Microscopic and Macroscopic Effect
% ------------- The G2G channel from CU-C to FC ------------- 
mCF = 1; OmgCF = 0.279; bCF = 0.251;
Z = random('Nakagami',mCF,OmgCF,[1,trials]);
gCF = random('Rician',Z,sqrt(bCF),[1,trials]);
hCF2= gCF.^2 .* LCF; 

% ------------- The G2A channels from CU-C to UAV ------------- 
mCU = 2;
gCULoS = sqrt(gamrnd(mCU,1/mCU,[1,trials]));
gCUNLoS= sqrt(exprnd(1,[1,trials]));
hCU2= gCULoS.^2 .* LCULoS .*(BCU)...
    + gCUNLoS.^2.* LCUNLoS.*(1-BCU);

% ------------- The G2A channels from CU-E to UAV -------------
mEU = 5;
gEULoS = sqrt(gamrnd(mEU,1/mEU,[1,trials]));
gEUNLoS= sqrt(exprnd(1,[1,trials]));
hEU2= gEULoS.^2 .* LEULoS .*(BEU)...
    + gEUNLoS.^2.* LEUNLoS.*(1-BEU);

% ------------- The G2A channels from UAV to FC ------------- 
mUF = 3;
gUFLoS = sqrt(gamrnd(mUF,1/mUF,[1,trials]));
gUFNLoS= sqrt(exprnd(1,[1,trials]));
hUF2= gUFLoS.^2 .* LUFLoS .*(BUF)...
    + gUFNLoS.^2.* LUFNLoS.*(1-BUF);

%% Residual Interference Normalized Power
% xi = db2pow(-10);
xi = db2pow(-30);
thCU2 = exprnd(xi,[1,trials]).*LCU;
thCF2 = exprnd(xi,[1,trials]).*LCF;
thUF2 = exprnd(xi,[1,trials]).*LUF;
%% Achievable Rates
R_F_xC1 = @(PC) (1/2)*log2( 1+PC*hCF2/noisePow );
R_ADMUC_U_xC1= @(PC,PE) (1/2)*log2( 1+(PC*hCU2)./(PE*hEU2+noisePow) );
R_ADMUC_U_xE = @(PC,PE) (1/2)*log2( 1+(PE*hEU2)./(PC*thCU2+noisePow));
R_ADMUE_U_xE = @(PC,PE) (1/2)*log2( 1+(PE*hEU2)./(PC*hCU2+noisePow) );

R_ADMFC_F_xC2= @(PC,PU) (1/2)*log2( 1+(PC*hCF2)./(PU*hUF2+noisePow) );
R_ADMFC_F_hxE= @(PC,PU) (1/2)*log2( 1+(PU*hUF2)./(PC*thCF2+noisePow));
R_ADMFE_F_hxE= @(PC,PU) (1/2)*log2( 1+(PU*hUF2)./(PC*hCF2+noisePow) );
R_ADMFE_F_xC2= @(PC,PU) (1/2)*log2( 1+(PC*hCF2)./(PU*thUF2+noisePow));

R_F_xC2 = @(PC) (1/2)*log2( 1+PC*hCF2/noisePow );