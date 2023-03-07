% Paper: Adaptive Decoding Mechanisms for UAV-enabled Double-Uplink Coordinated NOMA
% Authors: Thanh Luan Nguyen, Georges Kaddoum, Tri Nhu Do, Daniel Benevides da Costa, Zygmunt J. Haas
% Journal: IEEE Transactions on Vehicular Technology
% DOI: 10.1109/TVT.2023.3255001

%% 

clear all
N = 1e2;   % [sample] 

blue1  =[0, 0.4470, 0.7410];
orange1=[0.8500, 0.3250, 0.0980];
yellow1=[0.9290, 0.6940, 0.1250];
purple1=[0.4940, 0.1840, 0.5560];
green1 =[0.4660, 0.6740, 0.1880];
cyan1  =[0.3010, 0.7450, 0.9330];
red1   =[0.6350, 0.0780, 0.1840];
%%
R = 10;    % [m] Movement Radius of UEE
Rd = R/10; % [m] Maximum Deviation from the projection of UAV to UEE

pos_C = cell2mat(struct2cell(load('pos_C.mat')));
pos_E = cell2mat(struct2cell(load('pos_E.mat')));
pos_U = cell2mat(struct2cell(load('pos_U.mat'))); 

XF0 = 0; YF0 = 0; ZF0 = 0;
XC0 = 0; YC0 = 0; ZC0 = 0;

XE0 = 10.2947; YE0 = 11.7037; ZE0 = 0; 
XU0 = 10.2947; YU0 = 11.7037; H = 12.77; ZU0 = H;
D0 = sqrt(XE0^2+YE0^2); T0 = atan(YU0/XU0);

XF = 0; YF = 0; ZF = 0;
XC = pos_C(:,1); YC = pos_C(:,2); ZC = 0;
XE = XE0+pos_E(:,1); YE = YE0+pos_E(:,2); ZE = 0;
XU = XU0+pos_U(:,1); YU = YU0+pos_U(:,2); ZU = H;

%% Plot Captured Location
NN = 1:20:2500;

t = tiledlayout(1,2);

nexttile;
plot3(XF,YF,ZF,'^k','MarkerFaceColor','b'); hold on;
plot3(XE(NN),YE(NN),ZE*ones(length(NN),1),'.-','Color',blue1); hold on;
plot3(XC(NN),YC(NN),ZC*ones(length(NN),1),'.-','Color',red1); hold on;
plot3(XU(NN),YU(NN),ZU*ones(length(NN),1),'.-',...
    'Color',green1,'LineWidth',0.5,'MarkerSize',4); hold on;
grid on;
box on;
leg = legend('FC','UE-E','UE-C','UAV','Location','Best',...
    'Orientation','Horizontal','NumColumns',2);
leg.ItemTokenSize = [10,5];

xlabel('$X$ [m]','Interpreter','LaTex');
ylabel('$Y$ [m]','Interpreter','LaTex');
zlabel('$Z$ [m]','Interpreter','LaTex');
title('Captured Position','Interpreter','LaTex');

nexttile;
plot3(XF,YF,ZF,'^k','MarkerFaceColor','b'); hold on;
plot3(XE(NN),YE(NN),ZE*ones(length(NN),1),'.-','Color',blue1); hold on;
plot3(XC(NN),YC(NN),ZC*ones(length(NN),1),'.-','Color',red1); hold on;
plot3(XU(NN),YU(NN),ZU*ones(length(NN),1),'.-',...
    'Color',green1,'LineWidth',0.5,'MarkerSize',4); hold on;
grid on;
box on;
leg = legend('FC','UE-E','UE-C','UAV','Location','Best',...
    'Orientation','Horizontal','NumColumns',2);
leg.ItemTokenSize = [10,5];

xlabel('$X$ [m]','Interpreter','LaTex');
ylabel('$Y$ [m]','Interpreter','LaTex');
zlabel('$Z$ [m]','Interpreter','LaTex');
title('Projection to Oxy','Interpreter','LaTex');
view(2)

t.TileSpacing = 'none';
t.Padding = 'none';