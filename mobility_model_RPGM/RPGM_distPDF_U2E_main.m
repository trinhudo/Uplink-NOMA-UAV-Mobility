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

%% Lemma 6: Distance from C to F
dCF = sqrt((XC-XF).^2 + (YC-YF).^2 + (ZC-ZF).^2);

P3 = @(RR,r) 3*(RR^2-r.^2)/(257*pi).*(189-44*r.^2/RR^2-18*r.^4/RR^4); CP = @(RR) RR^4;

f_dCF = @(r) (2*pi*r.*P3(R,r)/CP(R)).*(r<R);
%% Lemma 7: Distance from E to F
dEF = sqrt((XE-XF).^2 + (YE-YF).^2 + (ZE-ZF).^2);

C_1 = @(D,RR) 33*D^4 + 54*D^2*RR^2 - 214*RR^4;
C_2 = @(D,RR) 114*D^2 + 54*RR^2;
C_3 = @(D,RR) 18*D^6 + 26*D^4*RR^2 - 233*D^2*RR^4 + 189*RR^6;
C_4 = @(D,RR) 162*D^2 + 26*RR^2;
C_5 = @(D,RR) 162*D^4 + 104*D^2*RR^2 - 233*RR^4;

f_1 = @(D,RR,r) 6*r/(257*pi*RR^8).*( -sqrt( (r.^2-(D-RR)^2).*((D+RR)^2-r.^2) )...
    .* ( C_1(D,RR) + 33*r.^4 + C_2(D,RR)*r.^2 ) + ( C_3(D,RR) + 18*r.^6 ...
        + C_4(D,RR)*r.^4 + C_5(D,RR)*r.^2 ).*acos((D^2+r.^2-RR^2)./(2*D*r)) );
f_2 = @(D,RR,r) 6*r/(257*RR^8).*( C_3(D,RR) + 18*r.^6 + C_4(D,RR)*r.^4 + C_5(D,RR)*r.^2 );

f_dEF= @(r) (f_1(D0,R,r).*(abs(D0-R) <= r).*(r <= D0+R)...
    + f_2(D0,R,r).*(r <= R-D0).*(D0 < R)).*( (max(D0-R,0)<r).*(r<D0+R) );
%% Lemma 8: Distance from E to U
dEU = sqrt((XE-XU).^2 + (YE-YU).^2 + (ZE-ZU).^2);

f_dEU= @(r) (2*pi*r.*P3(Rd,sqrt(r.^2-H^2))/CP(Rd)).*(H<r).*(r<sqrt(Rd^2+H^2));
%% Lemma 9: Distance from U to F
dUF = sqrt((XU-XF).^2 + (YU-YF).^2 + (ZU-ZF).^2);

f_dPUF= @(r) f_dEF(r);
f_dUF = @(r) (r./sqrt(r.^2-H^2).*f_dPUF(sqrt(r.^2-H^2)))...
    .* (sqrt((D0-R)^2+H^2)<r).*(r<sqrt((D0+R)^2+H^2));
%% Lemma 10: Distance from C to U
dCU = sqrt((XU-XC).^2 + (YU-YC).^2 + (ZC-ZU).^2);

C0 = 29.177*1e-6;
P0 = @(r) C0*(4-r.^2).^4.*(376 - 101*r.^2 + 48*r.^4 - 11*r.^6 + r.^8);

tilde_P0 = @(r,t) P0( sqrt(r.^2+D0^2-2*D0*r.*cos(t-T0))/R );

g_1 = @(r) integral(@(t) tilde_P0(r,t)...
        .* (T0-acos(-(4*R^2-r.^2-D0^2)./(2*r*D0)) < t)...
        .* (T0+acos(-(4*R^2-r.^2-D0^2)./(2*r*D0)) > t),-pi,pi,'ArrayValued',true);
g_2 = @(r) integral(@(t) tilde_P0(r,t),-pi,pi,'ArrayValued',true);

f_dCU = @(r) r./(2*pi*R^2).*g_1(sqrt(r.^2-H^2)).*(sqrt((D0-2*R)^2+H^2) <= r).*(r <= sqrt((D0+2*R)^2+H^2))...
           + r./(2*pi*R^2).*g_2(sqrt(r.^2-H^2)).*(r < sqrt((max(2*R-D0,0))^2+H^2));
%% Plot Captured Location
NN = 1:20:2500;

t = tiledlayout(3,3);

nexttile(5,[1 2]);
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
title('Captured Location','Interpreter','LaTex');

nexttile(8,[1 2]);
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
%% Plot PDF
NN = 100; L = 15;
% ====================== Distance from C to F =============================
xCF = linspace(0,R,NN);
yCF = f_dCF(xCF);

nexttile(1);
histogram(dCF,'normalization','pdf','NumBins',L); hold on; 
plot(xCF,yCF,'LineWidth',2);
leg = legend('Sim','Ana');
leg.ItemTokenSize = [10,5];
xlabel('$r$','Interpreter','LaTex');
ylabel('PDF');
title('PDF of $d_{\sf \bf CF}$','Interpreter','LaTex');

% ====================== Distance from E to F =============================
xEF = linspace(max(D0-R,0),D0+R,NN);
yEF = f_dEF(xEF);

nexttile(2);
histogram(dEF,'normalization','pdf','NumBins',L); hold on; 
plot(xEF,yEF,'LineWidth',2);
leg = legend('Sim','Ana');
leg.ItemTokenSize = [10,5];
xlabel('$r$','Interpreter','LaTex');
ylabel('PDF');
title('PDF of $d_{\sf \bf EF}$','Interpreter','LaTex');
% ====================== Distance from E to U =============================
xEU = linspace(H,sqrt(Rd^2+H^2),NN);
yEU = f_dEU(xEU);

nexttile(3);
histogram(dEU,'normalization','pdf','NumBins',L); hold on; 
plot(xEU,yEU,'LineWidth',2);
leg = legend('Sim','Ana');
leg.ItemTokenSize = [10,5];
xlabel('$r$','Interpreter','LaTex');
ylabel('PDF');
title('PDF of $d_{\sf \bf EU}$','Interpreter','LaTex');
% ====================== Distance from U to F =============================
xUF = linspace(sqrt((D0-R)^2+H^2),sqrt((D0+R)^2+H^2),NN);
yUF = f_dUF(xUF);

nexttile(4);
histogram(dUF,'normalization','pdf','NumBins',L); hold on; 
plot(xUF,yUF,'LineWidth',2);
leg = legend('Sim','Ana');
leg.ItemTokenSize = [10,5];
xlabel('$r$','Interpreter','LaTex');
ylabel('PDF');
title('PDF of $d_{\sf \bf UF}$','Interpreter','LaTex');
% ====================== Distance from C to U =============================
xCU = linspace(H,sqrt((D0+2*R)^2+H^2),NN);
yCU = f_dCU(xCU);

nexttile(7);
histogram(dCU,'normalization','pdf','NumBins',L); hold on; 
plot(xCU,yCU,'LineWidth',2);
xlabel('$r$','Interpreter','LaTex');
ylabel('PDF');
title('PDF of $d_{\sf \bf UC}$','Interpreter','LaTex');

leg = legend('Sim','Ana');
leg.ItemTokenSize = [10,5];
t.TileSpacing = 'none';
t.Padding = 'none';