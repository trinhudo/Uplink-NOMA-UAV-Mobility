% Paper: Adaptive Decoding Mechanisms for UAV-enabled Double-Uplink Coordinated NOMA
% Authors: Thanh Luan Nguyen, Georges Kaddoum, Tri Nhu Do, Daniel Benevides da Costa, Zygmunt J. Haas
% Journal: IEEE Transactions on Vehicular Technology
% DOI: 10.1109/TVT.2023.3255001

% Generalized Random Waypoint Mobility (GRWP) model on a line
%   1. All nodes are mobile (p_s = 0, no static nodes)
%   2. The pause time is set to zero.

%% Parameters

clear all

N = 1e4;

R = 10;  % [m] Movement Radius of UEE
Rd = R/10;     % [m] Maximum Deviation from the projection of UAV to UEE 
VE_min = 0.1;  % [m/sec] minimum velocity
VE_max = 0.2;  % [m/sec] maximum
VC_min = 0.1;  % [m/sec] minimum velocity
VC_max = 0.2;  % [m/sec] maximum
VU_min = 0.01; % [m/sec] minimum velocity
VU_max = 0.05; % [m/sec] maximum
T = R*0.025;    % [sec] sampling period

%% Initialize Simulation Part 1
% node (source) location
radial_WE_1 = R*sqrt(rand); % Radial-Coord of Source Waypoint
angula_WE_1 = 2*pi*rand;    % Angula-Coord of Source Waypoint
WE_1 = [radial_WE_1*cos(angula_WE_1) radial_WE_1*sin(angula_WE_1)]; % Cartesian Coordinates of Source Waypoint

% waypoint (destination) location
radial_WE_2 = R*sqrt(rand); % Radial-Coord of Destination Waypoint
angula_WE_2 = 2*pi*rand;    % Angula-Coord of Destination Waypoint
WE_2 = [radial_WE_2*cos(angula_WE_2) radial_WE_2*sin(angula_WE_2)]; % Cartesian Coordinates of Destination Waypoint

% distance between WE_1 and WE_2: Length of (WE_1,WE_2) Leg
DE_W12 = norm(WE_1-WE_2);

% direction of movement from WE_1 to WE_2
angula_E_12 = atan2(WE_2(2)-WE_1(2),WE_2(1)-WE_1(1));
dir_E_12 = [cos(angula_E_12) sin(angula_E_12)];

% total distance traveled from WE_1 to WE_2
LE_12 = 0;

% movement velocity on (WE_1,WE_2) Leg
VE_12 = (VE_max-VE_min)*rand+VE_min;

% movement velocity vector from WE_1 to WE_2
vecVE_12 = VE_12*dir_E_12;

%% Initialize Simulation Part 2
% UAV's initial location: uniformly distributed around WE_1;
radial_WU_1 = Rd*sqrt(rand);
angula_WU_1 = 2*pi*rand;
WU_1 = [radial_WU_1*cos(angula_WU_1) radial_WU_1*sin(angula_WU_1)];

% UAV's destination location: uniformly distributed around WE_1;
radial_WU_2 = Rd; 
angula_WU_2 = 2*pi*rand;
offset_WU_2 = [radial_WU_2*cos(angula_WU_2)...
               radial_WU_2*sin(angula_WU_2)];
WU_2 = offset_WU_2;

% distance between WU_1 and WU_2: Length of (WU_1,WU_2) Leg
DU_W12 = norm(WU_1-WU_2);

% movement direction of UAV
angula_U_12 = atan2(WU_2(2)-WU_1(2),WU_2(1)-WU_1(1));           
dir_U_12 = [cos(angula_U_12) sin(angula_U_12)];

% total distance traveled from WE_1 to WE_2
LU_12 = 0;

% movement velocity of UAV
VU_12 = (VU_max-VU_min)*rand+VU_min;

% movement velocity vector of UAV
vecVU_12 = VU_12*dir_U_12;

%% Initialize Simulation Part 3
radial_WC_1 = R*sqrt(rand); % Radial-Coord of Source Waypoint
angula_WC_1 = 2*pi*rand;    % Angula-Coord of Source Waypoint
WC_1 = [radial_WC_1*cos(angula_WC_1) radial_WC_1*sin(angula_WC_1)]; % Cartesian Coordinates of Source Waypoint

% waypoint (destination) location
radial_WC_2 = R*sqrt(rand); % Radial-Coord of Destination Waypoint
angula_WC_2 = 2*pi*rand;    % Angula-Coord of Destination Waypoint
WC_2 = [radial_WC_2*cos(angula_WC_2) radial_WC_2*sin(angula_WC_2)]; % Cartesian Coordinates of Destination Waypoint

% distance between WC_1 and WC_2: Length of (WC_1,WC_2) Leg
DC_W12 = norm(WC_1-WC_2);

% direction of movement from WC_1 to WC_2
angula_C_12 = atan2(WC_2(2)-WC_1(2),WC_2(1)-WC_1(1));
dir_C_12 = [cos(angula_C_12) sin(angula_C_12)];

% total distance traveled from WC_1 to WC_2
LC_12 = 0;

% movement velocity on (WC_1,WC_2) Leg
VC_12 = (VC_max-VC_min)*rand+VC_min;

% movement velocity vector from WC_1 to WC_2
vecVC_12 = VC_12*dir_C_12;

%% position data collection
pos_C = WC_1;
pos_E = WE_1;
pos_U = WU_1 + WE_1;
pos_R = WU_1;
pos_W = [WE_1;WE_2];

%% =========================================================================
while(1)
    % =========================================================================
    if (LE_12+VE_12*T >= DE_W12) % If next captured Location in (WE_2,WE_3) Leg
        
        % ------------ UEE ------------------------------------------------
        % location of next destination waypoint
        radial_WE_3 = R*sqrt(rand); % Radial-Coord of Next Destination Waypoint
        angula_WE_3 = 2*pi*rand;    % Angula-Coord of Next Destination Waypoint
        WE_3 = [radial_WE_3*cos(angula_WE_3) radial_WE_3*sin(angula_WE_3)];
        
        % distance between WE_2 and WE_3: Length of (WE_2,WE_3) Leg
        DE_W23 = norm(WE_3-WE_2);
        
        % direction of movement from WE_2 to WE_3
        angula_E_23 = atan2(WE_3(2)-WE_2(2),WE_3(1)-WE_2(1));
        dir_E_23 = [cos(angula_E_23) sin(angula_E_23)];
        
        % movement velocity on (WE_2,WE_3) Leg
        VE_23 = (VE_max-VE_min)*rand+VE_min;
        
        % movement velocity vector from WE_2 to WE_3
        vecVE_23 = VE_23*dir_E_23;
        
        % total distance traveled from WE_2 to WE_3
        LE_23 = (T-(DE_W12-LE_12)/VE_12)*VE_23;
        
        % ------------- update position -----------------------------------
        pos_E = [pos_E; pos_E(end,:) + (DE_W12-LE_12)*dir_E_12 + LE_23*dir_E_23];
        
        % ------------ update parameters ----------------------------------
        WE_2 = WE_3;
        LE_12 = LE_23;
        VE_12 = VE_23;
        DE_W12= DE_W23;
        vecVE_12 = vecVE_23;
        dir_E_12 = dir_E_23;
        
        pos_W = [pos_W ; WE_3];
    else % The next Captured Location is still in the Current Leg
        pos_E = [pos_E; pos_E(end,:) + vecVE_12*T];
        
        % Calculate the total traveled distance
        LE_12 = LE_12 + VE_12*T;
    end
    % =========================================================================
    if (LU_12+VU_12*T >= DU_W12) % If next captured Location in (WE_2,WE_3) Leg
        % ------------ UAV ------------------------------------------------
        % UAV's next destination: uniformly distributed around WE_3;
        radial_WU_3 = Rd*sqrt(rand); % Radial-Coord of Destination Waypoint
        angula_WU_3 = 2*pi*rand;     % Angula-Coord of Destination Waypoint
        offset_WU_3 = [radial_WU_3*cos(angula_WU_3)...
                       radial_WU_3*sin(angula_WU_3)];
        WU_3 = offset_WU_3; % Cartesian Coordinates of Destination Waypoint

        % distance between WE_2 and WE_3: Length of (WE_2,WE_3) Leg
        DU_W23 = norm(WU_3-WU_2);

        % movement direction of UAV
        angula_U_23 = atan2(WU_3(2)-WU_2(2),WU_3(1)-WU_2(1));
        dir_U_23 = [cos(angula_U_23) sin(angula_U_23)];
        
        % movement velocity of UAV
        VU_23 = (VU_max-VU_min)*rand+VU_min;

        % movement velocity vector of UAV
        vecVU_23 = VU_23*dir_U_23;

        % total distance traveled from WE_2 to WE_3
        LU_23 = (T-(DU_W12-LU_12)/VU_12)*VU_23;
        
        % ------------- update position -----------------------------------
        pos_R = [pos_R; pos_R(end,:) + (DU_W12-LU_12)*dir_U_12 + LU_23*dir_U_23];
        pos_U = [pos_U; pos_R(end,:) + pos_E(end,:)];

        % ------------ update parameters ----------------------------------
        WU_2 = WU_3;
        LU_12 = LU_23;
        VU_12 = VU_23;
        DU_W12= DU_W23;
        vecVU_12 = vecVU_23;
        dir_U_12 = dir_U_23;
    else % The next Captured Location is still in the Current Leg
        pos_R = [pos_R; pos_R(end,:) + vecVU_12*T];
        pos_U = [pos_U; pos_R(end,:) + pos_E(end,:)];
        % Calculate the total traveled distance
        LU_12 = LU_12 + VU_12*T;
    end
    % =========================================================================
    if (LC_12+VC_12*T >= DC_W12) % If next captured Location in (WC_2,WC_3) Leg
        
        % ------------ UCC ------------------------------------------------
        % location of next destination waypoint
        radial_WC_3 = R*sqrt(rand); % Radial-Coord of Next Destination Waypoint
        angula_WC_3 = 2*pi*rand;    % Angula-Coord of Next Destination Waypoint
        WC_3 = [radial_WC_3*cos(angula_WC_3) radial_WC_3*sin(angula_WC_3)];
        
        % distance between WC_2 and WC_3: Length of (WC_2,WC_3) Leg
        DC_W23 = norm(WC_3-WC_2);
        
        % direction of movement from WC_2 to WC_3
        angula_C_23 = atan2(WC_3(2)-WC_2(2),WC_3(1)-WC_2(1));
        dir_C_23 = [cos(angula_C_23) sin(angula_C_23)];
        
        % movement velocity on (WC_2,WC_3) Leg
        VC_23 = (VC_max-VC_min)*rand+VC_min;
        
        % movement velocity vector from WC_2 to WC_3
        vecVC_23 = VC_23*dir_C_23;
        
        % total distance traveled from WC_2 to WC_3
        LC_23 = (T-(DC_W12-LC_12)/VC_12)*VC_23;
        
        % ------------- update position -----------------------------------
        pos_C = [pos_C; pos_C(end,:) + (DC_W12-LC_12)*dir_C_12 + LC_23*dir_C_23];
        
        % ------------ update parameters ----------------------------------
        WC_2 = WC_3;
        LC_12 = LC_23;
        VC_12 = VC_23;
        DC_W12= DC_W23;
        vecVC_12 = vecVC_23;
        dir_C_12 = dir_C_23;
    else % The next Captured Location is still in the Current Leg
        pos_C = [pos_C; pos_C(end,:) + vecVC_12*T];
        
        % Calculate the total traveled distance
        LC_12 = LC_12 + VC_12*T;
    end
    %% Stop if Collected Enough Samples of E moves from W1 to W2
    if (length(pos_E) >= N || length(pos_U) >= N || length(pos_C) >= N)
        break;
    end
end
%%
XE0 = pos_E(:,1);
YE0 = pos_E(:,2);
XC0 = pos_C(:,1);
YC0 = pos_C(:,2);
XU0 = pos_U(:,1);
YU0 = pos_U(:,2);

% save('pos_C.mat','pos_C');
% save('pos_E.mat','pos_E');
% save('pos_U.mat','pos_U');

NN = 1:10:N;
figure;
plot(YC0(NN), XC0(NN), '.-r'); hold on;
plot(-18.85+YE0(NN), -13.49+XE0(NN), '.-b'); hold on;
plot(-18.85+YU0(NN), -13.49+XU0(NN), '.-k'); hold on;
plot(-18.85+pos_W(1:5,2), -13.49+pos_W(1:5,1),'or'); hold on;

% figure;
% ecdf(dE0); hold on;
% ecdf(dU0);