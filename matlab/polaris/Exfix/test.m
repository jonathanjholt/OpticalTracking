clear all
close all

%% State whether specimen is left or right knee

Right = true; % True for right knee, false for left knee
% hz = 20;
hz = 60;

angle = @(u_,v_) acosd(dot(u_,v_)/(norm(u_,2)*norm(v_,2))); %define a function to calculate the angle between two vectors
ucross=@(u_,v_) cross(u_,v_)/norm(cross(u_,v_),2); %define function to find unit cross product
uvector=@(a,b) (b-a)/norm(b-a,2); %define a function to find a unit vector from a to b

factor=1;

%% Import marker and probe tip locations for inital digitising position

% for the polaris coordinate systems, Y is width of view (parallel to floor), X is height, Z is depth of view 

% Points are later used to define body specific coordinate systems,
% calculated from anatomical landmarks

% Probe points are from stationary condition, only the probe
% columns should differ between data files.

% Find folder path for knee tracking data files:

folderPath='C:\Users\sande\Documents\Imperial\matlab\Files\Exfix_test';

% TIBIA
%data saved as 3 matrices within XYZdataT

% Lesser Trochancter circle 1
fileName=[folderPath,'\medial_less-troch.csv'];
cal = xlsread(fileName);
XYZdataP(1).points = cal(:,all(~isnan(cal)));

% Lesser Trochancter circle 2
fileName=[folderPath,'\posterolat_less-troch.csv'];
cal = xlsread(fileName);
XYZdataP(2).points = cal(:,all(~isnan(cal)));

% Lesser Trochancter circle 3
fileName=[folderPath,'\anterolat_less-troch.csv'];
cal = xlsread(fileName);
XYZdataP(3).points = cal(:,all(~isnan(cal)));

% Medial femoral epicondyle
fileName=[folderPath,'\med_fem_epicondyle.csv'];
cal = xlsread(fileName);
XYZdataP(4).points = cal(:,all(~isnan(cal)));

% Lateral femoral epicondyle
fileName=[folderPath,'\lat_fem_epicondyle.csv'];
cal = xlsread(fileName);
XYZdataP(5).points = cal(:,all(~isnan(cal)));

%Proximal point close to fracture (to be used as origin)
fileName=[folderPath,'\prox_fragment.csv'];
cal = xlsread(fileName);
XYZdataP(6).points = cal(:,all(~isnan(cal)));

%Distal point close to fracture (to be used as origin)
fileName=[folderPath,'\dist_fragment.csv'];
cal = xlsread(fileName);
XYZdataP(7).points = cal(:,all(~isnan(cal)));

% Find position of digitised tibial landmarks by 'finding' the probe
% columns in the respect digitised data files and extracting the points of
% the tip of the probe


% Medial Lesser Trochanter
less_troch_med=mean(XYZdataP(1).points(:,48:50))/factor;
Error.Less_troch_med=mean(XYZdataP(1).points(:,51))/factor;

% Anterolateral Lesser Trochanter
less_troch_antlat=mean(XYZdataP(2).points(:,48:50))/factor;
Error.Less_troch_antlat=mean(XYZdataP(2).points(:,51))/factor;

% Posterolater Lesser Trochanter
less_troch_postlat=mean(XYZdataP(3).points(:,48:50))/factor;
Error.Less_troch_postlat=mean(XYZdataP(3).points(:,51))/factor;

% Medial femoral epicondyle
fem_med=mean(XYZdataP(4).points(:,48:50))/factor;
Error.fem_med=mean(XYZdataP(4).points(:,51))/factor;

%Lateral femoral epicondyle
fem_lat=mean(XYZdataP(5).points(:,48:50))/factor;
Error.fem_lat=mean(XYZdataP(5).points(:,51))/factor;

%Proximal fracture
proximal=mean(XYZdataP(6).points(:,8:10))/factor;
Error.proximal=mean(XYZdataP(6).points(:,11))/factor;

%Distal fracture
distal=mean(XYZdataP(7).points(:,48:50))/factor;
Error.distal=mean(XYZdataP(7).points(:,51))/factor;



%%%%%%% END OF IMPORT OF ANATOMIC REFERENCE POINTS %%%%%%%%%%%%%%

%% Import position of tracking markers at inital digitisation position

% Same process as digitising the bony landmarks but instead of the probe
% points, the trackers themselves are used.

% NB: Pins can be named whatever is deemed fit. If something general, like
% 'Pin1' is used, then be sure to keep a record of which bone this tracker
% relates to for each experiment.

% As the knee should not have moved during the digitisation process, any of
% the XYZdata sets may be used.

% Find position and quaternions, using XYZdataT:
%CHANGE ONCE WE HAVE FULL DATA
% Tibia Tracker
ProxXYZ=mean(XYZdataP(1).points(:,28:30))/factor; % Prox T probe rows 28:30
ProxQ=mean(XYZdataP(1).points(:,24:27));
ProxRxRyRz=quaternion2euler(ProxQ);
Error.ProxTrac=mean(XYZdataP(1).points(:,31))/factor;

% Femur Tracker
DistXYZ=mean(XYZdataP(1).points(:,8:10))/factor; % Dist Y probe rows 8:10
DistQ=mean(XYZdataP(1).points(:,4:7));
DistRxRyRz=quaternion2euler(DistQ);
Error.FemTrac=mean(XYZdataP(1).points(:,11))/factor;

%% Define coordinate systems for each bone using digitised points

[less_troch_centre,rad,v1n,v2nb]=circlefit3d(less_troch_med,less_troch_antlat,less_troch_postlat);
% For tibia
[gTt0, originT, APAxisT, PDAxisT, MLAxisT] = defineBodyFixedFrameSpencer(less_troch_centre,fem_med,fem_lat,proximal, Right);
grt0 = [originT, 1]';%r is a point, T is a frame of reference.  This point is the location of the origin of the body fixed coordinate system

% For femur
[gTf0, originF, APAxisF, PDAxisF, MLAxisF] = defineBodyFixedFrameSpencer(less_troch_centre,fem_med,fem_lat,distal, Right);
grf0 = [originF, 1]';

H_0 = originT - originF;
coord=[less_troch_med;less_troch_antlat;less_troch_postlat;less_troch_centre;fem_med;fem_lat;distal;proximal]
scatter3(coord(:,1),coord(:,2),coord(:,3))
xlabel('x')
ylabel('y')
zlabel('z')
% xlim([0 300])
% ylim([0 300])
% zlim([1500 2000])