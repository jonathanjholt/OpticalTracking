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

folderPath='C:\Users\sande\Documents\Imperial\Documents\Exfix data\Cadaver testing ex fix all results-shared folder\femurSB01_left_20201120\Motion tracking';

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
proximal=mean(XYZdataP(6).points(:,48:50))/factor;
Error.proximal=mean(XYZdataP(6).points(:,51))/factor;

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

circLocs = [less_troch_med;less_troch_antlat;less_troch_postlat];
%for proximal end
[less_troch_centre, circleNormal, radius] = CircFit3D(circLocs);% For tibia
[gTt0, originT, APAxisT, PDAxisT, MLAxisT] = defineBodyFixedFrameFemur(less_troch_centre,fem_med,fem_lat,proximal, Right);
grt0 = [originT, 1]';%r is a point, T is a frame of reference.  This point is the location of the origin of the body fixed coordinate system

% For distal end
[gTf0, originF, APAxisF, PDAxisF, MLAxisF] = defineBodyFixedFrameFemur(less_troch_centre,fem_med,fem_lat,distal, Right);
grf0 = [originF, 1]';

H_0 = originT - originF;
% xvalues=[testi(1),testj(1),testk(1),];
% yvalues=[testi(2),testj(2),testk(2)];
% zvalues=[testi(3),testj(3),testk(3)];
%  figure
%  scatter3([less_troch_centre(1),fem_med(1),fem_lat(1)],[less_troch_centre(2),fem_med(2),fem_lat(2)],[less_troch_centre(3),fem_med(3),fem_lat(3)])
%  hold on 
% %  plot3([originT(1),testi(1)],[originT(2),testi(2)],[originT(3),testi(3)]);
% %   hold on
% %   plot3([originT(1),testj(1)],[originT(2),testj(2)],[originT(3),testj(3)]);
% %   hold on
%   plot3([originT(1),testk(1)],[originT(2),testk(2)],[originT(3),testk(3)]);
%   hold on
% %   plot3([originF(1),testi2(1)],[originF(2),testi2(2)],[originF(3),testi2(3)]);
% %   hold on
% %   plot3([originF(1),testj2(1)],[originF(2),testj2(2)],[originF(3),testj2(3)]);
% %   hold on
%   plot3([originF(1),testk2(1)],[originF(2),testk2(2)],[originF(3),testk2(3)]);
% %   scatter3(testi(1),testi(2),testi(3))
%   xlabel("X Label");
%   ylabel("Y Label");
%   zlabel("Z Label");
%   axis([-200 200 -200 200 -2000 -1600]);
%   pause(0.05);


%% Define fixed frame for each of the trackers in global coordinates

% Tibia tracker in global coordinates
[gTtt0]=defineTrackerFixedFrame_Exfix(ProxRxRyRz,ProxXYZ);

% Femur tracker in global coordinates
[gTft0]=defineTrackerFixedFrame_Exfix(DistRxRyRz,DistXYZ);


%% Generate constant transforms between the bony and tracker fixed coordinate systems

% This step generates a transform between each of the bone to their
% relevant trackers (rigid body assumed)

% For Tibia:

% For Tibia tracker to tibia frame
ttTtc=gTtt0\gTt0;

% For Tibia tracker to tibia origin
ttrtc=gTtt0\grt0;

% For Femur:

% For Femur tracker to femur frame
ftTfc=gTft0\gTf0;

% For Femur tracker to femur origin
ftrfc=gTft0\grf0;



%% Load tracked points

folderPath='C:\Users\sande\Documents\Imperial\Documents\Exfix data\Cadaver testing ex fix all results-shared folder\tibia&fibulaSB04_20201119\Motion Tracking';
fileName=[folderPath,'\imperial_ex_fix_TEST3.csv'];

cal = xlsread(fileName);
TrackedData.rawpoints = cal(:,all(~isnan(cal)));

TrackedData.points=TrackedData.rawpoints(~logical(TrackedData.rawpoints(:,11)== 0 | TrackedData.rawpoints(:,31)==0),:); % deletes any rows where trackers go missing (error appears to be 0)

tibiaXYZs=TrackedData.points(:,28:30)/factor;% tibia T probe rows 28:30
tibiaQs=TrackedData.points(:,24:27);
    for i = 1:size(tibiaQs,1)
        tibiaRxRyRzs(i,:)=quaternion2euler(tibiaQs(i,:));
    end

femurXYZs=TrackedData.points(:,8:10)/factor;% femur Y probe rows 8:10
femurQs=TrackedData.points(:,4:7);
    for i = 1:size(femurQs,1)
        femurRxRyRzs(i,:)=quaternion2euler(femurQs(i,:));
    end

Error.TibTracTrial=mean(TrackedData.points(:,31))/factor; 
Error.FemTracTrial=mean(TrackedData.points(:,11))/factor;
    
% Create matrices of tracker marker position and rotations in time

% Femur tracker:
[gTfti]=findTrackerFixedFrames_Exfix(femurRxRyRzs,femurXYZs);%how the femur tracker moves with time
% Tibia tracker:
[gTtti]=findTrackerFixedFrames_Exfix(tibiaRxRyRzs,tibiaXYZs);%how the tibia tracker moves with time



%% Calculate transformation matricies from body fixed to global fixed frames and motion relative to initial position


 for i = 1:size(gTfti,1)
     
     % Tibia Body fixed relative to global
     gTti{i,1} = gTtti{i,1}*ttTtc;
     % Position vectors of tibia origin in global
     grti{i,1} = gTtti{i,1}*ttrtc;
                                                 
     % Femur Body fixed relative to global
     gTfi{i,1} = gTfti{i,1}*ftTfc;
     % Position vectors of femur origin in global
     grfi{i,1} = gTfti{i,1}*ftrfc;
     
     % Calc motion of tibia bone relative to the femur
     fTt{i,1} = gTfi{i,1}\gTti{i,1};

     % Tibial origin point in femoral frame of reference
     frti{i,1} = gTfi{i,1}\(grti{i,1} - grfi{i,1});
     
        %5/6/19 calculates the position vectors of tibia and femur origin
        %relative to itself (in original fixed position)
        f0rfi{i,1}=gTf0\grfi{i,1};   
        t0rti{i,1}=gTt0\grti{i,1};
%         %5/6/19 calculates the tracker motion relative to itself (in original fixed
%         %position)
%         ft0Tfti{i}=gTft0\gTfti{i};
%         tt0Ttti{i}=gTtt0\gTtti{i};
        %5/6/19 calculates the body fixed motion relative to itself (in original fixed
        %position)
        f0Tfi{i,1}=gTf0\gTfi{i,1};
        t0Tti{i,1}=gTt0\gTti{i,1}; 
        %5/6/19 calculates the motion of the tibia bone relative to the femur,
        %using this transformation
        f0Tt{i,1} = f0Tfi{i,1}\t0Tti{i,1};
       

     %Begin richard new version (v5)
     ollie_I_(:,i)=gTfi{i,1}(1:3,1);%Femoral X axis unit vector, Grood and Suntay definition
     ollie_J_(:,i)=gTfi{i,1}(1:3,2);%Femoral Y axis unit vector, Grood and Suntay definition
     ollie_K_(:,i)=gTfi{i,1}(1:3,3);%Femoral Z axis unit vector, Grood and Suntay definition
                    
     ollie_i_(:,i)=gTti{i,1}(1:3,1);%Tibial x axis unit vector, Grood and Suntay definition
     ollie_j_(:,i)=gTti{i,1}(1:3,2);%Tibial y axis unit vector, Grood and Suntay definition
     ollie_k_(:,i)=gTti{i,1}(1:3,3);%Tibial z axis unit vector, Grood and Suntay definition
     
     ollie_e1_(:,i)=ollie_I_(:,i);%Femoral X axis in global reference frame, Grood and Suntay definition
     ollie_e3_(:,i)=ollie_k_(:,i);%Tibial z axis in global reference frame, Grood and Suntay definition
     ollie_e2_(:,i)=ucross(ollie_e3_(:,i),ollie_e1_(:,i));%Floating axis in global reference frame, Grood and Suntay definition
     
     %extract rotation and translaton matricies
     [ollie_anglestf(i,:),ollie_XYZtf(i,:)] = rotationsAndTranslations_v2(fTt{i}, Right);
     
          ollie_beta(i,1)=acosd(dot(ollie_I_(:,i),ollie_k_(:,i)));
        if Right
            ollie_IntExt(i,1)=asind(dot(-ollie_e2_(:,i),ollie_i_(:,i)));
            ollie_VarVal(i,1)=90-ollie_beta(i,1);
        else
            ollie_IntExt(i,1)=asind(dot(ollie_e2_(:,i),ollie_i_(:,i)));
            ollie_VarVal(i,1)=-(90-ollie_beta(i,1));
        end
        
        ollie_H_(i,:)=[grti{i,1}(1:3)-grfi{i,1}(1:3)]';%translation vector of from femoral origin to tibial origin (in global coordinate frame)
         if Right
            ollie_MedLat(i,1)=dot(ollie_H_(i,:),ollie_e1_(:,i));%projected onto the medial lateral axis e1
        else
            ollie_MedLat(i,1)=dot(ollie_H_(i,:),-ollie_e1_(:,i));
        end
     ollie_AntPost(i,1)=dot(ollie_H_(i,:),ollie_e2_(:,i));%projected onto the anterior posterior axis e2
     ollie_ProxDist(i,1)=-dot(ollie_H_(i,:),ollie_e3_(:,i));%projected onto the compression distraction axis e3, minus sign to make distraction +ve
        
     %5/6/19
      
        I_(:,i)=f0Tfi{i,1}(1:3,1);%Femoral X axis unit vector, Grood and Suntay definition
        J_(:,i)=f0Tfi{i,1}(1:3,2);%Femoral Y axis unit vector, Grood and Suntay definition
        K_(:,i)=f0Tfi{i,1}(1:3,3);%Femoral Z axis unit vector, Grood and Suntay definition
                    
        i_(:,i)=t0Tti{i,1}(1:3,1);%Tibial x axis unit vector, Grood and Suntay definition
        j_(:,i)=t0Tti{i,1}(1:3,2);%Tibial y axis unit vector, Grood and Suntay definition
        k_(:,i)=t0Tti{i,1}(1:3,3);%Tibial z axis unit vector, Grood and Suntay definition
                    
     e1_(:,i)=I_(:,i);%Femoral X axis in global reference frame, Grood and Suntay definition
     e3_(:,i)=k_(:,i);%Tibial z axis in global reference frame, Grood and Suntay definition
     e2_(:,i)=ucross(e3_(:,i),e1_(:,i));%Floating axis in global reference frame, Grood and Suntay definition
        
     %extract rotation and translaton matricies
        %5/6/19
       [anglestf(i,:),XYZtf(i,:)] = rotationsAndTranslations_v2(f0Tt{i}, Right);
     
     beta(i,1)=acosd(dot(I_(:,i),k_(:,i)));
        if Right
            IntExt(i,1)=asind(dot(-e2_(:,i),i_(:,i)));
            VarVal(i,1)=90-beta(i,1);
        else
            IntExt(i,1)=asind(dot(e2_(:,i),i_(:,i)));
            VarVal(i,1)=-(90-beta(i,1));
        end


            
%      %5/6/19
        H_(i,:)=[t0rti{i,1}(1:3)-f0rfi{i,1}(1:3)]';
        
%        FlexExt(i,1)=asind(-dot(e2_(:,i),K_(:,i)));
%        FlexExt(i,1)=acosd(dot(J_(:,1),e2_(:,i)));

        if Right
            MedLat(i,1)=dot(H_(i,:),e1_(:,i));%projected onto the medial lateral axis e1
        else
            MedLat(i,1)=dot(H_(i,:),-e1_(:,i));
        end
     AntPost(i,1)=dot(H_(i,:),e2_(:,i));%projected onto the anterior posterior axis e2

 end
 
%  figure
% for i=1:numel(grfi)
%   plot3(xval(i,1:2), yval(i,1:2), zval(i,1:2), ".");
%   xlabel("X Label");
%   ylabel("Y Label");
%   zlabel("Z Label");
%   axis([-200 0 -200 0 -1900 -1600]);
%   pause(0.05);
%   clf;
% end

 MinFlex= min(anglestf(:,1));
 MaxFlex = max(anglestf(:,1));
 RangeFlex = MaxFlex-MinFlex;
 FlexExt = anglestf(:,1)- anglestf(1,1);
 
 ollie_MinFlex = min(ollie_anglestf(:,1));
 ollie_MaxFlex = max(ollie_anglestf(:,1));
 ollie_RangeFlex = ollie_MaxFlex-ollie_MinFlex;
 ollie_FlexExt = ollie_anglestf(:,1)- ollie_anglestf(1,1);

%%Find peaks and throughs

time=zeros(numel(ollie_ProxDist),1);
for i=1:numel(ollie_ProxDist)
    time(i)=0.05*i;
end

[ProxDist_max,ProxDistmax_locs]=findpeaks(ollie_ProxDist,time,'MinPeakDistance',1.2);

[ProxDist_min,ProxDistmin_locs]=findpeaks(-ollie_ProxDist,time,'MinPeakDistance',1.2);


%remove beginning and ending with no loading (removing noise)
noise = find(abs(ProxDist_max-ollie_ProxDist(1))<0.25)
ProxDist_max(noise) = [];ProxDistmax_locs(noise)=[];
noise = find(abs(-ProxDist_min-ollie_ProxDist(1))<0.3);
ProxDist_min(noise) = [];ProxDistmin_locs(noise)=[];
ProxDist_min=-ProxDist_min;
disp(size(ProxDist_min));

%remove local peaks and troughs
wrong_peaks = find(ProxDist_max<ollie_ProxDist(1));
ProxDist_max(wrong_peaks) = [];ProxDistmax_locs(wrong_peaks)=[];
wrong_troughs = find(ProxDist_min>ollie_ProxDist(1));
ProxDist_min(wrong_troughs) = [];ProxDistmin_locs(wrong_troughs)=[];


%remove first 10 cycles
ProxDist_max=ProxDist_max(11:end-1);ProxDistmax_locs=ProxDistmax_locs(11:end-1);
ProxDist_min=ProxDist_min(11:end-1);ProxDistmin_locs=ProxDistmin_locs(11:end-1);

%Individual amplitudes for each cycle


Amplitude_complete=ProxDist_max-ProxDist_min;

% plot peaks and troughs

figure('Name','Movement Peaks and Troughs')


plot(time,ollie_ProxDist,ProxDistmax_locs,ProxDist_max,'o',ProxDistmin_locs,ProxDist_min,'x')

title ('Fracture distance')

xlabel ('Time (s)')

ylabel('Distance(mm)')



% mean of last 89/90 cycles

MeanpeakProxDist = mean(ProxDist_max);

MeanthroughProxDist = mean(ProxDist_min);
 
Amplitude=MeanpeakProxDist-MeanthroughProxDist;

disp("The amplitude of the fracture movement = "+ num2str(Amplitude))

%% Load cell sync code to be added here

%Using modified keydata_v2 code find first peak for Keydata+LC data and
%make new synced data by cutting off values before start point
                 
%% Plot key results
% close all

%
% time = 0:(1/hz):((length(anglestf)-1)/hz);
% figure
% subplot(2,3,1)
% plot(time,ollie_anglestf(:,1))
% ylabel('Flexion-Extension')
% xlabel('Time / s')
% title('Relative Rotations')
% subplot(2,3,2)
% plot(time,ollie_anglestf(:,3))
% ylabel('Internal-External')
% xlabel('Time / s')
% subplot(2,3,3)
% plot(time,ollie_anglestf(:,2))
% ylabel('Varus-Valgus')
% xlabel('Time / s')
% subplot(2,3,4)
% plot(time,ollie_AntPost(:,1))
% ylabel('Anterior-Posterior')
% xlabel('Time / s')
% title('Relative Translations')
% subplot(2,3,5)
% plot(time,ollie_ProxDist(:,1))
% ylabel('Proximal-Distal')
% xlabel('Time / s')
% subplot(2,3,6)
% plot(time,ollie_MedLat(:,1))
% ylabel('Medial-Lateral')
% xlabel('Time / s')                    
                   
% figure
% subplot(3,1,1)
% plot(FlexExt,VarVal)
% ylabel('Varus-Valgus')
% xlabel('Flexion')
% subplot(3,1,2)
% plot(FlexExt,IntExt)
% ylabel('Internal-External')
% xlabel('Flexion')
% subplot(3,1,3)
% plot(FlexExt,AntPost)
% ylabel('Anterior-Posterior')
% xlabel('Flexion')

%% Getting Data from Specific Flexion Angles

% [nextpeakstart,rowpos,FlexExtValues,VarValValues,IntExtValues,AntPosValues] = KeyData_v2(VarVal,IntExt,-FlexExt,AntPost);
% [ollie_nextpeakstart,ollie_rowpos,ollie_FlexExtValues,ollie_VarValValues,ollie_IntExtValues,ollie_AntPosValues] = KeyData_v2(ollie_VarVal,ollie_IntExt,-ollie_FlexExt,ollie_AntPost);

% [ollie_nextpeakstart,ollie_rowpos,ollie_FlexExtValues,ollie_VarValValues,ollie_IntExtValues,ollie_AntPosValues,ollie_ProxDistValues,ollie_MedLatValues] = KeyData_v3(ollie_VarVal,ollie_IntExt,-ollie_FlexExt,ollie_AntPost,ollie_ProxDist,ollie_MedLat);


%% 30 deg specific position
% 
% Static30deg.AP=mean(AntPost);
% Static30deg.IE=mean(IntExt);
% Static30deg.VV=mean(VarVal);
% Static30deg.FE=mean(FlexExt);
% Static30deg.ML=mean(MedLat);
% Static30deg.PD=mean(ProxDist);

%% Injury special

% figure
% subplot(2,3,1)
% plot(time,FlexExt)
% ylabel('Flexion-Extension')
% xlabel('Time / s')
% title('Relative Rotations')
% subplot(2,3,2)
% plot(time,IntExt)
% ylabel('Internal-External')
% xlabel('Time / s')
% subplot(2,3,3)
% plot(time,VarVal)
% ylabel('Varus-Valgus')
% xlabel('Time / s')
% subplot(2,3,4)
% plot(time,AntPost)
% ylabel('Anterior-Posterior')
% xlabel('Time / s')
% title('Relative Translations')
% subplot(2,3,5)
% plot(time,ProxDist)
% ylabel('Proximal-Distal')
% xlabel('Time / s')
% subplot(2,3,6)
% plot(time,MedLat)
% ylabel('Medial-Lateral')
% xlabel('Time / s')
% 
% Injury.AP.start=AntPost(1,:);
% Injury.AP.end=AntPost(length(anglestf),:);
% Injury.AP.max=max(AntPost);
% Injury.AP.min=min(AntPost);
% Injury.PD.start=ProxDist(1,:);
% Injury.PD.end=ProxDist(length(anglestf),:);
% Injury.PD.max=max(ProxDist);
% Injury.PD.min=min(ProxDist);
% Injury.VV.start=VarVal(1,:);
% Injury.VV.end=VarVal(length(anglestf),:);
% Injury.VV.max=max(VarVal);
% Injury.VV.min=min(VarVal);
% Injury.ML.start=MedLat(1,:);
% Injury.ML.end=MedLat(length(anglestf),:);
% Injury.ML.max=max(MedLat);
% Injury.ML.min=min(MedLat);
% Injury.IE.start=IntExt(1,:);
% Injury.IE.end=IntExt(length(anglestf),:);
% Injury.IE.max=max(IntExt);
% Injury.IE.min=min(IntExt);
% Injury.FE.start=anglestf(1,1);
% Injury.FE.end=anglestf(length(anglestf),1);
% Injury.FE.max=max(anglestf(:,1));
% Injury.FE.min=min(anglestf(:,1));
% 
% ollie_Injury.AP.start=ollie_AntPost(1,:);
% ollie_Injury.AP.end=ollie_AntPost(length(ollie_anglestf),:);
% ollie_Injury.AP.max=max(ollie_AntPost);
% ollie_Injury.AP.min=min(ollie_AntPost);
% ollie_Injury.PD.start=ollie_ProxDist(1,:);
% ollie_Injury.PD.end=ollie_ProxDist(length(ollie_anglestf),:);
% ollie_Injury.PD.max=max(ollie_ProxDist);
% ollie_Injury.PD.min=min(ollie_ProxDist);
% ollie_Injury.VV.start=ollie_VarVal(1,:);
% ollie_Injury.VV.end=ollie_VarVal(length(ollie_anglestf),:);
% ollie_Injury.VV.max=max(ollie_VarVal);
% ollie_Injury.VV.min=min(ollie_VarVal);
% ollie_Injury.ML.start=ollie_MedLat(1,:);
% ollie_Injury.ML.end=ollie_MedLat(length(ollie_anglestf),:);
% ollie_Injury.ML.max=max(ollie_MedLat);
% ollie_Injury.ML.min=min(ollie_MedLat);
% ollie_Injury.IE.start=ollie_IntExt(1,:);
% ollie_Injury.IE.end=ollie_IntExt(length(ollie_anglestf),:);
% ollie_Injury.IE.max=max(ollie_IntExt);
% ollie_Injury.IE.min=min(ollie_IntExt);
% ollie_Injury.FE.start=ollie_anglestf(1,1);
% ollie_Injury.FE.end=ollie_anglestf(length(ollie_anglestf),1);
% ollie_Injury.FE.max=max(ollie_anglestf(:,1));
% ollie_Injury.FE.min=min(ollie_anglestf(:,1));
% 
% figure
% subplot(2,3,1)
% plot(time,ollie_FlexExt)
% ylabel('Flexion-Extension')
% xlabel('Time / s')
% title('Relative Rotations')
% subplot(2,3,2)
% plot(time,ollie_IntExt)
% ylabel('Internal-External')
% xlabel('Time / s')
% subplot(2,3,3)
% plot(time,ollie_VarVal)
% ylabel('Varus-Valgus')
% xlabel('Time / s')
% subplot(2,3,4)
% plot(time,ollie_AntPost)
% ylabel('Anterior-Posterior')
% xlabel('Time / s')
% title('Relative Translations')
% subplot(2,3,5)
% plot(time,ollie_ProxDist)
% ylabel('Proximal-Distal')
% xlabel('Time / s')
% subplot(2,3,6)
% plot(time,ollie_MedLat)
% ylabel('Medial-Lateral')
% xlabel('Time / s')



