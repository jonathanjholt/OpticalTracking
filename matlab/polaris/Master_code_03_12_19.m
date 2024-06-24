% This file uses a modified version of Vegatestcode25_11_19 (KeydataOutput) to
% automatically output different loading conditions of the knee. Only
% outputs Keydata, for full graphs and raw data use original
% Vegatestcode25_11_19 script

folderPath='C:\Users\sholthof\Desktop\Imperial\Petra MCL work\Pilot study\Pilot study\ACL tension kinematics testing\Pilot 28112023\ACL stuff\'; %folder where kinematics data is stored
Digitise_file='C:\Users\sholthof\Desktop\Imperial\Petra MCL work\Pilot study\Pilot study\ACL tension kinematics testing\Pilot 28112023\Digitise'; %folder where digitise data is stored
LCPath='C:\Users\sande\Documents\Imperial\matlab\Files\Lat_cut_study\LW42\Load cell data\Men recon\'; %folder for LC data
OutputFile='C:\Users\sholthof\Desktop\Imperial\Petra MCL work\Pilot study\Pilot study\MCL_template.xlsx'; %excel file to output data to
LC_on=false; %true means load cell is connected, false is testing without LC
coeff=327; %values from calibration, converts mV/V to N
offset=-1;
injury=1;%Stage of testing, used to define which sheet to write to (e.g. if Intact is the first sheet, use injury=1 for intact knee data)
Right=false; %true for right knee, false for left knee
Total_kneestates=8; %total amount of tested kneestates

% KeydataOutput(Digitise_file,[folderPath,'Neutral.csv'],'Neutral.csv',[LCPath,'Neutral.log'],coeff,offset,LC_on,OutputFile,injury,Right,Total_kneestates);
% clear KeyDataOutput
% close all
% KeydataOutput(Digitise_file,[folderPath,'Anterior.csv'],'Anterior.csv',[LCPath,'Anterior.log'],coeff,offset,LC_on,OutputFile,injury,Right,Total_kneestates);
% clear KeyDataOutput
% close all
% KeydataOutput(Digitise_file,[folderPath,'Anterior_autograph.csv'],'Posterior.csv',[LCPath,'Posterior.log'],coeff,offset,LC_on,OutputFile,injury,Right,Total_kneestates);
% clear KeyDataOutput
% close all
% KeydataOutput(Digitise_file,[folderPath,'Anterior_SL.csv'],'External.csv',[LCPath,'External.log'],coeff,offset,LC_on,OutputFile,injury,Right,Total_kneestates);
% clear KeyDataOutput
% close all
% KeydataOutput(Digitise_file,[folderPath,'Anterior_intact.csv'],'Internal.csv',[LCPath,'Internal.log'],coeff,offset,LC_on,OutputFile,injury,Right,Total_kneestates);
% clear KeyDataOutput
% close all
% KeydataOutput(Digitise_file,[folderPath,'Anterior_injury.csv'],'Varus.csv',[LCPath,'Varus.log'],coeff,offset,LC_on,OutputFile,injury,Right,Total_kneestates);
% clear KeyDataOutput
% close all
KeydataOutput(Digitise_file,[folderPath,'Anterior_rope.csv'],'Valgus.csv',[LCPath,'Valgus.log'],coeff,offset,LC_on,OutputFile,injury,Right,Total_kneestates);
clear KeyDataOutput
close all
KeydataOutput(Digitise_file,[folderPath,'Anterior_tendon.csv'],'Anterior+IR.csv',[LCPath,'External+ant.log'],coeff,offset,LC_on,OutputFile,injury,Right,Total_kneestates);
clear KeyDataOutput
close all
% KeydataOutput(Digitise_file,[folderPath,'Pivot shift.csv'],'Pivot shift.csv',OutputFile,injury,Right,Total_kneestates);
% clear KeyDataOutput
% close all
% % KeydataOutput(Digitise_file,[folderPath,'Filler.csv'],'Internal+ant.csv',[LCPath,'Internal+ant.log'],coeff,offset,LC_on,OutputFile,injury,Right,Total_kneestates);
% % clear KeyDataOutput
% % close all