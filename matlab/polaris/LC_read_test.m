folderPath='C:\Users\sande\Documents\Imperial\matlab\Files\Lat_cut_study\LW41\Load cell data\ACL recon\'
fileName=[folderPath,'Internal.log']
t=importdata(fileName);
loadcelldata=t.data(:,1);
% loadcelldata_test=loadcelldata(:,1)
offset=-1;
coeff=327;
force=coeff*loadcelldata+offset;