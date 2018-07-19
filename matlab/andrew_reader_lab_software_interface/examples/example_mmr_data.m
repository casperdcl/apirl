%% EXAMPLE MLEM MARTIN PROJECTOR (ANY SPAN)
clear
close all
apirlPath = 'F:\workspace\apirl-code\trunk\';
addpath([apirlPath 'matlab' filesep 'andrew_reader_lab_software_interface']);
set_framework_environment(apirlPath);
%% WITH INTERFACE
%PETData = PETDataClass('/media/mab15/DATA/PatientData/FDG/Raw_PET/'); % prompts a window to locate rawdata file
%% TEST LIST-MODE
PETData = PETDataClass('F:\Downloads\F18CylinderPostUpgradeRawData\SmallExport');
ncf = PETData.NCF();
PETData.uncompress(PETData.Data.emission.n);
PETData.Reconstruct(1, 100, 0, 1);
PETData.Reconstruct(1, 100, 1, 1);

