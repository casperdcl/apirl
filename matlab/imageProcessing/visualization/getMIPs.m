%  *********************************************************************
%  Reconstruction Framework for Siemens Biograph mMR.  
%  Autor: Martín Belzunce. Kings College London.
%  Fecha de Creación: 11/03/2015
%  *********************************************************************
%  This function return the mips after rotating the image in the z axis.
function [mipTransverse, mipCoronal, mipSagital] = getMIPs(volume, angle)

volume = imrotate3(volume, angle, [0 0 1], 'crop');
% Get the Maximum Intensity ProjectionsL
% Transverse XY:
mipTransverse = max(volume,[],3);
% Coronal XZ:
mipCoronal = max(volume,[],2);
mipCoronal = permute(mipCoronal, [3 1 2]);
% Sagital YZ:
mipSagital = max(volume,[],1);
mipSagital = permute(mipSagital, [2 3 1]);





