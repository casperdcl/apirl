%% EXAMPLE MLEM MARTIN PROJECTOR (ANY SPAN)
clear all
close all
set_framework_environment();
% set_framework_environment(basePath, binaryPath);
%% DATA PATHS
sinogramFilename = '/media/mab15/DATA_BACKUP/Scans/PatientData/FDG_Patient_01/e7/PETSinoPlusUmap-Converted/PETSinoPlusUmap-00/PETSinoPlusUmap-00-sino-uncomp.s.hdr';
attenuationMap_filename = '/media/mab15/DATA_BACKUP/Scans/PatientData/FDG_Patient_01/e7/PETSinoPlusUmap-Converted/PETSinoPlusUmap-00/PETSinoPlusUmap-00-umap.v.hdr';
attenuationMapHardware_filename = '/media/mab15/DATA_BACKUP/Scans/PatientData/FDG_Patient_01/e7/PETSinoPlusUmap-Converted/PETSinoPlusUmap-00/PETSinoPlusUmap-00-umap-hardware.v.hdr';
scatterBinaryFilename = '/media/mab15/DATA_BACKUP/Scans/PatientData/FDG_Patient_01/e7/PETSinoPlusUmap-Converted/PETSinoPlusUmap-00/rawdata_sino/scatter_estim2d_000000.s';
t1DicomPath = '/media/mab15/DATA_BACKUP/Scans/PatientData/FDG_Patient_01/T1/';
%% INIT CLASS GPET
PET.scanner = 'mMR';
PET.method =  'otf_siddon_gpu';
PET.PSF.type = 'none';
PET.radialBinTrim = 0;
PET.Geom = '';
PET.method_for_randoms = 'from_ML_singles_matlab'; %'from_e7_binary_interfile';
PET.method_for_scatter = 'from_e7_binary_interfile';
% To change span:
PET.sinogram_size.span = 11; % Any span, 0 for multislice 2d, -1 for 2d.
PET.nSubsets = 1;
PET = classGpet(PET);

%% READ AND PREPARE DATA
% EMISSION SINOGRAM
[sinogram, delayedSinogram, structSizeSino3d] = interfileReadSino(sinogramFilename);
if PET.sinogram_size.span > 1
    sino_span = PET.apply_axial_compression_from_span1(sinogram);
elseif PET.sinogram_size.span == 0
    % Get the direct sinograms:
    sino_span = sinogram(:,:,1:PET.sinogram_size.nRings);
    delayedSinogram = delayedSinogram(:,:,1:PET.sinogram_size.nRings);
elseif PET.sinogram_size.span == -1
    % Get the central slice:
    ring = 32;
    sino_span = sinogram(:,:,ring);
    delayedSinogram = delayedSinogram(:,:,ring);
end
% ATTENUATION MAP
[attenuationMapHuman refMuMap] = interfileReadSiemensImage(attenuationMap_filename);
[attenuationMapHardware refMuMap] = interfileReadSiemensImage(attenuationMapHardware_filename);
attenuationMap = attenuationMapHuman + attenuationMapHardware;
if PET.sinogram_size.span == -1
    % Get only
    slice = 2*ring; % atten map has 127 slices and the direct rings 64.
    attenuationMap = attenuationMap(:,:,slice);
    attenuationMapHuman = attenuationMapHuman(:,:,slice);
end
% SCATTER
% The function it uses
% MR IMAGE
% for the mr image is important to have the correct bed position, that is
% set from any interfile header of the scan, for example the attenuation
% map:
PET.setBedPosition(attenuationMap_filename);
[MrInPet, refMr] = PET.getMrInPetImageSpace(t1DicomPath);
%% NORM
ncfs = PET.NCF(); % time-invariant.
%ncfs = PET.NCF('PETSinoPlusUmap-norm.n'); % time-variant.
%% ACFs
acfs = PET.ACF(attenuationMap, refMuMap);
acfs_human = PET.ACF(attenuationMapHuman, refMuMap);
%% RANDOMS
% if uses the randomsBinaryFilename = '/media/mab15/DATA/PatientData/FDG/PETSinoPlusUmap-Converted/PETSinoPlusUmap-00/temp/smoothed_rand_00.s';
% randoms = PET.R(randomsBinaryFilename);
randoms = PET.R(delayedSinogram);
%% SCATTER
% for 3d or multislice 2d:
if PET.sinogram_size.span >=0
    scatter = PET.S(scatterBinaryFilename, sino_span, ncfs, acfs_human, randoms); % Needs all that parameters to scale it.
else
    % For 2d select the sinogram ring:
    scatter = PET.S(scatterBinaryFilename, sino_span, ncfs, acfs_human, randoms, ring);
end
% This scatter is already normalzied
%% SENSITIVITY IMAGE
anf = acfs .* ncfs;
anf(anf~=0) = 1./anf(anf~=0);
sensImage = PET.Sensitivity(anf);
%% OP-OSEM
% additive term:
additive = (randoms + scatter).*ncfs.*acfs; % (randoms +scatter)./(afs*nfs) = (randoms+scatter)+
%additive = zeros(size(additive));
initial_image = PET.ones();
recon = PET.OPOSEM(sino_span, additive, sensImage, initial_image, ceil(60/PET.nSubsets));
% %%
% figure(1);
% recon = PET.ones();
% for i = 1 : 60
%     recon = recon.*PET.PT(sino_span./(PET.P(recon)+ additive + 1e-10))./(sensImage+1e-10);
%     imagesc(recon), axis equal,  drawnow;
% end