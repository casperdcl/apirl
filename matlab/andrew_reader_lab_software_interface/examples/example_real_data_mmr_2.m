%% EXAMPLE MLEM MARTIN PROJECTOR (ANY SPAN)
clear all
close all
%% CONFIGURE PATHS
% Check what OS I am running on:
if(strcmp(computer(), 'GLNXA64'))
    os = 'linux';
    pathBar = '/';
    sepEnvironment = ':';
elseif(strcmp(computer(), 'PCWIN') || strcmp(computer(), 'PCWIN64'))
    os = 'windows';
    pathBar = '\';
    sepEnvironment = ';';
else
    disp('OS not compatible');
    return;
end

% APIRL PATH
apirlPath = 'E:\apirl-code\trunk\';
addpath(genpath([apirlPath pathBar 'matlab']));
setenv('PATH', [getenv('PATH') sepEnvironment apirlPath pathBar 'build' pathBar 'bin']);
setenv('LD_LIBRARY_PATH', [getenv('LD_LIBRARY_PATH') sepEnvironment apirlPath pathBar 'build' pathBar 'bin']);
%% DATA PATHS
root = '\\Bioeng139-pc\pet-m\FDG_Patient_02';
sinogramFilename = [root,'\e7\data-Converted\data-00\data-00-sino-uncom_00.s.hdr'];
attenuationMap_filename = [root,'\e7\data-Converted\data-00\data-00-umap.v.hdr'];
attenuationMapHardware_filename = [root,'\e7\data-Converted\data-00\data-00-umap-hardware.v.hdr' ];
normalizationFilename = [root,'\e7\data-Converted\data-norm.n.hdr'];
scatterBinaryFilename = [root,'\e7\data-Converted\data-00\rawdata_sino\scatter_estim2d_000000.s'];
randomsBinaryFilename = [root,'\e7\data-Converted\data-00\rawdata_sino\smoothed_rand_00.s'];
% acfBinaryFilename = [root,'\e7\data-Converted\data-00\rawdata_sino\acf_00.a'];
% acfHumanBinaryFilename = [root,'\e7\data-Converted\data-00\rawdata_sino\acf_second_00.a'];
% ncfBinaryFilname = [root,'\e7\data-Converted\data-00\rawdata_sino\norm3d_00.a'];
t1DicomPath = [root,'\MPRAGE_image\'];


%% INIT CLASS GPET
param.scanner = 'mMR';
param.method =  'otf_siddon_cpu';
param.PSF.type = 'none';
param.method_for_randoms = 'from_ML_singles_matlab';% ; 'from_e7_binary_interfile'
param.method_for_scatter = 'from_e7_binary_interfile';
% To change span:
param.sinogram_size.span = 11; % Any span, 0 for multislice 2d, -1 for 2d.
param.nSubsets = 1;
PET = classGpet(param);
%%%% READ AND PREPARE DATA
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
ncfs = PET.NCF(normalizationFilename); % time-invariant.
%ncfs = PET.NCF('PETSinoPlusUmap-norm.n'); % time-variant.
%% ACFs
acfs = PET.ACF(attenuationMap, refMuMap);
acfs_human = PET.ACF(attenuationMapHuman, refMuMap);
%% RANDOMS
randoms = PET.R(randomsBinaryFilename);
% randoms = PET.R(delayedSinogram);
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

reconFiltered = PET.Gauss3DFilter(recon,4);
mlem_BQML = PET.BQML(reconFiltered,sinogramFilename,normalizationFilename);
mlem_SUV= PET.SUV(reconFiltered,sinogramFilename,normalizationFilename);

