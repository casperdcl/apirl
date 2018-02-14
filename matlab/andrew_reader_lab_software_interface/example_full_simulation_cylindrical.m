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
apirlPath = '/home/mab15/workspace/apirl-code/trunk/';
addpath(genpath([apirlPath pathBar 'matlab']));
setenv('PATH', [getenv('PATH') sepEnvironment apirlPath pathBar 'build' pathBar 'bin']);
setenv('LD_LIBRARY_PATH', [getenv('LD_LIBRARY_PATH') sepEnvironment apirlPath pathBar 'build' pathBar 'bin']);
%% INIT CLASS GPET
PET.scanner = 'cylindrical';
PET.method =  'otf_siddon_gpu';
PET.PSF.type = 'none';
PET.radialBinTrim = 0;
PET.Geom = '';
PET.random_algorithm = 'from_ML_singles_matlab';
PET = classGpet(PET);
%% SIMULATE A BRAIN PHANTOM WITH ATTENUATION, NORMALIZATION, RANDOMS AND SCATTER
%[sinogram, delayedSinogram, structSizeSino3d] = interfileReadSino('E:\PatientData\FDG\PETSinoPlusUmap-Converted\PETSinoPlusUmap-00\PETSinoPlusUmap-00-sino-uncomp.s.hdr');
load BrainMultiMaps_mMR.mat;

tAct = permute(MultiMaps_Ref.PET, [2 1 3]);
tAct = tAct(end:-1:1,:,:);
tMu = permute(MultiMaps_Ref.uMap, [2 1 3]);
tMu = tMu(end:-1:1,:,:);
pixelSize_mm = [2.08625 2.08625 2.03125];
xLimits = [-size(tAct,2)/2*pixelSize_mm(2) size(tAct,2)/2*pixelSize_mm(2)];
yLimits = [-size(tAct,1)/2*pixelSize_mm(1) size(tAct,1)/2*pixelSize_mm(1)];
zLimits = [-size(tAct,3)/2*pixelSize_mm(3) size(tAct,3)/2*pixelSize_mm(3)];
refAct = imref3d(size(tAct),xLimits,yLimits,zLimits);
refAt  = imref3d(size(tMu),xLimits,yLimits,zLimits);

% Change the image sie, to the one of the phantom:
PET.init_image_properties(refAct);
%% USE THIS FOR 2D
% % Change sinogram size:
% param.sinogram_size.nRadialBins = 300;
% param.sinogram_size.nAnglesBins = 300;
% param.nSubsets = 1;
% param.sinogram_size.span = -1;    % Span 0 is multi slice 2d.
% param.sinogram_size.nRings = 1; % for span 1, 1 ring.
% param.image_size.matrixSize = [refAct.ImageSize(1:2) 1];
% % iMAGE NEEDS TO BE 2D ALSO:
% tAct = tAct(:,:,60);
% tMu = tMu(:,:,60);

% %% USE THIS FOR MULTI-SLICE 2D
% % Change sinogram size:
% param.sinogram_size.nRadialBins = 300;
% param.sinogram_size.nAnglesBins = 300;
% param.nSubsets = 1;
% param.sinogram_size.span = 0;    % Span 0 is multi slice 2d.
% param.sinogram_size.nRings = refAct.ImageSize(3); % For multi slice 2d, the number of planes needs to be equal to the number of singorams.
%% USE THIS FOR 3D
% Change sinogram size:
param.sinogram_size.nRadialBins = 300;
param.sinogram_size.nAnglesBins = 300;
param.nSubsets = 10;
param.sinogram_size.span = 11;    % Span 0 is multi slice 2d.
param.sinogram_size.nRings = refAct.ImageSize(3); % For multi slice 2d, the number of planes needs to be equal to the number of singorams.
%%

PET.Revise(param);

% Counts to simulate:
counts = 1e9;
randomsFraction = 0.1;
scatterFraction = 0.35;
truesFraction = 1 - randomsFraction - scatterFraction;

% Geometrical projection:
y = PET.P(tAct); % for any other span

% Multiplicative correction factors:
acf= PET.ACF(tMu, refAct);
% Convert into factors:
af = acf;
af(af~=0) = 1./ af(af~=0);
% Introduce poission noise:
y = y.*af;
scale_factor = counts*truesFraction/sum(y(:));
y_poisson = poissrnd(y.*scale_factor);

% Additive factors:
r = PET.R(counts*randomsFraction);
% Poisson distribution:
r = poissrnd(r);

scatterFraction = 0.35;
counts_scatter = counts*scatterFraction;
s_withoutNorm = PET.S(y);
scale_factor_scatter = counts_scatter/sum(s_withoutNorm(:));
s_withoutNorm = s_withoutNorm .* scale_factor_scatter;
% noise for the scatter:
s = poissrnd(s_withoutNorm);
% Add randoms and scatter@
simulatedSinogram = y_poisson + s + r;


%% RECONSTRUCT
sensImage = PET.Sensitivity(af);
recon = PET.ones();
recon = PET.OPOSEM(simulatedSinogram,s+r, sensImage,recon, ceil(60/PET.nSubsets));

