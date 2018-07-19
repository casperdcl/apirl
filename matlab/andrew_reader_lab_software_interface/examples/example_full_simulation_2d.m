%% EXAMPLE MLEM MARTIN PROJECTOR (ANY SPAN)
clear all
close all
set_framework_environment();
% set_framework_environment(basePath, binaryPath);
%% INIT CLASS GPET
PET.scanner = 'mMR';
PET.method =  'otf_siddon_cpu';
PET.PSF.type = 'none';
PET.radialBinTrim = 0;
PET.Geom = '';
PET.sinogram_size.span = -1;
PET.random_algorithm = 'from_ML_singles_matlab';
PET = classGpet(PET);
%% SIMULATE A BRAIN PHANTOM WITH ATTENUATION, NORMALIZATION, RANDOMS AND SCATTER
%[sinogram, delayedSinogram, structSizeSino3d] = interfileReadSino('/media/mab15/DATA_BACKUP/Scans/PatientData/FDG_Patient_01/e7/PETSinoPlusUmap-Converted/PETSinoPlusUmap-00/PETSinoPlusUmap-00-sino-uncomp.s.hdr');
load BrainMultiMaps_mMR.mat;
%delayedSinogram_2d = delayedSinogram(:,:,60);
% Get only one slice of the phantom:
tAct = permute(MultiMaps_Ref.PET(:,:,60), [2 1 3]);
tAct = tAct(end:-1:1,:);
tMu = permute(MultiMaps_Ref.uMap(:,:,60), [2 1 3]);
tMu = tMu(end:-1:1,:);
% T1 image:
T1 = permute(MultiMaps_Ref.T1(:,:,60), [2 1 3]);
T1 = T1(end:-1:1,:,:);
% T1 image:
T2 = permute(MultiMaps_Ref.T2(:,:,60), [2 1 3]);
T2 = T2(end:-1:1,:,:);
pixelSize_mm = [2.08625 2.08625 2.03125];
xLimits = [-size(tAct,2)/2*pixelSize_mm(2) size(tAct,2)/2*pixelSize_mm(2)];
yLimits = [-size(tAct,1)/2*pixelSize_mm(1) size(tAct,1)/2*pixelSize_mm(1)];
zLimits = [-size(tAct,3)/2*pixelSize_mm(3) size(tAct,3)/2*pixelSize_mm(3)];
refAct = imref3d([size(tAct) 1],xLimits,yLimits,zLimits);
refAt  = imref3d([size(tMu) 1],xLimits,yLimits,zLimits);

% Change the image sie, to the one of the phantom:
PET.init_image_properties(refAct);

% Counts to simulate:
counts = 1e8;
randomsFraction = 0.1;
scatterFraction = 0.35;
truesFraction = 1 - randomsFraction - scatterFraction;

% Geometrical projection:
y = PET.P(tAct);

% Multiplicative correction factors:
ncf = PET.NCF;
acf= PET.ACF(tMu, refAct);
% Convert into factors:
n = ncf; a = acf;
n(n~=0) = 1./ n(n~=0); a(a~=0) = 1./ a(a~=0);
% Introduce poission noise:
y = y.*n.*a;
%%
scale_factor = counts*truesFraction/sum(y(:));
y_poisson = poissrnd(y.*scale_factor);

% Additive factors:
r = PET.R(counts*randomsFraction);
%r = PET.R(delayedSinogram_2d);  % Without a delayed sinograms, just
scale_factor_randoms = counts*randomsFraction./sum(r(:));
% Poisson distribution:
r = poissrnd(r.*scale_factor_randoms);

scatterFraction = 0.35;
counts_scatter = counts*scatterFraction;
s_withoutNorm = PET.S(y);
scale_factor_scatter = counts_scatter/sum(s_withoutNorm(:));
s_withoutNorm = s_withoutNorm .* scale_factor_scatter;
% noise for the scatter:
s = poissrnd(s_withoutNorm.*n);
% Add randoms and scatter@
simulatedSinogram = y_poisson + s + r;