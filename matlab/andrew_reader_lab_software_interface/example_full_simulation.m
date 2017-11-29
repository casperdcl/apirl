%% EXAMPLE MLEM MARTIN PROJECTOR (ANY SPAN)
clear all
close all
set_framework_environment();
% set_framework_environment(basePath, binaryPath);
%% INIT CLASS GPET
PET.scanner = 'mMR';
PET.method =  'otf_siddon_cpu';
PET.PSF.type = 'none';
%PET.method =  'otf_siddon_gpu';
%PET.PSF.type = 'shift-invar';
%PET.PSF.Width = 4.3;
%PET.nSubsets = 1;
PET.radialBinTrim = 0;
PET.Geom = '';
PET.random_algorithm = 'from_ML_singles_matlab';
PET.sinogram_size.span = 11;
%PET.verbosity = 1;
PET = classGpet(PET);
%% SIMULATE A BRAIN PHANTOM WITH ATTENUATION, NORMALIZATION, RANDOMS AND SCATTER
%[sinogram, delayedSinogram, structSizeSino3d] = interfileReadSino('E:\PatientData\FDG\PETSinoPlusUmap-Converted\PETSinoPlusUmap-00\PETSinoPlusUmap-00-sino-uncomp.s.hdr');
load BrainMultiMaps_mMR.mat;
%%
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
% T1 image:
T1 = permute(MultiMaps_Ref.T1, [2 1 3]);
T1 = T1(end:-1:1,:,:);

gaps = PET.gaps;

% Change the image size, to the one of the phantom:
PET.init_image_properties(refAct);

% Change the span size:
span = 11
numRings = 64;
maxRingDifference = 60;
PET.init_sinogram_size(span, numRings, maxRingDifference);

% Counts to simulate:
counts = 1e9;
randomsFraction = 0.1;
scatterFraction = 0.35;
truesFraction = 1 - randomsFraction - scatterFraction;

% Geometrical projection:
y = PET.P(tAct);

% Multiplicative correction factors:
ncf = PET.NCF;
acf = PET.ACF(tMu, refAct);
% Convert into factors:
n = ncf; a = acf;
n(n~=0) = 1./ n(n~=0); a(a~=0) = 1./ a(a~=0);
% Introduce poission noise:
y = y.*n.*a;
scale_factor = counts*truesFraction/sum(y(:));
y_poisson = poissrnd(y.*scale_factor);

% Additive factors:
r = PET.R(counts*randomsFraction); 
%r = PET.R(delayedSinogram);  % Without a delayed sinograms, just
scale_factor_randoms = counts*randomsFraction./sum(r(:));
% Poisson distribution:
r = poissrnd(r.*scale_factor_randoms);

counts_scatter = counts*scatterFraction;
s_withoutNorm = PET.S(y);
scale_factor_scatter = counts_scatter/sum(s_withoutNorm(:));
s_withoutNorm = s_withoutNorm .* scale_factor_scatter;
% noise for the scatter:
s = poissrnd(s_withoutNorm.*n);
% Add randoms and scatter@
simulatedSinogram = y_poisson + s + r;

%% SENSITIVITY IMAGE
anf = acf .* ncf;
anf(anf~=0) = 1./anf(anf~=0);
sensImage = PET.Sensitivity(anf);
%% OP-OSEM
% additive term:
additive = (r + s).*ncf.*acf; % (randoms +scatter)./(afs*nfs) = (randoms+scatter)+
recon = PET.ones();
recon = PET.OPOSEM(simulatedSinogram, additive, ...
  sensImage, recon, ceil(60/PET.nSubsets));
