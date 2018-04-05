function reconAPIRL = MLEM(fsSigma, fsPet, fsT1, varargin)
%% Simulate reconstruction of specified random-structured brainweb subject
%{
# Arguments:
fsSigma = 1, fsPet = 1.5, fsT1 = 0.75
fsSigma = 1, fsPet = 2, fsT1 = 1
fsSigma = 0, fsPet = 0.4, fsT1 = 0.2
# Options:
subj = 'subject_04'
counts = 500e6
numTumours = -1
# Note that counts < 1 will perform noise-free reconstruction
%}

subj = 'subject_04';
counts = 500e6;
numTumours = -1;
if nargin > 3, subj = varargin{1}; end
if nargin > 4, counts = varargin{2}; end
if nargin > 5, numTumours = varargin{3}; end

%% EXAMPLE MLEM MARTIN PROJECTOR (ANY SPAN)
%clear all, close all
set_framework_environment();
% set_framework_environment(basePath, binaryPath);

%% SIMULATE A BRAIN PHANTOM WITH ATTENUATION, NORMALIZATION, RANDOMS AND SCATTER
%[sinogram, delayedSinogram, structSizeSino3d] = interfileReadSino('E:\PatientData\FDG\PETSinoPlusUmap-Converted\PETSinoPlusUmap-00\PETSinoPlusUmap-00-sino-uncomp.s.hdr');
%load BrainMultiMaps_mMR.mat;
if strfind(subj, 'AD_')
  addpath('output/')
  MultiMaps_Ref = load(['reconAPIRL_real-' subj '.mat'], ['reconAPIRL']);
  MultiMaps_Ref = MultiMaps_Ref.reconAPIRL;
  MultiMaps_Ref.uMap = reshape(fread(fopen(['data-LM-00-umap-' subj '.v']), 'single'), [344 344 127]);
  MultiMaps_Ref.PET = MultiMaps_Ref.PET_psf_sharp;
  MultiMaps_Ref.PET = permute(MultiMaps_Ref.PET, [2 1 3]);
  MultiMaps_Ref.T1 = permute(MultiMaps_Ref.T1, [2 1 3]);
else
  addpath('brainweb.raws/')
  MultiMaps_Ref = readPetMr(fsSigma, fsPet, fsT1, subj);
end
%%
tAct = permute(MultiMaps_Ref.PET, [2 1 3]);
tAct = tAct(end:-1:1,:,:);
if numTumours
  if strfind(subj, 'AD_')
    tIntensity = max(max(max(medfilt3(tAct, [1 1 1] .* 3)))) / max(tAct(:))
    wmm = size(tAct, 1)*MultiMaps_Ref.voxelSize_mm(1);  % image width/[mm]
    tAct = permute(addTumours(permute(tAct, [3 2 1]), wmm, 5, tIntensity, numTumours), [3 2 1]);
  else
    tAct = permute(addTumours(permute(tAct, [3 2 1]), ...
                              344*2.08626, 5, 1.5, numTumours), [3 2 1]);
  end
end
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

% Change the span size:
span = 11
numRings = 64;
maxRingDifference = 60;

% Counts to simulate:
randomsFraction = 0.2  %.1
scatterFraction = 0.2  %.35
truesFraction = 1 - randomsFraction - scatterFraction;

noPSFpsf = 2.5;
psfPSF = 4.5;
nitersMlem = 100
nitersPsf = 300

%% INIT CLASS GPET
PET.scanner = 'mMR';
PET.method =  'otf_siddon_gpu';
PET.PSF.type = 'shift-invar';  % none
PET.PSF.Width = noPSFpsf;
PET.nSubsets = 1;
PET.radialBinTrim = 0;
PET.Geom = '';
PET.random_algorithm = 'from_ML_singles_matlab';
PET.sinogram_size.span = 11;
%PET.verbosity = 1;
PET.tempPath = '/dev/shm/cc16/temp/mlem';

%%
%parpool('local', 6);
use_gpus=[1 2 4] + 1;
reconMLEM = cell(6, 1);
%load('output/reconMLEM.mat')
%load('output/reconMLEMPSF.mat')
parfor i=0:5
noise_realisation=1 + floor(i / 2);
% N.B.: keep gpu const for each noise_realisation
gpu = use_gpus(1 + mod(noise_realisation - 1, length(use_gpus)));
gpuDevice(gpu);

%%
PETmlem = PET;
PETmlem.tempPath = [PET.tempPath sprintf('%d%s%.3g%.3g%.3g%.3g', i, subj, counts, fsSigma, fsPet, fsT1) '/'];
PETpsf = PETmlem;
PETpsf.PSF.Width = psfPSF;
PETm = classGpet(PETmlem);
PETp = classGpet(PETpsf);

PETm = init_recon(PETm, refAct, span, numRings, maxRingDifference);
PETp = init_recon(PETp, refAct, span, numRings, maxRingDifference);

% Geometrical projection:
%options.PSF.Width = 0;
%PET.Revise(options);
%tActRD = PET.Gauss3DFilter(tAct, PET.PSF.Width);
y = PETp.P(tAct);
[ncf, acf, n, y, y_poisson] = poisson_recon(...
  PETp, tMu, refAct, y, counts, truesFraction);

if mod(i, 2)
recon = do_recon(PETp, nitersPsf, y, y_poisson, n, ncf, acf, ...
  tMu, refAct, counts, truesFraction, randomsFraction, scatterFraction);
else
recon = do_recon(PETm, nitersMlem, y, y_poisson, n, ncf, acf, ...
  tMu, refAct, counts, truesFraction, randomsFraction, scatterFraction);
end
rmdir(PETmlem.tempPath, 's');
disp(noise_realisation);

%reconMLEM{noise_realisation} = recon;
reconMLEM{i + 1} = recon;

end  % noise_realisation

reconMLEM = permute(reshape(cell2mat(reconMLEM), ...
  [344, 6, 344, 127]), [2 1 3 4]);
reconAPIRL.PET = tAct;
reconAPIRL.T1 = T1;
reconAPIRL.mlem = reconMLEM(1:2:end, :,:,:);
reconAPIRL.mlem_psf = reconMLEM(2:2:end, :,:,:);
reconAPIRL.voxelSize_mm = pixelSize_mm;
reconAPIRL.psf_mm = [noPSFpsf, psfPSF];
save(['output/reconAPIRL_brainweb' ...
      sprintf('_%s-S_%.3g-NP_%.3g-NT1_%.3g-C_%.3g_t%d', ...
              subj, fsSigma, fsPet, fsT1, counts, -numTumours) ...
      '.mat'], 'reconAPIRL', '-v7.3')
%save('output/reconMLEMPSF.mat', 'reconMLEMPSF', '-v7.3')

end  % function MLEM


function PET = init_recon(PET, refAct, span, numRings, maxRingDifference)
gaps = PET.gaps;
% Change the image size, to the one of the phantom:
PET.init_image_properties(refAct);
PET.init_sinogram_size(span, numRings, maxRingDifference);
end  % function init_recon


function [ncf, acf, n, y, y_poisson] = poisson_recon(...
  PET, tMu, refAct, y, counts, truesFraction)

% Multiplicative correction factors:
ncf = PET.NCF;
acf = PET.ACF(tMu, refAct);

% Convert into factors:
n = ncf; a = acf;
n(n~=0) = 1./ n(n~=0); a(a~=0) = 1./ a(a~=0);

% Introduce poission noise:
y = y.*n.*a;
scale_factor = abs(counts)*truesFraction/sum(y(:));
y_poisson = y.*scale_factor;
if counts > 0
  y_poisson = poissrnd(y.*scale_factor);
end

end  % poisson_recon


function recon = do_recon(PET, niters, y, y_poisson, ...
  n, ncf, acf, tMu, refAct, counts, truesFraction, randomsFraction, scatterFraction)

% Additive factors:
r = PET.R(abs(counts)*randomsFraction);
%r = PET.R(delayedSinogram);  % Without a delayed sinograms, just
scale_factor_randoms = abs(counts)*randomsFraction./sum(r(:));
% Poisson distribution:
%r = poissrnd(r.*scale_factor_randoms);

counts_scatter = abs(counts)*scatterFraction;
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
  sensImage, recon, ceil(niters/PET.nSubsets));
%recon = PET.OPMLEM(simulatedSinogram, additive, sensImage, recon, niters);
end  % function do_recon
