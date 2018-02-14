%% EXAMPLE MLEM RADON MATLAB
clear all
close all
load BrainWeb.mat

% load g_truth
p.scanner = '2D_radon'; %'mMR'
p.method = 'otf_matlab'; %, 'otf_Gpu', 'pre-computed_matlab'
p.PSF.type = 'shift-invar';
p.PSF.Width = 4; %mm
%p.nSubsets = 1;
p.nIter = 30;

PET = classGpet(p);

y = PET.P(g_truth{1});
counts = 1e7/sum(y(:));

y = poissrnd(y.*counts)./counts;
Sensi = PET.PT(1);
% MLEM
x = PET.ones;
figure
for i = 1:p.nIter

    x = x./(Sensi).*PET.PT(y./(PET.P(x)+eps));
    drawnow, imagesc(x), axis image, colormap gray
end
%% EXAMPLE OSEM RADON MATLAB
load BrainWeb.mat

% load g_truth
p.scanner = '2D_radon'; %'mMR'
p.method = 'otf_matlab'; %, 'otf_Gpu', 'pre-computed_matlab'
p.PSF.type = 'shift-invar';
p.PSF.Width = 4; %mm
p.nSubsets = 10;
p.nIter = 3;
PET = classGpet(p);

y = PET.P(g_truth{1});
counts = 1e7/sum(y(:));

y = poissrnd(y.*counts)./counts;

% OSEM
x = PET.ones;
figure
for i = 1:PET.nIter
    for s= 1:PET.nSubsets
        ys=y(:,PET.sinogram_size.subsets(:,s));
        x = x./(PET.PT(1)).*PET.PT(ys./(PET.P(x,s)+eps),s);
    end

    drawnow, imagesc(x), axis image, colormap gray
end
%
% CHANGING SUBSET CONFIGURATION
PET.set_subsets(5);
x = PET.ones;
figure
for i = 1:PET.nIter
    for s= 1:PET.nSubsets
        ys=y(:,PET.sinogram_size.subsets(:,s));
        x = x./(PET.PT(1)).*PET.PT(ys./(PET.P(x,s)+eps),s);
    end

    drawnow, imagesc(x), axis image, colormap gray
end
%% EXAMPLE MLEM MARTIN PROJECTOR (DEFAULT SPAN 11)
clear all
close all
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
% CUDA PATH
cudaPath = '/usr/local/cuda/';
setenv('PATH', [getenv('PATH') sepEnvironment cudaPath pathBar 'bin']);
setenv('LD_LIBRARY_PATH', [getenv('LD_LIBRARY_PATH') sepEnvironment cudaPath pathBar 'lib64']);
% APIRL PATH
apirlPath = '/home/mab15/workspace/apirl-code/trunk/';
addpath(genpath([apirlPath pathBar 'matlab']));
setenv('PATH', [getenv('PATH') sepEnvironment apirlPath pathBar 'build' pathBar 'bin']);
setenv('LD_LIBRARY_PATH', [getenv('LD_LIBRARY_PATH') sepEnvironment apirlPath pathBar 'build' pathBar 'bin']);

fullFilename = '/home/mab15/workspace/KCL/Biograph_mMr/Mediciones/BRAIN_PETMR/SINOGRAMS/PET_ACQ_68_20150610155347_ima_AC_000_000.v.hdr';
[g_truth, refImage, bedPosition_mm, info]  = interfileReadSiemensImage(fullFilename);

objGpet.scanner = 'mMR';
                objGpet.method =  'otf_siddon_gpu';
                objGpet.PSF.type = 'none';
                objGpet.radialBinTrim = 0;
                objGpet.Geom = '';
PET = classGpet(objGpet);
% Change the image sie, to the one of the phantom:
init_image_properties(PET, refImage);

% EXAMPLE NCFS
%NCF = PET.NCF;
% Simulate sinogram
y = PET.P(g_truth);
counts = 1e8/sum(y(:));

y = poissrnd(y.*counts)./counts;
Sensi = PET.PT(1);

% MLEM
x = PET.ones;
sliceToShow = 50;
figure
for i = 1:PET.nIter

    x = x./(Sensi).*PET.PT(y./(PET.P(x)+eps));
    drawnow, imagesc(x(:,:,sliceToShow)), axis image, colormap gray
end
%% EXAMPLE MLEM MARTIN PROJECTOR (ANY SPAN)
clear all
close all
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
% CUDA PATH
cudaPath = '/usr/local/cuda/';
setenv('PATH', [getenv('PATH') sepEnvironment cudaPath pathBar 'bin']);
setenv('LD_LIBRARY_PATH', [getenv('LD_LIBRARY_PATH') sepEnvironment cudaPath pathBar 'lib64']);
% APIRL PATH
apirlPath = '/workspaces/Martin/apirl-code/trunk/';
addpath(genpath([apirlPath pathBar 'matlab']));
setenv('PATH', [getenv('PATH') sepEnvironment apirlPath pathBar 'build' pathBar 'bin']);
setenv('LD_LIBRARY_PATH', [getenv('LD_LIBRARY_PATH') sepEnvironment apirlPath pathBar 'build' pathBar 'bin']);

[g_truth, refImage]  = interfileRead('phantom.h33');

PET = classGpet();
% Change the image sie, to the one of the phantom:
init_image_properties(PET, refImage);
% Change the span size:
span = 51;
numRings = 64;
maxRingDifference = 60;

init_sinogram_size(PET, 51, 64, 60);

y = PET.P(g_truth);
counts = 1e8/sum(y(:));

y = poissrnd(y.*counts)./counts;
Sensi = PET.PT(1);
%
% MLEM
x = PET.ones;
sliceToShow = 50;
figure
for i = 1:PET.nIter

    x = x./(Sensi).*PET.PT(y./(PET.P(x)+eps));
    x(Sensi==0) = 0;
    drawnow, imagesc(x(:,:,sliceToShow)), axis image, colormap gray
    pause;
end

%% EXAMPLE MLEM MARTIN PROJECTOR 2D
clear all
close all

load BrainWeb.mat
% Resize image:

g_truth{1} = imresize(g_truth{1}, [344 344]);
g_truth{1}(g_truth{1}<0) = 0;
% load g_truth
p.scanner = '2D_mMR'; %'mMR'
p.method = 'otf_siddon_cpu'; %, 'otf_Gpu', 'pre-computed_matlab'
p.PSF.type = 'none';
p.PSF.Width = 4; %mm
%p.nSubsets = 1;
p.nIter = 30;

PET = classGpet(p);

y = PET.P(g_truth{1});
counts = 1e7/sum(y(:));

y = poissrnd(y.*counts)./counts;
Sensi = PET.PT(1);
% MLEM
x = PET.ones;
figure
for i = 1:p.nIter

    x = x./(Sensi).*PET.PT(y./(PET.P(x)+eps));
    x(Sensi==0) = 0;
    drawnow, imagesc(x), axis image, colormap gray
end
