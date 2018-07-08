%  *********************************************************************
%  Visualization tools.  
%  Author: Martín Belzunce. Kings College London.
%  Date: 05/12/2017
%  *********************************************************************

% Examples of how to use the function CreateAnimationComparingTwoMethods.

% Firs for individual type of animation and then for a sequence.
%% LOAD DATA
% Mlem data was saved in a cell array with the 3d image for every
% iteration:
load('/data/Results/FDG_11/MLEM/mlem_resampled.mat');
% Mr guided data was saved into an structure that contains the field
% "images" that is a cell array with the reconstructed images for every
% iteration, so images are in nonlocal_lange_bowsher_mr_voxels.images{i}.
% Also the images were trimmed before saving, therefore they are smaller
% than mlem.
load('/data/Results/FDG_11/nonlocal_lange_bowsher_psf/delta4.00e-02_B40_regparam2.23e+03/nonlocal_lange_bowsher_mr_voxels.mat');
%% GENERATE ANIMATION WITH FUNCTION
% For the animation both final images needs to have the same size, in this
% case that can be sorted out by chosing the rows and cols that need to be
% visualized.
% First go only through iterations
frame_time = 0.1;
opts0.slice = 80;
opts0.rows1 = 200:480; % rows and cols for image 1 (MLEM), chose values so your image is centred in the brain and also needs to match the size of the second method.
opts0.cols1 = 200:480;
opts0.rows2 = 1 : size(nonlocal_lange_bowsher_mr_voxels.images{1},1); % Second image was trimmed before saving to rows and cols 200:480, so now I use all the rows and cols available.
opts0.cols2 = 1 : size(nonlocal_lange_bowsher_mr_voxels.images{1},2);
opts0.cellArrayElements = [1:19 20:2:28 30:4:70 80:10:150 160:20:300]; % Iterations I want to show
opts0.title = '                           Positron Emission Tomography (PET)';
opts0.labelMethod1 = '                 Standard PET';
opts0.labelMethod2 = '              MRI-assisted PET';
opts0.fontName = 'arial';
opts0.fontSize = 60;
scale =1 ;
outputSize = [1080 1920];
[frames, outputFilename] = CreateAnimationComparingTwoMethods(mlem_resampled, nonlocal_lange_bowsher_mr_voxels.images, [0], frame_time, gray, scale, outputSize, sprintf('/data/Results/FDG_11/OnlyIterationsSlice_%d.gif', opts0.slice), opts0);
%% Only through slices
opts1.slices = 25:100;
opts1.rows1 = 200:480;
opts1.cols1 = 200:480;
opts1.rows2 = 1 : size(nonlocal_lange_bowsher_mr_voxels.images{1},1);
opts1.cols2 = 1 : size(nonlocal_lange_bowsher_mr_voxels.images{1},2);
opts1.cellArrayElement = 300; % Iteration I want to show
opts1.title = '                           Positron Emission Tomography (PET)';
opts1.labelMethod1 = '                 Standard PET';
opts1.labelMethod2 = '              MRI-assisted PET';
opts1.fontName = 'arial';
opts1.fontSize = 60;
scale =1 ;
[frames, outputFilename] = CreateAnimationComparingTwoMethods(mlem_resampled, nonlocal_lange_bowsher_mr_voxels.images, [1], frame_time, gray, scale, outputSize, sprintf('/data/Results/FDG_11/OnlySlicesIteration_%d.gif', opts1.cellArrayElement), opts1);
%% Only through MIPs
opts2.angles = [0:45:357];% 0:3:357];
opts2.rows1 = 200:480;
opts2.cols1 = 200:480;
opts2.rows2 = 1 : size(nonlocal_lange_bowsher_mr_voxels.images{1},1);
opts2.cols2 = 1 : size(nonlocal_lange_bowsher_mr_voxels.images{1},2);
opts2.rows2 = 200:480;
opts2.cols2 = 200:480;
opts2.cellArrayElement = 300; % Iteration I want to show
opts2.title = '                           Positron Emission Tomography (PET)';
opts2.labelMethod1 = '                 Standard PET';
opts2.labelMethod2 = '              MRI-assisted PET';
opts2.fontName = 'arial';
opts2.fontSize = 60;
scale = 1;
[frames, outputFilename] = CreateAnimationComparingTwoMethods(mlem_resampled, mlem_resampled, [2], frame_time, gray, scale, outputSize, 'RotatingMips.gif', opts2);
%% Create a Sequence
opts1.slices = [opts0.slice:100 25:opts0.slice-1]; % Change the slices to start from the same slice as the iteration update
opts2.transition = 6;
[frames, outputFilename] = CreateAnimationComparingTwoMethods(mlem_resampled, nonlocal_lange_bowsher_mr_voxels.images, [0 1 2], frame_time, gray, scale, outputSize, '/data/Results/FDG_11/FullSequence.gif', opts0, opts1, opts2);
v = VideoWriter('/data/Results/FDG_11/FullSequence.avi', 'Uncompressed AVI');
v.FrameRate = 10;
open(v);
for i = 1:size(frames,3);writeVideo(v,frames(:,:,i));end
close(v)

v = VideoWriter('/data/Results/FDG_11/FullSequence_2.avi', 'Uncompressed AVI');
v.FrameRate = 20;
open(v);
for i = 1:size(frames,3);writeVideo(v,frames(:,:,i));end
close(v)
% 
% v = VideoWriter('/data/Results/FDG_11/FullSequence.mpg', 'Motion JPEG 2000');
% v.FrameRate = 10;
% open(v);
% for i = 1:size(frames,3);writeVideo(v,frames(:,:,i));end
% close(v)

[frames, outputFilename] = CreateAnimationComparingTwoMethods(mlem_resampled, nonlocal_lange_bowsher_mr_voxels.images, [0 1 2], hot, scale, outputSize, '/data/Results/FDG_11/FullSequence_hot.gif', opts0, opts1, opts2);
gray = colormap(gray);
invgray = gray(end:-1:1,:);
[frames, outputFilename] = CreateAnimationComparingTwoMethods(mlem_resampled, nonlocal_lange_bowsher_mr_voxels.images, [0 1 2], invgray, scale, outputSize, '/data/Results/FDG_11/FullSequence_inv.gif', opts0, opts1, opts2);

frame_time = 0.05;
[frames, outputFilename] = CreateAnimationComparingTwoMethods(mlem_resampled, nonlocal_lange_bowsher_mr_voxels.images, [0 1 2], frame_time, gray, scale, outputSize, '/data/Results/FDG_11/FullSequence2.gif', opts0, opts1, opts2);
