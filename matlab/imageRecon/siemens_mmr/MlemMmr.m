%  *********************************************************************
%  Reconstruction Framework for Siemens Biograph mMR.  
%  Autor: Martín Belzunce. Kings College London.
%  Fecha de Creación: 06/05/2015
%  *********************************************************************
%  This function reconstructs a mmr sinogram that can be in apirl interfile
%  or siemens interfile. It receives the
%  header of the sinogram, a normalization file, an attenuation map and the outputh
%  path for the results. It returns a volume.
% The filename for the attenuation map must include only the path and the
% first parte of the name as is stored by the mMR. For example:
% attMapBaseFilename = 'path/PET_ACQ_190_20150220152253', where the
% filenames are: 
%   - 'path/PET_ACQ_190_20150220152253_umap_human_00.v.hdr'
%   - 'path/PET_ACQ_190_20150220152253_umap_hardware_00.v.hdr'
%  It receives the span to be used in the reconstruction. The sinogram
%  It receives two optional parameters
%  - pixelSize_mm: that is a vector with 3 elements [pixelSizeX_mm
%  pixelSizeY_mm pixelSizeZ_mm].
%  - numIterations: number of iterations.
%
%  The span of the reconstruction is the same of the sinogram.
%
% Examples:
%   volume = MlemMmr(sinogramFilename, normFilename, attMapBaseFilename, outputPath, pixelSize_mm, numInputIterations, optUseGpu)
function [volume overall_ncf_3d acfsSinogram randoms scatter] = MlemMmr(sinogramFilename, span, normFilename, attMapHumanFilename, attMapHardwareFilename, correctRandoms, correctScatter, outputPath, pixelSize_mm, numIterations, saveInterval, useGpu, stirMatlabPath, removeTempFiles)

mkdir(outputPath);
% Check what OS I am running on:
if(strcmp(computer(), 'GLNXA64'))
    os = 'linux';
    pathBar = '/';
elseif(strcmp(computer(), 'PCWIN') || strcmp(computer(), 'PCWIN64'))
    os = 'windows';
    pathBar = '\';
else
    disp('OS not compatible');
    return;
end
% Check if we have received pixel size:
if nargin ~= 14
    [volume overall_ncf_3d acfsSinogram randoms scatter] = MlemMmr(sinogramFilename, span, normFilename, attMapHumanFilename, attMapHardwareFilename, correctRandoms, correctScatter, outputPath, pixelSize_mm, numIterations, saveInterval, useGpu, stirMatlabPath, removeTempFiles)
end

imageSize_mm = [600 600 257.96875];
imageSize_pixels = ceil(imageSize_mm./pixelSize_mm);

%% READING THE SINOGRAMS
disp('Read input sinogram...');
% Read the sinograms:
if isstr(sinogramFilename)
    [sinograms, delayedSinograms, structSizeSino3dSpan1] = interfileReadSino(sinogramFilename);
    if structSizeSino3dSpan1.span == 1
        % Convert to span:
        [sinograms, structSizeSino3d] = convertSinogramToSpan(sinograms, structSizeSino3dSpan1, span);  
    else
        structSizeSino3d = structSizeSino3dSpan1;
        warning('The span parameter will be ignored because the input sinogram is not span 1.');
    end
else
    % binary sinogram (span-1?) or span span
    if size(sinogramFilename) == [344 252 4084]
        structSizeSino3d = getSizeSino3dFromSpan(344, 252, 64, ...
            296, 256, 1, 60);
        sinograms = sinogramFilename;
    else
        structSizeSino3d = getSizeSino3dFromSpan(344, 252, 64, ...
            296, 256, span, 60);
        if size(sinogramFilename) == [344 252 sum(structSizeSino3d.sinogramsPerSegment)]
            sinograms = sinogramFilename;
        else
            error('Invalid sinogram matrix size. It must be span 1 or "span".');
        end
    end
end
sinogramFilename = [outputPath pathBar 'sinogram'];
% Write the input sinogram:
interfileWriteSino(single(sinograms), sinogramFilename, structSizeSino3d);
%% CREATE INITIAL ESTIMATE FOR RECONSTRUCTION
disp('Creating inital image...');
% Inititial estimate:
initialEstimate = ones(imageSize_pixels, 'single');
filenameInitialEstimate = [outputPath pathBar 'initialEstimate'];
interfilewrite(initialEstimate, filenameInitialEstimate, pixelSize_mm);
%% NORMALIZATION FACTORS
if isstr(normFilename)
    disp('Computing the normalization correction factors...');
    % ncf:
    [overall_ncf_3d, scanner_time_invariant_ncf_3d, scanner_time_variant_ncf_3d, acquisition_dependant_ncf_3d, used_xtal_efficiencies, used_deadtimefactors, used_axial_factors] = ...
       create_norm_files_mmr(normFilename, [], [], [], [], structSizeSino3d.span);
    % invert for nf:
    overall_nf_3d = overall_ncf_3d;
    overall_nf_3d(overall_ncf_3d ~= 0) = 1./overall_nf_3d(overall_ncf_3d ~= 0);
else
    if ~isempty(normFilename)   % if empty no norm is used, if not used the matrix in normFilename proving is of the same size of the sinogram.
        if size(normFilename) ~= size(sinograms)
            error('The size of the normalization correction factors is incorrect.')
        end
        disp('Using the normalization correction factors received as a parameter...');
        overall_ncf_3d = normFilename;
        clear normFilename;
        % invert for nf:
        overall_nf_3d = overall_ncf_3d;
        overall_nf_3d(overall_ncf_3d ~= 0) = 1./overall_nf_3d(overall_ncf_3d ~= 0);
    else
        overall_ncf_3d = ones(size(sinograms));
        overall_nf_3d = overall_ncf_3d;
    end
end
   
%% ATTENUATION MAP
if isstr(attMapHumanFilename)
    % Check if its attenuation from siemens or a post processed image:
    if ~strcmp(attMapHumanFilename(end-3:end),'.h33')
        disp('Computing the attenuation correction factors from mMR mu maps...');
        % Read the attenuation map and compute the acfs.
        [attenMap_human, refAttenMapHum, bedPosition_mm, info]  = interfileReadSiemensImage(attMapHumanFilename);
        if isstr(attMapHardwareFilename)
            [attenMap_hardware, refAttenMapHard, bedPosition_mm, info]  = interfileReadSiemensImage(attMapHardwareFilename);
        end
        imageSizeAtten_mm = [refAttenMapHum.PixelExtentInWorldY refAttenMapHum.PixelExtentInWorldX refAttenMapHum.PixelExtentInWorldZ];
        % Compose both images:
        attenMap = attenMap_hardware + attenMap_human;
        % I need to translate because siemens uses an slightly displaced
        % center (taken from dicom images, the first pixel is -359.8493 ,-356.8832 
        %displacement_mm = [-1.5 -imageSizeAtten_mm(2)*size(attenMap_human,1)/2+356.8832 0];
        %[attenMap, Rtranslated] = imtranslate(attenMap, refAttenMapHum, displacement_mm,'OutputView','same');
        
    else
        disp('Computing the attenuation correction factors from post processed APIRL mu maps...');
        attenMap_human = interfileRead(attMapHumanFilename);
        attenMap = attenMap_human;
        infoAtten = interfileinfo(attMapHumanFilename); 
        imageSizeAtten_mm = [infoAtten.ScalingFactorMmPixel1 infoAtten.ScalingFactorMmPixel2 infoAtten.ScalingFactorMmPixel3];
    end   

    % Create ACFs of a computed phatoms with the linear attenuation
    % coefficients:
    acfFilename = ['acfsSinogram'];
    acfsSinogram = createACFsFromImage(attenMap, imageSizeAtten_mm, outputPath, acfFilename, structSizeSino3d, 0, useGpu);
    % After the projection read the acfs:
    acfFilename = [outputPath acfFilename];
else
    if isempty(attMapHumanFilename)
        acfFilename = '';
        acfsSinogram = ones(size(sinograms));
    elseif size(attMapBaseFilename) == size(sinograms)
        % I got the acfs:
        acfsSinogram = attMapHumanFilename;
    end
end
%% GENERATE AND SAVE ATTENUATION AND NORMALIZATION FACTORS AND CORECCTION FACTORS
disp('Generating the ANF sinogram...');
% Save:
outputSinogramName = [outputPath 'NF'];
interfileWriteSino(single(overall_nf_3d), outputSinogramName, structSizeSino3d);

% Compose with acfs:
atteNormFactors = overall_nf_3d;
atteNormFactors(acfsSinogram ~= 0) = overall_nf_3d(acfsSinogram ~= 0) ./acfsSinogram(acfsSinogram ~= 0);
anfFilename = [outputPath 'ANF'];
interfileWriteSino(single(atteNormFactors), anfFilename, structSizeSino3d);
clear atteNormFactors;
%clear normFactorsSpan11;
%% RANDOMS CORRECTION
randoms = zeros(size(sinograms));
if numel(size(correctRandoms)) == numel(size(sinograms))
    if size(correctRandoms) == size(sinograms)
        % The input is the random estimate:
        randoms = correctRandoms;
    else
        numIterations = 5;
        [randoms, singles] = estimateRandomsFromDelayeds(delayedSinograms, structSizeSino3dSpan1, numIterations, structSizeSino3d.span);
    end
else
    % If not, we expect a 1 to estimate randoms or a 0 to not cprrect for
    % them:
    if(correctRandoms)
        % Stir computes randoms that are already normalized:
        numIterations = 5;
        [randoms, singles] = estimateRandomsFromDelayeds(delayedSinograms, structSizeSino3dSpan1, numIterations, structSizeSino3d.span);
    end
end
%% SCATTER CORRECTION
% Check if we hve the scatter estimate as a parameter. If not we compute it
% at the end:
scatter = zeros(size(sinograms));
if numel(size(correctScatter)) == numel(size(sinograms))
    if size(correctScatter) == size(sinograms)
        scatter = correctScatter;
    end
end
%% ADDITIVE SINOGRAM
additive = (randoms .* overall_ncf_3d  + scatter).* acfsSinogram; % (randoms +scatter.*norm)./(attenFactors*nprmFactrs)
% write additive singoram:
additiveFilename = [outputPath 'additive'];
interfileWriteSino(single(additive), additiveFilename, structSizeSino3d);
%% GENERATE MLEM RECONSTRUCTION FILES FOR APIRL
if useGpu == 0
    disp('Mlem reconstruction...');
    saveIntermediate = 0;
    outputFilenamePrefix = [outputPath 'reconImage'];
    filenameMlemConfig = [outputPath 'mlem.par'];
    CreateMlemConfigFileForMmr(filenameMlemConfig, [sinogramFilename '.h33'], [filenameInitialEstimate '.h33'], outputFilenamePrefix, numIterations, [],...
        saveInterval, saveIntermediate, [anfFilename '.h33'], [additiveFilename '.h33']);
    % Execute APIRL: 
    status = system(['MLEM ' filenameMlemConfig]) 

else
    
    disp('cuMlem reconstruction...');
    saveIntermediate = 0;
    outputFilenamePrefix = [outputPath 'reconImage'];
    filenameMlemConfig = [outputPath 'cumlem.par'];
    CreateCuMlemConfigFileForMmr(filenameMlemConfig, [sinogramFilename '.h33'], [filenameInitialEstimate '.h33'], outputFilenamePrefix, numIterations, [],...
        saveInterval, saveIntermediate, [anfFilename '.h33'], [additiveFilename '.h33'], 0, 576, 576, 512);
    % Execute APIRL: 
    status = system(['cuMLEM ' '"' filenameMlemConfig '"'])
end
%% READ RESULTS
% Read interfile reconstructed image:
volume = interfileRead([outputFilenamePrefix '_final.h33']);
%% REMOVE TEMPORARY FILES
if removeTempFiles
    delete([outputPath '*sensitivity*']);
    delete([outputPath 'acfs*']);
    delete([outputPath 'additive*']);
    delete([outputPath 'NF*']);
end
%% IF NEEDS TO CORRECT SCATTER, ESTIMATE IT AFTER RECONSTRUCTION:
if (numel(correctScatter) == 1) 
    if (correctScatter == 1)
        thresholdForTail = 1.01;
        stirScriptsPath = [stirMatlabPath pathBar 'scripts'];
        % The scatter needs the image but also the acf to scale, and in the case of
        % the mr is better if this acf include the human?
        if isstr(attMapHumanFilename)
            acfFilename = 'acfsOnlyHuman';
            [attenMap_human, refAttenMapHum, bedPosition_mm, info]  = interfileReadSiemensImage(attMapHumanFilename);
            acfsOnlyHuman = createACFsFromImage(attenMap_human, imageSizeAtten_mm, outputPath, acfFilename, structSizeSino3d, 0, useGpu);
            acfFilename = [outputPath acfFilename];
        else
            acfsOnlyHuman = acfsSinogram;
        end

        % The emission sinogram needs to be normalized and corrected for randoms:
        outputPathScatter = [outputPath pathBar 'Scatter1' pathBar];
        if ~isdir(outputPathScatter)
            mkdir(outputPathScatter);
        end
        [scatter_1, structSizeSino, mask] = estimateScatterWithStir(volume, attenMap, pixelSize_mm, sinograms, randoms, overall_ncf_3d, acfsOnlyHuman, structSizeSino3d, outputPathScatter, stirScriptsPath, thresholdForTail);

        % Reconstruct again with the scatter:
        interfileWriteSino(single(scatter_1), [outputPathScatter 'scatter'], structSizeSino3d);
        % Remove files:
        if removeTempFiles
            delete([outputPath 'acfs*']);
            delete([outputPathScatter 'tail*']);
            delete([outputPathScatter 'acfs*']);
            delete([outputPathScatter 'emission*']);
            delete([outputPathScatter 'scatter.s']);
            delete([outputPathScatter 'scatter.hs']);
        end
        
        % Normalize the scatter:
        normScatter = scatter_1 .* overall_nf_3d;
        % Plot profiles to test:
        profileSinogram = sum(sinograms(:,126,:),3);
        profileRandoms = sum(randoms(:,126,:),3);
        profileNormScatter = sum(normScatter(:,126,:),3);
        figure;
        subplot(1,3,1);
        title('Scatter Iteration 1');
        plot([profileSinogram profileRandoms profileNormScatter (profileRandoms+profileNormScatter)]);
        legend('Sinogram', 'Randoms', 'Scatter', 'Randoms+Scatter');

        outputPathScatter = [outputPath pathBar 'Scatter2' pathBar];
        if ~isdir(outputPathScatter)
            mkdir(outputPathScatter);
        end
        
        % New dditive sinogram:
        additive = (randoms .* overall_ncf_3d  + scatter_1).* acfsSinogram; % (randoms +scatter.*norm)./(attenFactors*nprmFactrs)
        % write additive singoram:
        additiveFilename = [outputPathScatter 'additive'];
        interfileWriteSino(single(additive), additiveFilename, structSizeSino3d);

        if useGpu == 0
            disp('Mlem reconstruction...');
            saveIntermediate = 0;
            outputFilenamePrefix = [outputPathScatter 'reconImage'];
            filenameMlemConfig = [outputPathScatter 'mlem.par'];
            CreateMlemConfigFileForMmr(filenameMlemConfig, [sinogramFilename '.h33'], [filenameInitialEstimate '.h33'], outputFilenamePrefix, numIterations, [],...
                saveInterval, saveIntermediate, [anfFilename '.h33'], [additiveFilename '.h33']);
            % Execute APIRL: 
            status = system(['OSEM ' filenameMlemConfig]) 
        else
            disp('cuMlem reconstruction...');
            saveIntermediate = 0;
            outputFilenamePrefix = [outputPathScatter 'reconImage'];
            filenameMlemConfig = [outputPathScatter 'cuosem.par'];
            CreateCuMlemConfigFileForMmr(filenameMlemConfig, [sinogramFilename '.h33'], [filenameInitialEstimate '.h33'], outputFilenamePrefix, numIterations, [],...
                saveInterval, saveIntermediate, [anfFilename '.h33'], [additiveFilename '.h33'], 0, 576, 576, 512);
            % Execute APIRL: 
            status = system(['cuMLEM ' filenameMlemConfig]) 
        end
        volume = interfileRead([outputFilenamePrefix '_final.h33']);
        % Remove files:
        if removeTempFiles
            delete([outputPathScatter 'tail*']);
            delete([outputPathScatter 'acfs*']);
            delete([outputPathScatter 'emission*']);
            delete([outputPathScatter 'scatter.s']);
            delete([outputPathScatter 'scatter.hs']);
            delete([outputPathScatter '*sensitivity*']);
            delete([outputPathScatter 'additive*']);
        end
        
        [scatter_2, structSizeSino, mask] = estimateScatterWithStir(volume, attenMap, pixelSize_mm, sinograms, randoms, overall_ncf_3d, acfsOnlyHuman, structSizeSino3d, outputPathScatter, stirScriptsPath, thresholdForTail);
        interfileWriteSino(single(scatter_2), [outputPathScatter 'scatter'], structSizeSino3d);
        % Normalize the scatter:
        normScatter = scatter_2 .* overall_nf_3d;
        % Plot profiles to test:
        profileSinogram = sum(sinograms(:,126,:),3);
        profileRandoms = sum(randoms(:,126,:),3);
        profileNormScatter = sum(normScatter(:,126,:),3);
        subplot(1,3,2);
        title('Scatter Iteration 2');
        plot([profileSinogram profileRandoms profileNormScatter (profileRandoms+profileNormScatter)]);
        legend('Sinogram', 'Randoms', 'Scatter', 'Randoms+Scatter');
        
        outputPathScatter = [outputPath pathBar 'ScatterFinal' pathBar];
        if ~isdir(outputPathScatter)
            mkdir(outputPathScatter);
        end
        scatter = (scatter_1+scatter_2)./2;
        interfileWriteSino(single(scatter), [outputPathScatter 'scatter'], structSizeSino3d);
        % Normalize the scatter:
        normScatter = scatter .* overall_nf_3d;
        % Plot profiles to test:
        profileSinogram = sum(sinograms(:,126,:),3);
        profileRandoms = sum(randoms(:,126,:),3);
        profileNormScatter = sum(normScatter(:,126,:),3);
        subplot(1,3,3);
        title('Final Scatter');
        plot([profileSinogram profileRandoms profileNormScatter (profileRandoms+profileNormScatter)]);
        legend('Sinogram', 'Randoms', 'Scatter', 'Randoms+Scatter');
        
        % New dditive sinogram:
        additive = (randoms .* overall_ncf_3d  + scatter).* acfsSinogram; % (randoms +scatter.*norm)./(attenFactors*nprmFactrs)
        % write additive singoram:
        additiveFilename = [outputPathScatter 'additive'];
        interfileWriteSino(single(additive), additiveFilename, structSizeSino3d);
        
        if useGpu == 0
            disp('Mlem reconstruction...');
            saveIntermediate = 0;
            outputFilenamePrefix = [outputPathScatter 'reconImage'];
            filenameMlemConfig = [outputPathScatter 'mlem.par'];
            CreateOsemConfigFileForMmr(filenameMlemConfig, [sinogramFilename '.h33'], [filenameInitialEstimate '.h33'], outputFilenamePrefix, numIterations, [],...
                saveInterval, saveIntermediate, [anfFilename '.h33'], [additiveFilename '.h33']);
            % Execute APIRL: 
            status = system(['MLEM ' filenameMlemConfig]) 
        else
            disp('cuMlem reconstruction...');
            saveIntermediate = 0;
            outputFilenamePrefix = [outputPathScatter 'reconImage'];
            filenameMlemConfig = [outputPathScatter 'cuomlem.par'];
            CreateCuMlemConfigFileForMmr(filenameMlemConfig, [sinogramFilename '.h33'], [filenameInitialEstimate '.h33'], outputFilenamePrefix, numIterations, [],...
                saveInterval, saveIntermediate, [anfFilename '.h33'], [additiveFilename '.h33'], 0, 576, 576, 512);
            % Execute APIRL: 
            status = system(['cuMLEM ' filenameMlemConfig]) 
        end
        volume = interfileRead([outputFilenamePrefix '_final.h33']);
    end
end

if removeTempFiles
    if (numel(correctScatter) == 1) 
        if (correctScatter == 1)
            delete([outputPathScatter 'acfs*']);
            delete([outputPathScatter 'emission*']);
            delete([outputPathScatter 'scatter.s']);
            delete([outputPathScatter 'scatter.hs']);
            delete([outputPathScatter '*sensitivity*']);
            delete([outputPathScatter 'additive*']);
            delete([outputPathScatter 'tail*']);
        end
    end
    delete([outputPath 'sinogram*']);
    delete([outputPath 'ANF*']);
end