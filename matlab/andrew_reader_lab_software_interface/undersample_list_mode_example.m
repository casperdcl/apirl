%% EXAMPLE MLEM MARTIN PROJECTOR (ANY SPAN)
clear all; close all
set_framework_environment();
%% DATA PATHS

rootDir = '/home/cc16/Downloads/AD_patients';
patient = '3'

dataPath = [rootDir '/P0' patient '/e7/data-Converted/data-LM-00/'];
outputPath = [rootDir '/P0' patient '/'];
nameDataSet = 'data-LM-00';
listModeFilename = [dataPath nameDataSet '.l'];
%sinogramFilename = [dataPath '/Scan2-00-sino-uncomp.s.hdr'];
attenuationMap_filename = [dataPath nameDataSet '-umap.v.hdr']; % This might be the one overwritten using mumap_registered_2 of reconstruct_without_scatter_and_register_ct
attenuationMapHardware_filename = [dataPath nameDataSet '-umap-hardware.v.hdr'];
%normFilename = [dataPath '../Scan2-norm.n'];
t1Path = [rootDir '/P0' patient '/T1/'];

%% INIT CLASS GPET
paramPET.scanner = 'mMR';
paramPET.method =  'otf_siddon_gpu';
paramPET.PSF.type = 'none';
paramPET.radialBinTrim = 0;
paramPET.Geom = '';
paramPET.method_for_normalization = 'from_e7_binary_interfile';
paramPET.method_for_randoms = 'from_e7_binary_interfile';%'from_ML_singles_matlab'; %
paramPET.method_for_scatter = 'from_e7_binary_interfile';
% To change span:
paramPET.sinogram_size.span = 11; % Any span, 0 for multislice 2d, -1 for 2d.
paramPET.nSubsets = 1;
% With PSF
paramPET.PSF.type = 'shift-invar';
paramPET.PSF.Width = 2.5;
PET = classGpet(paramPET);
paramPET.PSF.type = 'shift-invar';
paramPET.PSF.Width = 4.5;
PET_psf = classGpet(paramPET);
%% PET DATA CLASS TO UNDERSAMPLE
chunk_size_events = 1e7;
countReductionFractors = [0.901 0.201]; % It's not very intuitive how was implemented, this number is not the fraction of counts but the fraction to remove.
numRealizations = [3 1];
PETData = PETDataClass([dataPath]);
for i = 1 : numel(countReductionFractors)
    [totalPrompts,totalRandoms, totalWords, outFileHdr, output_listmode_file] = PETData.undersample_mMR_listmode_data(listModeFilename,countReductionFractors(i),chunk_size_events, numRealizations(i)); % The last factor is num realizations per call
    for j = 1 : numRealizations(i)
        pathThisRealization = [dataPath sprintf('count_reduction_factor_%.2f_r%d', countReductionFractors(i), j)];
        if ~isdir(pathThisRealization)
            mkdir(pathThisRealization);
        end
    end
end
%% NOW PREPARE THE DATA FOR CALLING E7 TOOLS
cd(dataPath);
for i = 1 : numel(countReductionFractors)
    for j = 1 : numRealizations(i)
        cd(dataPath);
        % Only if we don't need to rerun that:
        outFileHdr{j} = [sprintf('./data-LM-00-%.1f-r%d.hdr', (1-countReductionFractors(i))*100, j)];
        output_listmode_file{j} = [sprintf('./data-LM-00-%.1f-r%d.l', (1-countReductionFractors(i))*100, j)];
        pathThisRealization = [sprintf('count_reduction_factor_%.2f_r%d', countReductionFractors(i), j)];
        movefile(outFileHdr{j}, pathThisRealization);
        movefile(output_listmode_file{j}, pathThisRealization);
        % Now mumap:
        copyfile('data-LM-00-umap*', pathThisRealization);
        copyfile('../data-norm*', pathThisRealization);
        cd(pathThisRealization);
        try
            PETData = PETDataClass('./');
        catch
            disp('Keep running...');
        end
        PETData = PETDataClass('./');
        %PETData.Reconstruct(1, 100, 0, 0);
        %PETData.Reconstruct(1, 100, 1, 0);
        PETData.NCF();
    end
end
%% RECONSTRUCT
cd(dataPath);
% T1 image:
[T1InPetFov, refImagePetFov] = PET.getMrInPetImageSpace(t1Path);
patient_ad.T1 = T1InPetFov;
for i = 1 : numel(countReductionFractors)
    patient_ad.counts(i).undersampling_rate_perc = (1-countReductionFractors(i))*100;
    for j = 1 : numRealizations(i)
        pathThisRealization = [sprintf('count_reduction_factor_%.2f_r%d/', countReductionFractors(i), j)];
        sinogramFilename = [dataPath pathThisRealization nameDataSet sprintf('-%.1f-r%d-sino-0.s.hdr', (1-countReductionFractors(i))*100, j)];
        scatterBinaryFilename = [dataPath pathThisRealization '/rawdata_sino_00/scatter_estim2d_000000.s'];
        randomsBinaryFilename = [dataPath pathThisRealization '/rawdata_sino_00/smoothed_rand_00.s'];
        normBinaryFilename = [dataPath pathThisRealization '/rawdata_sino_00/norm3d_00.a'];
        %% READ AND PREPARE DATA
        % EMISSION SINOGRAM
        [sinogram, delayedSinogram, structSizeSino3d, info] = interfileReadSino(sinogramFilename);
        sino_span = PET.apply_axial_compression_from_span1(sinogram);
        % ATTENUATION MAP
        [attenuationMapHuman refMuMap] = interfileReadSiemensImage(attenuationMap_filename);
        [attenuationMapHardware refMuMap] = interfileReadSiemensImage(attenuationMapHardware_filename);
        attenuationMap = attenuationMapHuman + attenuationMapHardware;
        % MR IMAGE
        % for the mr image is important to have the correct bed position, that is
        % set from any interfile header of the scan, for example the attenuation
        % map:
        PET.setBedPosition(attenuationMap_filename);
        %% NORM
        ncfs = PET.NCF(normBinaryFilename);
        nf = ncfs;
        nf(nf~=0) = 1./nf(nf~=0);
        %% RANDOMS
        randoms = PET.R(randomsBinaryFilename);
        %% ACFs
        acfs = PET.ACF(attenuationMap, refMuMap);
        acfs_human = PET.ACF(attenuationMapHuman, refMuMap);
        %% SCATTER
        scatter = PET.S(scatterBinaryFilename, sino_span, ncfs, acfs_human, randoms); % Needs all that parameters to scale it.
        %% SENSITIVITY IMAGE
        anf = acfs .* ncfs;
        anf(anf~=0) = 1./anf(anf~=0);
        %% OP-MLEM
        reconPath = [outputPath 'OPMLEM_' sprintf('count_reduction_factor_%.2f_r%d/', countReductionFractors(i), j)];
        numIterations = 301;
        saveInterval = 100;
        % additive term:
        additive = randoms + scatter.*nf; % this is the correct way to do it, but AN needs to be included in the projector and backprojector.
        initial_image = PET.ones();
        sensImage = PET.Sensitivity(anf);
        opmlem = PET.OPMLEMsaveIter(sino_span, anf, additive, sensImage, initial_image, numIterations, reconPath, saveInterval);
        reconPath = [outputPath 'OPMLEM_psf_' sprintf('count_reduction_factor_%.2f_r%d/', countReductionFractors(i), j)];
        sensImage = PET_psf.Sensitivity(anf);
        opmlem_psf = PET_psf.OPMLEMsaveIter(sino_span, anf, additive, sensImage, initial_image, numIterations, reconPath, saveInterval);
        patient_ad.counts(i).realization(j).mlem_iter_100 = opmlem{1};
        patient_ad.counts(i).realization(j).mlem_iter_300 = opmlem{3};
        patient_ad.counts(i).realization(j).mlem_psf_iter_100 = opmlem_psf{1};
        patient_ad.counts(i).realization(j).mlem_psf_iter_300 = opmlem_psf{3};
        patient_ad.counts(i).realization(j).counts = sum(sino_span(:));
        save([outputPath 'casper_patient_ad' patient '.mat'], 'patient_ad', '-v7');
    end
end
