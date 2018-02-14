
    set_framework_environment();


%% Data Path
count_levels = 100;
patientNumber = 'P02';
ImgInfo.image_tag = 'AD-P2'; % needs to be changed for each patient
imgSaveFolder = ['\\bioeng202-pc\PET-M\AD_patients\' patientNumber '\e7\data-Converted\data-LM-00\sino_rawdata_' num2str(count_levels) '\temp\' ]; % needs to be specified

dataPath = ['\\bioeng202-pc\PET-M\AD_patients\' patientNumber '\e7\data-Converted\data-LM-00\'];
rawsinoName = ['\\bioeng202-pc\PET-M\AD_patients\' patientNumber '\e7\data-Converted\data-LM-00\sino_rawdata_\' num2str(count_levels)];
sino_hdr_names = ['\\bioeng202-pc\PET-M\AD_patients\' patientNumber '\e7\data-Converted\data-LM-00\data-LM-00-sino-' num2str(count_levels) '-0.s.hdr'];
norm_hdr_dir = ['\\bioeng202-pc\PET-M\AD_patients\' patientNumber '\e7\data-Converted\data-norm.n.hdr'];

%%
jsReconPath = '\\bioeng202-pc\\PET-M\e7_tools\C-JSRecon12-14-OCT-2015\IF2Dicom.js';
IF2DICOM_txt = 'Run-05-data-LM-00-IF2Dicom.txt';
%%
opt.tempPath ='';
opt.method =  'otf_siddon_gpu';
opt.nSubsets = 1;
opt.nIter = 100;
opt.PSF.Width = 0;
PET = classGpet(opt);

% some image info
ImgInfo.ReconMethod = 'OP-MLEM';
ImgInfo.matrixSize = PET.image_size.matrixSize;
ImgInfo.voxelSize_mm = PET.image_size.voxelSize_mm;
ImgInfo.nIter = PET.nIter;
ImgInfo.nSubsets = PET.nSubsets;

%% get Prompts, RS , AN
if 0 %generate them from e7 tool sinograms
    data = PETDataClass(rawsinoName,'get_sino_rawdata');
    Prompts = data.Prompts;
    RS = data.RS;
    AN = data.AN;
    save([rawsinoName 'rawdata_matlab.mat'], 'Prompts','RS','AN','-v7.3')
else % load from a matlab file
    load([rawsinoName 'rawdata_matlab.mat'], 'Prompts','RS','AN');
end


SensImg = PET.Sensitivity(AN);

mlem = PET.OPMLEM(Prompts,RS, SensImg,PET.ones,PET.nIter);


if ~exist(imgSaveFolder,'dir'), mkdir(imgSaveFolder); end

ImageInterFileName = [imgSaveFolder ImgInfo.image_tag '-' ImgInfo.ReconMethod '-' num2str(count_levels) '%.v'];


[Imgn,scaleFactor]= PET.BQML(mlem,sino_hdr_names,norm_hdr_dir);

ImgInfo.MinMaxImgCount = [min(Imgn(:)),max(Imgn(:))];
ImgInfo.scaleFactor = scaleFactor;
% write the interfile of the un-filtered image
wrireSiemensImageIntefileHdrs(ImageInterFileName,sino_hdr_names,norm_hdr_dir,ImgInfo)
writeBinaryFile(permute(Imgn,[2,1,3]),ImageInterFileName)

% post-filter with 3 mm kernel, write interfiles and convert into dicom
Img = PET.Gauss3DFilter(Imgn,3);
ImgInfo.MinMaxImgCount = [min(Img(:)),max(Img(:))];

ImageInterFileName = [imgSaveFolder ImgInfo.image_tag '-' ImgInfo.ReconMethod '-' num2str(count_levels) '%-3mm.v'];
wrireSiemensImageIntefileHdrs(ImageInterFileName,sino_hdr_names,norm_hdr_dir,ImgInfo)
writeBinaryFile(permute(Img,[2,1,3]),ImageInterFileName)

% populate DICOM files
[Path] = fileparts(ImageInterFileName);
copyfile([dataPath IF2DICOM_txt],Path);
IF2DICOM_txt_i = [Path PET.bar IF2DICOM_txt];
command = ['cscript ' jsReconPath ' "' [ImageInterFileName '.hdr'] '" "' IF2DICOM_txt_i '"'];
[status,~] = system(command);
if status, error('Dicom2IF failed\n'); end
delete(IF2DICOM_txt_i);

