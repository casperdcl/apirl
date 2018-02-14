%% Generate sinograms based on the example_real_data_mmr_2.m
run('example_real_data_mmr_2.m')

%% Build Prior's object
opt.ImageSize= PET.image_size.matrixSize;
opt.sWindowSize = 5; % search window size
opt.lWindowSize = 1; % patch-size
opt.imCropFactor = 2;
G = PriorsClass(opt);

filteredMR = PET.Gauss3DFilter(MrInPet,2.5);
%% Bowsher
opt.prior = 'tikhonov';
opt.weight_method = 'bowsher';
opt.B = 50;
opt.nl_weights = G.W_Bowsher(filteredMR,opt.B);
opt.tolerance = 1e-4;
opt.save = 0;
opt.nIter = 150;
opt.display = 1;
opt.lambda = 5e4;
opt.message = 'Bowsher-B50';
opt.save_i = ['D:\Multi_parametricPET_MR\Results\ClinicalData\FDG_P2_TEST\' opt.message '\'];

BW = G.MAP_OSEM(PET,sino_span,additive, sensImage, opt,initial_image);
%% Local joint entropy
opt.prior = 'approx_joint_entropy';
opt.sigma_a = 15;
opt.sigma_f = 0.2;
opt.je_weights = G.W_JointEntropy(filteredMR,opt.sigma_a);
opt.tolerance = 1e-4;
opt.save = 0;
opt.nIter = 150;
opt.display = 1;
opt.lambda = 5e4;
opt.message = 'JointEntropy-a15-f0.2';
opt.save_i = ['D:\Multi_parametricPET_MR\Results\ClinicalData\FDG_P2_TEST\' opt.message '\'];
JE = G.MAP_OSEM(PET,sino_span,additive, sensImage, opt,initial_image);


