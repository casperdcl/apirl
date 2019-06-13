function mmaps = readPetMr(sigma, petNoise, t1Noise, varargin)
% @param sigma  : float
% @param petNoise  : float
% @param t1Noise  : float
% @return mmaps  : struct [PET, uMap, T1]

prefix = 'subject_04';
if nargin > 3, prefix=varargin{1}; end
suffix = '_mMR.mat';

data = load([prefix sprintf('_sigma%.3g_noise%.3g', sigma, petNoise) suffix]);
mmaps.PET = data.MultiMaps_Ref.PET;

data = load([prefix sprintf('_sigma%.3g_noise%.3g', sigma, t1Noise) suffix]);
mmaps.uMap = data.MultiMaps_Ref.uMap;
mmaps.T1 = data.MultiMaps_Ref.T1;

end  % readPetMr
