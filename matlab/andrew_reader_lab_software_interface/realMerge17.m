%% NOTE: need to run on pat = AD_1, 2, ...
pre = 'output/reconAPIRL_real_nosharp'
pat = 'AD_2'
postC = '-C_'
post = '_t0.mat'
%%

% full (0.7) count `mlem(_psf)_fc`
d = load([pre '_' pat postC '3.01e+08' post], ['reconAPIRL']);
mlem_fc = squeeze(d.reconAPIRL.mlem(1, :, :, :));
mlem_psf_fc = squeeze(d.reconAPIRL.mlem_psf(1, :, :, :));

% append `fc` to low (0.1) count
d = load([pre '_' pat postC '4.3e+07' post], ['reconAPIRL']);
reconAPIRL = d.reconAPIRL;
reconAPIRL.PET_fc = mlem_fc;
reconAPIRL.PET_psf_fc = mlem_psf_fc;

% mirror
reconAPIRL.PET = reconAPIRL.PET(end:-1:1, :, :);
reconAPIRL.T1 = reconAPIRL.T1(end:-1:1, :, :);
reconAPIRL.mlem = reconAPIRL.mlem(:, end:-1:1, :, :);
reconAPIRL.mlem_psf = reconAPIRL.mlem_psf(:, end:-1:1, :, :);
reconAPIRL.PET_fc = reconAPIRL.PET_fc(end:-1:1, :, :);
reconAPIRL.PET_psf_fc = reconAPIRL.PET_psf_fc(end:-1:1, :, :);
%reconAPIRL.PET_psf = reconAPIRL.PET_psf(end:-1:1, :, :);

figure(1);
reconShow(reconAPIRL);
input('continue? ');

% NOTE: save to `pre-` not `pre_`
save([pre '-' pat postC '4.3e+07' post], '-v7.3', 'reconAPIRL');
