pre = 'output/reconAPIRL_real_nosharp_AD_2-C_'
post = '_t0.mat'
%%

% full (0.7) count `mlem(_psf)_fc`
d = load([pre '3.01e+08' post], ['reconAPIRL']);
mlem_fc = squeeze(d.reconAPIRL.mlem(1, :, :, :));
mlem_psf_fc = squeeze(d.reconAPIRL.mlem_psf(1, :, :, :));

% append `fc` to low (0.1) count
d = load([pre '4.3e+07' post], ['reconAPIRL']);
reconAPIRL = d.reconAPIRL;
reconAPIRL.mlem_fc = mlem_fc;
reconAPIRL.mlem_psf_fc = mlem_psf_fc;
save([pre '4.3e+07' post], '-v7.3', 'reconAPIRL');
