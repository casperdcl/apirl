%%
PAT='1'; pat=load(['output/casper_patient_ad' PAT '.mat'], ['patient_ad' PAT]);
%%
pat=getfield(pat, ['patient_ad' PAT]);
if PAT == '1'
  tmp = pat.T1(:,:,end-6:end);
  pat.T1(:,:,8:end) = pat.T1(:,:,1:end-7);
  pat.T1(:,:,1:7) = tmp;
end
reconAPIRL.T1=single(pat.T1);
reconAPIRL.PET=single(pat.counts(2).realization(1).mlem_iter_300);
reconAPIRL.PET_psf=single(pat.counts(2).realization(1).mlem_psf_iter_300);
mlem=zeros([3 344 344 127], 'single');
psf=zeros([3 344 344 127], 'single');
for i=1:3
  mlem(i,:,:,:)=single(pat.counts(1).realization(i).mlem_iter_100);
  psf(i, :,:,:)=single(pat.counts(1).realization(i).mlem_psf_iter_300);
end
reconAPIRL.mlem=mlem;
reconAPIRL.mlem_psf=psf;
reconAPIRL.psf_mm=4.5;
reconAPIRL.voxelSize_mm=[2.08625 2.08625 2.03125];
%%
input('continue? ');
save(['output/reconAPIRL_real-AD_' PAT '.mat'], 'reconAPIRL', '-v7.3');
