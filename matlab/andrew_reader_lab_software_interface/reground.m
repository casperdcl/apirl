%%
tVersion = -1;
PAT = num2str(-tVersion);
%%
fname = ['output/reconAPIRL_real-AD_' PAT '.mat'];
uMap = permute(reshape(fread(fopen(['data-LM-00-umap-AD_' PAT '.v']), 'single'), [344 344 127]), [2 1 3]);
reconAPIRL = load(fname);
reconAPIRL = reconAPIRL.reconAPIRL
%%
sigmas = ones([1, 3]) * reconAPIRL.psf_mm ./ (reconAPIRL.voxelSize_mm * 2*sqrt(2*log(2)));
sieves = imgaussfilt3(reconAPIRL.PET_psf, sigmas);

im3d = sharper(sieves);

%{
tIntensity = max(max(max(medfilt3(im3d, [1 1 1] .* 3)))) / max(im3d(:))
wmm = size(im3d, 1)*reconAPIRL.voxelSize_mm(1);  % image width/[mm]
im3d = permute(addTumours(permute(im3d, [3 2 1]), wmm, 5, tIntensity, tVersion), [3 2 1]);
%}

%%
subplot(121), imshow(im3d(:,:,64), [])
subplot(122), imshow(uMap(:,:,64), [])

reconAPIRL.PET_psf_sharp = im3d;
reconAPIRL.uMap = uMap;
save(fname, '-v7.3', 'reconAPIRL');
%%
x = 110; X = 128;
y = 128; Y = 128;
z = 64;

im = im3d(x:end-X,y:end-Y,z);
figure, imshow(im, [])
