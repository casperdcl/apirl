function ha = reconShow(reconAPIRL, varargin)

z = 64; if nargin>1, z = varargin{1}; end

% fieldnames of images
names = fieldnames(reconAPIRL);
for i=length(names):-1:1
  dims = ndims(getfield(reconAPIRL, names{i}));
  if (dims == 3) || (dims == 4)
    i = i - 1;
  else
    names(i) = [];
  end
end

% subplot dimensions
dimX = ceil(sqrt(length(names)));
dimY = dimX;
while dimX * (dimY - 1) >= length(names),
  dimY = dimY - 1;
end

[ha, pos] = tight_subplot(dimY, dimX, [.02 .01],[.01 .02],[.01 .01]);
%axes(ha(1)), imshow(1-reconAPIRL.PET(110:end-120,120:end-120,64), []), title('(a) MLEM')

for i=1:length(names)
im = getfield(reconAPIRL, names{i});

%names{i}, ndims(im), size(im)
if ndims(im) == 3, im = im(:, :, z);
elseif ndims(im) == 4, im = im(1, :, :, z);
else continue; end

axes(ha(i)), imshow(squeeze(im), []), title(names{i});
end  % for

linkaxes;
end  % function reconShow
