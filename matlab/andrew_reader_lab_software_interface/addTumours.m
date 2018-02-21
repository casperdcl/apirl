function im3d = addTumours(im3d, wmm, varargin)
% im3d  : 3-dimesional matrix.
% wmm  : `im3d` width in mm.
%% VARARGS
% diam  : minimum tumour diameter [default: 6mm].
% maxScale  : max(im3d(:)) * `maxScale` is tumour intensity [default: 1.5].
% num  : number of tumours (implies const diameter specified) [default: 0].
% sigma  : minimum blur for tumours [default: 0mm].

%% defaults
dmmSmall = 5;  % ~5mm tumour diameter
maxScale = 1.5;  % max(im3d(:)) * `maxScale` is tumour intensity
numT = 0;
blur = 0;  % 0mm
%% varargs
if nargin>2, dmmSmall=varargin{1}; end
if nargin>3, maxScale=varargin{2}; end
if nargin>4, numT=varargin{3}; end
if nargin>5, blur=varargin{4}; end

%% main
rSmall = dmmSmall / (2 * wmm);
sigma = blur / wmm;

% coords: depth, height, width
[d, h, w] = size(im3d);
m = max(im3d(:));
if ~m, m = 1; end

% NB: tumours specified by:
% fractional positions d/h/w, radius/[w], value/[max], sigma/[w]
pad = 1.5;  % place tumours in circle of diameter = 2 / (2+pad)
if numT > 0
  z = 0.5;
  for rt=rand([2 numT])
    x, y = pol2cart(rt(1) *2*pi, rt(2));
    % centre about 0.5
    y = (y + 1) / 2;
    x = (x + 1) / 2;
    im3d = max(im3d, tumour(d, h, w, m, ...
        [z, y, x, rSmall, maxScale, sigma]), 1/pad);
  end
else
  if numT == -2
  tBig = [0.5, 0.45, 0.4, rSmall*2.0, maxScale, sigma];
  tSmall = [0.5, 0.475, 0.55, rSmall, maxScale, sigma];
  tBigBlur = [0.5, 0.625, 0.55, tBig(4), maxScale, sigma + tBig(4)/3.0];
  tBiggest = [0.5, 0.625, 0.45, rSmall*3.0, maxScale, sigma];
  else
  tBig = [0.5, 0.45, 0.43, rSmall*2.0, maxScale, sigma];
  tSmall = [0.5, 0.53, 0.57, rSmall, maxScale, sigma];
  tBigBlur = [0.5, 0.625, 0.43, tBig(4), maxScale, sigma + tBig(4)/3.0];
  tBiggest = [0.5, 0.5, 0.46, rSmall*3.0, maxScale, sigma];
  end  % numT

  im3d = max(im3d, tumour(d, h, w, m, tBig, 1/pad));
  im3d = max(im3d, tumour(d, h, w, m, tSmall, 1/pad));
  im3d = max(im3d, tumour(d, h, w, m, tBigBlur, 1/pad));
  im3d = max(im3d, tumour(d, h, w, m, tBiggest, 1/pad));
end  % numT
im3d = im3d ./ max(im3d(:));
end  % addTumours


function im3d = tumour(d, h, w, m, t, varargin)
%% Arguments:
% image depth/height/width, maxscale, tumour params
%% Options:
% scale (default 1) towards centre for tumour z/y/x

scale = 1; if nargin>5, scale=varargin{1}; end
[phi, theta, r] = cart2sph(t(3) - 0.5, t(2) - 0.5, t(1) - 0.5);
[x, y, z] = sph2cart(phi, theta, r * scale);
z = z + 0.5;
y = y + 0.5;
x = x + 0.5;

T = [z, y, x, t(4:end)] .* [d h w w m w];

% [D, H, W, R, V, S] = T;
S = T(6);

[X, Y, Z] = meshgrid(1:w, 1:h, 1:d);
im3d = ((X-T(3)).^2 + (Y-T(2)).^2 + (Z-T(1)).^2 <= T(4)^2) .* T(5);
if S > 0
  im3d = imgaussfilt3(im3d, S);
end

im3d = permute(im3d, [3 2 1]);

end  % tumour

