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
if numT
  pad = 1.5;  % place tumours in circle of diameter = 2 / (2+pad)
  for rt=rand([2 numT])
    y = (sin(rt(2) * 2 * pi) * rt(1) + 1) / (2+pad) + pad/(2*(2+pad));
    x = (cos(rt(2) * 2 * pi) * rt(1) + 1) / (2+pad) + pad/(2*(2+pad));
    im3d = max(im3d, tumour(d, h, w, m, ...
        [0.5, y, x, rSmall, maxScale, sigma]));
  end
else
  tBig = [0.5, 0.45, 0.4, rSmall*2.0, maxScale, sigma];
  tSmall = [0.5, 0.475, 0.55, rSmall, maxScale, sigma];
  tBigBlur = [0.5, 0.625, 0.55, tBig(4), maxScale, sigma + tBig(4)/3.0];
  tBiggest = [0.5, 0.625, 0.45, rSmall*3.0, maxScale, sigma];
  %tBig = [0.5, 0.45, 0.43, rSmall*2.0, maxScale, sigma];
  %tSmall = [0.5, 0.53, 0.57, rSmall, maxScale, sigma];
  %tBigBlur = [0.5, 0.625, 0.43, tBig(4), maxScale, sigma + tBig(4)/3.0];
  %tBiggest = [0.5, 0.5, 0.46, rSmall*3.0, maxScale, sigma];

  im3d = max(im3d, tumour(d, h, w, m, tBig));
  im3d = max(im3d, tumour(d, h, w, m, tSmall));
  im3d = max(im3d, tumour(d, h, w, m, tBigBlur));
  im3d = max(im3d, tumour(d, h, w, m, tBiggest));
end  % numT
im3d = im3d ./ max(im3d(:));
end  % addTumours


function im3d = tumour(d, h, w, m, t)

T = t .* [d h w w m w];
% [D, H, W, R, V, S] = T;
S = T(6);

[X, Y, Z] = meshgrid(1:w, 1:h, 1:d);
im3d = ((X-T(3)).^2 + (Y-T(2)).^2 + (Z-T(1)).^2 <= T(4)^2) .* T(5);
if S > 0
  im3d = imgaussfilt3(im3d, S);
end

im3d = permute(im3d, [3 2 1]);

end  % tumour
