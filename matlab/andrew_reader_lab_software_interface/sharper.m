function im = sharper(im)

im = (im - mean(im(:))) / std(im(:));

if ndims(im) == 3
mfilt = @(im, sigma) medfilt3(im, sigma);
imsh = @(im, varargin) imsharpen3(im);
else
mfilt = @(im, sigma) medfilt2(im, sigma);
imsh = @(im, varargin) imsharpen(im, varargin{:});
end

%im = mfilt(im, [1, 1] * 2);
%im = (im .^ 1.3);
%im = imgaussfilt(im, 1);
%im = mfilt(im, [1, 1] .* 3);
%im = mfilt(im, [1, 1] .* 3);
%im = mfilt(im, [1, 1] .* 3);

im = imsh(im, 'Threshold',0.8);
%im = imsh(im, 'Threshold',0.8);
%im = imsh(im, 'Threshold',0.8);
im = (im - min(im(:))) .^ 1.3;

%im = mfilt(im, [1, 1] * 3);
%im = (im - min(im(:))) .^ 1.3;

%im = imsh(im, 'Radius',4);
%im = mfilt(im, [1, 1] * 3);
%im = imsh(im, 'Radius',4);
%im = (im - min(im(:))) .^ 1.2;
%im = imsh(im, 'Radius',4);
%im = (im - min(im(:))) .^ 1.2;

end  % segment


function im = imsharpen3(im, varargin)

blur = imgaussfilt3(im, [1 1 1] .* 0.6);  % radius
%blur = blur .* (abs(blur - im) >= abs(im .* 0.3));  % threshold
im = im + (im - blur) .* 8;  % amount

end  % imsharpen3
