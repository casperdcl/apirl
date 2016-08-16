function [BW,maskedImage] = segmentSmallColdSpot(X)
%segmentImage Segment image using auto-generated code from imageSegmenter App
%  [BW,MASKEDIMAGE] = segmentImage(X) segments image X using auto-generated
%  code from the imageSegmenter App. The final segmentation is returned in
%  BW, and a masked image is returned in MASKEDIMAGE.

% Auto-generated by imageSegmenter app on 16-Aug-2016
%----------------------------------------------------


% Normalize input data to range in [0,1].
Xmin = min(X(:));
Xmax = max(X(:));
X = (X - Xmin) ./ (Xmax - Xmin);

% Create empty mask.
BW = false(size(X));

% Flood fill
row = 11;
column = 23;
tolerance = 5.000000e-02;
addedRegion = grayconnected(X, row, column, tolerance);
BW = BW | addedRegion;

% Create masked image.
maskedImage = X;
maskedImage(~BW) = 0;