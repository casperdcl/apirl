%  *********************************************************************
%  Reconstruction Framework for Siemens Biograph mMR.  
%  Autor: Martín Belzunce. Kings College London.
%  Fecha de Creación: 01/07/2015
%  *********************************************************************
%  This function resamples an image to get the pixel in a new coordinates.
%  Uses a spline interpolation. It receives as parameters the image, a 3-D
%  spatial referencing object (imref3d) that contains the spatial reference
%  of the input image, and another one of the the spatial reference for the
%  coordinate system wanted for the output image.
%  It return the resampledImage in the new coorindate system and the new
%  refimage object.
%
%  Example:
%   [resampledImage, refResampledImage] = ImageResample(image, refImage, newRefImage)
function [resampledImage, refResampledImage] = ImageResample(image, refImage, newRefImage, interpolation)

if nargin < 4
    interpolation = 'linear';
end

if numel(refImage.ImageSize) == 3
    % Coordiantes of the input image:
    xCoordIn_mm = refImage.XWorldLimits(1) + refImage.PixelExtentInWorldX/2: refImage.PixelExtentInWorldX : refImage.XWorldLimits(end);
    yCoordIn_mm = refImage.YWorldLimits(1) + refImage.PixelExtentInWorldY/2 : refImage.PixelExtentInWorldY : refImage.YWorldLimits(end);
    zCoordIn_mm = refImage.ZWorldLimits(1) + refImage.PixelExtentInWorldZ/2: refImage.PixelExtentInWorldZ : refImage.ZWorldLimits(end);
    % grid:
    [xIn_mm, yIn_mm, zIn_mm] = meshgrid(xCoordIn_mm, yCoordIn_mm, zCoordIn_mm);

    % New coordinates for the image:
    xCoordOut_mm = newRefImage.XWorldLimits(1) + newRefImage.PixelExtentInWorldX/2: newRefImage.PixelExtentInWorldX : newRefImage.XWorldLimits(end);
    yCoordOut_mm = newRefImage.YWorldLimits(1) + newRefImage.PixelExtentInWorldY/2 : newRefImage.PixelExtentInWorldY : newRefImage.YWorldLimits(end);
    zCoordOut_mm = newRefImage.ZWorldLimits(1) + newRefImage.PixelExtentInWorldZ/2: newRefImage.PixelExtentInWorldZ : newRefImage.ZWorldLimits(end);
    % grid:
    [xOut_mm, yOut_mm, zOut_mm] = meshgrid(xCoordOut_mm, yCoordOut_mm, zCoordOut_mm);

    % Interpolate the image to the new coordinate system:
    resampledImage = interp3(xIn_mm,yIn_mm,zIn_mm,double(image),xOut_mm,yOut_mm,zOut_mm, interpolation,0);
    resampledImage(isnan(resampledImage)) = 0;
    refResampledImage = newRefImage;

else
    % Coordiantes of the input image:
    xCoordIn_mm = refImage.XWorldLimits(1) + refImage.PixelExtentInWorldX/2: refImage.PixelExtentInWorldX : refImage.XWorldLimits(end);
    yCoordIn_mm = refImage.YWorldLimits(1) + refImage.PixelExtentInWorldY/2 : refImage.PixelExtentInWorldY : refImage.YWorldLimits(end);
    % grid:
    [xIn_mm, yIn_mm] = meshgrid(xCoordIn_mm, yCoordIn_mm);

    % New coordinates for the image:
    xCoordOut_mm = newRefImage.XWorldLimits(1) + newRefImage.PixelExtentInWorldX/2: newRefImage.PixelExtentInWorldX : newRefImage.XWorldLimits(end);
    yCoordOut_mm = newRefImage.YWorldLimits(1) + newRefImage.PixelExtentInWorldY/2 : newRefImage.PixelExtentInWorldY : newRefImage.YWorldLimits(end);
    % grid:
    [xOut_mm, yOut_mm] = meshgrid(xCoordOut_mm, yCoordOut_mm);

    % Interpolate the image to the new coordinate system:
    resampledImage = interp2(xIn_mm,yIn_mm,double(image),xOut_mm,yOut_mm, interpolation,0);
    resampledImage(isnan(resampledImage)) = 0;
    refResampledImage = newRefImage;
end

