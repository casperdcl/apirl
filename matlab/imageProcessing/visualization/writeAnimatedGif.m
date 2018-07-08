%  *********************************************************************
%  Visualization tools.  
%  Author: Martín Belzunce. Kings College London.
%  Date: 10/02/2016
%  *********************************************************************

% Function that writes an animated gif from a 3d volume:
function filename = writeAnimatedGif(volume, filename, frame_time, colormap,adjustCotnrast)

if nargin < 3
    frame_time = 0.05;
end
if nargin < 4
    numColors = 255;
    colormap = repmat([0 : 1/numColors : 1]',1,3);
else
    numColors = size(colormap,1)-1;
end
if nargin < 5
    adjustCotnrast = 1; % Use the whole range
end
if ndims(volume) == 4
    % indexed volume
    maxImage = max(max(max(volume))).*adjustCotnrast;
    for n = 1 : size(volume,4)
        %imind = uint8(round(double(volume(:,:,slices(n))./maxImage.*numColors)));
        %imind = double(volume(:,:,slices(n))./maxImage);
        % Write to the GIF File 
        imwrite(double(volume(:,:,:,n)),filename,'gif','WriteMode','append','DelayTime', frame_time); 
    end
else
    slices = find(mean(mean(volume))>0);
    maxImage = max(max(max(volume))).*adjustCotnrast;
    for n = 1 : numel(slices)
        imind = uint8(round(double(volume(:,:,slices(n))./maxImage.*numColors)));
        %imind = double(volume(:,:,slices(n))./maxImage);
        % Write to the GIF File 
        if n == 1 
          imwrite(imind,colormap,filename,'gif', 'Loopcount',inf); 
        else 
          imwrite(imind,colormap,filename,'gif','WriteMode','append','DelayTime', frame_time); 
        end 
    end
end