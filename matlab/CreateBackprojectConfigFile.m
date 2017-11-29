%  *********************************************************************
%  Reconstruction Framework for Siemens Biograph mMR.  
%  Autor: Martín Belzunce. Kings College London.
%  Fecha de Creación: 15/06/2015
%  *********************************************************************
% Function that creates a config file to project an image with apirl.

% Backproject Parameters :=
% ; Archivo de configuración de reconstrucción MLEM para paper en Parallel Computing.
% input file := ../../inputData/Project_Nema_Siddon.h33
% input type := Sinogram3D
% ; Projectors:
% backprojector := Siddon
% output image := ../../../../images3D/constantImage_FullFov_256_256_81.h33
% output filename := Project_PointSources_
% END :=

% To backproject only a subset, add the following parameters:
% number of subsets := 21
% subset index := 5
%
% scanner types: 'cylindrical', 'mMR'
% scanner_properties: for cylindrical scanner_properties.radius_mm,
% scanner_properties.axialFov_mm, scanner_properties.radiusFov_mm

function CreateBackprojectConfigFile(configfilename, inputFile, outputSample, scanner, scanner_properties, outputFilename, numberOfSubsets, subsetIndex, useGpu, numSamples, numAxialSamples)
if nargin == 6
    useGpu = 0;
    numSamples = 1;
    numAxialSamples = 1;
elseif nargin == 7
    numSamples = 1;
    numAxialSamples = 1;
elseif nargin == 8
    numAxialSamples = 1;
end

% Handle the number of subsets:
if isempty(numberOfSubsets)
    numberOfSubsets = 0;
end
if(isempty(subsetIndex))
    subsetIndex = 0;
end

% Get info of the sinogram size to choose the type of sinogram:
[info, structSizeSino] = getInfoFromInterfile(inputFile);

if strcmp(scanner, 'mMR')
    if isfield(structSizeSino, 'sinogramsPerSegment')
        if numel(structSizeSino.sinogramsPerSegment) == 1      % If its span a field a 3d sinogram with axial compression of 121 (1 segment)
            if ~isfield(structSizeSino, 'span')
                if structSizeSino.numZ == 1
                    sinogramType = 'Sinogram2DinSiemensMmr';
                else
                    sinogramType = 'Sinograms2DinSiemensMmr';
                end
            elseif structSizeSino.span >= 1
                sinogramType = 'Sinogram3DSiemensMmr';
            else
                if structSizeSino.numZ == 1
                    sinogramType = 'Sinogram2DinSiemensMmr';
                else
                    sinogramType = 'Sinograms2DinSiemensMmr';
                end
            end
        else
            sinogramType = 'Sinogram3DSiemensMmr';
        end
    else
        if structSizeSino.numZ == 1
            sinogramType = 'Sinogram2DinSiemensMmr';
        else
            sinogramType = 'Sinograms2DinSiemensMmr';
        end
    end
elseif strcmp(scanner, 'cylindrical')
    if isfield(structSizeSino, 'sinogramsPerSegment')
        if numel(structSizeSino.sinogramsPerSegment) == 1
            if structSizeSino.numZ == 1
                sinogramType = 'Sinogram2D';
            else
                sinogramType = 'Sinograms2D';
            end
        else
            sinogramType = 'Sinogram3D';
        end
    else
        if structSizeSino.numZ == 1
            sinogramType = 'Sinogram2D';
        else
            sinogramType = 'Sinograms2D';
        end
    end
end

% Primero genero el archivo de encabezado.
fid = fopen(configfilename, 'w');
if(fid == -1)
    fprintf('No se pudo crear el archivo %s.', configfilename);
end
% Ahora debo ir escribiendo los campos. Algunos son fijos, y otros
% dependerán de la imagen:
fprintf(fid,'Backproject Parameters :=\n');
fprintf(fid,'input type := %s\n', sinogramType);
if strcmp(scanner, 'cylindrical')
    fprintf(fid,'cylindrical pet radius (in mm) := %f\n', scanner_properties.radius_mm);
    fprintf(fid,'axial fov (in mm) := %f\n', scanner_properties.axialFov_mm);
    fprintf(fid,'radius fov (in mm) := %f\n', scanner_properties.radialFov_mm);
end
if useGpu == 0
    fprintf(fid,'backprojector := Siddon\n');
else
    fprintf(fid,'backprojector := CuSiddonProjector\n');
    fprintf(fid,'backprojector block size := {576,1,1}\n');
    %fprintf(fid,'gpu id := 0\n');
    g = gpuDevice;
    fprintf(fid,['gpu id := ' num2str(g.Index - 1) '\n']);
end
fprintf(fid,'siddon number of samples on the detector := %d\n', numSamples);
fprintf(fid,'siddon number of axial samples on the detector := %d\n', numAxialSamples);
% If we have subsets put it in the config file:
if numberOfSubsets ~= 0
    fprintf(fid,'number of subsets := %d\n', numberOfSubsets);
    fprintf(fid,'subset index := %d\n', subsetIndex-1); % Subset index is 0 base-index in apirl, while 1 in matlab.
end
fprintf(fid,'input file := %s\n', inputFile);
fprintf(fid,'output image := %s\n', outputSample);
fprintf(fid,'output filename := %s\n', outputFilename);

fclose(fid);
