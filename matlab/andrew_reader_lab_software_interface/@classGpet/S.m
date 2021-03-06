% *********************************************************************
% Reconstruction Framework for Siemens Biograph mMR.  
% class: Gpet
% Authors: Martin Belzunce, Abolfazl Mehranian. Kings College London.
% Date: 08/02/2016
% *********************************************************************
% Computes the scatter using one of the availables methods.
function s=S(varargin)
    objGpet = varargin{1};
    h = fspecial('gaussian',30 ,15);
    s = [];
    if ~strcmpi(objGpet.scanner,'mMR')&& ~strcmpi(objGpet.scanner,'cylindrical')
        error('NCFs are only available for mMR and cylindrical scanners.');
    end
    if nargin == 2 % Simple simulation, smooth a sinogram.
        % filename or a sinogram to smooth?
        if strcmpi(objGpet.method_for_scatter,'from_e7_binary_interfile') && ischar(varargin{2})
            filename = varargin{2};
            fid = fopen(filename, 'r');
            if fid == -1
                error('The scatter binary file was not found.');
            end
            s = fread(fid,[], 'single');
            if numel(s) ~= [344 252 127]
                error('Unexpected file size for %s', filename);
            end
            s = reshape(s, [344 252 127]);
            % Convert into the sinogram size:
            s = objGpet.iSSRB(s);
        else
            s = zeros(size(varargin{2}));
            % 1 D filter for each projection:
            for i = 1 : size(s,3)
                s(:,:,i) = imfilter(varargin{2}(:,:,i), h, 'same', 'circular');
            end
            s(isnan(s)) = 0;
        end
        % Fit tail?
    elseif nargin == 3 % Sss simulation, need activty image and attenuation.
        if strcmpi(objGpet.method_for_scatter,'e7_tools')
            % Call e7 tools:
        elseif strcmpi(objGpet.method_for_scatter,'from_e7_binary_interfile')
            error('Invalid number of parameters');
        elseif strcmpi(objGpet.scatter_algorithm,'sss_stir')
            % Call stir:
        end
    elseif nargin == 6 % S(objGpet,scatter_3D, emission_sinogram, ncf, acf, randoms)
        if strcmpi(objGpet.method_for_scatter,'from_e7_binary_interfile') && ischar(varargin{2})
            filename = varargin{2};
            fid = fopen(filename, 'r');
            if fid == -1
                error('The scatter binary file was not found.');
            end
            s = fread(fid, inf, 'single');
            if numel(s) ~= [344*252*127]
                error('Unexpected file size for %s', filename);
            end
            s = reshape(s, [344 252 127]);
            if objGpet.sinogram_size.span >= 1
                % Convert into the sinogram size:
                s = objGpet.iSSRB(s);
                % Scale it:
                %s = objGpet.scatter_scaling(s, varargin{3}, varargin{4}, varargin{5}, varargin{6});
            elseif objGpet.sinogram_size.span == 0
                rings = 1 : objGpet.sinogram_size.nRings;
                planes = 1 : objGpet.sinogram_size.nRings/(127+1) : objGpet.sinogram_size.nRings;
                [X1,Y1,Z1] = meshgrid(1:252,1:344,planes);
                [X2,Y2,Z2] = meshgrid(1:252,1:344,rings);
                s = interp3(X1,Y1,Z1,s,X2,Y2,Z2);
                s = objGpet.scatter_scaling(s, varargin{3}, varargin{4}, varargin{5}, varargin{6});
            end
        end
    elseif nargin == 7 % S(objGpet,scatter_3D, emission_sinogram, ncf, acf, randoms, ring)
        if strcmpi(objGpet.method_for_scatter,'from_e7_binary_interfile') && ischar(varargin{2})
            filename = varargin{2};
            fid = fopen(filename, 'r');
            if fid == -1
                error('The scatter binary file was not found.');
            end
            s = fread(fid, inf, 'single');
            if numel(s) ~= [344*252*127]
                error('Unexpected file size for %s', filename);
            end
            s = reshape(s, [344 252 127]);
            if objGpet.sinogram_size.span == -1
                plane = varargin{7}*(127+1)/ 64;
                s = s(:,:,plane);
                s = objGpet.scatter_scaling(s, varargin{3}, varargin{4}, varargin{5}, varargin{6});
            end
        end
    end
 end