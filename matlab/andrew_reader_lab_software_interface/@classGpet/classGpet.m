% *********************************************************************
% Reconstruction Framework for Siemens Biograph mMR.
% class: Gpet
% Authors: Martin Belzunce, Abolfazl Mehranian. Kings College London.
% Date: 08/02/2016
% *********************************************************************
% Definition of the class classGpet to manage projector/backprojectors and
% other tools used in the image reconstruction of Siemens Biograph mMR
% data.

classdef classGpet < handle
    properties (SetAccess = private)
        % Type of scanner. Options: '2D_radon', 'mMR', 'cylindrical'
        scanner
        % Properties of the scanner. This is an structure with additional
        % parameters specific for each scanner, for example for cylindrical
        % scanner it has the radius_mm propery.
        scanner_properties
        % Sinogram size:
        sinogram_size     % Struct with the size of the sinogram, it can be 2d, multislice2d and 3d.
        % Image size:
        image_size
        % Bed position (introduces an offset in the axial axis):
        bed_position_mm
        % Ref structure:
        ref_image
        % Ref structure for the native image space for the scanner:
        ref_native_image
        % Projector/Backrpojector. Options:
        % 'pre-computed_matlab','otf_matlab', 'otf_siddon_cpu','otf_siddon_gpu'
        method
        % PSF. Options: 'shift-invar', 'shift-var', 'none'
        PSF
        % Number of iterations:
        nIter
        % Number of subsets:
        nSubsets
        % Radial Bin trimming
        radialBinTrim
        % Temporary files path:
        tempPath
        % Flag to delete the temp folder when we're finished (SAM ELLIS,
        % 24/08/2016)
        deleteTemp
        % Asymmertic pre-computed sysetm matrix
        Geom
        % Method to estimate randoms:
        method_for_randoms
        % Method to estimate scatter:
        method_for_scatter
        % Method for normalization:
        method_for_normalization    % 'cbn_expansion', 'from_e7_binary_interfile'
        % Operative system.
        os
        % bar for the paths.
        bar
        % number of Siddon rays per LOR (Transverse direction)
        nRays
        % number of Siddon rays per LOR in the axial direction
        % (axialDirection)
        nAxialRays
        % Verbosity level: 0, 1, 2. (0 silent)
        verbosity
        % Prioir for map reconstruction.
        Prior
    end

    methods
        % Constructors:
        function objGpet = classGpet(varargin) % Default options: (). From file (filename)
            if(strcmp(computer(), 'GLNXA64')) % If linux, call wine64.
                objGpet.os = 'linux';
                objGpet.bar = '/';
            else
                objGpet.os = 'windows';
                objGpet.bar = '\';
            end

            if (nargin == 0)
                objGpet.scanner = 'mMR';
            else
                if isfield(varargin{1},'scanner')
                    objGpet.scanner = varargin{1}.scanner;
                elseif length(varargin)>=2 && ischar(varargin{1}) % pair name,values
                    i = find(cellfun(@any,strfind(varargin,'scanner')));
                    objGpet.scanner = varargin{i+1};
                else
                    objGpet.scanner = 'mMR';
                end
            end

            objGpet.initScanner();
            objGpet.bed_position_mm = 0;
            % Init the native image space:
            objGpet.init_ref_image();
            objGpet.ref_native_image = objGpet.ref_image;
            objGpet.nRays = 1;
            objGpet.nAxialRays = 1;
            objGpet.method =  'otf_siddon_cpu';
            objGpet.PSF.type = 'shift-invar';
            objGpet.PSF.Width = 4; %mm
            objGpet.radialBinTrim = 0;
            objGpet.Geom = '';
            % SAM ELLIS EDIT (25/08/2016): changed the '.' to pwd so that
            % when deleting temp files, don't need to be in the original
            % directory
            objGpet.tempPath = [pwd objGpet.bar 'temp' objGpet.bar];
            objGpet.deleteTemp = false; % SAM ELLIS 24/08/2016
            objGpet.verbosity = 0;
            objGpet.method_for_randoms = 'from_ML_singles_matlab';
            objGpet.method_for_normalization = 'cbn_expansion';
            user_opts = []; % empty options if there weren't any field related with priors.
            if nargin == 1
                % Read configuration from file or from struct:
                if isstruct(varargin{1})
                    if ~isfield(varargin{1}, 'scanner') || ~isfield(varargin{1}, 'method')
                        disp('Configuration for''otf_siddon_cpu''');
                    end
                    Revise(objGpet,varargin{1});
                    % Options for the prior:
                    if isfield(varargin{1}, 'Prior')
                        user_opts = varargin{1}.Prior;
                    end
                elseif ischar(varargin{1}) % Read configuration from file:
                    objGpet.readConfigFromFile(varargin{1})
                end
            else
                % Read configuration from name/value pairs:
                objGpet = varargin_pair(objGpet, varargin);
            end
            if ~isdir(objGpet.tempPath)
                mkdir(objGpet.tempPath)
            end

            if (objGpet.sinogram_size.span ~=11) || (objGpet.sinogram_size.maxRingDifference ~=60 ) % default values
                init_sinogram_size(objGpet, objGpet.sinogram_size.span, objGpet.sinogram_size.nRings, objGpet.sinogram_size.maxRingDifference);
            end
            if objGpet.sinogram_size.span == 0
                objGpet.image_size.voxelSize_mm(3) = objGpet.image_size.voxelSize_mm(3)*objGpet.image_size.matrixSize(3)/objGpet.sinogram_size.nRings;
                objGpet.image_size.matrixSize(3)=objGpet.sinogram_size.nRings;
                warning('Overriding the image size to the number of rings to work in multislice 2d. Now the image size is %dx%dx%d and voxel size is %fx%fx%f.', ...
                    objGpet.image_size.matrixSize(1), objGpet.image_size.matrixSize(2), objGpet.image_size.matrixSize(3),...
                    objGpet.image_size.voxelSize_mm(1), objGpet.image_size.voxelSize_mm(2), objGpet.image_size.voxelSize_mm(3));
            elseif objGpet.sinogram_size.span == -1
                % Keep the voxel sizeobjGpet.image_size.voxelSize_mm(3)
                objGpet.image_size.matrixSize(3)=objGpet.sinogram_size.nRings;

                % SAM ELLIS EDIT: removed this warning, since span = -1
                % gives nRings = 1 as desired, so no warning necessary
                %                 warning('Overriding the image size to the number of rings to work in multislice 2d. Now the image size is %dx%dx%d and voxel size is %fx%fx%f.', ...
                %                     objGpet.image_size.matrixSize(1), objGpet.image_size.matrixSize(2), objGpet.image_size.matrixSize(3),...
                %                     objGpet.image_size.voxelSize_mm(1), objGpet.image_size.voxelSize_mm(2), objGpet.image_size.voxelSize_mm(3));
            end
            objGpet.init_ref_image();

            if ~isempty(user_opts) % User opts should have been loaded a few lines before when we check the input parameters.
                p.ImageSize = objGpet.image_size.matrixSize;
                p.imCropFactor = [3,3,0];
                p.sWindowSize = 5;
                p.lWindowSize = 1;
                if isstruct(varargin{1}) % have to generalize it
                    p = getFiledsFromUsersOpt(p,user_opts);
                end
                objGpet.Prior = PriorsClass(p);
            end
        end

        function objGpet = readConfigFromFile(objGpet, strFilename)
            error('todo: read configuration from file');
        end

        % Function that intializes the scanner:
        function initScanner(objGpet)
            if strcmpi(objGpet.scanner,'2D_radon')
                disp('Configuration for ''2D_radon'' scanner');
                objGpet.G_2D_radon_setup();
            elseif strcmpi(objGpet.scanner,'mMR')
                disp('Configuration for ''mMR'' scanner');
                objGpet.G_mMR_setup();
            elseif strcmpi(objGpet.scanner,'cylindrical')
                disp('Configuration for ''cylindrical'' scanner');
                objGpet.G_cylindrical_setup();
            else
                error('unkown scanner')
            end
        end



        function G_2D_radon_setup(objGpet)
            % Default parameter, only if it havent been loaded by the
            % config previously:
            if isempty(objGpet.image_size)
                objGpet.image_size.matrixSize = [512, 512, 1];
                objGpet.image_size.voxelSize_mm = [1, 1, 1];
            end
            if isempty(objGpet.sinogram_size)
                x = radon(ones(objGpet.image_size.matrixSize(1:2)),0:179);
                x = size(x);
                objGpet.sinogram_size.nRadialBins = x(1);
                objGpet.sinogram_size.nAnglesBins = x(2);
                objGpet.sinogram_size.nSinogramPlanes = 1;
                % EDIT: SAM ELLIS - 23/08/2016
                % NEED TO ADD OTHER DEFAULTS FOR 2D RADON
                objGpet.sinogram_size.span = -1;
                objGpet.sinogram_size.nRings = 1;
                objGpet.sinogram_size.maxRingDifference = 0;
            end
            if isempty(objGpet.nSubsets)
                objGpet.nSubsets = 1;
            end
            if isempty(objGpet.nIter)
                objGpet.nIter = 40;
            end
            objGpet.sinogram_size.matrixSize = [objGpet.sinogram_size.nRadialBins objGpet.sinogram_size.nAnglesBins objGpet.sinogram_size.nSinogramPlanes];
            objGpet.osem_subsets(objGpet.nSubsets, objGpet.sinogram_size.nAnglesBins);
        end

        function G_mMR_setup(objGpet)
            % Default parameter, only if it havent been loaded by the
            % config previously:
            if isempty(objGpet.sinogram_size)
                objGpet.sinogram_size.nRadialBins = 344;
                objGpet.sinogram_size.nAnglesBins = 252;
                objGpet.sinogram_size.nRings = 64;
                objGpet.sinogram_size.nSinogramPlanes = 837;
                objGpet.sinogram_size.maxRingDifference = 60;
                objGpet.sinogram_size.nPlanesPerSeg = [127   115   115    93    93    71    71    49    49    27    27];
                objGpet.sinogram_size.span = 11;
                objGpet.sinogram_size.nSeg = 11;
                objGpet.scanner_properties.radius_mm = 328;
                objGpet.scanner_properties.sinogramDepthOfInteraction_mm = 6.7;
                objGpet.scanner_properties.LORDepthOfInteraction_mm = 9.6;
                objGpet.scanner_properties.planeSep_mm = 2.03125;
                objGpet.scanner_properties.nCrystalsPerRing = 504;
                objGpet.scanner_properties.binSize_mm = pi*objGpet.scanner_properties.radius_mm/objGpet.scanner_properties.nCrystalsPerRing;%2.0445;

            else

            end
            if isempty(objGpet.nSubsets)
                objGpet.nSubsets = 21;
            end
            if isempty(objGpet.nIter)
                objGpet.nIter = 3;
            end

            objGpet.sinogram_size.matrixSize = [objGpet.sinogram_size.nRadialBins objGpet.sinogram_size.nAnglesBins objGpet.sinogram_size.nSinogramPlanes];
            objGpet.image_size.matrixSize =[344, 344, 127];
            objGpet.image_size.voxelSize_mm = [2.08626 2.08626 2.03125];

            objGpet.osem_subsets(objGpet.nSubsets, objGpet.sinogram_size.nAnglesBins);
        end

        function G_cylindrical_setup(objGpet)
            % Default parameter, only if it havent been loaded by the
            % config previously. By default use the same size of mmr, additionally needs the scanner radius:

            if isempty(objGpet.sinogram_size)
                objGpet.scanner_properties.radius_mm = 328;
                objGpet.sinogram_size.nRadialBins = 344;
                objGpet.sinogram_size.nAnglesBins = 252;
                objGpet.sinogram_size.nRings = 64;
                objGpet.sinogram_size.nSinogramPlanes = 837;
                objGpet.sinogram_size.maxRingDifference = 60;
                objGpet.sinogram_size.nPlanesPerSeg = [127   115   115    93    93    71    71    49    49    27    27];
                objGpet.sinogram_size.span = 11;
                objGpet.sinogram_size.nSeg = 11;
                %add maxDiff
                %                 objGpet.radialBinTrim = 0;

            end

            if ~isfield(objGpet.scanner_properties, 'radius_mm')
                objGpet.scanner_properties.radius_mm = 370;
            end

            if isempty(objGpet.nSubsets)
                objGpet.nSubsets = 21;
            end
            if isempty(objGpet.nIter)
                objGpet.nIter = 3;
            end

            objGpet.sinogram_size.matrixSize = [objGpet.sinogram_size.nRadialBins objGpet.sinogram_size.nAnglesBins objGpet.sinogram_size.nSinogramPlanes];
            objGpet.image_size.matrixSize =[344, 344, 127];
            objGpet.image_size.voxelSize_mm = [2.08626 2.08626 2.03125];

            objGpet.osem_subsets(objGpet.nSubsets, objGpet.sinogram_size.nAnglesBins);

            % Init the scanner parameters depending on the desired image.
            % The FOV used covers the whole image:
            objGpet.scanner_properties.radialFov_mm = max(objGpet.image_size.matrixSize(2)*objGpet.image_size.voxelSize_mm(2),objGpet.image_size.matrixSize(1)*objGpet.image_size.voxelSize_mm(1))./2;
            objGpet.scanner_properties.axialFov_mm = objGpet.image_size.matrixSize(3)*objGpet.image_size.voxelSize_mm(3);
        end

        function set_subsets(objGpet, numSubsets)
            % Update number of iteration to keep the subset%iterations
            % constant:
            objGpet.nIter = objGpet.nIter*objGpet.nSubsets/numSubsets;
            objGpet.nSubsets = numSubsets;

            objGpet.osem_subsets(objGpet.nSubsets, objGpet.sinogram_size.nAnglesBins);
        end

        function objGpet=init_image_properties(objGpet, refImage)
            objGpet.image_size.matrixSize = refImage.ImageSize;

            if refImage.ImageSize(3) >1
                objGpet.image_size.voxelSize_mm = [refImage.PixelExtentInWorldY refImage.PixelExtentInWorldX refImage.PixelExtentInWorldZ];
            else
                objGpet.image_size.voxelSize_mm = [refImage.PixelExtentInWorldY refImage.PixelExtentInWorldX];
            end
        end

        function objGpet =init_ref_image(objGpet)
            refImage = objGpet.init_ref_structure(objGpet.image_size.matrixSize, objGpet.image_size.voxelSize_mm);
            objGpet.ref_image = refImage;
        end

        function refImage = init_ref_structure(objGpet, matrixSize, voxelSize_mm)
            origin_mm = [-voxelSize_mm(2)*matrixSize(2)/2 -voxelSize_mm(1)*matrixSize(1)/2 ...
                -voxelSize_mm(3)*matrixSize(3)/2];
            XWorldLimits= [origin_mm(2) origin_mm(2)+voxelSize_mm(2)*matrixSize(2)];
            YWorldLimits= [origin_mm(1) origin_mm(1)+voxelSize_mm(1)*matrixSize(1)];
            ZWorldLimits= [-voxelSize_mm(3)*matrixSize(3)/2 voxelSize_mm(3)*matrixSize(3)/2] + [objGpet.bed_position_mm objGpet.bed_position_mm];
            refImage = imref3d(matrixSize, XWorldLimits, YWorldLimits, ZWorldLimits);

        end

        % This functions get the field of the struct used in Apirl and
        % converts to the struct sinogram_size used in this framework.
        function structSizeSino = get_sinogram_size_for_apirl(objPETRawData)
            structSizeSino.numR = objPETRawData.sinogram_size.nRadialBins;
            structSizeSino.numTheta = objPETRawData.sinogram_size.nAnglesBins;
            structSizeSino.numSinogramPlanes = objPETRawData.sinogram_size.nSinogramPlanes;
            structSizeSino.span = objPETRawData.sinogram_size.span;
            structSizeSino.numZ = objPETRawData.sinogram_size.nRings;
            structSizeSino.rFov_mm = 0; % not used
            structSizeSino.zFov_mm = 0; % not used
            if structSizeSino.span > 0
                structSizeSino.numSegments = objPETRawData.sinogram_size.nSeg;
                structSizeSino.sinogramsPerSegment = objPETRawData.sinogram_size.nPlanesPerSeg;
                structSizeSino.minRingDiff = objPETRawData.sinogram_size.minRingDiffs;
                structSizeSino.maxRingDiff = objPETRawData.sinogram_size.maxRingDiffs;
                structSizeSino.maxAbsRingDiff = objPETRawData.sinogram_size.maxRingDifference;
                structSizeSino.numPlanesMashed = objPETRawData.sinogram_size.numPlanesMashed;
            end
        end

    end

    % Methods in a separate file:
    methods (Access = private)

        ii = bit_reverse(objGpet, mm);
        lambda = Project_preComp(objGpet,X,g,Angles,RadialBins,dir);
        g = init_precomputed_G (objGpet);

        function objGpet =  varargin_pair(objGpet, varargs)
            % example PET = classGpet('scanner','mMR',
            % 'nSubsets',14,'PSF.Width', 2);

            npair = floor(length(varargs) / 2);
            if 2*npair ~= length(varargs),
                error('need names and values in pairs');
            end
            args = {varargs{1:2:end}};
            vals = {varargs{2:2:end}};

            prop = properties(objGpet);

            for ii=1:npair
                arg = args{ii};
                val = vals{ii};

                if ~ischar(arg)
                    error('unknown option of class %s', class(arg))
                end

                if sum(strcmpi(prop, arg )) > 0
                    objGpet.(arg) = val;
                    continue;
                end

                % for pairs with sub-fields such as : ('PSF.Width', 3)
                idot = strfind(arg, '.');
                if ~isempty(idot)
                    arg1 = arg([1:(idot-1)]);
                    arg2 = arg([(idot+1):end]);

                    if ~(sum(strcmpi(prop, arg1 )) > 0)
                        error('unknown property')
                    end
                    subProp = fieldnames(objGpet.(arg1));
                    if ~(sum(strcmpi(subProp, arg2 )) > 0)
                        error('unknown sub-property')
                    end
                    s = struct('type', {'.', '.'}, 'subs', ...
                        { arg([1:(idot-1)]), arg([(idot+1):end]) });

                    objGpet = subsasgn(objGpet, s, val);
                end

            end

        end

        function nPlanePerSeg = Planes_Seg (objGpet)

            a = (objGpet.sinogram_size.maxRingDifference -(objGpet.sinogram_size.span+1)/2);
            b = floor(a/objGpet.sinogram_size.span)+ floor(objGpet.sinogram_size.maxRingDifference/a);
            if b > 0
                nseg = 2*b+1;
            else
                nseg = 1;
            end

            a = ones(objGpet.sinogram_size.nRings);
            minRingDiff = [0,(objGpet.sinogram_size.span+1)/2:objGpet.sinogram_size.span:objGpet.sinogram_size.maxRingDifference];

            s = zeros(1,(nseg+1)/2);
            for j = 1:(nseg+1)/2
                s(j) = length(diag(a,minRingDiff(j)));
            end

            if objGpet.sinogram_size.span>1, s = 2*s-1; end

            nPlanePerSeg = zeros(1,nseg);
            nPlanePerSeg(1) = s(1);
            nPlanePerSeg(2:2:end) = s(2:end);
            nPlanePerSeg(3:2:end) = s(2:end);

        end

        function Scatter3D = iSSRB(objGpet,Scatter2D)

            nPlanePerSeg = Planes_Seg(objGpet);

            mo = cumsum(nPlanePerSeg)';
            no =[[1;mo(1:end-1)+1],mo];

            Scatter3D = zeros(objGpet.sinogram_size.matrixSize,'single');
            % For span 1 this method fails, we use the other available in
            % the library for any span conversion:
            if objGpet.sinogram_size.span > 1

                Scatter3D(:,:,no(1,1):no(1,2),:) = Scatter2D;

                for i = 2:2:length(nPlanePerSeg)

                    delta = (nPlanePerSeg(1)- nPlanePerSeg(i))/2;
                    indx = nPlanePerSeg(1) - delta;

                    Scatter3D (:,:,no(i,1):no(i,2),:) = Scatter2D(:,:,delta+1:indx,:);
                    Scatter3D (:,:,no(i+1,1):no(i+1,2),:) = Scatter2D(:,:,delta+1:indx,:);
                end
            else
                % issrb is for an input sinogram span of span 121
                structSizeSino3d = getSizeSino3dFromSpan(objGpet.sinogram_size.nRadialBins, objGpet.sinogram_size.nAnglesBins, objGpet.sinogram_size.nRings, ...
                    0, 0, 121, objGpet.sinogram_size.maxRingDifference);
                [Scatter3D, structSizeSino3dSpanN] = convertSinogramToSpan(Scatter2D, structSizeSino3d,  objGpet.sinogram_size.span);
            end
        end

        function scatter_3D = scatter_scaling(objGpet,scatter_3D, emission_sinogram, ncf, acf, randoms)

            gaps = ncf==0;
            scatter_3D =  scatter_3D./ncf;
            scatter_3D(gaps)=0;

            Trues = emission_sinogram - randoms;
            Trues(Trues<0) = 0;
            for i = 1: size(acf,3)
                acf_i = acf(:,:,i);
                %                 mask = acf_i <= min(acf(:));
                mask = acf_i <1.03;
                scatter_3D_tail = scatter_3D(:,:,i).*mask;
                Trues_tail = Trues(:,:,i).*mask;

                scale_factor = sum(Trues_tail(:))./sum(scatter_3D_tail(:));
                scatter_3D(:,:,i) = scatter_3D(:,:,i)/scale_factor;

            end

        end

    end
    methods (Access = public)
        % Now vec division is public, it can be useful externally too:
        % SAM ELLIS EDIT (18/07/2016): new method to allow easy division without
        % using a small additive term to avoid div. by zero.
        function c = vecDivision(objGpet,a,b)
            % element-by-element division of two vectors, a./b, BUT avoiding
            % division by 0
            c = a;
            c(b~=0) = a(b~=0)./b(b~=0);
            c(b==0) = 0;
        end

        % Project:
        m = P(objGpet, x,subset_i, localNumSubsets);
        % Backproject:
        x = PT(objGpet,m, subset_i, localNumSubsets);
        % Normalization correction factors:
        [n, n_ti, n_tv, gaps] = NCF(varargin);
        % Attenuation correction factors:
        a=ACF(varargin);
        % Randoms:
        r=R(varargin);
        % Scatter
        s=S(varargin);
        %
        osem_subsets(objGpet, nsub,nAngles);

        function gaps = gaps(objGpet)
            if strcmp(objGpet.scanner,'mMR')
                if objGpet.sinogram_size.span > 0
                    crystalMasks = ones(504,64);
                    crystalMasks(9:9:end,:) = 0;
                    gaps = createSinogram3dFromDetectorsEfficency(crystalMasks, objGpet.get_sinogram_size_for_apirl(), 0);
                else
                    crystalMasks = ones(504);
                    crystalMasks(9:9:end,:) = 0;
                    gaps = createSinogram2dFromDetectorsEfficency(crystalMasks, objGpet.get_sinogram_size_for_apirl(), 2, 0);
                end
            else
                gaps = ones(objGpet.sinogram_size.matrixSize);
            end
        end

        function x = ones(objGpet)
            if objGpet.sinogram_size.span == -1
                %x = ones([objGpet.image_size.matrixSize(1) objGpet.image_size.matrixSize(2)]);
                [x0,y0] = meshgrid(-objGpet.image_size.matrixSize(1)/2+1:objGpet.image_size.matrixSize(2)/2);
                x = (x0.^2+y0.^2)<(objGpet.image_size.matrixSize(1)/2.5)^2;
            else
                [x0,y0] = meshgrid(-objGpet.image_size.matrixSize(1)/2+1:objGpet.image_size.matrixSize(2)/2);
                x = (x0.^2+y0.^2)<(objGpet.image_size.matrixSize(1)/2.5)^2;
                x = repmat(x,[1,1,objGpet.image_size.matrixSize(3)]);
            end
            x = single(x);
        end

        function x = zeros(objGpet)

            x = zeros(objGpet.image_size.matrixSize);
        end

        init_sinogram_size(objGpet, inSpan, numRings, maxRingDifference);

        function sino_size = get_sinogram_size(objGpet)
            sino_size = objGpet.sinogram_size;
        end

        function setBedPosition(objGpet, siemensInterfile)
            [auxImage, auxRef, objGpet.bed_position_mm] = interfileReadSiemensImage(siemensInterfile);
            objGpet.init_ref_image();
        end

        % Function that maps MR into PET image space:
        function [MrInPet, refResampledImage] = getMrInPetImageSpace(objGpet, pathMrDicom)
            % Read dicom image:
            [imageMr, refMrImage, affineMatrix, dicomInfo] = ReadDicomImage(pathMrDicom, '', 1);
            % Convert into the nes image space
            [MrInPet, refResampledImage] = ImageResample(imageMr, refMrImage, objGpet.ref_image);
        end

        % Function that reads MR and  configures PET to reconstruct in that
        % space. The first input parameter is the path of the dicom image
        % and the second (optional) is logical value that says if the
        % object update the image size to refImagePetFov.
        function [imageMr, refImageMr, MrInPetFov, refImagePetFov] = getMrInNativeImageSpace(objGpet, pathMrDicom, updateImageSize)
            if nargin == 2
                updateImageSize = 0;
            end
            % Read dicom image:
            [imageMr, refImageMr, affineMatrix, dicomInfo] = ReadDicomImage(pathMrDicom, '', 1);
            % Update the ref_image to this pixel size:
            newVoxelSize= [refImageMr.PixelExtentInWorldY refImageMr.PixelExtentInWorldX refImageMr.PixelExtentInWorldZ];
            ratio = objGpet.image_size.voxelSize_mm ./ newVoxelSize;
            newVoxelSize_mm = newVoxelSize;
            newMatrixSize = round(objGpet.image_size.matrixSize .* ratio);
            refImagePetFov = objGpet.init_ref_structure(newMatrixSize, newVoxelSize_mm);
            % Update the image size if required:
            if updateImageSize
                objGpet.image_size.voxelSize_mm = newVoxelSize;
                objGpet.image_size.matrixSize = round(objGpet.image_size.matrixSize .* ratio);
                objGpet.init_ref_image();
                %refImagePet = objGpet.ref_image;
            end
            % And finally complete the MR image into PET FOV:
            [MrInPetFov, refImagePetFov] = ImageResample(imageMr, refImageMr, refImagePetFov);
        end

        % function that maps an MR in reference space to PET Fov with a FOV
        % trimming for storage reasons
        function [MrInPetFov, refImagePetFov] = mapMrRefSpaceToPetFov(objGpet,imageMr,refImageMr,FovReductionFactor)
            if nargin ==3
                voxelSize_mm = objGpet.image_size.voxelSize_mm;
                matrixSize = objGpet.image_size.matrixSize;
            else % for the case of trimming the FOV of MrInPetFov images
                FovReductionFactor = max(2.5, FovReductionFactor);
                I = floor(objGpet.image_size.matrixSize(1)/FovReductionFactor);
                J = floor(objGpet.image_size.matrixSize(2)/FovReductionFactor);
                matrixSize = [length((I:(objGpet.image_size.matrixSize(1)-I-1))+1),length((J:(objGpet.image_size.matrixSize(2)-J-1))+1),objGpet.image_size.matrixSize(3)];
                voxelSize_mm = [objGpet.ref_native_image.PixelExtentInWorldX, objGpet.ref_native_image.PixelExtentInWorldY, objGpet.ref_native_image.PixelExtentInWorldZ];

                origin_mm = [-voxelSize_mm(2)*matrixSize(2)/2 -voxelSize_mm(1)*matrixSize(1)/2 ...
                    -voxelSize_mm(3)*matrixSize(3)/2];
                XWorldLimits= [origin_mm(2) origin_mm(2)+voxelSize_mm(2)*matrixSize(2)];
                YWorldLimits= [origin_mm(1) origin_mm(1)+voxelSize_mm(1)*matrixSize(1)];
                ZWorldLimits= objGpet.ref_native_image.ZWorldLimits;
                refImagePet = imref3d(matrixSize, XWorldLimits, YWorldLimits, ZWorldLimits);

                voxelSize_mm = [refImagePet.PixelExtentInWorldX, refImagePet.PixelExtentInWorldY, refImagePet.PixelExtentInWorldZ];
                matrixSize = refImagePet.ImageSize;
            end

            % Update the ref_image to this pixel size:
            newVoxelSize= [refImageMr.PixelExtentInWorldY refImageMr.PixelExtentInWorldX refImageMr.PixelExtentInWorldZ];
            ratio = voxelSize_mm ./ newVoxelSize;
            newVoxelSize_mm = newVoxelSize;
            newMatrixSize = round(matrixSize .* ratio);

            if nargin==3
                refImagePetFov = objGpet.init_ref_structure(newMatrixSize, newVoxelSize_mm);
            else

                origin_mm = [-newVoxelSize_mm(2)*newMatrixSize(2)/2 -newVoxelSize_mm(1)*newMatrixSize(1)/2 ...
                    -newVoxelSize_mm(3)*newMatrixSize(3)/2];
                XWorldLimits= [origin_mm(2) origin_mm(2)+newVoxelSize_mm(2)*newMatrixSize(2)];
                YWorldLimits= [origin_mm(1) origin_mm(1)+newVoxelSize_mm(1)*newMatrixSize(1)];
                ZWorldLimits= refImagePet.ZWorldLimits;
                refImagePetFov = imref3d(newMatrixSize, XWorldLimits, YWorldLimits, ZWorldLimits);
            end



            [MrInPetFov, refImageMrFov] = ImageResample(imageMr, refImageMr, refImagePetFov);
        end

        % Converts sinogram
        function [sino_compressed, structSizeSino3dSpanN] = apply_axial_compression_from_span1(objGpet, sinogram)
            structSizeSino3d = getSizeSino3dFromSpan(objGpet.sinogram_size.nRadialBins, objGpet.sinogram_size.nAnglesBins, objGpet.sinogram_size.nRings, ...
                0, 0, 1, objGpet.sinogram_size.maxRingDifference);
            if numel(sinogram) ~= (objGpet.sinogram_size.nRadialBins*objGpet.sinogram_size.nAnglesBins*sum(structSizeSino3d.sinogramsPerSegment))
                error('The input sinogram needs to be span 1.');
            end
            [sino_compressed, structSizeSino3dSpanN] = convertSinogramToSpan(sinogram, structSizeSino3d,  objGpet.sinogram_size.span);
        end

        function mask = get_fov_maks(objGpet, radiusFov_mm)
            if nargin == 1
                radiusFov_mm = 596/2;
            end
            mask = objGpet.ones;
            radiusFov_pixels = radiusFov_mm/objGpet.image_size.voxelSize_mm(1);
            [X,Y,Z] = meshgrid(1:objGpet.image_size.matrixSize(2),1:objGpet.image_size.matrixSize(1),1:objGpet.image_size.matrixSize(3));
            indicesOutMask = ((X-objGpet.image_size.matrixSize(2)/2).^2+(Y-objGpet.image_size.matrixSize(1)/2).^2) > radiusFov_pixels.^2;
            mask(indicesOutMask) = 0;
        end

        function senseImg = Sensitivity(objGpet, AN)
            % SAM ELLIS EDIT: IF AN IS DOUBLE, THEN LET THE SENSEIMG BE
            % DOUBLE TOO
            classAN = whos('AN');
            classAN = classAN.class;

            senseImg = zeros([objGpet.image_size.matrixSize, objGpet.nSubsets],classAN) ;
            % SAM ELLIS EDIT (23/08/2016): only show messages if using more
            % than one subset
            if objGpet.nSubsets > 1
                for n = 1:objGpet.nSubsets
                    fprintf('%d, ',n);
                    senseImg(:,:,:,n) = objGpet.PT(AN,n);
                end
            else
                senseImg = objGpet.PT(AN);
            end
            fprintf('Done.\n');
        end

        function display(objGpet) %#ok<DISPLAY>
            disp(objGpet)
            methods(objGpet)
        end

        function objGpet = Revise(objGpet,opt)
            % to revise the properties of a given object without
            % re-instantiation
            vfields = fieldnames(opt);
            prop = properties(objGpet);
            for i = 1:length(vfields)
                % if is an struct, update only the received fields in that
                % structure, for example for sinogram_size:
                if isstruct(opt.(vfields{i}))
                    field = vfields{i};
                    if sum(strcmpi(prop, field )) > 0 && ~strcmpi('Prior', field ) % prior is intialized in other part
                        subprop = fieldnames(objGpet.(field));
                        secondary_fields = fieldnames(opt.(field));
                        for j = 1 : length(secondary_fields)
                            subfield = secondary_fields{j};
                            if sum(strcmpi(subprop, subfield )) > 0
                                objGpet.(field).(secondary_fields{j}) = opt.(field).(secondary_fields{j});
                            end
                        end
                    end
                else
                    % It isnt
                    field = vfields{i};
                    if sum(strcmpi(prop, field )) > 0
                        objGpet.(field) = opt.(field);
                    end
                end
            end
            objGpet.init_sinogram_size(objGpet.sinogram_size.span, objGpet.sinogram_size.nRings, objGpet.sinogram_size.maxRingDifference);
            if isfield(opt,'nSubsets')
                objGpet.set_subsets(opt.nSubsets)
            end

            % SAM ELLIS EDIT (23/08/2016): for 2D radon, need to update
            % sinogram size based on new image size
            if any(strcmpi(vfields,'image_size')) && strcmpi(objGpet.scanner,'2D_radon')
                x = radon(ones(objGpet.image_size.matrixSize(1:2)),0:179);
                x = size(x);
                objGpet.sinogram_size.nRadialBins = x(1);
                objGpet.sinogram_size.nAnglesBins = x(2);
                objGpet.sinogram_size.nSinogramPlanes = 1;
                objGpet.sinogram_size.span = -1;
                objGpet.sinogram_size.nRings = 1;
                objGpet.sinogram_size.maxRingDifference = 0;
            end
        end

        function Img = OPMLEM(objGpet,Prompts,RS, SensImg,Img, nIter)
            for i = 1:nIter
                % SAM ELLIS EDIT (18/07/2016): replaced vector divisions by vecDivision
                Img = Img.*objGpet.vecDivision(objGpet.PT(objGpet.vecDivision(Prompts,objGpet.P(Img)+ RS)),SensImg);
                Img = max(0,Img);
            end
        end

        function image_iter = OPMLEMsaveIter(objGpet,Prompts, AN, RS, SensImg, initialEstimate, nIter, outputPath, saveInterval)
            if ~isdir(outputPath)
                mkdir(outputPath);
            end
            image = initialEstimate;
            k=1;
            for i = 1:nIter
                % SAM ELLIS EDIT (18/07/2016): replaced vector divisions by vecDivision
                image = image.*objGpet.vecDivision(objGpet.PT(AN.*objGpet.vecDivision(Prompts,AN.*objGpet.P(image)+ RS)),SensImg);
                image = max(0,image);
                if rem(i,saveInterval) == 0 %
                    interfilewrite(single(image), [outputPath 'mlem_iter_' num2str(i)], objGpet.image_size.voxelSize_mm); % i use i instead of i+1 because i=1 is the inital estimate
                    image_iter{k} = image;
                    k = k + 1;
                end
            end
        end

        % Opmlem with downsample
        function [image, image_ds] = OPMLEM_DS(objGpet, prompts, anf, additive, initialEstimate, numIterations, outputPath, saveInterval)
            if nargin >=7
                if ~isdir(outputPath)
                    mkdir(outputPath);
                end
            end
            % Create grids for downsample and high sample:
            x_lowres = objGpet.ref_native_image.XWorldLimits+objGpet.ref_native_image.PixelExtentInWorldX(1)/2 : objGpet.ref_native_image.PixelExtentInWorldX : objGpet.ref_native_image.XWorldLimits(2)-objGpet.ref_native_image.PixelExtentInWorldX/2;
            y_lowres = objGpet.ref_native_image.YWorldLimits+objGpet.ref_native_image.PixelExtentInWorldY(1)/2 : objGpet.ref_native_image.PixelExtentInWorldY : objGpet.ref_native_image.YWorldLimits(2)-objGpet.ref_native_image.PixelExtentInWorldY/2;
            z_lowres = objGpet.ref_native_image.ZWorldLimits+objGpet.bed_position_mm+objGpet.ref_native_image.PixelExtentInWorldZ(1)/2 : objGpet.ref_native_image.PixelExtentInWorldZ : objGpet.ref_native_image.ZWorldLimits(2)+objGpet.bed_position_mm-objGpet.ref_native_image.PixelExtentInWorldZ/2;
            x_highres = objGpet.ref_image.XWorldLimits+objGpet.ref_image.PixelExtentInWorldX(1)/2 : objGpet.ref_image.PixelExtentInWorldX : objGpet.ref_image.XWorldLimits(2)-objGpet.ref_image.PixelExtentInWorldX/2;
            y_highres = objGpet.ref_image.YWorldLimits+objGpet.ref_image.PixelExtentInWorldY(1)/2 : objGpet.ref_image.PixelExtentInWorldY : objGpet.ref_image.YWorldLimits(2)-objGpet.ref_image.PixelExtentInWorldY/2;
            z_highres = objGpet.ref_image.ZWorldLimits+objGpet.ref_image.PixelExtentInWorldZ(1)/2 : objGpet.ref_image.PixelExtentInWorldZ : objGpet.ref_image.ZWorldLimits(2)-objGpet.ref_image.PixelExtentInWorldZ/2;
            [X_lowres,Y_lowres,Z_lowres] = meshgrid(x_lowres, y_lowres, z_lowres);
            [X_highres,Y_highres,Z_highres] = meshgrid(x_highres, y_highres, z_highres);
            % Masks:
            mask_highres = (X_highres.^2+Y_highres.^2)<(objGpet.image_size.matrixSize(1)*objGpet.image_size.voxelSize_mm(1)/2.5)^2;
            mask = (X_lowres.^2+Y_lowres.^2)<(objGpet.ref_native_image.ImageSize(1)*objGpet.ref_native_image.PixelExtentInWorldX(1)/2.5)^2;
            % PET ds:
            paramPET.scanner = objGpet.scanner;
            paramPET.method =  objGpet.method;
            paramPET.PSF.type = objGpet.PSF.type;
            paramPET.PSF.Width = objGpet.PSF.Width;
            paramPET.sinogram_size.span = objGpet.sinogram_size.span;
            paramPET.nSubsets = objGpet.nSubsets;
            paramPET.verbosity = 0;
            paramPET.tempPath = objGpet.tempPath;
            PET_lowres = classGpet(paramPET);

            % The sensitivity image needs to include the interpolation
            % matrix:
            sensImage = PET_lowres.Sensitivity(anf);
            sensImg_highres = interp3(X_lowres, Y_lowres, Z_lowres, sensImage, X_highres, Y_highres, Z_highres, 'linear', 0); %imresize(sensImage, PET_highres.image_size.matrixSize, 'bicubic'); % High resolution image
            %sensImg_highres = interp3_lin(X_lowres, Y_lowres, Z_lowres, sensImage, X_highres, Y_highres, Z_highres); % this is my function to have matched projector/backprojector
            image = initialEstimate;
            k = 1;
            for i = 1:numIterations
                % Projection:
                low_res_image = interp3(X_highres, Y_highres, Z_highres, image, X_lowres, Y_lowres, Z_lowres, 'linear', 0); %imresize(opmlem{end}, PET_highres.image_size.matrixSize, 'bicubic');
                %low_res_image = interp3_lin(X_highres, Y_highres, Z_highres, image, X_lowres, Y_lowres, Z_lowres); % this is my function to have matched projector/backprojector
                projected = anf.*PET_lowres.P(low_res_image) + additive;
                % Backproject:
                backprojected_image = PET_lowres.PT(anf.*PET_lowres.vecDivision(prompts, projected)).*mask;
                % transpose of interpolation (high sample)
                backprojected_image_highres = interp3(X_lowres, Y_lowres, Z_lowres, backprojected_image, X_highres, Y_highres, Z_highres, 'linear', 0).*mask_highres; %imresize(opmlem{end}, PET_highres.image_size.matrixSize, 'bicubic');
                %backprojected_image_highres = interp3_lin(X_lowres, Y_lowres, Z_lowres, backprojected_image, X_highres, Y_highres, Z_highres).*mask_highres;
                % Update image
                image = image.*PET_lowres.vecDivision(backprojected_image_highres, sensImg_highres);
                image = max(0,image);
                if nargin>=7
                    if rem(i-1,saveInterval) == 0 % -1 to save the first iteration
                        image_ds{k} = image;
                        interfilewrite(single(image_ds{k}), [outputPath 'mlem_ds_iter_' num2str(i)], [objGpet.ref_image.PixelExtentInWorldX objGpet.ref_image.PixelExtentInWorldY objGpet.ref_image.PixelExtentInWorldZ]); % i use i instead of i+1 because i=1 is the inital estimate
                        k = k + 1;
                    end
                else
                    image_ds = image;
                end
            end
        end

        function [senseImg] = Sensitivity_DS(objGpet, anf)
        % Create grids for downsample and high sample:
            x_lowres = objGpet.ref_native_image.XWorldLimits+objGpet.ref_native_image.PixelExtentInWorldX(1)/2 : objGpet.ref_native_image.PixelExtentInWorldX : objGpet.ref_native_image.XWorldLimits(2)-objGpet.ref_native_image.PixelExtentInWorldX/2;
            y_lowres = objGpet.ref_native_image.YWorldLimits+objGpet.ref_native_image.PixelExtentInWorldY(1)/2 : objGpet.ref_native_image.PixelExtentInWorldY : objGpet.ref_native_image.YWorldLimits(2)-objGpet.ref_native_image.PixelExtentInWorldY/2;
            z_lowres = objGpet.ref_native_image.ZWorldLimits+objGpet.bed_position_mm+objGpet.ref_native_image.PixelExtentInWorldZ(1)/2 : objGpet.ref_native_image.PixelExtentInWorldZ : objGpet.ref_native_image.ZWorldLimits(2)+objGpet.bed_position_mm-objGpet.ref_native_image.PixelExtentInWorldZ/2;
            x_highres = objGpet.ref_image.XWorldLimits+objGpet.ref_image.PixelExtentInWorldX(1)/2 : objGpet.ref_image.PixelExtentInWorldX : objGpet.ref_image.XWorldLimits(2)-objGpet.ref_image.PixelExtentInWorldX/2;
            y_highres = objGpet.ref_image.YWorldLimits+objGpet.ref_image.PixelExtentInWorldY(1)/2 : objGpet.ref_image.PixelExtentInWorldY : objGpet.ref_image.YWorldLimits(2)-objGpet.ref_image.PixelExtentInWorldY/2;
            z_highres = objGpet.ref_image.ZWorldLimits+objGpet.ref_image.PixelExtentInWorldZ(1)/2 : objGpet.ref_image.PixelExtentInWorldZ : objGpet.ref_image.ZWorldLimits(2)-objGpet.ref_image.PixelExtentInWorldZ/2;
            [X_lowres,Y_lowres,Z_lowres] = meshgrid(x_lowres, y_lowres, z_lowres);
            [X_highres,Y_highres,Z_highres] = meshgrid(x_highres, y_highres, z_highres);
            % PET ds:
            paramPET.scanner = objGpet.scanner;
            paramPET.method =  objGpet.method;
            paramPET.PSF.type = objGpet.PSF.type;
            paramPET.PSF.Width = objGpet.PSF.Width;
            paramPET.sinogram_size.span = objGpet.sinogram_size.span;
            paramPET.nSubsets = objGpet.nSubsets;
            paramPET.verbosity = 0;
            paramPET.tempPath = objGpet.tempPath;
            PET_lowres = classGpet(paramPET);

            % The sensitivity image needs to include the interpolation
            % matrix:
            senseImg = PET_lowres.Sensitivity(anf);
            senseImg = interp3(X_lowres, Y_lowres, Z_lowres, senseImg, X_highres, Y_highres, Z_highres, 'linear', 0); %imresize(sensImage, PET_highres.image_size.matrixSize, 'bicubic'); % High resolution image
        end

        function image_iters = OPOSEMsaveIter(objGpet,Prompts, AN, RS, SensImg, initialEstimate, nIter, outputPath, saveInterval)
            k=1;
            image = initialEstimate;
            for i = 1:nIter
                for j = 1:objGpet.nSubsets
                    % SAM ELLIS EDIT (18/07/2016): replaced vector divisions by vecDivision
                    image = image.*objGpet.vecDivision(objGpet.PT(AN.*objGpet.vecDivision(Prompts,AN.*objGpet.P(image,j)+ RS),j),SensImg(:,:,:,j));
                    image = max(0,image);
                    if rem(k-1,saveInterval) == 0 % -1 to save the first iteration
                        interfilewrite(single(image), [outputPath 'oposem_subiter_' num2str(k)], objGpet.image_size.voxelSize_mm); % i use i instead of i+1 because i=1 is the inital estimate
                        image_iters{k} = image;
                        k = k + 1;
                    end
                end
            end
        end

        function Img = OPOSEM(objGpet,Prompts,RS, SensImg,Img, nIter, arg)
            opt.display = 0;
            opt = getFiledsFromUsersOpt(opt,arg);
            if opt.display, fprintf('\n'); t = tic; end
            for i = 1:nIter
                for j = 1:objGpet.nSubsets
                    % SAM ELLIS EDIT (18/07/2016): replaced vector divisions by vecDivision
                    Img = Img.*objGpet.vecDivision(objGpet.PT(objGpet.vecDivision(Prompts,objGpet.P(Img,j)+ RS),j),SensImg(:,:,:,j));
                    % Img = Img.*objGpet.PT(Prompts./(objGpet.P(Img,j)+ RS + 1e-5),j)./(SensImg(:,:,:,j)+1e-5);
                    Img = max(0,Img);
                end
                if opt.display, fprintf('Iter: %4d/%4d\r', i, nIter); end
            end
            if opt.display, fprintf('\n'); toc(t), end
        end

        % SAM ELLIS EDIT: new sensitivity image for basis function
        % coefficients
        function alphSens = basisSensitivity(objGpet,AN,basisMat)

            % define function to multiply image by matrix and reshape
            % (assuming numel(alph) ==  numel(Img))
            multiplyByBases = @(alphIn,basisMat) reshape(basisMat*double(alphIn(:)),size(alphIn));

            % get sensitivity image, using input att.*norm factors
            SensImg = Sensitivity(objGpet, AN);

            % calcuate new sens image in terms of the alpha coefficients
            alphSens = multiplyByBases(SensImg,basisMat');
        end

        % SAM ELLIS EDIT: new iterative method using a basis function matrix (full matrix required)
        function [alph,Img] = BFEM(objGpet,basisMat,Prompts,RS, AN, alphSens, alph, nIter)
            % assume the basisFcns are in a matrix for now..

            % define function to multiply image by matrix and reshape
            % (assuming numel(alph) ==  numel(Img))
            multiplyByBases = @(alphIn,basisMat) reshape(basisMat*double(alphIn(:)),size(alphIn));

            for ii = 1:nIter

                % forward model:
                FP_guess = AN.*objGpet.P(multiplyByBases(alph,basisMat)) + RS;

                % ratio with measured data and backproject (including
                % transpose basis function multiplication)
                corr_fact = multiplyByBases(objGpet.PT(AN.*objGpet.vecDivision(Prompts,FP_guess)),basisMat');

                % divide by sens image and multiply by old alpha values
                alph = alph.*objGpet.vecDivision(corr_fact,alphSens);
            end

            Img = multiplyByBases(alph,basisMat);

        end

        % SAM ELLIS EDIT (24/08/2016): new method to delete the temp folder
        % when we're done
        function delete(objGpet)
            if (objGpet.deleteTemp ~= false) && (objGpet.deleteTemp ~= 0)
                % if the flag is not false, then it is true, and delete the
                % temp folder and its contents if it exists
                if exist(objGpet.tempPath,'dir') == 7
                    rmdir(objGpet.tempPath,'s');
                end
            else
            end
        end

        function Img = MAPEM(objGpet,Prompts,RS, SensImg,Img, nIter,arg)
            % First check if the Prior class has already been intialized
            % default parameters
            opt.PetOptimizationMethod = 'DePierro';%'OSL'
            opt.PetPriorType = 'Quadratic'; %'Bowsher' 'JointBurgEntropy'
            opt.PetRegularizationParameter = 1;
            opt.PetPreCompWeights = 1;
            opt.TVsmoothingParameter = 0.1;
            opt.LangeDeltaParameter = 1;
            opt.BowsherB = 70;
            opt.PriorMrImage =[];
            opt.MrSigma = 0.1; % JBE
            opt.PetSigma  = 10; %JBE
            opt.display = 0;
            % check if the object already was initilized:
            if isempty(objGpet.Prior)
                opt.ImageSize = objGpet.image_size.matrixSize;
                opt.imCropFactor = [4,4,0];
                opt.sWindowSize = 5;
                opt.lWindowSize = 1;
                opt = getFiledsFromUsersOpt(opt,arg);
                objGpet.Prior = PriorsClass(opt);
            else
                % whos_prior = whos('objGpet.Prior');
                % if strcmp(whos_prior.class, 'PriorsClass')
                % Only initilize the parameters for this reconstruction
                opt = getFiledsFromUsersOpt(opt,arg);
            end

            if opt.display, figure; end

            if strcmpi(opt.PetPriorType,'Quadratic')
                if opt.PetPreCompWeights == 1
                    W0 = 1/objGpet.Prior.nS;
                else
                    W0 = opt.PetPreCompWeights./repmat(sum(opt.PetPreCompWeights,2),[1,objGpet.Prior.nS]); %for weighted quadratic
                end
            elseif strcmpi(opt.PetPriorType,'Bowsher')
                if opt.PetPreCompWeights == 1
                    fprintf('calculating Bowsher weighting coeffcients\n')
                    if isempty(opt.PriorMrImage), error('MR image should be provided in opt.PriorMrImage\n'); end
                    W0 = objGpet.Prior.W_Bowsher(opt.PriorMrImage,opt.BowsherB);
                    W0 = W0./repmat(sum(W0,2),[1,objGpet.Prior.nS]);
                else % better to precomputed Bowsher weighting coeffcients into opt.PetPreCompWeights
                    W0 = opt.PetPreCompWeights;
                end
            elseif strcmpi(opt.PetPriorType,'JointBurgEntropy')
                if opt.PetPreCompWeights == 1
                    fprintf('calculating MR-based Gaussian weighting coeffcients\n');
                    if isempty(opt.PriorMrImage), error('MR image should be provided in opt.PriorMrImage\n'); end
                    W0 = objGpet.Prior.W_JointEntropy(opt.PriorMrImage,opt.MrSigma);
                else % better to precomputed MR weighting coeffcients into opt.PetPreCompWeights
                    W0 = opt.PetPreCompWeights;
                end
            end
           if opt.display, fprintf('Prior: %s, Method: %s\n',opt.PetPriorType,opt.PetOptimizationMethod); end
            for i = 1:nIter
                if opt.display, fprintf('Iteration: %d\n',i); end
                if strcmpi(opt.PetOptimizationMethod,'DePierro')
                    xn = Img;
                    x_em = xn.*objGpet.vecDivision(objGpet.PT(objGpet.vecDivision(Prompts,objGpet.P(xn)+ RS)),SensImg);

                    if strcmpi(opt.PetPriorType,'JointBurgEntropy')
                        W0 = objGpet.Prior.W_JointEntropy(xn,opt.PetSigma).*W0;
                        W0 = W0./repmat(sum(W0,2),[1,objGpet.Prior.nS]);
                    end
                    W = objGpet.Prior.Wd.*W0;
                    wj = objGpet.Prior.UndoImCrop(reshape(sum(W,2),objGpet.Prior.CropedImageSize));
                    B = SensImg - opt.PetRegularizationParameter / 2*objGpet.Prior.UndoImCrop(reshape(sum(W.*objGpet.Prior.GraphDivCrop(xn),2),objGpet.Prior.CropedImageSize));
                    Img = 2*x_em.*SensImg./(B + sqrt(B.^2+4*opt.PetRegularizationParameter.*SensImg.*x_em.*wj + 1e-5));
                    Img = max(0,Img);
                elseif strcmpi(opt.PetOptimizationMethod,'OSL')
                    xn = Img;
                    if strcmpi(opt.PetPriorType,'JointBurgEntropy')
                        W0 = objGpet.Prior.W_JointEntropy(xn,opt.PetSigma).*W0;
                        W0 = W0./repmat(sum(W0,2),[1,objGpet.Prior.nS]);
                    elseif strcmpi(opt.PetPriorType,'sTV')
                        imgGrad = objGpet.Prior.GraphGradCrop(xn);
                        Norm = repmat(sqrt(sum(abs(imgGrad).^2,2)+ opt.TVsmoothingParameter),[1,objGpet.Prior.nS]);
                        W0 = 1./(Norm+eps);
                    elseif strcmpi(opt.PetPriorType,'lange')
                        imgGrad = objGpet.Prior.GraphGradCrop(xn);
                        Norm = repmat(sqrt(sum(objGpet.Prior.Wd.*imgGrad.^2,2) + opt.TVsmoothingParameter),[1,objGpet.Prior.nS]);
                        W0 = zeros(size(Norm));
                        W0(Norm~=0) = 1./Norm(Norm~=0).*(1 + opt.LangeDeltaParameter./(opt.LangeDeltaParameter + Norm(Norm~=0)));
                    end
                    W = objGpet.Prior.Wd.*W0;
                    dP = -2* sum(W.*objGpet.Prior.GraphGradCrop(xn),2);
                    dP = opt.PetRegularizationParameter*objGpet.Prior.UndoImCrop(reshape(dP,objGpet.Prior.CropedImageSize));
                    Img = xn.*objGpet.vecDivision(objGpet.PT(objGpet.vecDivision(Prompts,objGpet.P(xn)+ RS)),SensImg + dP + 1e-5);
                    Img = max(0,Img);
                end
                if opt.display
                    if objGpet.image_size.matrixSize(3)==1
                        imshow(Img,[]);
                    else %if 3D recon, use opt.display =x, where x is a transvese slice
                        imshow(Img(:,:,opt.display),[]);
                    end
                end
            end
        end

        function Img = MAPEM2(objGpet,Prompts,RS, SensImg,Img, nIter,arg)
            % First check if the Prior class has already been intialized
            % default parameters
            opt.PetOptimizationMethod = 'DePierro';%'OSL'
            opt.PetPriorType = 'Quadratic'; %'TV' 'Lange'
            opt.PetSimilarityKernel = 'local';  %'Bowsher' 'JointBurgEntropy'
            opt.PriorImplementation = 'matlab'; % 2 options: 'matlab', 'mex-cuda'
            opt.PetRegularizationParameter = 1;
            opt.PetPreCompWeights = 1;
            opt.TVsmoothingParameter = 0.1;
            opt.LangeDeltaParameter = 1;
            opt.BowsherB = 70;
            opt.PriorImage =[];
            opt.MrSigma = 0.1; % JBE
            opt.PetSigma  = 10; %JBE
            opt.display = 0;
            opt.save = 0;
            opt.outputPath = [];
            opt.saveInterval = 0;
            opt.ImageSize = objGpet.image_size.matrixSize;
            opt.imCropFactor = [4,4,0];
            opt.sWindowSize = 5;
            opt.lWindowSize = 1;
            opt = getFiledsFromUsersOpt(opt,arg);
            objGpet.Prior = PriorsClass(opt);

            if opt.display, figure; end


            if opt.display, fprintf('Prior: %s, Method: %s\n',opt.PetPriorType,opt.PetOptimizationMethod); end
            for i = 1:nIter
                if opt.display, fprintf('Iteration: %d\n',i); end
                if strcmpi(opt.PetOptimizationMethod,'DePierro')
                    xn = Img;
                    x_em = xn.*objGpet.vecDivision(objGpet.PT(objGpet.vecDivision(Prompts,objGpet.P(xn)+ RS)),SensImg);

                    if strcmpi(opt.PetPriorType,'JointBurgEntropy')
                        W0 = objGpet.Prior.W_JointEntropy(xn,opt.PetSigma).*W0;
                        W0 = W0./repmat(sum(W0,2),[1,objGpet.Prior.nS]);
                    end
                    W = objGpet.Prior.Wd.*W0;
                    wj = objGpet.Prior.UndoImCrop(reshape(sum(W,2),objGpet.Prior.CropedImageSize));
                    B = SensImg - opt.PetRegularizationParameter / 2*objGpet.Prior.UndoImCrop(reshape(sum(W.*objGpet.Prior.GraphDivCrop(xn),2),objGpet.Prior.CropedImageSize));
                    Img = 2*x_em.*SensImg./(B + sqrt(B.^2+4*opt.PetRegularizationParameter.*SensImg.*x_em.*wj + 1e-5));
                    Img = max(0,Img);
                elseif strcmpi(opt.PetOptimizationMethod,'OSL')
                    xn = Img;
                    dP = opt.PetRegularizationParameter*objGpet.Prior.dPrior(xn,opt);
                    Img = xn.*objGpet.vecDivision(objGpet.PT(objGpet.vecDivision(Prompts,objGpet.P(xn)+ RS)),SensImg + dP + 1e-5);
                    Img = max(0,Img);
                end
                if opt.display
                    if objGpet.image_size.matrixSize(3)==1
                        imshow(Img,[]);
                    else %if 3D recon, use opt.display =x, where x is a transvese slice
                        imshow(Img(:,:,opt.display),[]);
                    end
                end
                if opt.save
                    if rem(i,opt.saveInterval) == 0 %
                        interfilewrite(single(Img), [opt.outputPath 'map_iter_' num2str(i)], [objGpet.image_size.matrixSize]); % i use i instead of i+1 because i=1 is the inital estimate
                    end
                end
            end
        end

        gf3d = Gauss3DFilter (objGpet, data, fwhm);
        [Img,totalScaleFactor, info] = BQML(objGpet,Img,sinogramInterFileFilename,normalizationInterFileFilename);
        Img = SUV(objGpet,Img,sinogramInterFileFilename,normalizationInterFileFilename);


        function image_ds = MAPEM_DS(objGpet,prompts, anf, additive, Img, nIter,arg, outputPath, saveInterval)
            % First initilization of the parameters needed to downsample in
            % system matrix:
            x_lowres = objGpet.ref_native_image.XWorldLimits+objGpet.ref_native_image.PixelExtentInWorldX(1)/2 : objGpet.ref_native_image.PixelExtentInWorldX : objGpet.ref_native_image.XWorldLimits(2)-objGpet.ref_native_image.PixelExtentInWorldX/2;
            y_lowres = objGpet.ref_native_image.YWorldLimits+objGpet.ref_native_image.PixelExtentInWorldY(1)/2 : objGpet.ref_native_image.PixelExtentInWorldY : objGpet.ref_native_image.YWorldLimits(2)-objGpet.ref_native_image.PixelExtentInWorldY/2;
            z_lowres = objGpet.ref_native_image.ZWorldLimits+objGpet.bed_position_mm+objGpet.ref_native_image.PixelExtentInWorldZ(1)/2 : objGpet.ref_native_image.PixelExtentInWorldZ : objGpet.ref_native_image.ZWorldLimits(2)+objGpet.bed_position_mm-objGpet.ref_native_image.PixelExtentInWorldZ/2;
            x_highres = objGpet.ref_image.XWorldLimits+objGpet.ref_image.PixelExtentInWorldX(1)/2 : objGpet.ref_image.PixelExtentInWorldX : objGpet.ref_image.XWorldLimits(2)-objGpet.ref_image.PixelExtentInWorldX/2;
            y_highres = objGpet.ref_image.YWorldLimits+objGpet.ref_image.PixelExtentInWorldY(1)/2 : objGpet.ref_image.PixelExtentInWorldY : objGpet.ref_image.YWorldLimits(2)-objGpet.ref_image.PixelExtentInWorldY/2;
            z_highres = objGpet.ref_image.ZWorldLimits+objGpet.ref_image.PixelExtentInWorldZ(1)/2 : objGpet.ref_image.PixelExtentInWorldZ : objGpet.ref_image.ZWorldLimits(2)-objGpet.ref_image.PixelExtentInWorldZ/2;
            [X_lowres,Y_lowres,Z_lowres] = meshgrid(x_lowres, y_lowres, z_lowres);
            [X_highres,Y_highres,Z_highres] = meshgrid(x_highres, y_highres, z_highres);
            % Masks:
            mask_highres = (X_highres.^2+Y_highres.^2)<(objGpet.image_size.matrixSize(1)*objGpet.image_size.voxelSize_mm(1)/2.5)^2;
            mask = (X_lowres.^2+Y_lowres.^2)<(objGpet.ref_native_image.ImageSize(1)*objGpet.ref_native_image.PixelExtentInWorldX(1)/2.5)^2;
            % PET ds:
            paramPET.scanner = objGpet.scanner;
            paramPET.method =  objGpet.method;
            paramPET.PSF.type = objGpet.PSF.type;
            paramPET.PSF.Width = objGpet.PSF.Width;
            paramPET.sinogram_size.span = objGpet.sinogram_size.span;
            paramPET.nSubsets = objGpet.nSubsets;
            paramPET.verbosity = 0;
            paramPET.tempPath = objGpet.tempPath;
            PET_lowres = classGpet(paramPET);

            % default parameters
            opt.PetOptimizationMethod = 'DePierro';%'OSL'
            opt.PetPriorType = 'Quadratic'; %'Bowsher' 'JointBurgEntropy'
            opt.PetRegularizationParameter = 1;
            opt.PetPreCompWeights = 1;
            opt.BowsherB = 70;
            opt.PriorMrImage =[];
            opt.TVsmoothingParameter = 0.1;
            opt.MrSigma = 0.1; % JBE
            opt.PetSigma  = 10; %JBE
            opt.display = 0;
            opt.save = 0;
            opt.outputPath = objGpet.tempPath;
            opt.saveInterval = 10;
            % check if the object already was initilized:
            if isempty(objGpet.Prior)
                opt.ImageSize = objGpet.image_size.matrixSize;
                opt.imCropFactor = [4,4,0];
                opt.sWindowSize = 5;
                opt.lWindowSize = 1;
                opt = getFiledsFromUsersOpt(opt,arg);
                objGpet.Prior = PriorsClass(opt);
            else
                % whos_prior = whos('objGpet.Prior');
                % if strcmp(whos_prior.class, 'PriorsClass')
                % Only initilize the parameters for this reconstruction
                opt = getFiledsFromUsersOpt(opt,arg);
            end

            if opt.display, figure; end

            if strcmpi(opt.PetPriorType,'Quadratic')
                if opt.PetPreCompWeights == 1
                    W0 = 1/objGpet.Prior.nS;
                else
                    W0 = opt.PetPreCompWeights./repmat(sum(opt.PetPreCompWeights,2),[1,objGpet.Prior.nS]); %for weighted quadratic
                end
            elseif strcmpi(opt.PetPriorType,'Bowsher')
                if opt.PetPreCompWeights == 1
                    fprintf('calculating Bowsher weighting coeffcients\n')
                    if isempty(opt.PriorMrImage), error('MR image should be provided in opt.PriorMrImage\n'); end
                    W0 = objGpet.Prior.W_Bowsher(opt.PriorMrImage,opt.BowsherB);
                    W0 = W0./repmat(sum(W0,2),[1,objGpet.Prior.nS]);
                else % better to precomputed Bowsher weighting coeffcients into opt.PetPreCompWeights
                    W0 = opt.PetPreCompWeights;
                end
            elseif strcmpi(opt.PetPriorType,'JointBurgEntropy')
                if opt.PetPreCompWeights == 1
                    fprintf('calculating MR-based Gaussian weighting coeffcients\n');
                    if isempty(opt.PriorMrImage), error('MR image should be provided in opt.PriorMrImage\n'); end
                    W0 = objGpet.Prior.W_JointEntropy(opt.PriorMrImage,opt.MrSigma);
                else % better to precomputed MR weighting coeffcients into opt.PetPreCompWeights
                    W0 = opt.PetPreCompWeights;
                end
            end

            % SENS IMAGE
            sensImage = PET_lowres.Sensitivity(anf);
            sensImg_highres = interp3(X_lowres, Y_lowres, Z_lowres, sensImage, X_highres, Y_highres, Z_highres, 'linear', 0);
            fprintf('Prior: %s, Method: %s\n',opt.PetPriorType,opt.PetOptimizationMethod);
            k = 1;
            for i = 1:nIter
                if opt.display,fprintf('Iteration: %d\n',i); end
                if strcmpi(opt.PetOptimizationMethod,'DePierro')
                    xn = Img;
                    % Projection:
                    low_res_image = interp3(X_highres, Y_highres, Z_highres, xn, X_lowres, Y_lowres, Z_lowres, 'linear', 0); %imresize(opmlem{end}, PET_highres.image_size.matrixSize, 'bicubic');
                    projected = anf.*PET_lowres.P(low_res_image) + additive;
                    % Backproject:
                    backprojected_image = PET_lowres.PT(anf.*objGpet.vecDivision(prompts, projected)).*mask;
                    % transpose of interpolation (high sample)
                    backprojected_image_highres = interp3(X_lowres, Y_lowres, Z_lowres, backprojected_image, X_highres, Y_highres, Z_highres, 'linear', 0).*mask_highres; %imresize(opmlem{end}, PET_highres.image_size.matrixSize, 'bicubic');
                    % Update image
                    x_em = xn.*objGpet.vecDivision(backprojected_image_highres, sensImg_highres);

                    if strcmpi(opt.PetPriorType,'JointBurgEntropy')
                        W0 = objGpet.Prior.W_JointEntropy(xn,opt.PetSigma).*W0;
                        W0 = W0./repmat(sum(W0,2),[1,objGpet.Prior.nS]);
                    end
                    W = objGpet.Prior.Wd.*W0;
                    wj = objGpet.Prior.UndoImCrop(reshape(sum(W,2),objGpet.Prior.CropedImageSize));
                    B = sensImg_highres - opt.PetRegularizationParameter / 2*objGpet.Prior.UndoImCrop(reshape(sum(W.*objGpet.Prior.GraphDivCrop(xn),2),objGpet.Prior.CropedImageSize));
                    Img = 2*x_em.*sensImg_highres./(B + sqrt(B.^2+4*opt.PetRegularizationParameter.*sensImg_highres.*x_em.*wj + 1e-5));
                    Img = max(0,Img);
                elseif strcmpi(opt.PetOptimizationMethod,'OSL')
                    xn = Img;
                    if strcmpi(opt.PetPriorType,'JointBurgEntropy')
                        W0 = objGpet.Prior.W_JointEntropy(xn,opt.PetSigma).*W0;
                        W0 = W0./repmat(sum(W0,2),[1,objGpet.Prior.nS]);
                        elseif strcmpi(opt.PetPriorType,'sTV')
                        imgGrad = objGpet.Prior.GraphGradCrop(xn);
                        Norm = repmat(sqrt(sum(abs(imgGrad).^2,2)+ opt.TVsmoothingParameter),[1,objGpet.Prior.nS]);
                        W0 = 1./(Norm+eps);
                    end
                    W = objGpet.Prior.Wd.*W0;
                    dP = -2* sum(W.*objGpet.Prior.GraphGradCrop(xn),2);
                    dP = opt.PetRegularizationParameter*objGpet.Prior.UndoImCrop(reshape(dP,objGpet.Prior.CropedImageSize));
                    Img = xn.*objGpet.vecDivision(objGpet.PT(objGpet.vecDivision(Prompts,objGpet.P(xn)+ additive)),SensImg + dP + 1e-5);
                    Img = max(0,Img);
                end
                if opt.display
                    if objGpet.image_size.matrixSize(3)==1
                        imshow(Img,[]);
                    else %if 3D recon, use opt.display =x, where x is a transvese slice
                        imshow(Img(:,:,opt.display),[]);
                    end
                end
                if opt.save
                    if rem(i-1,opt.saveInterval) == 0 % -1 to save the first iteration
                    image_ds{k} = Img;
                    interfilewrite(single(image_ds{k}), [outputPath 'map_ds_iter_' num2str(i)], [objGpet.ref_image.PixelExtentInWorldX objGpet.ref_image.PixelExtentInWorldY objGpet.ref_image.PixelExtentInWorldZ]); % i use i instead of i+1 because i=1 is the inital estimate
                    k = k + 1;
                    end
                else
                    image_ds = Img;
                end
            end
        end

        function image_ds = MAPEM2_DS(objGpet,prompts, anf, additive, Img, nIter,arg, outputPath)
            % First initilization of the parameters needed to downsample in
            % system matrix, that can be done in gpu or cpu:
            x_lowres = objGpet.ref_native_image.XWorldLimits+objGpet.ref_native_image.PixelExtentInWorldX(1)/2 : objGpet.ref_native_image.PixelExtentInWorldX : objGpet.ref_native_image.XWorldLimits(2)-objGpet.ref_native_image.PixelExtentInWorldX/2;
            y_lowres = objGpet.ref_native_image.YWorldLimits+objGpet.ref_native_image.PixelExtentInWorldY(1)/2 : objGpet.ref_native_image.PixelExtentInWorldY : objGpet.ref_native_image.YWorldLimits(2)-objGpet.ref_native_image.PixelExtentInWorldY/2;
            z_lowres = objGpet.ref_native_image.ZWorldLimits+objGpet.bed_position_mm+objGpet.ref_native_image.PixelExtentInWorldZ(1)/2 : objGpet.ref_native_image.PixelExtentInWorldZ : objGpet.ref_native_image.ZWorldLimits(2)+objGpet.bed_position_mm-objGpet.ref_native_image.PixelExtentInWorldZ/2;
            x_highres = objGpet.ref_image.XWorldLimits+objGpet.ref_image.PixelExtentInWorldX(1)/2 : objGpet.ref_image.PixelExtentInWorldX : objGpet.ref_image.XWorldLimits(2)-objGpet.ref_image.PixelExtentInWorldX/2;
            y_highres = objGpet.ref_image.YWorldLimits+objGpet.ref_image.PixelExtentInWorldY(1)/2 : objGpet.ref_image.PixelExtentInWorldY : objGpet.ref_image.YWorldLimits(2)-objGpet.ref_image.PixelExtentInWorldY/2;
            z_highres = objGpet.ref_image.ZWorldLimits+objGpet.ref_image.PixelExtentInWorldZ(1)/2 : objGpet.ref_image.PixelExtentInWorldZ : objGpet.ref_image.ZWorldLimits(2)-objGpet.ref_image.PixelExtentInWorldZ/2;
            [X_lowres,Y_lowres,Z_lowres] = meshgrid(single(x_lowres), single(y_lowres), single(z_lowres));
            [X_highres,Y_highres,Z_highres] = meshgrid(single(x_highres), single(y_highres), single(z_highres));
            % If running a gpu, create a copy in gpu:
%             if strcmp(objGpet.method, 'otf_siddon_gpu')
%                 gpuX_lowres = gpuArray(X_lowres); gpuY_lowres = gpuArray(Y_lowres); gpuZ_lowres = gpuArray(Z_lowres);
%                 gpuX_highres = gpuArray(X_highres); gpuY_highres = gpuArray(Y_highres); gpuZ_highres = gpuArray(Z_highres);
%             end
            % Masks:
            mask_highres = (X_highres.^2+Y_highres.^2)<(objGpet.image_size.matrixSize(1)*objGpet.image_size.voxelSize_mm(1)/2.5)^2;
            mask = (X_lowres.^2+Y_lowres.^2)<(objGpet.ref_native_image.ImageSize(1)*objGpet.ref_native_image.PixelExtentInWorldX(1)/2.5)^2;
            % PET ds:
            paramPET.scanner = objGpet.scanner;
            paramPET.method =  objGpet.method;
            paramPET.PSF.type = objGpet.PSF.type;
            paramPET.PSF.Width = objGpet.PSF.Width;
            paramPET.sinogram_size.span = objGpet.sinogram_size.span;
            paramPET.nSubsets = objGpet.nSubsets;
            paramPET.verbosity = 0;
            paramPET.tempPath = objGpet.tempPath;
            PET_lowres = classGpet(paramPET);

            % default parameters
            opt.PetOptimizationMethod = 'DePierro';%'OSL'
            opt.PetPriorType = 'Quadratic'; %'Bowsher' 'JointBurgEntropy'
            opt.PetSimilarityKernel = 'local';
            opt.PriorImplementation = 'matlab';
            opt.PetRegularizationParameter = 1;
            opt.PetPreCompWeights = 1;
            opt.BowsherB = 70;
            opt.PriorImage = [];
            opt.LangeDeltaParameter = [];
            opt.TVsmoothingParameter = 0.1;
            opt.MrSigma = 0.1; % JBE
            opt.PetSigma  = 10; %JBE
            opt.display = 0;
            opt.save = 0;
            opt.outputPath = objGpet.tempPath;
            opt.saveInterval = 10;

            opt.ImageSize = objGpet.image_size.matrixSize;
            opt.imCropFactor = [4,4,0];
            opt.sWindowSize = 5;
            opt.lWindowSize = 1;
            opt = getFiledsFromUsersOpt(opt,arg);
            objGpet.Prior = PriorsClass(opt);

            if opt.display, figure; end

            % SENS IMAGE
            sensImage = PET_lowres.Sensitivity(anf);
            %if strcmp(objGpet.method, 'otf_siddon_gpu')
            %    sensImg_highres = gather(interp3(gpuX_lowres, gpuY_lowres, gpuZ_lowres, gpuArray(sensImage), gpuX_highres, gpuY_highres, gpuZ_highres, 'linear', 0));
            %else
                sensImg_highres = interp3(X_lowres, Y_lowres, Z_lowres, sensImage, X_highres, Y_highres, Z_highres, 'linear', 0);
            %end
            fprintf('Prior: %s, Method: %s\n',objGpet.Prior.PriorType,opt.PetOptimizationMethod);
            k = 1;
            if opt.display, fprintf('Prior: %s, Method: %s\n',objGpet.Prior.PriorType,opt.PetOptimizationMethod); end

            for i = 1:nIter
                tic
                if opt.display,fprintf('Iteration: %d\n',i); end
                xn = Img;
                % Projection:
                %if strcmp(objGpet.method, 'otf_siddon_gpu')
                %    % if using gpu, run the interpolation also in gpu
                %    gpuXn = gpuArray(xn);
                %    low_res_image = gather(interp3(gpuX_highres, gpuY_highres, gpuZ_highres, gpuXn, gpuX_lowres, gpuY_lowres, gpuZ_lowres, 'linear', 0)); %imresize(opmlem{end}, PET_highres.image_size.matrixSize, 'bicubic');
                %    clear gpuXn
                %else
                    low_res_image = interp3(X_highres, Y_highres, Z_highres, xn, X_lowres, Y_lowres, Z_lowres, 'linear', 0); %imresize(opmlem{end}, PET_highres.image_size.matrixSize, 'bicubic');
                %end
                projected = anf.*PET_lowres.P(low_res_image) + additive;
                % Backproject:
                backprojected_image = PET_lowres.PT(anf.*objGpet.vecDivision(prompts, projected)).*mask;
                % transpose of interpolation (high sample)
                %if strcmp(objGpet.method, 'otf_siddon_gpu')
                %    gpuBackprojected_image = gpuArray(backprojected_image);
                %    backprojected_image_highres = gather(interp3(gpuX_lowres, gpuY_lowres, gpuZ_lowres, gpuBackprojected_image, gpuX_highres, gpuY_highres, gpuZ_highres, 'linear', 0)).*mask_highres; %imresize(opmlem{end}, PET_highres.image_size.matrixSize, 'bicubic');
                %    clear gpuBackprojected_image
                %else
                    backprojected_image_highres = gather(interp3(X_lowres, Y_lowres, Z_lowres, backprojected_image, X_highres, Y_highres, Z_highres, 'linear', 0)).*mask_highres; %imresize(opmlem{end}, PET_highres.image_size.matrixSize, 'bicubic');
                %end
                if strcmpi(opt.PetOptimizationMethod,'DePierro')
                    % Update image
                    x_em = xn.*objGpet.vecDivision(backprojected_image_highres, sensImg_highres);

                    if strcmpi(opt.PetPriorType,'JointBurgEntropy')
                        W0 = objGpet.Prior.W_JointEntropy(xn,opt.PetSigma).*W0;
                        W0 = W0./repmat(sum(W0,2),[1,objGpet.Prior.nS]);
                    end
                    W = objGpet.Prior.Wd.*W0;
                    wj = objGpet.Prior.UndoImCrop(reshape(sum(W,2),objGpet.Prior.CropedImageSize));
                    B = sensImg_highres - opt.PetRegularizationParameter / 2*objGpet.Prior.UndoImCrop(reshape(sum(W.*objGpet.Prior.GraphDivCrop(xn),2),objGpet.Prior.CropedImageSize));
                    Img = 2*x_em.*sensImg_highres./(B + sqrt(B.^2+4*opt.PetRegularizationParameter.*sensImg_highres.*x_em.*wj + 1e-10));
                    Img = max(0,Img);
                elseif strcmpi(opt.PetOptimizationMethod,'OSL')
                    dP = opt.PetRegularizationParameter*objGpet.Prior.dPrior(xn,opt);
                    Img = xn.*objGpet.vecDivision(backprojected_image_highres, sensImg_highres + dP + 1e-10);
                    Img = max(0,Img);
                end
                if opt.display
                    if objGpet.image_size.matrixSize(3)==1
                        imshow(Img,[]);
                    else %if 3D recon, use opt.display =x, where x is a transvese slice
                        imshow(Img(:,:,opt.display),[]);
                    end
                end
                if opt.save
                    if rem(i-1,opt.saveInterval) == 0 % -1 to save the first iteration
                        image_ds{k} = Img;
                        interfilewrite(single(image_ds{k}), [outputPath 'map_ds_iter_' num2str(i)], [objGpet.ref_image.PixelExtentInWorldX objGpet.ref_image.PixelExtentInWorldY objGpet.ref_image.PixelExtentInWorldZ]); % i use i instead of i+1 because i=1 is the inital estimate
                        k = k + 1;
                    end
                else
                    image_ds = Img;
                end
                timeIter = toc;
                if objGpet.verbosity>0,fprintf('Iteration time: %f sec\n',timeIter); end
            end
        end

        % Opmlem with downsample
        function image_mac = OPMLEM_MAC(objGpet, prompts, anf, additive, initialEstimate, numIterations, outputPath, saveInterval)
            if ~isdir(outputPath)
                mkdir(outputPath);
            end
            % The projector is only span 1:
            paramPET.scanner = objGpet.scanner;
            paramPET.method =  objGpet.method;
            paramPET.PSF.type = objGpet.PSF.type;
            paramPET.sinogram_size.span = 1;
            paramPET.nSubsets = objGpet.nSubsets;
            paramPET.verbosity = 1;
            paramPET.image_size = objGpet.image_size;
            PET_span1 = classGpet(paramPET);

            % Check anf size:
            if (size(anf) == objGpet.sinogram_size.matrixSize)
                % MAC-2: axial compression in the span of the input
                % sinogram.
                mac_method = 2;
            elseif (size(anf) == PET_span1.sinogram_size.matrixSize)
                % MAC-1
                mac_method = 1;
            else
                error('The anf matrix needs to be span 1 or of the same span as the classGpet object.');
            end
            % Compute the sensitivity:
            if mac_method == 1
                %[anf_spanN, structSizeSino3dSpanN] = objGpet.apply_axial_compression_from_span1(anf);
                sensImage = PET_span1.Sensitivity(anf);
            elseif mac_method == 2
                % Expand to span 1:
                [anf_span1, structSizeSino3dSpan1] = convertSinogramToSpan(anf, objGpet.get_sinogram_size_for_apirl(),  1);
                sensImage = PET_span1.Sensitivity(anf_span1);
            end
            image = initialEstimate;
            k = 1;
            for i = 1:numIterations
                % Project span 1 sinogram:
                projected = PET_span1.P(image);
                if mac_method == 1
                    projected = projected .* anf;
                    % comprss into span N:
                    [projected_spanN, structSizeSino3dSpanN] = objGpet.apply_axial_compression_from_span1(projected);
                elseif mac_method == 2
                    % First compress:
                    [projected_spanN, structSizeSino3dSpanN] = objGpet.apply_axial_compression_from_span1(projected);
                    % Normalize compression:
                    for s = 1 : size(projected_spanN,3)
                        projected_spanN(:,:,s) = projected_spanN(:,:,s)./objGpet.sinogram_size.numPlanesMashed(s);
                    end
                    % Now multiply by anf:
                    projected_spanN = projected_spanN .* anf;
                end
                % Add additive term:
                projected_spanN = projected_spanN + additive;

                % Expand correction sinogram to span 1:
                correctionSinogram = objGpet.vecDivision(prompts, projected_spanN);
                if mac_method == 1
                    % Firs go to span 1:
                    [correctionSinogram, structSizeSino3dSpan1] = convertSinogramToSpan(correctionSinogram, objGpet.get_sinogram_size_for_apirl(),  1);
                    % Then multiply by anf:
                    correctionSinogram = correctionSinogram .* anf;
                elseif mac_method == 2
                    % First multiply by anf:
                    correctionSinogram = correctionSinogram .* anf;
                    % Expand to span 1:
                    [correctionSinogram, structSizeSino3dSpan1] = convertSinogramToSpan(correctionSinogram, objGpet.get_sinogram_size_for_apirl(),  1);
                end
                % Backproject with span 1:
                backprojected_image = PET_span1.PT(correctionSinogram);
                % Update image
                image = image.*objGpet.vecDivision(backprojected_image, sensImage);
                image = max(0,image);
                if rem(i-1,saveInterval) == 0 % -1 to save the first iteration
                    image_mac{k} = image;
                    interfilewrite(single(image_mac{k}), [outputPath 'mlem_ds_iter_' num2str(i)], [objGpet.ref_image.PixelExtentInWorldX objGpet.ref_image.PixelExtentInWorldY objGpet.ref_image.PixelExtentInWorldZ]); % i use i instead of i+1 because i=1 is the inital estimate
                    k = k + 1;
                end
            end
        end

    end

end


