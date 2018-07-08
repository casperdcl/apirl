
function out = ParseInterfileDataFile(directory,report)

% Parse the interfile data folder of SIEMSE mMR scanner
% output: number of list-mode data, sinograms, normalization files and
% their interfile headers plus the number of unclassified files


if nargin==1, report = 0; end

nCompressedSinogramFiles = 0;
nUncompressedSinogramFiles = 0;
nListModeFiles = 0;
nNormFiles = 0;
nHumanUmaps = 0;
nHardwareUmaps = 0;
nHardUmapMhdrs = 0;
nHumanUmapMhdrs = 0;
nCompressedSinogramMhdrs = 0;
nUncompressedSinogramMhdrs = 0;
%             Unclassified = 0;
[CompressedSinogramHdrs,CompressedSinogramMhdrs,UncompressedSinogramHdrs,UncompressedSinogramMhdrs,...
    ListModeHdrs, NormFileHdrs, HumanUmapHdrs, HumanUmapMhdrs,HardwareUmapHdrs,HardwareUmapMhdrs] = deal(struct('hdr',[]));
if(strcmp(computer(), 'GLNXA64'))
    bar = '/';
else
    bar = '\';
end
% if directory already ends with a bar
% if strcmp(directory(end),bar) 
%     bar = '';
% end

directory = [directory bar];
listing = dir(directory);
% Dir return the string sorted but not sorts in a way that you get '1',
% '11', '2', etc, then the dynamic data is messed up, I correct that by
% resorting by length, the data is already ordered alphabetically for the
% same length:
lengthString(1) = 0; lengthString(2) = 0; % For '.' and '..'
for i = 3:length(listing)
    lengthString(i) = numel(listing(i).name);
end
[values, indices] = sort(lengthString);
listing = listing(indices);
for i = 3:length(listing)
    [~, name, ext] = fileparts(listing(i).name); % This name doesnt include the path.
    if ~isempty(name) && ~isempty(ext)
        if strcmpi(ext,'.l')
            nListModeFiles = nListModeFiles + 1;
            hdrFilename = [directory name '.hdr'];
            ListModeHdrs = getHeaderFiles(hdrFilename,ListModeHdrs,nListModeFiles);
        elseif strcmpi(ext,'.s')
            hdrFilename = [directory name '.s.hdr'];
            hdr = getInfoFromInterfile(hdrFilename);
            % exclude the non-emission sinograms, i.e. norm, scatter sinograms
            % with NumberOfBytesPerPixel==4
            if strcmpi(hdr.PetDataType,'emission') && strcmpi(hdr.DataFormat,'sinogram') && hdr.NumberOfBytesPerPixel==2
                if strcmpi(hdr.Compression,'on')
                    nCompressedSinogramFiles = nCompressedSinogramFiles + 1;
                    CompressedSinogramHdrs = getHeaderFiles(hdrFilename,CompressedSinogramHdrs,nCompressedSinogramFiles);
                else
                    nUncompressedSinogramFiles = nUncompressedSinogramFiles + 1;
                    UncompressedSinogramHdrs = getHeaderFiles(hdrFilename,UncompressedSinogramHdrs,nUncompressedSinogramFiles);
                end
            end
        elseif strcmpi(ext,'.n')
            hdrFilename = [directory name '.n.hdr'];
            nNormFiles = nNormFiles + 1;
            NormFileHdrs = getHeaderFiles(hdrFilename,NormFileHdrs,nNormFiles);
        elseif strcmpi(ext,'.v')
            hdrFilename = [directory name '.v.hdr'];
            if ~isempty(strfind(name,'umap'))
                if ~isempty(strfind(name,'hardware'))
                    nHardwareUmaps = nHardwareUmaps + 1;
                    HardwareUmapHdrs = getHeaderFiles(hdrFilename,HardwareUmapHdrs,nHardwareUmaps);
                else
                    nHumanUmaps = nHumanUmaps + 1;
                    HumanUmapHdrs = getHeaderFiles(hdrFilename,HumanUmapHdrs,nHumanUmaps);
                end
            end
        elseif strcmpi(ext,'.mhdr')
            MhdrFilename = [directory listing(i).name];
            Mhdr = getInfoFromMhdrInterfile(MhdrFilename);
            if strcmpi(Mhdr.DataDescription,'image')
                % check if the file is for umap(_human) or umap_hardware
                if ~isempty(strfind(Mhdr.NameOfDataFile,'umap'))
                    if ~isempty(strfind(Mhdr.NameOfDataFile,'hardware'))
                       nHardUmapMhdrs = nHardUmapMhdrs + 1; 
                       HardwareUmapMhdrs(nHardUmapMhdrs).hdr = Mhdr;
                       HardwareUmapMhdrs(nHardUmapMhdrs).hdrFilename = MhdrFilename;
                    else
                        nHumanUmapMhdrs = nHumanUmapMhdrs + 1;
                        HumanUmapMhdrs(nHumanUmapMhdrs).hdr = Mhdr;
                        HumanUmapMhdrs(nHumanUmapMhdrs).hdrFilename = MhdrFilename;
                    end
                end
            elseif strcmpi(Mhdr.DataDescription,'sinogram')
                % exclude those of scatters, norms
                % NumberOfEmissionDataTypes ==1
                if Mhdr.NumberOfEmissionDataTypes==2                 
                    % check theat the .s is available, sometimes is only
                    % the mhdr when the list mode is used:
                    % check if the Mhdr.NameOfDataFile referes to compressed or uncompressed sinograms
                    if exist(Mhdr.NameOfDataFile)
                        hdr = getInfoFromInterfile(Mhdr.NameOfDataFile);
                        if strcmpi(hdr.Compression,'on')
                            nCompressedSinogramMhdrs = nCompressedSinogramMhdrs + 1;
                            CompressedSinogramMhdrs(nCompressedSinogramMhdrs).hdr = Mhdr;
                        else
                            nUncompressedSinogramMhdrs = nUncompressedSinogramMhdrs + 1;
                            UncompressedSinogramMhdrs(nUncompressedSinogramMhdrs).hdr = Mhdr;
                        end
                    end
                end
            end
        end
    else
        %         Unclassified = Unclassified + 1;
    end
end

out.CompressedSinogramHdrs = CompressedSinogramHdrs;
out.CompressedSinogramMhdrs = CompressedSinogramMhdrs;
out.nCompressedSinogramFiles = nCompressedSinogramFiles;

out.UncompressedSinogramHdrs = UncompressedSinogramHdrs;
out.UncompressedSinogramMhdrs = UncompressedSinogramMhdrs;
out.nUncompressedSinogramFiles = nUncompressedSinogramFiles;

out.ListModeHdrs = ListModeHdrs;
out.nListModeFiles = nListModeFiles;

out.NormFileHdrs = NormFileHdrs;
out.nNormFiles = nNormFiles;

out.HumanUmapHdrs = HumanUmapHdrs;
out.HumanUmapMhdrs = HumanUmapMhdrs;
out.nHumanUmaps = nHumanUmaps;

out.HardwareUmapHdrs = HardwareUmapHdrs;
out.HardwareUmapMhdrs = HardwareUmapMhdrs;
out.nHardwareUmaps = nHardwareUmaps;




if report
    
    fprintf('InterFiles: %d Compressed Sinograms, %d UnCompressed Sinograms,%d ListModeFiles, %d NormFiles, %d HumanUmaps, %d HardwareUmaps\n',...
        nCompressedSinogramFiles, nUncompressedSinogramFiles,nListModeFiles, nNormFiles,...
        nHumanUmaps, nHardwareUmaps);
 
end
end

function dataStruct = getHeaderFiles(hdrFilename,dataStruct,n)
if exist(hdrFilename,'file')
    dataStruct(n).hdr = getInfoFromInterfile(hdrFilename);
    dataStruct(n).hdrFilename = hdrFilename;
else
    error('could not find the header file: %s\n',hdrFilename);
end
end