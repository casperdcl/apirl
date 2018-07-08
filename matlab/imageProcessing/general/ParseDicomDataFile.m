

function out = ParseDicomDataFile(directory,report)

% Parse the DICOM data folder of SIEMSE mMR scanner
% output: number of list-mode data, sinograms, normalization files and
% their interfile headers plus the number of unclassified files

if nargin==1, report = 0; end
d = dir(directory);

[Sinograms,listModes,listModeLarges,NormFiles] = deal(struct('hdr',[]));
nSinograms = 0;
nListModes = 0;
nListModeLarges = 0;
nNormFiles = 0;
Unclassified = 0;
if(strcmp(computer(), 'GLNXA64'))
    bar = '/';
else
    bar = '\';
end
% if directory already ends with a bar
if strcmp(directory(end),bar) 
    bar = '';
end
    
for i= 3:length(d)
    
    filename = [directory bar d(i).name];
    [~,name,ext]=fileparts(filename);
    
    if ~isempty(name) && ~isempty(ext) && isdicom(filename) && ~strcmp(ext, '.s')  && ~strcmp(ext, '.l') % I sclude known extension of binary files because the first bytes can be consdered as headers of adicom
        dicomHdr = dicominfo(filename);
        if ~isempty(dicomHdr.Private_0029_1008)
            dataType = dicomHdr.Private_0029_1008;
        elseif ~isempty(dicomHdr.Private_0029_1108)
            dataType = dicomHdr.Private_0029_1108;
        else
            error('Could not find data type from SIEMENS private tags');
        end

        switch lower(dataType)
            case lower('MRPETSINO')
                nSinograms = nSinograms + 1;
                Sinograms(nSinograms).hdr = getInterfileHdrFromDicom(dicomHdr);
            case lower('MRPETLM')
                nListModes = nListModes + 1;
                listModes(nListModes).hdr = getInterfileHdrFromDicom(dicomHdr);
            case lower('MRPETLM_LARGE')
                nListModeLarges = nListModeLarges + 1;
                listModeLarges(nListModeLarges).hdr = getInterfileHdrFromDicom(dicomHdr);
                % And add the binary filename:
                bfFileName = [directory bar name '.bf'];
                if exist(bfFileName)
                    listModeLarges(nListModeLarges).bf = bfFileName;
                else
                    error('List mode large file format header found, but the binary file (.bf) is not present in the folder.' );
                end
            case lower('MRPETNORM')
                nNormFiles = nNormFiles + 1;
                NormFiles(nNormFiles).hdr = getInterfileHdrFromDicom(dicomHdr);
            otherwise
                Unclassified = Unclassified + 1;
        end
    else
        Unclassified = Unclassified + 1;
    end
end
    
out.SinogramHdrs = Sinograms;
out.nSinograms = nSinograms;

out.listModeHdrs = listModes;
out.nListModes = nListModes;

out.listModeLargeHdrs = listModeLarges;
out.nListModeLarges = nListModeLarges;

out.NormFileHdrs = NormFiles;
out.nNormFiles = nNormFiles;

out.Unclassified = Unclassified;
    
if report
    fprintf('DICOM Files: %d Sinograms, %d ListModeFiles, %d ListModeLargeFiles, %d NormFiles\n',nSinograms, nListModes, nListModeLarges, nNormFiles);
end