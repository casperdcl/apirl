function hdr = getInterfileHdrFromDicom(d)
% d: a Siemsne dicom header or a dicom filename

if(strcmp(computer(), 'GLNXA64'))
    os = 'linux';
    pathBar = '/';
elseif(strcmp(computer(), 'PCWIN') || strcmp(computer(), 'PCWIN64'))
    os = 'windows';
    pathBar = '\';
else
    disp('OS not compatible');
    return;
end

if ~isstruct(d), d = dicominfo(d); end

if isempty(d.Private_0029_1010)
    error('Interfile subheader was not found in SIEMENS Private tags')
end

[pathstr] = fileparts(d.Filename) ;
hdrFilename = [pathstr pathBar 'info.hdr'];
fid = fopen(hdrFilename,'w');
fprintf(fid,'%s\n',char(d.Private_0029_1010)');
fclose(fid);

hdr = getInfoFromInterfile(hdrFilename);
delete(hdrFilename);
end