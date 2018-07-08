
function mdhrFilename = writeInterfileMhdr(interfileName)

if(strcmp(computer(), 'GLNXA64'))
    pathBar = '/';
elseif(strcmp(computer(), 'PCWIN') || strcmp(computer(), 'PCWIN64'))
    pathBar = '\';
else
    disp('OS not compatible');
    return;
end

if ~exist(interfileName,'file')
    error('could not find the interfile header: %s\n',interfileName);
end

info = getInfoFromInterfile(interfileName);

if isfield(info,'DataType')
    dataType = info.DataType; % sinogram | image
elseif isfield(info,'PetDataType')
    if strcmpi(info.PetDataType,'emission')
        dataType = 'sinogram';
    elseif strcmpi(info.PetDataType,'image')
        dataType = 'image';
    else
        error('Unknown datatype');
    end
else
    error('could not find datatype');
end

[Path,name,ext] = fileparts(interfileName);
if ~isempty(ext)
    [~,name2] = fileparts(name);
end
mdhrFilename = [Path pathBar name2 '.mhdr'];
fid = fopen(mdhrFilename,'w');

fprintf(fid,'!INTERFILE:=\n');
fprintf(fid,'%%comment:=SMS-MI sinogram common attributes\n');
fprintf(fid,'!originating system:=2008\n');
fprintf(fid,'%%SMS-MI header name space:=sinogram main header\n');
fprintf(fid,'%%SMS-MI version number:=3.1\n');
fprintf(fid,'\n');
fprintf(fid,'!GENERAL DATA:=\n');
fprintf(fid,'data description:=%s\n',dataType);%sinogram | image
fprintf(fid,'exam type:=wholebody\n');
fprintf(fid,'%%study date (yyyy:mm:dd):=%s\n',info.StudyDateYyyyMmDd);
fprintf(fid,'%%study time (hh:mm:ss GMT+00:00):=%s\n',info.StudyTimeHhMmSsGmt0000);
fprintf(fid,'%%type of detector motion:=step and shoot\n'); %SinogramType
fprintf(fid,'\n');
fprintf(fid,'%%DATA MATRIX DESCRIPTION:=\n');
fprintf(fid,'number of time frames:=1\n');
fprintf(fid,'%%number of horizontal bed offsets:=1\n');
fprintf(fid,'number of time windows:=1\n');

if isfield(info,'NumberOfScanDataTypes') && info.NumberOfScanDataTypes==2
    fprintf(fid,'%%number of emission data types:=%d\n',info.NumberOfScanDataTypes);
    fprintf(fid,'%%emission data type description [1]:=%s\n',info.ScanDataTypeDescription1);
    fprintf(fid,'%%emission data type description [2]:=%s\n',info.ScanDataTypeDescription2);
    fprintf(fid,'%%number of transmission data types:=0\n');
    
end

fprintf(fid,'%%scan direction:=in\n');
fprintf(fid,'\n');
fprintf(fid,'%%DATA SET DESCRIPTION:=\n');
fprintf(fid,'!total number of data sets:=%d\n',info.TotalNumberOfDataSets);
fprintf(fid,'%%data set [1]:={');

if strcmpi(dataType,'image')
    fprintf(fid,'30,');
    fprintf(fid,'%s', [name,ext]);
    fprintf(fid,',UNKNOWN}\n');
else
    fprintf(fid,'0,');
    fprintf(fid,'%s', [name,ext]);
    fprintf(fid,',%s}\n',name);
end

fclose(fid);
