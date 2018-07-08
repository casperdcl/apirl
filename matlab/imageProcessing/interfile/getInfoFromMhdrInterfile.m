%  *********************************************************************
%  Reconstruction Framework for Siemens Biograph mMR.  
%  Autor: Martín Belzunce. Kings College London.
%  Fecha de Creación: 10/02/2015
%  *********************************************************************
%  function [info, structSizeSino]  = getInfoFromInterfile(filename)
% 
%  This function reads the header of an interfile sinograms and gets some
%  useful information from it. Such as the filename of the raw data, the
%  singles rates per bucket,..
%  It is based on the interfileinfo of matlab.
% Optionally it also returns an structure with the size of the sinogram.


function info = getInfoFromMhdrInterfile(filename)

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

info = [];
% check header file extension
if (isempty(filename) || ~ischar(filename))
    error('getInfoFromSiemensIntf:filenameNotChar', ...
          'Filename must be a character array.')
end
[fpath,name,ext] = fileparts(filename);
if isempty(ext)
    filename = [filename '.hdr'];
end
% Como los interfile yo los manejo siempre con el h33 e i33 siempre en el
% mismo path, pero que se puedan leer y crear desde otro path. Cuando el
% nombre del interfile a leer tiene un path además del nombre, dicho path
% también debo agregarselo al "data file name" que figura en el h33:
relativePath = '';
barras = strfind(filename, pathBar);
if ~isempty(barras)
    relativePath = filename(1 : barras(end));
end

% open file for parsing
fid = fopen(filename);
if fid == -1
    error('Images:interfileinfo:invalidFilename', ...
          'Can not find header file %s.', filename);
end

% initialize variables
bad_chars = '!()[]/-_%:+*';    % Added %,:, siemens interfile sues it.
dates = {'DateOfKeys' 'ProgramDate' 'PatientDob' 'StudyDate'};
times = {'StudyTime' 'ImageStartTime'};
found_header = 0;
found_end = 0;
line_num = 0;

% parse through the file
while (true)
    line_txt = fgetl(fid);
    % stop if no more lines
    if (line_txt == -1)
        break;
    end
    
    % Strip out comments.  Interfile v3.3 spec, paragraph I.4.H: "Key-value
    % pairs may have comments appended to them by preceding the comment with
    % a semicolon <;>.  Conversion programs can ignore all characters
    % including and following a semicolon  <;> to the end of line code.
    % Where no key is stated, for example when an ASCII line starts with a
    % semicolon, the whole line is treated as comment.
    line_txt = regexprep(line_txt, ';.*$', '');
        
    if (sum(isspace(line_txt))+sum(line_txt==0)) == length(line_txt)
        % Line is empty, skip to the next.
        continue;

    else
        line_num = line_num+1;
        % find index of separator and issue warning if not found
        sep_ind = strfind(line_txt, ':=');
        if (isempty(sep_ind))
            fclose(fid);
            % if no separator on first non-empty line, then not in INTERFILE format
            if isempty(info)
                error('Images:interfileinfo:invalidFile', ...
                      '%s is not a valid INTERFILE file.', filename);
                
            % if not on first non-empty line, then invalid expression
            else
                error('Images:interfileinfo:noSeparator', ...
                      'Invalid expression in line %s of %s.', num2str(line_num), filename);
            end
        
        else
            field_str_ind = 1;
            value_str_ind = sep_ind+2;
            field = '';
            
            % parse string to extract field
            while (true)
                [str, count, errmsg, nextindex] = sscanf(line_txt(field_str_ind:sep_ind-1), '%s', 1);
                % check for duplicate header
                if (strcmp(str, '!INTERFILE'))
                    if (found_header == 1)
                        fclose(fid);
                        error('Images:interfileinfo:duplicateHeader', ...
                              'Duplicate Interfile header in line %d of %s.', line_num, filename);
                        
                    else
                        found_header = 1;
                    end
                end
                
                % break if no field in rest of string
                if (count == 0)
                    break;
                end
                
                % concatenate strings to form field
                if (strcmp(str, 'ID'))
                    field = [field str];
                    
                else
                    str = lower(str);
                    i = 1;
                    
                    % remove illegal characters
                    while (i <= length(str))
                        k = strfind(bad_chars, str(i));
                        if (~isempty(k))
                            if (k >= 6)
                                str = [str(1:i-1) upper(str(i+1)) str(i+2:length(str))];

                            else
                                str = [str(1:i-1) str(i+1:length(str))];
                            end
                            % Added:
                            % If I take out one chrarachter I need to check
                            % that position again, if not when two non valid
                            % characters are found the second one is not detected:
                            i = i-1;
                        end

                        i = i+1;
                    end
                    if ~isempty(str) % Just in case a word is formed by all illegal characters
                        field = [field upper(str(1)) str(2:length(str))];
                    end
                end
                
                field_str_ind = field_str_ind+nextindex-1;
            end
            
            % remove extra spaces from beginning of value string
            for i = value_str_ind:length(line_txt)
                if (~isspace(line_txt(i)))
                    break;
                end
            end
            
            value = strcat(line_txt(i:length(line_txt)), '');
            if (strcmp(field, 'VersionOfKeys'))
                if (~strcmp(value, '3.3'))
                    fclose(fid);
                    err_id = 'Images:interfileinfo:unsupportedVersion';
                    err_msg = 'Unsupported version of keys detected.';
                    error(err_id, err_msg);
                end
            end
            
            if isempty(value)
                value = '';
            end
                
            x = str2double(value);
            if (~isnan(x)) % && (isempty(strfind(dates, field))) && (isempty(strfind(times, field))))
                value = x;
            end
            
            % close file if end-of-file marker encountered
            if (strcmp(field, 'EndOfInterfile'))
                found_end = 1;
                break;
                
            else
                % check for header
                if (found_header == 0)
                    fclose(fid);
                    err_id = 'Images:interfileinfo:noHeader';
                    err_msg = 'Interfile header not found.';
                    error(err_id, err_msg);

                % store field and value
                elseif (~strcmp(field, 'Interfile'))
                    if (isfield(info, field))
                        if (ischar(info.(field)))
                            info.(field) = {info.(field) value};
                            
                        elseif (iscell(info.(field)))
                            info.(field){length(info.(field))+1} = value;
                            
                        else
                            info.(field) = [info.(field) value];
                        end
                        
                    else
                        info.(field) = value;
                    end
                end
            end
        end
    end
end

% In the siemens interfile there is no enfofinterfile line.
% % check for end of file marker
% if (found_end == 0)
%     fclose(fid);
%     err_id = 'Images:interfileinfo:unexpectedEOF';
%     err_msg = 'Unexpected end of file.';
%     error(err_id, err_msg);
% end

% close file
fclose(fid);
% Add the relative path to the binary filename, if there is no path in the
% interfile:

% extract the hdr's name from
% %data set [1]:={1000000000,PET_ACQ_68_20150610155347_umap_hardware_00.v.hdr,UNKNOWN}
Hdr = info.DataSet1;
i = strfind(Hdr,',');
Hdr = Hdr(i(1)+1:i(2)-1);
% remove any lealing or trailing white space
Hdr = strtrim(Hdr);
% add path if it doesn't include it.
[Path,name,ext] = fileparts(Hdr);
if ~isdir(Path) || isempty(Path)
    Hdr = [fpath pathBar name ext];
end

info.NameOfDataFile = Hdr;
info.NameOfMhdrFile = filename;



