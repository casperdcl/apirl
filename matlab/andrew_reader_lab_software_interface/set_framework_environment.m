function set_framework_environment(apirlPath, apirlBinaryPath)

if nargin == 0
    % APIRL PATH
    apirlPath = '../../';
    apirlBinaryPath = apirlPath;
elseif nargin == 1
    apirlBinaryPath = apirlPath;
end

%% CONFIGURE PATHS
setenv('PATH', [getenv('PATH') pathsep ...
    apirlBinaryPath filesep 'build' filesep 'bin']);
setenv('LD_LIBRARY_PATH', [getenv('LD_LIBRARY_PATH') pathsep ...
    apirlBinaryPath filesep 'build' filesep 'bin']);
% Check what OS I am running on:
if(strcmp(computer(), 'GLNXA64'))
    os = 'linux';
elseif(strcmp(computer(), 'PCWIN') || strcmp(computer(), 'PCWIN64'))
    os = 'windows';
else
    disp('OS not compatible');
    return;
end
% Matlab path:
addpath(genpath([apirlPath filesep 'matlab']));
