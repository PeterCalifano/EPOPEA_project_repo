clc

% kernelpool = '..\EPOPEA_project_repo\MA\EPOPEA_metakernel.tm';
kernelpool = fullfile('..\EPOPEA_metakernel.tm'); % for Windows
% kernelpool = fullfile('..\EPOPEA_metakernelUNIX.tm'); % for UNIX/MAC

% kernelpool = '..\EPOPEA_metakernel.tm'

cspice_furnsh(kernelpool);

et0 = cspice_str2et('2022-03-10 11:00 UTC');
cspice_spkezr('EARTH', et0, 'ECLIPJ2000', 'NONE', 'SSB')

cspice_kclear();