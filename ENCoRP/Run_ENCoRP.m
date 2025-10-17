%
% Housekeeping
%
clear
close all
clc

%
% Check that adequate python version is available and loaded
%
pyver = pyenv;

try
    ver = split(pyver.Version,'.');
    ver = [str2double(ver{1}), str2double(ver{2})];
catch
    ver = [-1, -1];
end
if (pyver.Status == "Loaded" && ver(1) < 3) || (pyver.Status == "Loaded" && ver(1) == 3 && ver(2) < 8)
    added = false;
    disp('Your Python version is insufficient to run the backend of the interface.')
    disp('Please choose a version 3.8 or newer.')
    disp('To do so, restart MATLAB, then type pyenv(Version="3.9") into the command window.')
elseif (pyver.Status == "NotLoaded" && ((ver(1) < 3) || (ver(1) == 3 && ver(2) < 8))) || isnan(str2double(pyver.Version))
    added = true;
    suitableversion = false;
    for v1 = 3:3
        for v2 = 8:14
            try 
                pyenv(Version=[num2str(v1),'.',num2str(v2)])
                suitableversion = true;
                break
            end
        end
        if suitableversion
            break
        end
    end
    if ~suitableversion
        disp('Could not locate any suitable Python versions (checked 3.8-3.14).')
        added = false;
    end
    if added
        disp('A suitable version of Python was found and added to the environment.')
    end
else
    added = true;
end

%
% Check MATLAB versioning
%
matver = version;

start = false;

if contains(matver,'R2022b')
    start = true;
else
    disp('The GUI may fail at some juncture due to a difference in MATLAB version.')
    disp('It is recommended to use MATLAB 2022b.')
end

%
% Clean workspace from previous sessions
%
if ~(start && added)
    response = questdlg('ENCoRP failed to detect sufficient MATLAB and Python versions. Proceed anyway?','Question','Yes','No','No');
    if strcmp(response,'Yes')
        start = true;
        added = true;
    end
end

if start && added
    
    if exist('Data/VIOut.mat','file')
        delete('Data/VIOut.mat')
    end
    if exist('Data/ModelHandoff.mat','file')
        delete('Data/ModelHandoff.mat')
    end
    if exist('Radiative/GUIData.mat','file')
        delete('Radiative/GUIData.mat')
    end

    %
    % Clear cached python files
    %
    run('VI/reloadPy.m')

    %
    % Prompt path restoration to avoid conflicting file caches
    %
    Answer = questdlg('If you have more than one copy of this interface, restoring the default path may avoid conflicting file issues. Restore the default path?','Question','Yes','No','Yes');
    if strcmp(Answer,'Yes')
        disp('Restoring default path...')
        restoredefaultpath
        disp('Done.')
        pause(1)
    end
    clear
    clc

    %
    % Start ENCoRP
    %
    addpath('Interface/')
    run('Interface/ENCoRP')

end
