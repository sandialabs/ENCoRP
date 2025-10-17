
%
% function reloadPy()
%
% Summary: 
% MATLAB establishes Python classes the first time .py files are executed. 
% These classes are not inherently overwritten when changes are made to the 
% underlying .py files and the code re-run. In other words, MATLAB will use 
% the first version of each Python file regardless of later changes. This 
% function clears the MATLAB classes, permitting modified Python code to be
% executed without restarting MATLAB. This function must be executed after
% changes to the .py files have been saved and prior to MATLAB calling the 
% .py files.
%
% Inputs:
%   None
%
% Outputs:
%   None
%
% Example Execution:
%   >> reloadPy
%

function reloadPy()

    %
    % Disable warning
    %
    warning('off','MATLAB:ClassInstanceExists')

    %
    % Clear classes
    %
    clear classes

    %
    % Reload AssembleData.py
    %
    mod = py.importlib.import_module('AssembleData');
    py.importlib.reload(mod);
    
    %
    % Reload compVInodes_Functions.py
    %
    mod = py.importlib.import_module('compVInodes_Functions');
    py.importlib.reload(mod);

    %
    % Reload compVInodes_Nwire.py
    %
    mod = py.importlib.import_module('compVInodes_Nwire');
    py.importlib.reload(mod);

    % Note that Compute_VI.py is not reloaded as this is called directly.

end
