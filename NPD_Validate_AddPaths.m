function NPD_Validate_AddPaths()
if strmatch(getenv('computername'),'FREE')
    
    spmpath = 'C:\Users\timot\Documents\GitHub\spm12';
    preFix_MLAD = 'C:\Users\timot\Documents\Work\MATLAB ADDONS\';
    preFix_Git = 'C:\Users\timot\Documents\GitHub\';
elseif strmatch(getenv('computername'),'SFLAP-2')
    spmpath = 'C:\Users\Tim\Documents\spm12';
    preFix_MLAD = 'C:\Users\Tim\Documents\MATLAB_ADDONS\';
    preFix_Git = 'C:\Users\Tim\Documents\Work\GIT\';
elseif strmatch(getenv('computername'),'DESKTOP-94CEG1L-2')
    spmpath = 'C:\Users\timot\Documents\GitHub\spm12';
    preFix_MLAD = 'C:\Users\timot\Documents\Work\MATLAB ADDONS\';
    preFix_Git = 'C:\Users\timot\Documents\GitHub\';
    
end
addpath([preFix_MLAD 'linspecer'])
% addpath('C:\spm12'); spm eeg; close all
addpath(spmpath); spm eeg; close all
addpath([preFix_MLAD 'Neurospec\neurospec21'])
addpath([preFix_MLAD 'Neurospec\PartialDirected'])
addpath([preFix_MLAD 'TWtools'])
addpath([preFix_MLAD 'linspecer'])
addpath([preFix_Git 'NPD_Validate\For_Paper\Figures'])
addpath(genpath([preFix_Git 'NPD_Validate\main_functions']))

addpath([preFix_MLAD 'DrosteEffect-BrewerMap-221b913'])

