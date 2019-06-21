function NPD_Validate_AddPaths()
if strmatch(getenv('computername'),'FREE')
spmpath = 'C:\spm12';
% spmpath = 'C:\Users\timot\Documents\GitHub\spm12';
preFix_MLAD = 'C:\Users\twest\Documents\Work\MATLAB ADDONS\';
% preFix_MLAD = 'C:\Users\timot\Documents\Work\MATLAB ADDONS\';
preFix_Git = 'C:\Users\twest\Documents\Work\GitHub\';
% preFix_Git = 'C:\Users\timot\Documents\GitHub\';

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

