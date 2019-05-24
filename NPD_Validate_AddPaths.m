function NPD_Validate_AddPaths()
if strmatch(getenv('computername'),'FREE')
    addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\linspecer')
% addpath('C:\spm12'); spm eeg; close all
addpath('C:\spm12'); spm eeg; close all
addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\Neurospec\neurospec21')
addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\Neurospec\PartialDirected')
addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\TWtools')
addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\linspecer')
addpath('C:\Users\twest\Documents\Work\GitHub\NPD_Validate\For_Paper\Figures')
addpath('C:\Users\twest\Documents\Work\GitHub\NPD_Validate\main_functions\wrappers')
addpath('C:\Users\twest\Documents\Work\GitHub\NPD_Validate\main_functions\metrics')
else
    
    
    
addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\linspecer')
% addpath('C:\spm12'); spm eeg; close all
addpath('C:\Users\Tim\Documents\spm12'); spm eeg; close all
addpath('C:\Users\Tim\Documents\MATLAB_ADDONS\Neurospec\neurospec21')
addpath('C:\Users\Tim\Documents\MATLAB_ADDONS\Neurospec\PartialDirected')
addpath('C:\Users\Tim\Documents\MATLAB_ADDONS\TWtools')
addpath('C:\Users\Tim\Documents\MATLAB_ADDONS\linspecer')
addpath('C:\Users\Tim\Documents\Work\GIT\NPD_Validate\For_Paper\Figures')
addpath('C:\Users\Tim\Documents\Work\GIT\NPD_Validate\main_functions\wrappers')
addpath('C:\Users\Tim\Documents\Work\GIT\NPD_Validate\main_functions\metrics')
end