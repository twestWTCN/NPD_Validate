% addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\linspecer')
% addpath('C:\spm12'); spm eeg; close all
% addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\Neurospec\neurospec21')
% addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\Neurospec\PartialDirected')


clear; close all
load('C:\Users\twest\Documents\Work\PhD\NPD_Validate\For_Paper\NMM_validation\data\ABC_commondrive_NMM\ABC_NMM.mat')
% Set up data
for i = 1:size(xsimsave,2)
    data(i).trial{1} = xsimsave{i}{1};
    data(i).label = {'X1','X2','X3'};
    data(i).time{1} = linspace(0,size(data(i).trial{1},2).*R.IntP.dt,size(data(i).trial{1},2));
    data(i).fsample = 1/R.IntP.dt;
end

NMM_NPD_connectivity_effects(data,4,1,2,'NMM_Connectivity')


N = 3;
for ncov = 1:4
figure(2+(10*ncov))
for i = 1:N^2
subplot(N,N,i); ylim([0 1.2]);xlim([2 98])
end
subplot(N,N,1); ylim([0 0.5]);xlim([2 98])
subplot(N,N,5); ylim([0 0.5]);xlim([2 98])
subplot(N,N,9); ylim([0 0.5]);xlim([2 98])
legend({'Power','NPD-Zero','NPD','NPDx2'})
set(gcf,'Position',[1077         347         867         734])

figure(3+(10*ncov))
for i = 1:N^2
subplot(N,N,i); ylim([-0.1 0.1]);xlim([-75 75])
end
legend({'Power','NPD-Zero','NPD','NPDx2'})
set(gcf,'Position',[1077         347         867         734])
end