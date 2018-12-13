% addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\linspecer')
% addpath('C:\spm12'); spm eeg; close all
% addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\Neurospec\neurospec21')
% addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\Neurospec\PartialDirected')
% addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\TWtools')

%% Simulation 1
fname = 'Figure1';
N = 3; % # of nodes
MO = 3;% model order
C = repmat(repmat(0.5,N,N).*eye(N),1,1,MO);
C(:,:,2) = -C(:,:,2);
C(2,1,2) = 0.5;
C(3,1,3) = 0.5;
NCV     = eye(N).*0.3;
mvarconsim_npdver_F1(C,NCV,1)

figure(1)
set(gca, 'Color', 'None'); box off
set(gcf,'Position',[704   678   536   321])

figure(2)
subplot(N,N,1); ylim([0 0.2])
subplot(N,N,5); ylim([0 0.2])
subplot(N,N,9); ylim([0 0.2])
legend({'Power','Coherence','NPD','NPDx2','Granger'})
set(gcf,'Position',[1077         347         800         700])

