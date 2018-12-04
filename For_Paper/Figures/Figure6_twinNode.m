close all; clear all
addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\TWtools')
%% Simulation 6- Incomplete Signals for Conditioning
% Compute Partialization
fname = 'Figure6';
N = 4; % # of nodes
MO = 3;% model order
C = repmat(repmat(0.5,N,N).*eye(N),1,1,MO);
C(:,:,2) = -C(:,:,2);
C(:,:,3) = zeros(N);

C(2,1,3) = 0.5;
C(3,1,3) = 0.5;
C(4,2,3) = 0.5;
C(4,3,3) = 0.5;
% C(3,3,1:3) = [-0.5 0.5 0.5];

NCV     = eye(N).*0.3;

% NCV = nearestSPD(NCV);
% mvarconsim_npdver_F2b(C,NCV,25)
mvarconsim_npdver_F6(C,NCV,1)
N = N-1;
% figure(1)
figure(2)
for i = 1:N^2
    subplot(N,N,i); ylim([0 1.2])
end
subplot(N,N,1); ylim([0 0.2])
subplot(N,N,5); ylim([0 0.2])
subplot(N,N,9); ylim([0 0.2])

legend({'Power','NPD','NPDx2','Granger'})
set(gcf,'Position',[471         137        1066         979])

figure(3)
legend({'NPD','NPDx2'})
set(gcf,'Position',[471         137        1066         979])

%%
% addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\linspecer')
% addpath('C:\spm12'); spm eeg; close all
% addpath(genpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\Neurospec'))
% addpath(genpath('C:\Users\twest\Documents\Work\PhD\NPD_Validate'))