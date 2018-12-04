close all; clear all
addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\TWtools')
%% Simulation 5- Recurrent Connectivity
% Compute Partialization
fname = 'Figure5';
N = 4; % # of nodes
MO = 4;% model order
C = repmat(repmat(0,N,N).*eye(N),1,1,MO);
% Autospectra
C(1,1,1:4) = [0     0   0.35 -0.35];% "High Beta"
C(2,2,1:4) = [0.6 -0.6  0    0]; % "Low beta"
C(3,3,1:4) = [0     0   0    0];   % Null"
C(4,4,1:4) = [0     0   0    0];   % Null"

% C(3,3,1:3) = [-0.5 0.5 0.5];

C(2,1,3) = 0.4;
C(3,1,2) = 0.4;
C(3,2,2) = 0.4;
C(4,3,3) = 0.4;
C(1,4,3) = 0.4;

NCV     = eye(N).*0.2;
NCV(2,2) = 0.1;
NCV(3,3) = 0.1;
NCV(4,4) = 0.1;

% mvarconsim_npdver_F2b(C,NCV,25)
mvarconsim_npdver_F5(C,NCV,1)

% figure(1)
figure(2)
for i = 1:N^2
    subplot(N,N,i); ylim([0 0.65])
end
subplot(N,N,1); ylim([0 0.35])
subplot(N,N,6); ylim([0 0.35])
subplot(N,N,11); ylim([0 0.35])
subplot(N,N,16); ylim([0 0.35])

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