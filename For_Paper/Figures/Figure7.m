close all; clear all
addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\TWtools')
%% Simulation 7- Incomplete Signals for Conditioning
% Route X->Z->Y
fname = 'Figure7';
N = 3; % # of nodes
MO = 3;% model order
C = repmat(repmat(0.5,N,N).*eye(N),1,1,MO);
C(:,:,2) = -C(:,:,2);
C(:,:,3) = zeros(N);
C(2,1,3) = 0.3;
C(3,2,3) = 0.3;
% C(3,3,1:3) = [-0.5 0.5 0.5];
NCV     = eye(N).*0.3;

figure(1)
wrapper_Fig7_IncompleteCond(C,NCV,25)
ylim([0 0.5])
set(gcf,'Position',[680   678   483   393])
% Route X->Y->Z
C = repmat(repmat(0.5,N,N).*eye(N),1,1,MO);
C(:,:,2) = -C(:,:,2);
C(:,:,3) = zeros(N);
C(3,1,3) = 0.3;
C(2,3,3) = 0.3;
figure(2)
wrapper_Fig7_IncompleteCond(C,NCV,25)
ylim([0 0.5])
set(gcf,'Position',[680   678   483   393])

% Route X->Y->Z->X
C = repmat(repmat(0.5,N,N).*eye(N),1,1,MO);
C(:,:,2) = -C(:,:,2);
C(:,:,3) = zeros(N);
C(3,1,3) = 0.3;
C(2,3,3) = 0.3;
C(1,2,3) = 0.3;
figure(3)
wrapper_Fig7_IncompleteCond(C,NCV,25)
ylim([0 0.5])
set(gcf,'Position',[680   678   483   393])



%%
% addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\linspecer')
% addpath('C:\spm12'); spm eeg; close all
% addpath(genpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\Neurospec'))
% addpath(genpath('C:\Users\twest\Documents\Work\PhD\NPD_Validate'))