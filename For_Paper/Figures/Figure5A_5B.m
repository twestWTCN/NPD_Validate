close all; clear all
%% Simulation 2
% Effects of Signal Mixing
fname = 'Figure4';
N = 3; % # of nodes
MO = 3;% model order
C = repmat(repmat(0.5,N,N).*eye(N),1,1,MO);
C(:,:,2) = -C(:,:,2);
C(2,1,2) = 0.5;
C(3,1,3) = 0.5;
NCV     = eye(N).*0.3;

% 5A - FC vs SNR
% % mvarconsim_npdver_F2b(C,NCV,25)
figure(1)
mvarconsim_npdver_F3b(C,NCV,25)
% 

% 5A - FC vs SNR Asymmetry
figure(2)
C = repmat(repmat(0.5,N,N).*eye(N),1,1,MO);
C(:,:,2) = -C(:,:,2);
C(2,1,2) = 0.35;
C(1,2,2) = 0.35;
mvarconsim_npdver_F3b_asysm(C,NCV,25)


