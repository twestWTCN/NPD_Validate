close all; clear all
%% Simulation 5 (SNR and Asym SNR Sweep and Sig Mix Sweep)
% Effects of Signal Mixing
N = 3; % # of nodes
MO = 3;% model order
C = repmat(repmat(0.5,N,N).*eye(N),1,1,MO);
C(:,:,2) = -C(:,:,2);
C(2,1,2) = 0.5;
C(3,1,3) = 0.5;
NCV     = eye(N).*0.3;

% 5A - FC vs SNR
fname = 'Figure5A';
figure(1)
wrapper_Fig5A_SymSNR(C,NCV,25)

% 5A - FC vs SNR Asymmetry
fname = 'Figure5B';
figure(2)
C = repmat(repmat(0.5,N,N).*eye(N),1,1,MO);
C(:,:,2) = -C(:,:,2);
C(2,1,2) = 0.35;
C(1,2,2) = 0.35;
wrapper_Fig5B_AsymSNR(C,NCV,25)

% Figure 5C- Continuous sweep of sig mixing
fname = 'Figure5C';
figure(3)
N = 3; % # of nodes
MO = 3;% model order
C = repmat(repmat(0.5,N,N).*eye(N),1,1,MO);
C(:,:,2) = -C(:,:,2);
C(2,1,2) = 0.5;
C(3,1,3) = 0.5;
NCV     = eye(N).*0.3;
wrapper_Fig5C_contSigmMix(C,NCV,25)
