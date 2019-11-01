close all; clear
%% Simulation 5 (SNR and Asym SNR Sweep and Sig Mix Sweep)
% Effects of Signal Mixing
N = 2; % # of nodes
MO = 3;% model order
C = repmat(repmat(0.5,N,N).*eye(N),1,1,MO);
C(:,:,2) = -C(:,:,2);
C(2,1,2) = 0.35;
C(1,2,2) = 0.35;
NCV     = eye(N).*0.3;

% 5A - FC vs SNR
% figure(1)
% wrapper_Fig5A_SymSNR(C,NCV,25)

% 5A - FC vs SNR Asymmetry
figure(2)
wrapper_Fig5B_AsymSNR(C,NCV,25)

% Figure 5C- Continuous sweep of sig mixing
N = 3; % # of nodes
MO = 3;% model order
C = repmat(repmat(0.5,N,N).*eye(N),1,1,MO);
C(:,:,2) = -C(:,:,2);
C(2,1,2) = 0.5;
C(3,1,3) = 0.5;
NCV     = eye(N).*0.3;
figure(3)
wrapper_Fig5C_contSigmMix(C,NCV,25)
