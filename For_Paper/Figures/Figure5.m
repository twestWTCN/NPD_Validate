% NPD_Validate_AddPaths()
close all
%% This script will reproduce Figure 5 - 
% Investigating the effects of signal-to-noise ratios (SNR),
% SNR asymmetries, and instantaneous linear mixing upon functional
% connectivity measures: coherence, non-parametric directionality
% (NPD), and non-parametric Granger causality (NPG).

fresh = 0; % (1) Run permutation tests; (0) Load precomputed tables.
permrun = 1; % (1) FFT Shuffle; (2) Phase Randomize
% Setup the Simulation
N = 2; % # of nodes
MO = 3;% model order
C = repmat(repmat(0.5,N,N).*eye(N),1,1,MO);
C(:,:,2) = -C(:,:,2);
C(2,1,2) = 0.35;
C(1,2,2) = 0.35;
NCV     = eye(N).*0.3;

% These are the main routines:
% 5A - FC vs SNR
figure(1)
wrapper_Fig5A_SymSNR(C,NCV,25,permrun,0)
ylim([0 1])

% 5B - FC vs SNR Asymmetry
figure(2)
wrapper_Fig5B_AsymSNR(C,NCV,25,permrun,0)

% Figure 5C- Continuous sweep of sig mixing
% Setup the simulation
N = 3; % # of nodes
MO = 3;% model order
C = repmat(repmat(0.5,N,N).*eye(N),1,1,MO);
C(:,:,2) = -C(:,:,2);
C(2,1,2) = 0.5;
C(3,1,3) = 0.5;
NCV     = eye(N).*0.3;

figure(3)
wrapper_Fig5C_contSigmMix(C,NCV,25,permrun,1)
figure(3)
ylim([0 1.2])