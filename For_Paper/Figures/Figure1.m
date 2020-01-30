% NPD_Validate_AddPaths()
close all
%% This script will reproduce Figure 1 - 
% Three-node simulation of MVAR model to compare functional connectivity measures.
permrun = 0; % (1) FFT Shuffle; (2) Phase Randomize
% Setup the Simulation
fname = 'Figure1';
N = 3; % # of nodes
MO = 3;% model order
% Define Connectivity and Autonomous Oscillations
C = repmat(repmat(0.5,N,N).*eye(N),1,1,MO);
C(:,:,2) = -C(:,:,2);
C(2,1,2) = 0.5;
C(3,1,3) = 0.5;
NCV     = eye(N).*0.3;
fstord = 2;
sndord = [2 3];

% This is the main routine:
wrapper_Fig1_CommonInput(C,NCV,1,permrun)

% Now Configure plots:
figure(1)
set(gca, 'Color', 'None'); box off
set(gcf,'Position',[704   678   536   321])

figure(2)
subplot(N,N,1); ylim([0 0.175])
subplot(N,N,5); ylim([0 0.175])
subplot(N,N,9); ylim([0 0.175])
legend({'Power','Coherence','NPD','NPDx2','Granger'})
set(gcf,'Position',[1077         347         800         700])

