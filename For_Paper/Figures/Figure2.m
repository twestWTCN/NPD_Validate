% NPD_Validate_AddPaths()
close all
%% This script will reproduce Figure 2 - 
% Analysis of the effects of signal-to-noise ratio (SNR) upon estimators of
% directed functional connectivity.

permrun = 0; % (1) FFT Shuffle; (2) Phase Randomize
% Setup the Simulation
N = 3; % # of nodes
MO = 3;% model order 
C = repmat(repmat(0.5,N,N).*eye(N),1,1,MO);
C(:,:,2) = -C(:,:,2);
C(2,1,2) = 0.5;
C(3,1,3) = 0.5;
NCV     = eye(N).*0.3;
fstord = 2;
sndord = [2 3];

% This is the main routine:
wrapper_Fig2_SymSNR(C,NCV,3,permrun)

% Now Configure plots:
figure(2)
for i = 1:N^2
subplot(N,N,i); ylim([0 1])
end
subplot(N,N,1); ylim([0 0.25])
subplot(N,N,5); ylim([0 0.25])
subplot(N,N,9); ylim([0 0.25])
legend({'Power','Coherence','NPD','Granger'})
set(gcf,'Position',[1077         347         867         734])

figure(3)
set(gcf,'Position',[1077         347         867         734])
legend('\lambda = 0','\lambda = 0.3','\lambda = 0.6')


