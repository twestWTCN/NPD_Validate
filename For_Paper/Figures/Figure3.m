% NPD_Validate_AddPaths()
close all
%% This script will reproduce Figure 3 - 
% Analysis of the effects of unequal signal-to-noise ratios,
% measured as a difference of the SNRs between X and Y 
%(??SNR?_XY) upon symmetrical directed functional connectivity (dFC).

% Setup the Simulation
N = 3; % # of nodes
MO = 3;% model order
C = repmat(repmat(0.5,N,N).*eye(N),1,1,MO);
C(:,:,2) = -C(:,:,2);
C(2,1,2) = 0.35;
C(1,2,2) = 0.35;
NCV     = eye(N).*0.3;

% This is the main routine:
wrapper_Fig3_AsymSNR(C,NCV,3)

% Now Configure plots:
figure(2)
for i = 1:N^2
subplot(N,N,i); ylim([-0.5 0.5])
end
subplot(N,N,1); ylim([-0.25 0.2])
subplot(N,N,5); ylim([0 0.2])
subplot(N,N,9); ylim([0 0.2])
legend({'Power','NPD','Granger'})
set(gcf,'Position',[1077         347         867         734])
