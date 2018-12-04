close all; clear all
%% Simulation 2
% Effects of Symetric SNR
fname = 'Figure3';
N = 3; % # of nodes
MO = 3;% model order
C = repmat(repmat(0.5,N,N).*eye(N),1,1,MO);
C(:,:,2) = -C(:,:,2);
C(2,1,2) = 0.5;
C(3,1,3) = 0.5;
NCV     = eye(N).*0.3;
fstord = 2;
sndord = [2 3];
mvarconsim_npdver_F3(C,NCV,3)
%
figure(3)
set(gcf,'Position',[1077         347         867         734])
legend('\lambda = 0','\lambda = 0.3','\lambda = 0.6')
%
figure(2)
subplot(N,N,1); ylim([0 0.2])
subplot(N,N,5); ylim([0 0.2])
subplot(N,N,9); ylim([0 0.2])
legend({'Power','Coherence','NPD','Granger'})
set(gcf,'Position',[1077         347         867         734])

