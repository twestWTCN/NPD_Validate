close all; clear;
% 
%% Figure 4 - Effects of Signal Mixing
fname = 'Figure4';
N = 3; % # of nodes
MO = 3;% model order
C = repmat(repmat(0.5,N,N).*eye(N),1,1,MO);
C(:,:,2) = -C(:,:,2);
C(2,1,2) = 0.5;
C(3,1,3) = 0.5;
NCV     = eye(N).*0.3;
fstord = 2;
sndord = [2 3];
wrapper_Fig4_SigMix(C,NCV,3,fstord,sndord,fname)

figure(3)
set(gcf,'Position',[1077         347         867         734])
legend('\lambda = 0','\lambda = 0.3','\lambda = 0.6')

figure(2)
subplot(N,N,1); ylim([0 0.2])
subplot(N,N,5); ylim([0 0.2])
subplot(N,N,9); ylim([0 0.2])
legend({'Power','Zero-NPD','NPD','Granger'})
set(gcf,'Position',[1077         347         867         734])


% %% Simulation 3 - Increase the delay
% fname = 'Figure3';
% N = 3; % # of nodes
% MO = 3;% model order
% C1 = repmat(repmat(0.5,N,N).*eye(N),1,1,MO);
% C1(:,:,2) = -C1(:,:,2);
% 
% C2 = repmat(repmat(0,N,N).*eye(N),1,1,MO);
% C2(:,:,2) = -C2(:,:,2);
% C2(2,1,2) = 1;
% C2(3,1,3) = 1;
% 
% C = cat(3,C1,C2);
% NCV     = eye(N).*0.3;
% fstord = 2;
% sndord = [2 3];
% mvarconsim_npdver_F2(C,NCV,3,fstord,sndord,fname)
% 
% figure(3)
% set(gcf,'Position',[1077         347         867         734])
% legend('\lambda = 0','\lambda = 0.5','\lambda = 0.8')
% 
% figure(2)
% subplot(N,N,1); ylim([0 0.2])
% subplot(N,N,5); ylim([0 0.2])
% subplot(N,N,9); ylim([0 0.2])
% legend({'Power','Zero-NPD','NPD','Granger'})
% set(gcf,'Position',[1077         347         867         734])
