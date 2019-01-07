clear; close all
% Figure 3 - Asymmetric SNR
N = 3; % # of nodes
MO = 3;% model order

C = repmat(repmat(0.5,N,N).*eye(N),1,1,MO);
C(:,:,2) = -C(:,:,2);
C(2,1,2) = 0.35;
C(1,2,2) = 0.35;
NCV     = eye(N).*0.3;

wrapper_Fig3_AsymSNR(C,NCV,3)


figure(2)
for i = 1:N^2
subplot(N,N,i); ylim([0 0.6])
end
subplot(N,N,1); ylim([0 0.3])
subplot(N,N,5); ylim([0 0.3])
subplot(N,N,9); ylim([0 0.3])
legend({'Power','NPD','Granger'})
set(gcf,'Position',[1077         347         867         734])
