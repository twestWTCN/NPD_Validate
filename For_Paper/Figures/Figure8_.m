% NPD_Validate_AddPaths()
% addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\TWtools')
close all

%% This script will reproduce Figure 8 - The effects of incomplete
% signal observation upon estimation of directed functional connectivity:
% non-parametric Granger causality (NPG); non-parametric directionality
% (NPD); and NPD conditioned on reference signal Z (NPD(Z)). 

% Serial Routing
% Route X->Z->Y
N = 3; % # of nodes
MO = 3;% model order
C = repmat(repmat(0.5,N,N).*eye(N),1,1,MO);
C(:,:,2) = -C(:,:,2);
C(:,:,3) = zeros(N);
C(2,1,3) = 0.3;
C(3,2,3) = 0.3;
% C(3,3,1:3) = [-0.5 0.5 0.5];
NCV     = eye(N).*0.3;

figure(1)
wrapper_Fig8_IncompleteCond(C,NCV,25)
ylim([0 0.5])
set(gcf,'Position',[680   678   483   393])

% Feedforward Routing
% Route X->Y->Z
C = repmat(repmat(0.5,N,N).*eye(N),1,1,MO);
C(:,:,2) = -C(:,:,2);
C(:,:,3) = zeros(N);
C(3,1,3) = 0.3;
C(2,3,3) = 0.3;

figure(2)
wrapper_Fig8_IncompleteCond(C,NCV,25)
ylim([0 0.5])
set(gcf,'Position',[680   678   483   393])

% Reccurent Routing
% Route X->Y->Z->X
C = repmat(repmat(0.5,N,N).*eye(N),1,1,MO);
C(:,:,2) = -C(:,:,2);
C(:,:,3) = zeros(N);
C(3,1,3) = 0.3;
C(2,3,3) = 0.3;
C(1,2,3) = 0.3;

figure(3)
wrapper_Fig8_IncompleteCond(C,NCV,25)
ylim([0 0.5])
set(gcf,'Position',[680   678   483   393])

