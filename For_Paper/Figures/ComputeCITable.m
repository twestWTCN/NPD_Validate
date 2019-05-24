% NPD_Validate_AddPaths()
% % % mvarconsim_npdver(C,NCV,fstord,sndord,fname)
close all

%% Compute CI table
ntrials = 1000;
N = 3; % # of nodes
MO = 3;% model order
C = repmat(repmat(0.5,N,N).*eye(N),1,1,MO);
C(:,:,2) = -C(:,:,2);
NCV     = eye(N).*0.3;
wrapper_CI_bootstrap(C,NCV,ntrials)

