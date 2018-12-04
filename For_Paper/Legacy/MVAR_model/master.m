% addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\linspecer')
% addpath('C:\spm12'); spm eeg; close all
% addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\Neurospec\neurospec21')
% addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\Neurospec\PartialDirected')
% mvarconsim_npdver(C,NCV,fstord,sndord,fname)
close all
% Simulation 1
fname = '3node1con';
N = 3; % # of nodes
MO = 2;% model order
C = repmat(repmat(0.5,N,N).*eye(N),1,1,MO);
C(:,:,2) = -C(:,:,2);
C(3,1,1) = 0.5;
NCV     = eye(N).*0.3;
fstord = 2;
sndord = [2 3];
mvarconsim_npdver(C,NCV,fstord,sndord,fname)

close all
% Simulation 2
fname = '3node2con';
N = 3; % # of nodes
MO = 2;% model order
C = repmat(repmat(0.5,N,N).*eye(N),1,1,MO);
C(:,:,2) = -C(:,:,2);
C(1,2,1) = 0.5;
C(3,1,1) = 0.5;
NCV     = eye(N).*0.3;
fstord = 1;
sndord = [1 3];
mvarconsim_npdver(C,NCV,fstord,sndord,fname)

close all
% Simulation 4
fname = '4node_part_ring';
N = 4; % # of nodes
MO = 2;% model order
C = repmat(repmat(0.5,N,N).*eye(N),1,1,MO);
C(:,:,2) = -C(:,:,2);
C(2,1,1) = 0.5;
C(3,2,1) = 0.5;
C(4,3,1) = 0.5;

NCV     = eye(N).*0.3;
fstord = 2;
sndord = [2 3];
mvarconsim_npdver(C,NCV,fstord,sndord,fname)