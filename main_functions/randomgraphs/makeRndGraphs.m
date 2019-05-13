function [CMat,NCV] = makeRndGraphs(nc,nreps)

% Create a connectivity matrix of 1-4 connections
% possible connections
clist = combvec(1:3,1:3);
clist(:,diff(clist)==0) = []; % remove self
% Make template connection matrix
A = repmat(repmat(0.5,3,3).*eye(3),1,1,3);
A(:,:,2) = -A(:,:,1);
% noise covarianc matrix
NCV     = eye(3).*0.3;


for i = 1:nc
    for n = 1:nreps
        csel = randsample(1:size(clist,2),i);
        osel = randi(3,1,size(csel,1),i);
        CMat{i,n} = placeCon(A,clist(:,csel),osel,0.3);
        
%         cfg             = [];
%         cfg.ntrials     = 1;
%         cfg.triallength = 250;
%         cfg.fsample     = 200;
%         cfg.nsignal     = 3;
%         cfg.method      = 'ar';
%         cfg.params = CMat{i,n};
%         cfg.noisecov = NCV;
%         data{i,n} = ft_connectivitysimulation(cfg);
    end
end

% save('C:\Users\Tim\Documents\Work\GIT\NPD_Validate\main_functions\benchmark\data\data','data')
% save('C:\Users\Tim\Documents\Work\GIT\NPD_Validate\main_functions\benchmark\data\CMat','CMat')
