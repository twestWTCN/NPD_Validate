function [CMat,NCV] = makeRndGraphs(nc,nreps,sz)

% Create a connectivity matrix of 1-4 connections
% possible connections
clist = combvec(1:sz,1:sz);
clist(:,diff(clist)==0) = []; % remove self
% Make template connection matrix
A = repmat(repmat(0.5,sz,sz).*eye(sz),1,1,sz);
A(:,:,2) = -A(:,:,1);
% noise covarianc matrix
NCV     = eye(sz).*0.3;

cstrength = [0.3 0.3 0.3 0.3];
for i = 1:nc
    for n = 1:nreps
        csel = randsample(1:size(clist,2),i);
        osel = randi(sz,1,size(csel,1),i);
        CMat{i,n} = placeCon(A,clist(:,csel),osel,cstrength(i));
        
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
