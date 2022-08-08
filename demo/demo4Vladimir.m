NPD_Validate_AddPaths

% You need to add the neurospec paths and fieldtrip %
% Setup Connectivity Simulation
N = 3; % # of nodes
MO = 3;% model order
C = repmat(repmat(0.5,N,N).*eye(N),1,1,MO);
C(:,:,2) = -C(:,:,2);
C(2,1,2) = 0.5;
C(3,1,3) = 0.5;
NCV     = eye(N).*0.3;

[eigVals,flag] = checkMVARStability(C)

% Simular AR model
cfg             = [];
cfg.ntrials     = 1;
cfg.triallength = 250;
cfg.fsample     = 200;
cfg.nsignal     = N;
cfg.method      = 'ar';
cfg.params = C;
cfg.noisecov = NCV;
data              = ft_connectivitysimulation(cfg);

bstrp = 0; %option for bootstrapping
condind = 1; % index for conditioning
NPDord = 8; % 2^N (order for NPD)

[Hz lags npdspctrm npdspctrmZ npdspctrmW nscohspctrm npdcrcv npdcrcvZ npdcrcvW] = computeNPD(data,condind,NPDord,bstrp);
% npdspctrm is the main (unconditioned) NPD
% it is organized into cells {i=1:2,j=1:3} where i=1 is the metric and i=2
% is the confidence interval. j=1 is instantaneous; j=2 is foward; j=3
% backward. Each cell has an N x N x Hz matrix corresponding to the
% connectivity matrix. 
% npdspctrmZ: estimated conditioned on condind
% npdspctrmW: estimated conditioned on condind + 1
% nscohspctrm: Neurospec coherence estimate
% npdcrcv#: these are the time domain prewhitened cumulants.

% Create the coherence structure
coh.freq= Hz; coh.cohspctrm = nscohspctrm{1}; coh.ci = nscohspctrm{2};

% Plot the results
cmap = hsv(3);
plotNSCoherence(coh,cmap(1,:),N,1,'-',bstrp)
%     plotNPD_zero(Hz,npdspctrm,data,cmap(1,:),plotfig,linestyle)
% NPD
plotNPD(Hz,npdspctrm,data,cmap(3,:),1,'--',bstrp)
legend('NPD','Coherence')