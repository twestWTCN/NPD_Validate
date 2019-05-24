function [] = wrapper_CI_bootstrap(C,NCV,ntrials)
Nsig = size(C,1);

for trial = 1:ntrials
    cfg             = [];
    cfg.ntrials     = 1;
    cfg.triallength = 500;
    cfg.fsample     = 200;
    cfg.nsignal     = Nsig;
    cfg.method      = 'ar';
    cfg.params = C;
    cfg.noisecov = NCV;
    data              = ft_connectivitysimulation(cfg);    % Weighted mixture of signals
    
    %% Compute Spectra
    freq = computeSpectra(data,[],Nsig,0,[]);
    %% NPD
    [Hz lags npdspctrm npdspctrmZ npdspctrmW nscohspctrm npdcrcv] = computeNPD(data,1);
    coh.freq= Hz; coh.cohspctrm = nscohspctrm;  
    %% GRANGER
    [Hz granger grangerft] = computeGranger(freq,[],Nsig,0,[],0)
    
    
    
end


% cfg = [];
% cfg.viewmode = 'butterfly';  % you can also specify 'butterfly'
% ft_databrowser(cfg, data);
