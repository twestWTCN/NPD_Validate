function [] = wrapper_Fig2_SymSNR(C,NCV,NC)

NCvec = [0 0.75 3];
NCvec = sqrt(NCvec); % rescale to yield variance
nc_col_sc = [1 0.9 0.8];

cmap = linspecer(4);
lsstyles = {'-','-.',':'};

% Run the MVAR Simulation (fieldtrip implementation uses BSMART)
rng(1236)
Nsig = size(C,1);
cfg             = [];
cfg.ntrials     = 1;
cfg.triallength = 250;
cfg.fsample     = 200;
cfg.nsignal     = Nsig;
cfg.method      = 'ar';
cfg.params = C;
cfg.noisecov = NCV;
X              = ft_connectivitysimulation(cfg);
segOrd = 8; % 2^n length of segment used for FFT

for ncov = 1:NC
    if ncov == 1
        perm = 1;
        permtype = 2;
    else
        perm = 0;
        permtype = 0;
    end
    
    % Plotting parameters   
    linestyle = lsstyles{ncov};
    cmapn = cmap.*nc_col_sc(ncov);
    plotfig =1;
    
    %% Adjust Symmetrical SNR
    data = X;
    randproc = randn(size(data.trial{1}));
    for i = 1:size(randproc,1)
        s = X.trial{1}(i,:);
        s = (s-mean(s))./std(s);
        n = ((NCvec(ncov)*1).*randproc(i,:));
        y = s + n;
        snr = var(s)/var(n);
        snrbank(ncov,i) = snr;
        data.trial{1}(i,:) = y;
    end
    
    %% Compute Power with Fieldtrip (used for NPD and NPG computations)
    figure(2)
    datalength = (2^segOrd)./cfg.fsample;
    freq = computeSpectra(data,[0 0 0],Nsig,plotfig,linestyle,-1,datalength);
        
    %% Compute and Plot NPD
    figure(2)
    [Hz lags npdspctrm npdspctrmZ npdspctrmW nscohspctrm npdcrcv] = ft_computeNPD(freq,cfg.fsample,1,segOrd,perm,permtype);
    coh.freq= Hz; coh.cohspctrm = nscohspctrm{1}; coh.ci = nscohspctrm{2};
    
    % Plot Coherence (Neurospec version)
    plotNSCoherence(coh,cmapn(1,:),Nsig,plotfig,linestyle,perm)
    % Plot zero-lag
    %     plotNPD_zero(Hz,npdspctrm,data,cmap(1,:),plotfig,linestyle)
    % Plot NPD
    plotNPD(Hz,npdspctrm,data,cmapn(3,:),plotfig,linestyle,perm)
    % Plot Partialized NPD
    %     plotNPD(Hz,npdspctrmZ,data,cmapn(4,:),plotfig,linestyle)
    %     plotNPD(Hz,npdspctrmW,data,cmap(4,:),plotfig,linestyle)
    
    %% Compute and plot GRANGER
    figure(2)
    [Hz granger grangerft] = computeGranger(freq,0,perm,permtype)
    plotNPD(grangerft.freq,granger,freq,cmapn(2,:),1,linestyle,perm)
    %% NPD CORR
    figure(3)
    plotNPDXCorr(lags,npdcrcv,data,[0 0 0],plotfig,linestyle)
end
SNRDB = 10*log10(snrbank(:,1))
a= 1;
