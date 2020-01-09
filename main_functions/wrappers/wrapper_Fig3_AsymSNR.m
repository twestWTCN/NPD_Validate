function [] = wrapper_Fig3_AsymSNR(C,NCV,NC)
% Asymetric SNRs
% Setup std for the channels:
NCvec = [0.01  1      10
         0.01  0.01   0.01
         0.01  0.01   0.01];
NCvec = sqrt(NCvec); % convert to variance

nc_col_sc = [1 0.9 0.8]; % rescale colors
cmap = linspecer(4);
lsstyles = {'-','-.',':'};% linestyles

% Run the MVAR Simulation (fieldtrip implementation uses BSMART)
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
    
    %% Compute Signal Mixing
    data = X;
    randproc = randn(size(data.trial{1}));
    for i = 1:size(randproc,1)
        s = data.trial{1}(i,:);
        s = (s-mean(s))./std(s);
        n = ((NCvec(i,ncov)*1).*randproc(i,:));
        y = s + n;
        snr = var(s)/var(n);
        snrbank(ncov,i) = snr;
        data.trial{1}(i,:) = y;
    end
    
    %% Plot Example Trace
    figure(1)
    plot(data.time{1},data.trial{1}+[0 2.5 5]');
    xlim([100 105])
    xlabel('Time'); ylabel('Amplitude'); grid on
    legend({'X1','X2','X3'})
    
    %% Compute Power with Fieldtrip (used for NPD and NPG computations)
    figure(2)
    datalength = (2^segOrd)./cfg.fsample;
    freq = computeSpectra(data,[0 0 0],Nsig,plotfig,linestyle,-1,datalength);
    
    %% Coherence
    %     figure(2)
    %     computeCoherence(freq,cmap(1,:),Nsig,plotfig,linestyle)
    
    %% Compute and Plot NPD
    figure(2)
    [Hz lags npdspctrm npdspctrmZ npdspctrmW nscohspctrm npdcrcv] = ft_computeNPD(freq,cfg.fsample,1,segOrd,perm,permtype);
    coh.freq= Hz; coh.cohspctrm = nscohspctrm{1}; coh.ci = nscohspctrm{2};

    % Plot NPD differential 
    plotDiffNPD(Hz,npdspctrm,data,cmapn(3,:),plotfig,linestyle,perm)
    
    %% Compute and plot NPG
    figure(2)
    [Hz granger grangerft] = computeGranger(freq,0,perm,permtype)
    % Plot NPG differential 
    plotDiffNPD(grangerft.freq,granger,freq,cmapn(2,:),1,linestyle,perm)
end
 % Review the estimated SNRs
SNRDB = 10*log10(snrbank(:,1));
SNRBASE = 10*log10(snrbank(:,2));
SRAT = SNRDB-SNRBASE;

a= 1;
