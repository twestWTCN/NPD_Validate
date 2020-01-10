function [] = wrapper_Fig1_CommonInput(C,NCV,NC)

nc_col_sc = 1; % rescale colors for each

cmap = linspecer(4); % lineplot colors
lsstyles = {'-','-.',':'}; % linestyles

% Run the MVAR Simulation (fieldtrip implementation uses BSMART)
rng(8532)
Nsig = size(C,1);
cfg             = [];
cfg.ntrials     = 1;
cfg.triallength = 250; % 250s long
cfg.fsample     = 200; % 200 Hz
cfg.nsignal     = Nsig;
cfg.method      = 'ar';
cfg.params = C;
cfg.noisecov = NCV;
data = ft_connectivitysimulation(cfg);

segOrd = 8; % 2^n length of segment used for FFT
for ncov = 1:NC
    if ncov == 1
        perm = 1; % for first run- do permutation testing
        permtype = 2; % Permtype I is FFT segment shuffling; Permtype II is phase randomization
    else
        perm = 0;
        permtype = 0;
    end
    
    % Plotting parameters
    cmapn = cmap.*nc_col_sc(ncov);
    linestyle = lsstyles{ncov};
    plotfig =1;

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
    plotNPD(Hz,npdspctrmZ,data,cmapn(4,:),plotfig,linestyle,perm)
    %     plotNPD(Hz,npdspctrmW,data,cmap(4,:),plotfig,linestyle)
    
    %% Compute and plot GRANGER
    figure(2)
    computeGranger(freq,1,perm,permtype)
    plotNPD(grangerft.freq,granger,freq,cmapn(2,:),1,linestyle,perm)
    %% NPD CORR
    %     figure(3)
    %     plotNPDXCorr(lags,npdcrcv,data,[0 0 0],plotfig,linestyle)
    
    a =1;
end
% save('C:\Users\Tim\Documents\Work\GIT\NPD_Validate\precomp_CI_table\F1_CItab','npdCi','npGCci')
% load('C:\Users\Tim\Documents\Work\GIT\NPD_Validate\precomp_CI_table\F1A_CItab','npdCi','npGCci')