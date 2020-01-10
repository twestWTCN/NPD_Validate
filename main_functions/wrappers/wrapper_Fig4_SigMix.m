function [] = wrapper_Fig4_SigMix(C,NCV,NC)
% Signal Mixings
NCvec = [0 0.45 1.2];
nc_col_sc = [1 0.9 0.8];

cmap = linspecer(4);
lsstyles = {'-','-.',':'};

% Run the MVAR Simulation (fieldtrip implementation uses BSMART)
rng(2131)
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
    sigmix = repmat(NCvec(ncov)/(Nsig-1),Nsig,Nsig).*~eye(Nsig);
    sigmix = sigmix+eye(Nsig).*1; %(1-NCvec(ncov));
    data.trial{1} = (data.trial{1} -mean(data.trial{1},2))./std(data.trial{1},[],2);
    data.trial{1} = sigmix*data.trial{1};
    
    dumdata = randn(size(data.trial{1}));
    shvar(:,:,ncov) = corr((sigmix*dumdata)');

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
    % Plot NPD
    plotNPD(Hz,npdspctrm,data,cmapn(3,:),plotfig,linestyle,perm)
    
    %% Compute and plot GRANGER
    figure(2)
    [Hz granger grangerft] = computeGranger(freq,0,perm,permtype)
    plotNPD(grangerft.freq,granger,freq,cmapn(2,:),1,linestyle,perm)
   
    %% NPD CORR
    figure(3)
    plotNPDXCorr(lags,npdcrcv,data,[0 0 0],0,linestyle)
end

a = 1;
