function [] = wrapper_Fig9C_sigmix(X,NC)
% Simulates instantaneous signal mixing in empirical data
% SNRs

% Setup Mixings
NCvec = [0 0.075 0.15];
nc_col_sc = [1 0.9 0.8];

% plotting
cmap = linspecer(4);
lsstyles = {'-','-.',':'};

Nsig = 2;
segOrd = 8; % 2^n length of segment used for FFT
for ncov = 1:NC
    if ncov == 1
        perm = 1;
        permtype = 2;
    else
        perm = 0;
        permtype = 0;
    end
    
    % plotting parameters
    cmapn = cmap.*nc_col_sc(ncov);
    linestyle = lsstyles{ncov};
    plotfig =1;
    
    %% Alter Signal Mixing
    data = X;
    sigmix = repmat(NCvec(ncov)/(Nsig-1),Nsig,Nsig).*~eye(Nsig);
    sigmix = sigmix+eye(Nsig).*1; %(1-NCvec(ncov));
    data.trial{1} = (data.trial{1} -mean(data.trial{1},2))./std(data.trial{1},[],2);
    data.trial{1} = sigmix*data.trial{1};
    shvar(:,:,ncov) = corr(data.trial{1}');    

    %% Compute Power with Fieldtrip (used for NPD and NPG computations)
    datalength = (2^segOrd)./cfg.fsample;
    freq = computeSpectra(data,[0 0 0],Nsig,plotfig,linestyle,-1,datalength);
    
    %% NPD
    [Hz lags npdspctrm npdspctrmZ npdspctrmW nscohspctrm npdcrcv] = ft_computeNPD(freq,cfg.fsample,1,segOrd,perm,permtype);
    coh.freq= Hz; coh.cohspctrm = nscohspctrm{1}; coh.ci = nscohspctrm{2};
    
    %% GRANGER
    [gHz granger grangerft] = computeGranger(freq,0,perm,permtype);
    
    plotSTN_M2_data(Hz,npdspctrm,data,cmapn(3,:),1,linestyle,bstrp)
    plotSTN_M2_data(gHz,granger,data,cmapn(2,:),1,linestyle,bstrp)
end
