function [] = wrapper_Fig9A_symSNR(X,NC,permrun,fresh)
% Simulates symmetrical SNR changes in empirical data
% SNRs
% Setup SNRs to test
NCvec = [0 0.75 3];
NCvec = sqrt(NCvec); % convert to std
nc_col_sc = [1 0.9 0.8];
NCbase = 1;
% plotting
cmap = linspecer(4);
lsstyles = {'-','-.',':'};

Nsig = 2;
segOrd = 8; % 2^n length of segment used for FFT
for ncov = 1:NC
    if ncov == NCbase
        if fresh == 1
            perm = 1;
            permtype = permrun;
        elseif fresh == 0
            % Can switch off if you have run once to save the results see
            % L75-76
            perm = 0;
            permtype = 0;
        end
    else
        perm = 0;
        permtype = 0;
    end
    
    % plotting parameters
    linestyle = lsstyles{ncov};
    cmapn = cmap.*nc_col_sc(ncov);
    plotfig =0;
    
    %% Alter SNR of Data
    data = X;
    randproc = randn(size(data.trial{1}));
    for i = 1:size(randproc,1)
        s = data.trial{1}(i,:);
        s = (s-mean(s))./std(s);
        n = ((NCvec(ncov)*1).*randproc(i,:));
        y = s + n;
        snr = var(s)/var(n);
        snrbank(ncov,i) = snr;
        data.trial{1}(i,:) = y;
    end
    
    %% Compute Power with Fieldtrip (used for NPD and NPG computations)
    datalength = (2^segOrd)./data.fsample;
    freq = computeSpectra(data,[0 0 0],Nsig,plotfig,linestyle,2,datalength);
    
    %% NPD
    [Hz lags npdspctrm npdspctrmZ npdspctrmW nscohspctrm npdcrcv] = ft_computeNPD(freq,data.fsample,1,segOrd,perm,permtype);
    coh.freq= Hz; coh.cohspctrm = nscohspctrm{1}; coh.ci = nscohspctrm{2};
    
    %% GRANGER
    [gHz granger grangerft] = computeGranger(freq,0,perm,permtype);
    plotSTN_M2_data(Hz,npdspctrm,data,cmapn(3,:),1,linestyle,perm)
    plotSTN_M2_data(gHz,granger,data,cmapn(2,:),1,linestyle,perm)
end

