function [] = wrapper_Fig8_IncompleteCond(C,NCV,NC)
% Simulate and analyses different routing + partialized NPD

% Sweep across SNRs of hidden node
NCvec = logspace(log10(0.001),log10(2000),NC);
NCvec = sqrt(NCvec); % convert to std


nc_col_sc = [1 0.9 0.8];
cmap = linspecer(4);
lsstyles = {'-','-.',':'};

% Run the MVAR Simulation (fieldtrip implementation uses BSMART)
rng(95342)
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
%         perm = 1;
%         permtype = 2;
        % Can switch off if you have run once to save the results see
        % L75-76
        perm = 0;
        permtype = 0;        
    else
        perm = 0;
        permtype = 0;
    end
    
    % Plotting parameters
    plotfig =0;
    linestyle = lsstyles{1};
    cmapn = cmap.*nc_col_sc(1);
    
    %% Modify SNR of hidden node
    data = X;
    randproc = randn(size(data.trial{1}));
    for i = 2
        s = data.trial{1}(i,:);
        s = (s-mean(s))./std(s);
        n = ((NCvec(ncov)*1).*randproc(i,:));
        y = s + n;
        snr = var(s)/var(n);
        snrbank(ncov,i) = snr;
        data.trial{1}(i,:) = y;
    end
   
    %% Compute Power with Fieldtrip (used for NPD and NPG computations)
    datalength = (2^segOrd)./cfg.fsample;
    freq = computeSpectra(data,[0 0 0],Nsig,plotfig,linestyle,-1,datalength);
    pow = mean(abs(squeeze(freq.fourierspctrm(:,1,:))),1); % Convert to power
    npPow(ncov) = max(pow(freq.freq>42 & freq.freq<62));
    
    %% Compute NPD with Neurospec
    coh = [];
    [Hz lags npdspctrm npdspctrmZ npdspctrmW nscohspctrm npdcrcv] = ft_computeNPD(freq,cfg.fsample,2 ,segOrd,perm,permtype);
    coh.freq= Hz; coh.cohspctrm = nscohspctrm{1};
    nscoh(ncov) = max(nscohspctrm{1}(3,1,Hz>42 & Hz<62)); %-max(nscohspctrm(1,2,Hz>42 & Hz<62));
    npd(ncov) = max(npdspctrm{1,2}(3,1,Hz>42 & Hz<62)); %- max(npdspctrm{3}(2,1,Hz>42 & Hz<62));
    npdz(ncov) = max(npdspctrmZ{1,2}(3,1,Hz>42 & Hz<62)); %- max(npdspctrm{3}(2,1,Hz>42 & Hz<62));
    npdCi(ncov) = mean(npdspctrm{2,2}(3,1,:));
    
    %% Compute NPG
    [Hz granger grangerft] = computeGranger(freq,1,perm,permtype)
    npGCmv(ncov) = max(granger{1,2}(3,1,grangerft.freq>42 & grangerft.freq<62));
    npGCci(ncov) = mean(granger{2,2}(3,1,:));
end
% save('C:\Users\Tim\Documents\Work\GIT\NPD_Validate\precomp_CI_table\F7_CItab','npdCi','npGCci')
load('C:\Users\Tim\Documents\Work\GIT\NPD_Validate\precomp_CI_table\F7_CItab','npdCi','npGCci')


SNRDB = 10*log10(snrbank(:,2));
A = SNRDB;

% scatter(NCvec,npPow,40,cmapn(1,:),'filled')
% scatter(NCvec(1,:),nscoh,40,cmapn(1,:),'filled')
hold on
scatter(A,npd,40,cmapn(3,:),'Marker','+','LineWidth',3)
plot(A,repmat(npdCi(1),1,size(A,1)),'LineStyle','--','color',cmapn(3,:))

scatter(A,npGCmv,40,cmapn(2,:),'filled')
plot(A,repmat(npGCci(1),1,size(A,1)),'LineStyle','--','color',cmapn(2,:))
scatter(A,npdz,40,cmapn(4,:),'filled')
% scatter(A,npGCmv,40,cmapn(1,:),'filled')
grid on
xlabel('X2 SNR');ylabel('FC X_1 \rightarrow X_3')
legend({'NPD','Granger','NPDx2'},'Location','SouthEast')
title('Effects of Incomplete Signals for Conditioning')

