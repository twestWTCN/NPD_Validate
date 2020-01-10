function [] = wrapper_Fig5B_AsymSNR(C,NCV,NC)
% Sweep across asymmetrical SNRs
NCvec = repmat(0.001,2,NC);
NCvec(1,1:floor(NC/2)) = logspace(log10(1000),log10(0.01),floor(NC/2));
NCvec(2,ceil(NC/2)+1:end) = logspace(log10(0.01),log10(1000),floor(NC/2));
NCvec = sqrt(NCvec); % convert to std

cmap = linspecer(4);

% Run the MVAR Simulation (fieldtrip implementation uses BSMART)
rng(2321)
Nsig = size(C,1);
cfg             = [];
cfg.ntrials     = 1;
cfg.triallength = 500;
cfg.fsample     = 200;
cfg.nsignal     = Nsig;
cfg.method      = 'ar';
cfg.params = C;
cfg.noisecov = NCV;
X              = ft_connectivitysimulation(cfg);
segOrd = 8; % 2^n length of segment used for FFT

for ncov = 1:NC
    disp(ncov)
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
    linestyle = '-';
    cmapn = cmap;
    plotfig =0;
    
    %% Compute Asymmetrical SNR
    data = X;
    randproc = randn(size(data.trial{1}));
    for i = 1:size(randproc,1)
        s = data.trial{1}(i,:);
        s = (s-mean(s))./std(s);
        n  =((NCvec(i,ncov)*1).*randproc(i,:));
        y = s+n;
        snr = var(s)/var(n);
        disp(snr)
        snrbank(ncov,i) = snr;
        %         y = (y-mean(y))./std(y);
        data.trial{1}(i,:) = y;
    end
    
    %% Compute Power with Fieldtrip (used for NPD and NPG computations)
    datalength = (2^segOrd)./cfg.fsample;
    freq = computeSpectra(data,[0 0 0],Nsig,plotfig,linestyle,1,datalength);
    pow = mean(abs(squeeze(freq.fourierspctrm(:,1,:))),1);
    npPow(ncov) = max(pow(freq.freq>42 & freq.freq<62));
    
    %% Compute NPD with Neurospec
    coh = [];
    [Hz lags npdspctrm npdspctrmZ npdspctrmW nscohspctrm npdcrcv] = ft_computeNPD(freq,cfg.fsample,1,segOrd,perm,permtype);
    coh.freq= Hz; coh.cohspctrm = nscohspctrm{1};
    nscoh(ncov) = max(nscohspctrm{1}(1,2,Hz>42 & Hz<62))-max(nscohspctrm{1}(1,2,Hz>42 & Hz<62));
    npd(ncov) = max(npdspctrm{1,2}(1,2,Hz>42 & Hz<62))- max(npdspctrm{1,2}(2,1,Hz>42 & Hz<62));
    npdCi(ncov) = mean(npdspctrm{2,2}(1,2,:));
    
    %% Compute NPG
    [Hz granger grangerft] = computeGranger(freq,0,perm,permtype)
    npGC(ncov) = max(granger{1,2}(1,2,grangerft.freq>42 & grangerft.freq<62))-max(granger{1,2}(2,1,grangerft.freq>42 & grangerft.freq<62));
    npGCci(ncov) = mean(granger{2,2}(1,2,:));
end
% save('C:\Users\Tim\Documents\Work\GIT\NPD_Validate\precomp_CI_table\F5B_CItab','npdCi','npGCci')
load('C:\Users\Tim\Documents\Work\GIT\NPD_Validate\precomp_CI_table\F5B_CItab','npdCi','npGCci')

figure(2)
% snrbank
% scatter(NCvec,npPow,40,cmapn(1,:),'filled')
SNRDB = 10*log10(snrbank(:,1));
SNRBASE = 10*log10(snrbank(:,2));
A = SNRDB-SNRBASE;
scatter(A,nscoh,50,cmapn(1,:),'Marker','+','LineWidth',3);

hold on
scatter(A,npd,40,cmapn(3,:),'filled')
plot(A,repmat(npdCi(1),1,size(A,1)),'LineStyle','--','color',cmapn(3,:))
plot(A,-repmat(npdCi(1),1,size(A,1)),'LineStyle','--','color',cmapn(3,:))

scatter(A,npGC,40,cmapn(2,:),'filled')
plot(A,repmat(npGCci(1),1,size(A,1)),'LineStyle','--','color',cmapn(2,:))
plot(A,-repmat(npGCci(1),1,size(A,1)),'LineStyle','--','color',cmapn(2,:))

grid on
xlabel('SAsym');ylabel('FC Difference')
legend({'Coherence','NPD','Granger'})
title('FC vs SNR Asymetry')
xlim([-60 60])
% xlim([-0.5 1])