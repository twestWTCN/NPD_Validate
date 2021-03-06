function [] = wrapper_Fig5A_SymSNR(C,NCV,NC)
NCvec = logspace(log10(0.001),log10(1000),NC);
% NCvec = linspace(0.001,10,NC);
NCvec = sqrt(NCvec); % Convert to std

for i = 1:NC
    NCtits{i} = ['SigLeak = ' num2str(NCvec(i))];
end
rng(12312)
cmap = linspecer(4);
lsstyles = {'-','-.',':'};
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


for ncov = 1:NC
    if ncov == floor(NC/2)
        bstrp = 1;
    else
        bstrp = 0;
    end
    disp(ncov)
    data = X;
        plotfig =0;
    linestyle = '-';
    cmapn = cmap;
    ncv = NCvec(ncov);
    
    %% Compute SNR
    randproc = randn(size(data.trial{1}));
    for i = 1:size(randproc,1)
        s = X.trial{1}(i,:);
        s = (s-mean(s))./std(s);
        n = ((NCvec(ncov)*1).*randproc(i,:));
        %         si = filter(bfilt,afilt,s);
        %         ni = filter(bfilt,afilt,n);
        y = s + n;
        snr = var(s)/var(n);
        snrbank(ncov,i) = snr;
        %         snrbp = var(si)/var(ni);
        snrbp = computeBandLimSNR(s,n,[45 55],data);
        snrbpbank(ncov,i) = snrbp;
        data.trial{1}(i,:) = y;
    end
    
    %% Power
    freq = computeSpectra(data,[0 0 0],Nsig,plotfig,linestyle);
    pow = mean(abs(squeeze(freq.fourierspctrm(:,1,:))),1);
    npPow(ncov) = max(pow(freq.freq>42 & freq.freq<62));

    %% NPD
    coh = [];
    [Hz lags npdspctrm npdspctrmZ npdspctrmW nscohspctrm npdcrcv] = computeNPD(data,1,8,1,bstrp);
    coh.freq= Hz; coh.cohspctrm = nscohspctrm{1};
    nscoh(ncov) = max(nscohspctrm{1}(1,2,Hz>42 & Hz<62));
    npd(ncov) = max(npdspctrm{1,3}(1,2,Hz>42 & Hz<62));
    npdCi(ncov) = mean(npdspctrm{2,2}(1,2,:));
    % NS Coh
    plotNSCoherence(coh,cmapn(1,:),Nsig,plotfig,linestyle)
    % NPD
    plotNPD(Hz,npdspctrm,data,cmapn(3,:),plotfig,linestyle)

    %% GRANGER
    [Hz granger grangerft] = computeGranger(freq,cmapn(2,:),Nsig,plotfig,linestyle,1,bstrp);
    npGC(ncov) = max(granger{1,3}(1,2,grangerft.freq>42 & grangerft.freq<62));
    npGCci(ncov) = mean(granger{2,3}(1,2,:));
end
save('C:\Users\Tim\Documents\Work\GIT\NPD_Validate\precomp_CI_table\F5A_CItab','npdCi','npGCci')
load('C:\Users\Tim\Documents\Work\GIT\NPD_Validate\precomp_CI_table\F5A_CItab','npdCi','npGCci')

% SNRDB = 10*log10(snrbank(:,1));
SNRDB = 10*log10(snrbpbank(:,1));

A = SNRDB;
% scatter(NCvec,npPow,40,cmapn(1,:),'filled')
pa(1) = scatter(A,nscoh,50,cmapn(1,:),'Marker','+','LineWidth',3);;
hold on
% [param,stat]=sigm_fit(A,nscoh)
% hold on
pa(2) = scatter(A,npd,40,cmapn(3,:),'filled');
plot(A,repmat(npdCi(floor(NC/2)),1,size(A,1)),'LineStyle','--','color',cmapn(3,:))

% [param,stat]=sigm_fit(A,npd)
pa(3) = scatter(A,npGC,40,cmapn(2,:),'filled');
plot(A,repmat(npGCci(floor(NC/2)),1,size(A,1)),'LineStyle','--','color',cmapn(2,:))

% [param,stat]=sigm_fit(A,npGC); %,[],[.0043 0.5884 0.689 -3.76])
  grid on
xlabel('SNR_{dB}');ylabel('FC Magnitude')
legend(pa,{'Coherence','NPD','Granger'})
title('FC vs SNR')
xlim([-inf inf])