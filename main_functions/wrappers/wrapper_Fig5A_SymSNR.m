function [] = wrapper_Fig5A_SymSNR(C,NCV,NC)
NCvec = logspace(log10(0.005),log10(100),NC);
% NCvec = linspace(0.001,10,NC);
NCvec = sqrt(NCvec); % Convert to std

for i = 1:NC
    NCtits{i} = ['SigLeak = ' num2str(NCvec(i))];
end

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
    disp(ncov)
    data = X;
        plotfig =0;
    linestyle = '-';
    cmapn = cmap;
    ncv = NCvec(ncov);
    
    %% Compute SNR
    randproc = randn(size(data.trial{1}));
    for i = 1:size(randproc,1)
        s = data.trial{1}(i,:);
        s = (s-mean(s))./std(s);
        n = (ncv.*randproc(i,:));
        y = s + n;
        snr = var(s)/var(n);
        snrbank(ncov,i) = snr;
        data.trial{1}(i,:) = y;
    end
    
    %% Power
    freq = computeSpectra(data,[0 0 0],Nsig,plotfig,linestyle);
    pow = mean(abs(squeeze(freq.fourierspctrm(:,1,:))),1);
    npPow(ncov) = max(pow(freq.freq>42 & freq.freq<62));

    %% NPD
    coh = [];
    [Hz lags npdspctrm npdspctrmZ npdspctrmW nscohspctrm npdcrcv] = computeNPD(data,1);
    coh.freq= Hz; coh.cohspctrm = nscohspctrm{1};
    nscoh(ncov) = max(nscohspctrm{1}(2,1,Hz>42 & Hz<62));
    npd(ncov) = max(npdspctrm{3}(2,1,Hz>42 & Hz<62));
    
    % NS Coh
    plotNSCoherence(coh,cmapn(1,:),Nsig,plotfig,linestyle)
    % NPD
    plotNPD(Hz,npdspctrm,data,cmapn(3,:),plotfig,linestyle)

    %% GRANGER
    [Hz granger grangerft] = computeGranger(freq,cmapn(2,:),Nsig,plotfig,linestyle,0);
    npGC(ncov) = 0; %max(granger{1,3}(1,2,grangerft.freq>42 & grangerft.freq<62));
end

SNRDB = 10*log10(snrbank(:,1));

A = SNRDB;
% scatter(NCvec,npPow,40,cmapn(1,:),'filled')
scatter(A,nscoh,50,cmapn(1,:),'Marker','+','LineWidth',3);
hold on
% [param,stat]=sigm_fit(A,nscoh)
% hold on
scatter(A,npd,40,cmapn(3,:),'filled')
% [param,stat]=sigm_fit(A,npd)
scatter(A,npGC,40,cmapn(2,:),'filled')
% [param,stat]=sigm_fit(A,npGC); %,[],[.0043 0.5884 0.689 -3.76])
grid on
xlabel('SNR_{dB}');ylabel('FC Magnitude')
legend({'Coherence','NPD','Granger'})
title('FC vs SNR')
xlim([-inf inf])