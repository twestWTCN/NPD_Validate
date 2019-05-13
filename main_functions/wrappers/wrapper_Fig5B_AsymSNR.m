function [] = wrapper_Fig5B_AsymSNR(C,NCV,NC)

% NCvec = linspace(0.01,50,NC);
NCvec = repmat(0.001,2,NC);
NCvec(1,1:floor(NC/2)) = logspace(log10(1000),log10(0.01),floor(NC/2));
NCvec(2,ceil(NC/2)+1:end) = logspace(log10(0.01),log10(1000),floor(NC/2));
NCvec = sqrt(NCvec); % convert to std

nc_col_sc = [1 0.9 0.8];
for i = 1:NC
    NCtits{i} = ['SigLeak = ' num2str(NCvec(i))];
end
rng(2321)

cmap = linspecer(4);
lsstyles = {'-','-.',':'};
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

for ncov = 1:NC
    if ncov == 1
        bstrp = 0;
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
    
    %% Power
    freq = computeSpectra(data,[0 0 0],Nsig,plotfig,linestyle);
    pow = mean(abs(squeeze(freq.fourierspctrm(:,1,:))),1);
    npPow(ncov) = max(pow(freq.freq>42 & freq.freq<62));
    
    %% NPD
    coh = [];
    [Hz lags npdspctrm npdspctrmZ npdspctrmW nscohspctrm npdcrcv] = computeNPD(data,1,8,1,bstrp);
    coh.freq= Hz; coh.cohspctrm = nscohspctrm{1};
    nscoh(ncov) = max(nscohspctrm{1}(1,2,Hz>42 & Hz<62))-max(nscohspctrm{1}(1,2,Hz>42 & Hz<62));
    npd(ncov) = max(npdspctrm{1,2}(1,2,Hz>42 & Hz<62))- max(npdspctrm{1,2}(2,1,Hz>42 & Hz<62));
    npdCi(ncov) = mean(npdspctrm{2,2}(1,2,:));
    % NS Coh
    plotNSCoherence(coh,cmapn(1,:),Nsig,plotfig,linestyle)
    % NPD
    %     plotNPD(Hz,npdspctrm,data,cmapn(3,:),plotfig,linestyle)
    
    %% GRANGER
    [Hz granger grangerft] = computeGranger(freq,cmapn(2,:),Nsig,plotfig,linestyle,0,bstrp);
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