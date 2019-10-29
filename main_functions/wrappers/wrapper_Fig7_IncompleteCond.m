function [] = wrapper_Fig7_IncompleteCond(C,NCV,NC)

NCvec = logspace(log10(0.001),log10(2000),NC);
NCvec = sqrt(NCvec); % convert to std

% NCvec = [0 0.5 1.5];
nc_col_sc = [1 0.9 0.8];


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
    if ncov == 1
        bstrp = 0;
    else
        bstrp = 0;
    end
    rng(124321)
    plotfig =0;
    linestyle = lsstyles{1};
    cmapn = cmap.*nc_col_sc(1);
    ncv = NCvec(ncov);
    data = X;
    randproc = randn(size(data.trial{1}));
    for i = 2
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
    
    
    %% Plot Example Trace
    %     figure(1)
    %     plot(X.time{1},X.trial{1}+[0 2.5 5]');
    %     xlim([100 105])
    %     xlabel('Time'); ylabel('Amplitude'); grid on
    %     legend({'X1','X2','X3'})
    
    freq = computeSpectra(data,[0 0 0],Nsig,plotfig,linestyle);
    pow = mean(abs(squeeze(freq.fourierspctrm(:,1,:))),1);
    npPow(ncov) = max(pow(freq.freq>42 & freq.freq<62));
    
    %% NPD
    coh = [];
    [Hz lags npdspctrm npdspctrmZ npdspctrmW nscohspctrm npdcrcv] = computeNPD(data,2,8,1,bstrp);
    coh.freq= Hz; coh.cohspctrm = nscohspctrm{1};
    nscoh(ncov) = max(nscohspctrm{1}(3,1,Hz>42 & Hz<62)); %-max(nscohspctrm(1,2,Hz>42 & Hz<62));
    npd(ncov) = max(npdspctrm{1,2}(3,1,Hz>42 & Hz<62)); %- max(npdspctrm{3}(2,1,Hz>42 & Hz<62));
    npdz(ncov) = max(npdspctrmZ{1,2}(3,1,Hz>42 & Hz<62)); %- max(npdspctrm{3}(2,1,Hz>42 & Hz<62));
    npdCi(ncov) = mean(npdspctrm{2,2}(3,1,:));
    
    % NS Coh
    plotNSCoherence(coh,cmapn(1,:),Nsig,plotfig,linestyle)
    % NPD
    %     plotNPD(Hz,npdspctrm,data,cmapn(3,:),1,linestyle)
    %     plotNPD(Hz,npdspctrmW,data,cmapn(4,:),1,linestyle)
    %
    %% GRANGER
    %     [Hz granger grangerft] = computeGranger(freq,cmapn(2,:),Nsig,plotfig,linestyle,0);
    %     npGCpw(ncov) = max(granger{1,2}(3,1,grangerft.freq>42 & grangerft.freq<62));
    [Hz granger grangerft] = computeGranger(freq,cmapn(1,:),Nsig,plotfig,linestyle,1,bstrp);
    npGCmv(ncov) = max(granger{1,2}(3,1,grangerft.freq>42 & grangerft.freq<62));
    npGCci(ncov) = mean(granger{2,2}(3,1,:));
    %     %% NPD CORR
    %         figure(3)
    %         plotNPDXCorr(lags,npdcrcv,data,cmapn(3,:),0,linestyle)
    %         plotNPDXCorr(lags,npdcrcvW,data,cmapn(4,:),0,linestyle)
    a =1;
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

% cfg = [];
% cfg.viewmode = 'butterfly';  % you can also specify 'butterfly'
% ft_databrowser(cfg, data);
