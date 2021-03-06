function [] = mvarconsim_npdver_F6r(C,NCV,NC)

NCvec = logspace(log10(0.01),log10(10),NC);
% NCvec = [0 0.5 1.5];
nc_col_sc = [1 0.9 0.8];


cmap = linspecer(4);
lsstyles = {'-','-.',':'};
Nsig = size(C,1);
cfg             = [];
cfg.ntrials     = 1;
cfg.triallength = 500;
cfg.fsample     = 100;
cfg.nsignal     = Nsig;
cfg.method      = 'ar';
cfg.params = C;
cfg.noisecov = NCV;
data              = ft_connectivitysimulation(cfg);


for ncov = 1:NC
rng(124321)
        plotfig =0;
    linestyle = lsstyles{1};
    cmapn = cmap.*nc_col_sc(1);
    ncv = NCvec(ncov);
    X = data;
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


    %% Plot Example Trace
%     figure(1)
%     plot(X.time{1},X.trial{1}+[0 2.5 5]');
%     xlim([100 105])
%     xlabel('Time'); ylabel('Amplitude'); grid on
%     legend({'X1','X2','X3'})
    
    freq = computeSpectra(X,[0 0 0],Nsig,plotfig,linestyle);
    pow = mean(abs(squeeze(freq.fourierspctrm(:,1,:))),1);
    npPow(ncov) = max(pow(freq.freq>42 & freq.freq<62));

    %% NPD
    coh = [];
    [Hz lags npdspctrm npdspctrmZ npdspctrmW nscohspctrm npdcrcv] = computeNPD(X,1);
    coh.freq= Hz; coh.cohspctrm = nscohspctrm;
    nscoh(ncov) = max(nscohspctrm(3,1,Hz>16 & Hz<22)); %-max(nscohspctrm(1,2,Hz>42 & Hz<62));
    npd(ncov) = max(npdspctrm{2}(3,1,Hz>16 & Hz<22)); %- max(npdspctrm{3}(2,1,Hz>42 & Hz<62));
    npdz(ncov) = max(npdspctrmW{2}(3,1,Hz>16 & Hz<22)); %- max(npdspctrm{3}(2,1,Hz>42 & Hz<62));
    % NS Coh
    computeNSCoherence(coh,cmapn(1,:),Nsig,plotfig,linestyle)
    % NPD
%     plotNPD(Hz,npdspctrm,data,cmapn(3,:),1,linestyle)
%     plotNPD(Hz,npdspctrmW,data,cmapn(4,:),1,linestyle)
% 
    %% GRANGER
    [dum1 dum2 grangerft] = computeGranger(freq,cmapn(2,:),Nsig,plotfig,linestyle);
    npGC(ncov) = max(grangerft.grangerspctrm(1,3,grangerft.freq>16 & grangerft.freq<22)); %-max(granger.grangerspctrm(1,2,granger.freq>42 & granger.freq<62));
    
%     %% NPD CORR
%         figure(3)
%         plotNPDXCorr(lags,npdcrcv,data,cmapn(3,:),0,linestyle)
%         plotNPDXCorr(lags,npdcrcvW,data,cmapn(4,:),0,linestyle)
    a =1;
end


SNRDB = 10*log10(snrbank(:,2));
A = SNRDB;

% scatter(NCvec,npPow,40,cmapn(1,:),'filled')
% scatter(NCvec(1,:),nscoh,40,cmapn(1,:),'filled')
hold on
scatter(A,npd,40,cmapn(3,:),'Marker','+','LineWidth',3)
scatter(A,npGC,40,cmapn(2,:),'filled')
scatter(A,npdz,40,cmapn(4,:),'filled')
grid on
xlabel('X2 SNR');ylabel('FC X_1 \rightarrow X_3')
legend({'NPD','Granger','NPDx2'},'Location','SouthEast')
title('Effects of Incomplete Signals for Conditioning')

% cfg = [];
% cfg.viewmode = 'butterfly';  % you can also specify 'butterfly'
% ft_databrowser(cfg, data);
