function [] = mvarconsim_npdver_F3b_asysm(C,NCV,NC)

% NCvec = linspace(0,2,NC);
NCvec = logspace(log10(0.01),log10(10),NC);
NCvec(2:3,:) = repmat(0.01,2,NC);
nc_col_sc = [1 0.9 0.8];
for i = 1:NC
    NCtits{i} = ['SigLeak = ' num2str(NCvec(i))];
end

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
data              = ft_connectivitysimulation(cfg);


for ncov = 1:NC
    disp(ncov)
    X = data;
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
    freq = computeSpectra(X,[0 0 0],Nsig,plotfig,linestyle);
    pow = mean(abs(squeeze(freq.fourierspctrm(:,1,:))),1);
    npPow(ncov) = max(pow(freq.freq>42 & freq.freq<62));

    %% NPD
    coh = [];
    [Hz lags npdspctrm npdspctrmZ npdspctrmW nscohspctrm npdcrcv] = computeNPD(X,1);
    coh.freq= Hz; coh.cohspctrm = nscohspctrm;
    nscoh(ncov) = max(nscohspctrm(2,1,Hz>42 & Hz<62))-max(nscohspctrm(1,2,Hz>42 & Hz<62));
    npd(ncov) = max(npdspctrm{2}(2,1,Hz>42 & Hz<62))- max(npdspctrm{3}(2,1,Hz>42 & Hz<62));
    
    % NS Coh
    computeNSCoherence(coh,cmapn(1,:),Nsig,plotfig,linestyle)
    % NPD
%     plotNPD(Hz,npdspctrm,data,cmapn(3,:),plotfig,linestyle)

    %% GRANGER
    [Hz granger grangerft] = computeGranger(freq,cmapn(2,:),Nsig,plotfig,linestyle);
    npGC(ncov) = max(grangerft.grangerspctrm(2,1,grangerft.freq>42 & grangerft.freq<62))-max(grangerft.grangerspctrm(1,2,grangerft.freq>42 & grangerft.freq<62));
end

a =1;
figure(2)
% snrbank 
% scatter(NCvec,npPow,40,cmapn(1,:),'filled')
SNRDB = 10*log10(snrbank(:,1));
SNRBASE = 10*log10(snrbank(:,2));
A = SNRDB./SNRBASE;
scatter(A,nscoh,50,cmapn(1,:),'Marker','+','LineWidth',3);
% [param,stat]=sigm_fit(A,nscoh)
hold on
scatter(A,npd,40,cmapn(3,:),'filled')
% [param,stat]=sigm_fit(A,npd)
scatter(A,npGC,40,cmapn(2,:),'filled')
% [param,stat]=sigm_fit(A,npGC)

grid on
xlabel('SAsym');ylabel('FC Difference')
legend({'Coherence','NPD','Granger'})
title('FC vs SNR Asymetry')
xlim([-0.5 1])